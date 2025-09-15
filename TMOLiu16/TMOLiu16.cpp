/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                       Author: Jakub Krystufek                                     *
*                       Brno 2024                                                   *
*                                                                                   *
*                       Implementation of the SPDecolor method                      *
*                       Credits to Matlab implementation                            *
*                       https://github.com/yqx7150/SPDecolor                        *
*                                                                                   *
************************************************************************************/

/* --------------------------------------------------------------------------- *
 * --------------------------------------------------------------------------- */

#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include "TMOLiu16.h"

TMOLiu16::TMOLiu16()
{
   SetName(L"Liu16");
   SetDescription(L"Semiparametric Decolorization With Laplacian-Based Perceptual Quality Metric");

   sigmaParameter.SetName(L"Sigma");
   sigmaParameter.SetDescription(L"Sigma");
   sigmaParameter.SetDefault(0.005);
   sigmaParameter = 0.005;
   sigmaParameter.SetRange(0.0, 1.0);
   this->Register(sigmaParameter);
}

TMOLiu16::~TMOLiu16()
{
}

cv::Mat computeHorizontalGradient(const cv::Mat& img) {
   cv::Mat Rx(cv::Size(img.cols-1, img.rows), CV_64F, cv::Scalar(0));
   for (int i = 0; i<img.rows; i++) {
      for (int j = 0; j<img.cols-1; j++) {
         Rx.at<double>(i, j) = img.at<double>(i, j)-img.at<double>(i, j+1);
      }
   }
   return Rx;
}

cv::Mat computeVerticalGradient(const cv::Mat& img) {
   cv::Mat Ry(cv::Size(img.cols, img.rows-1), CV_64F, cv::Scalar(0));
   for (int i = 0; i<img.rows-1; i++) {
      for (int j = 0; j<img.cols; j++) {
         Ry.at<double>(i, j) = img.at<double>(i, j)-img.at<double>(i+1, j);
      }
   }
   return Ry;
}

cv::Mat weak_order(const cv::Mat& polyGrad) {
   double level = 0.05;

   cv::Mat Rg = polyGrad.col(0);
   cv::Mat Gg = polyGrad.col(2);
   cv::Mat Bg = polyGrad.col(5);

   // Calculate positive weak order
   cv::Mat posWeakOrder = cv::Mat::zeros(Rg.size(), Rg.type());

   for (int i = 0; i<Rg.rows; i++) {
      if (Rg.at<double>(i)>level&&Gg.at<double>(i)>level&&Bg.at<double>(i)>level) {
         posWeakOrder.at<double>(i) = 1;
      }
   }

   cv::Mat negWeakOrder = cv::Mat::zeros(Rg.size(), Rg.type());
   for (int i = 0; i<Rg.rows; i++) {
      if (Rg.at<double>(i)<-level&&Gg.at<double>(i)<-level&&Bg.at<double>(i)<-level) {
         negWeakOrder.at<double>(i) = 1;
      }
   }

   posWeakOrder.convertTo(posWeakOrder, CV_64F);
   negWeakOrder.convertTo(negWeakOrder, CV_64F);

   cv::Mat alf = posWeakOrder-negWeakOrder;

   return alf;
}

cv::Mat TMOLiu16::calculatePolyGrad(const cv::Mat& img, int order) {
   // Resize image to roughly constant pixel count
   double factor = 64.0/std::sqrt(img.rows*img.cols);

   cv::Mat resized;
   cv::resize(img, resized, cv::Size(), factor, factor, cv::INTER_NEAREST);

   // Split color channels of resized image
   std::vector<cv::Mat> channels(3);
   cv::split(resized, channels);

   cv::Mat ImR = channels[0];
   cv::Mat ImG = channels[1];
   cv::Mat ImB = channels[2];

   // Create feature vector matrix
   std::vector<cv::Mat> featuresMatrices = {
       ImR,             // R
       ImR.mul(ImR),    // R^2
       ImG,             // G
       ImG.mul(ImR),    // G*R
       ImG.mul(ImG),    // G^2
       ImB,             // B
       ImB.mul(ImR),    // B*R
       ImB.mul(ImG),    // B*G
       ImB.mul(ImB)     // B^2
   };

   // Transform featuresMatrices to a single row matrix
   cv::Mat features;
   for (int i = 0; i<featuresMatrices.size(); i++) {
      featuresMatrices[i] = featuresMatrices[i].reshape(0, 1);
      features.push_back(featuresMatrices[i]);
   }

   // Random permutation of the feature matrix rows
   cv::Mat featureMatrix = features.t();

   std::vector<int> indices(featureMatrix.rows);
   std::iota(indices.begin(), indices.end(), 0);
   cv::randShuffle(indices);

   // Compute initial polyGrad by subtraction
   cv::Mat polyGrad = cv::Mat::zeros(featureMatrix.size(), featureMatrix.type());
   for (int i = 0; i<featureMatrix.rows; i++) {
      featureMatrix.row(i).copyTo(polyGrad.row(indices[i]));
   }
   polyGrad = featureMatrix-polyGrad;

   // Split color channels of original image
   cv::split(img, channels);
   ImR = channels[0];
   ImG = channels[1];
   ImB = channels[2];

   cv::Mat P1;
   // Compute polynomial gradients for each combination
   for (int r = 0; r<=order; ++r) {
      for (int g = 0; g<=order; ++g) {
         for (int b = 0; b<=order; ++b) {
            if ((r+g+b)<=order&&(r+g+b)>0) {
               cv::Mat curIm = cv::Mat::ones(img.size(), CV_64F);

               cv::pow(ImR, r, curIm);

               // If channel is zero we can save some computation time
               if (g>0) {
                  cv::Mat temp;
                  cv::pow(ImG, g, temp);
                  curIm = curIm.mul(temp);
               }

               // If channel is zero we can save some computation time
               if (b>0) {
                  cv::Mat temp;
                  cv::pow(ImB, b, temp);
                  curIm = curIm.mul(temp);
               }

               // Compute gradients
               cv::Mat gradX = computeHorizontalGradient(curIm);
               cv::Mat gradY = computeVerticalGradient(curIm);

               gradX = gradX.reshape(0, 1).t(); // Flatten Rx to single column
               gradY = gradY.reshape(0, 1).t(); // Flatten Ry to single column

               cv::Mat tmp;
               // Merge gradients to single row
               cv::vconcat(gradX, gradY, tmp);

               // Add as new column to P1
               if (P1.empty()) {
                  P1 = tmp;
               }
               else {
                  cv::hconcat(P1, tmp, P1);
               }
            }
         }
      }
   }

   cv::vconcat(polyGrad, P1, polyGrad);

   return polyGrad;
}

cv::Mat calculateCg(const cv::Mat& polyGrad) {
   // Select columns
   std::vector<int> indices = { 0, 2, 5 };
   cv::Mat selectedRows = cv::Mat(polyGrad.rows, indices.size(), polyGrad.type());
   for (size_t i = 0; i<indices.size(); ++i) {
      polyGrad.col(indices[i]).copyTo(selectedRows.col(i));
   }

   // Square the values
   cv::Mat squared;
   cv::pow(selectedRows, 2, squared);

   // Sum the squares row-wise
   cv::Mat rowSums;
   cv::reduce(squared, rowSums, 1, cv::REDUCE_SUM, CV_64F);

   // Calculate the square root of the sums and normalize
   cv::Mat Cg;
   cv::sqrt(rowSums, Cg);
   Cg = Cg/1.41;

   return Cg;
}

cv::Mat wei_update_matrix(const cv::Mat& polyOriginal) {
   // Create a new matrix that excludes columns 0, 2, and 5 
   std::vector<int> includeRows = { 1,3,4,6,7,8 };

   cv::Mat poly(polyOriginal.rows, includeRows.size(), polyOriginal.type());
   for (size_t i = 0; i<includeRows.size(); ++i) {
      polyOriginal.col(includeRows[i]).copyTo(poly.col(i));
   }

   cv::Mat polyTransposed = poly.t();
   cv::Mat A = polyTransposed*poly;

   // Solve the linear system A * Mt = B
   cv::Mat Mt;
   cv::solve(A, polyTransposed, Mt, cv::DECOMP_SVD); // Using SVD decomposition to solve the system

   return Mt;
}

double energyCalcu(const cv::Mat& Cg, const cv::Mat& polyGrad, const cv::Mat& wei, double sigma) {
   cv::Mat polyGradWei = polyGrad*wei;

   // Calculate differences
   cv::Mat negDiff = polyGradWei-Cg;
   cv::Mat posDiff = polyGradWei+Cg;

   negDiff = negDiff.mul(negDiff)/sigma;
   posDiff = posDiff.mul(posDiff)/sigma;

   cv::exp(-negDiff, negDiff);
   cv::exp(-posDiff, posDiff);

   // Sum of the differences
   cv::Mat expSum = negDiff+posDiff;

   // Logarithmic transformation
   cv::Mat P;
   cv::log(expSum, P);
   P = -P;

   // Mean of P
   cv::Scalar meanValue = cv::mean(P);

   return meanValue[0]; 
}

cv::Mat grayImContruct(const std::vector<double>& wei, const cv::Mat& Im, int order) {
   cv::Mat channels[3];
   cv::split(Im, channels);

   cv::Mat ImR = channels[0];
   cv::Mat ImG = channels[1];
   cv::Mat ImB = channels[2];

   int n = ImR.rows;
   int m = ImR.cols;
   cv::Mat grayIm = cv::Mat::zeros(n, m, CV_64F);

   int kk = 0;
   for (int r = 0; r<=order; ++r) {
      for (int g = 0; g<=order; ++g) {
         for (int b = 0; b<=order; ++b) {
            if ((r+g+b)<=order&&(r+g+b)>0) {
               cv::Mat curIm = cv::Mat::ones(Im.size(), CV_64F);

               cv::pow(ImR, r, curIm);

               if (g>0) {
                  cv::Mat temp;
                  cv::pow(ImG, g, temp);
                  curIm = curIm.mul(temp);
               }
               if (b>0) {
                  cv::Mat temp;
                  cv::pow(ImB, b, temp);
                  curIm = curIm.mul(temp);
               }
               grayIm = grayIm+wei[kk++]*curIm;
            }
         }
      }
   }

   // Normalize grayIm to range [0, 1]
   double minVal, maxVal;
   cv::minMaxLoc(grayIm, &minVal, &maxVal);
   grayIm = (grayIm-minVal)/(maxVal-minVal);

   return grayIm;
}


int TMOLiu16::Transform()
{
   // Parameters and initial setup
   double pre_E = std::numeric_limits<double>::infinity();
   double E = 0;
   double tol = std::pow(10, -4);
   int maxIter = 15;

   pSrc->Convert(TMO_RGB, false);
   pDst->Convert(TMO_RGB, false);

   double* pSourceData = pSrc->GetData();
   double* pDestinationData = pDst->GetData();

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();

   cv::Mat sourceImg(height, width, CV_64FC3, pSourceData);

   cv::Mat polyGrad = calculatePolyGrad(sourceImg, 2);

   cv::Mat Cg = calculateCg(polyGrad);

   cv::Mat filteredPolyGrad(0, polyGrad.cols, polyGrad.type());
   cv::Mat filteredCg(0, Cg.cols, Cg.type());

   auto thr = 0.05;
   for (int i = 0; i<Cg.rows; i++) {
      if (Cg.at<double>(i)>=thr) {
         filteredPolyGrad.push_back(polyGrad.row(i));
         filteredCg.push_back(Cg.row(i));
      }
   }

   cv::Mat alf = weak_order(filteredPolyGrad);

   cv::Mat Mt = wei_update_matrix(filteredPolyGrad);

   cv::Mat wei = cv::Mat::zeros(9, 1, CV_64F);
   std::vector<int> weiIndex = { 0, 1, 0, 1, 1, 0, 1, 1, 1 };
   std::vector<double> values = { 0.1140, 0.5870, 0.2989 };

   // Reverse the vector 'values'
   std::reverse(values.begin(), values.end());

   // Assign reversed values to 'wei' based on 'weiIndex'
   for (size_t i = 0, j = 0; i<weiIndex.size(); ++i) {
      if (weiIndex[i]==0) {
         if (j<values.size()) {
            wei.at<double>(i, 0) = values[j++];
         }
      }
   }

   // Selecting columns 1, 3, 6 from polyGrad
   std::vector<int> indices = { 0, 2, 5 };
   cv::Mat selectedColumns = cv::Mat(filteredPolyGrad.rows, indices.size(), filteredPolyGrad.type());
   for (size_t i = 0; i<indices.size(); ++i) {
      filteredPolyGrad.col(indices[i]).copyTo(selectedColumns.col(i));
   }

   // Select weights
   cv::Mat selectedWeights = cv::Mat(indices.size(), 1, wei.type());
   for (size_t i = 0; i<indices.size(); ++i) {
      selectedWeights.at<double>(i, 0) = wei.at<double>(indices[i], 0);
   }

   cv::Mat weightedSum = selectedColumns*selectedWeights;

   cv::Mat Cg1 = filteredCg-weightedSum;
   cv::Mat Cg2 = filteredCg+weightedSum;

   int iterCount = 0;
   while (std::abs(E-pre_E)>tol) {
      iterCount++;
      pre_E = E;

      weightedSum = filteredPolyGrad*wei;

      // Compute differences and sums
      cv::Mat diff = weightedSum-filteredCg;
      cv::Mat summation = weightedSum+filteredCg;

      diff = diff.mul(diff);
      summation = summation.mul(summation);

      double scale = -0.5/(sigmaParameter*sigmaParameter);
      diff *= scale;
      summation *= scale;

      cv::exp(diff, diff);
      cv::exp(summation, summation);

      // Combine with alf
      cv::Mat G_pos = (1.0+alf).mul(0.5*diff);
      cv::Mat G_neg = (1.0-alf).mul(0.5*summation);

      cv::Mat ExpSum = G_pos+G_neg;
      cv::Mat ExpTerm;

      cv::Mat numerator = G_pos.mul(Cg1)-G_neg.mul(Cg2);

      cv::Mat zeros = ExpSum==0;
      zeros.convertTo(zeros, CV_64F);

      // To avoid division by zero
      cv::Mat denominator = ExpSum+zeros;  
      cv::divide(numerator, denominator, ExpTerm); 

      cv::Mat weiTemp = Mt*ExpTerm;

      int j = 0;
      for (size_t i = 0; i<weiIndex.size(); ++i) {
         if (weiIndex[i]==1) {
            if (j<weiTemp.rows) {
               wei.at<double>(i, 0) = weiTemp.at<double>(j++, 0);
            }
         }
      }

      E = energyCalcu(filteredCg, filteredPolyGrad, wei, sigmaParameter);

      if (iterCount>maxIter) {
         break;
      }
   }

   // Construct the grayscale image
   cv::Mat grayIm = grayImContruct(wei, sourceImg, 2);

   // Save the grayscale image to the destination data
   int N = width*height;
   for (int i = 0; i<N; i++) {
      pDestinationData[i*3] = grayIm.at<double>(i/width, i%width);
      pDestinationData[i*3+1] = grayIm.at<double>(i/width, i%width);
      pDestinationData[i*3+2] = grayIm.at<double>(i/width, i%width);
   }

   return 0;
}
