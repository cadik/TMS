/* --------------------------------------------------------------------------- *
 * TMOLu14.cpp: implementation of the TMOLu14 class.   *
 * --------------------------------------------------------------------------- */

#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm> // For std::shuffle
#include <random>    // For std::default_random_engine
#include <chrono>    // For std::chrono::system_clock
#include "TMOLu14.h"


 /* --------------------------------------------------------------------------- *
  * Constructor serves for describing a technique and input parameters          *
  * --------------------------------------------------------------------------- */
TMOLu14::TMOLu14()
{
   SetName(L"TMOLu14");
   SetDescription(L"Semiparametric Decolorization With Laplacian-Based Perceptual Quality Metric");

   std::cout<<"init"<<std::endl;

   sigmaParameter.SetName(L"Sigma");
   sigmaParameter.SetDescription(L"ParameterDescription");
   sigmaParameter.SetDefault(0.005);
   sigmaParameter = 0.005;
   sigmaParameter.SetRange(0.0, 1.0);
   this->Register(sigmaParameter);
}

TMOLu14::~TMOLu14()
{
}

cv::Mat computeHorizontalGradient(const cv::Mat& img) {
   cv::Mat Rx(cv::Size(img.cols-1, img.rows), CV_64F, cv::Scalar(0)); // Initialize with zeros
   for (int i = 0; i<img.rows; i++) {
      for (int j = 0; j<img.cols-1; j++) {
         Rx.at<double>(i, j) = img.at<double>(i, j)-img.at<double>(i, j+1);
      }
   }
   return Rx;
}

cv::Mat computeVerticalGradient(const cv::Mat& img) {
   cv::Mat Ry(cv::Size(img.cols, img.rows-1), CV_64F, cv::Scalar(0)); // Initialize with zeros
   for (int i = 0; i<img.rows-1; i++) {
      for (int j = 0; j<img.cols; j++) {
         Ry.at<double>(i, j) = img.at<double>(i, j)-img.at<double>(i+1, j);
      }
   }
   return Ry;
}


void printMat(const cv::Mat& mat, const char* name = "mat") {
   std::cout<<name<<": ["<<mat.rows<<","<<mat.cols<<"]"<<std::endl;
   //   Print first 10 values of gradX
   auto startX = 0;
   auto startY = 0;
   auto endX = startX+9;
   auto endY = startY+9;

   for (int j = startX; j<endX; j++) {
      std::cout<<"[";
      for (int i = startY; i<endY; i++) {
         if (j>=mat.rows||i>=mat.cols)
         {
            break;
         }

         std::cout<<mat.at<double>(j, i)<<",";
      }
      std::cout<<"]"<<std::endl;
   }
}

cv::Mat weak_order(const cv::Mat& polyGrad) {
   double level = 0.05;

   // Extract specific columns for Rg, Gg, Bg
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

   // Convert logical arrays to double
   posWeakOrder.convertTo(posWeakOrder, CV_64F);
   negWeakOrder.convertTo(negWeakOrder, CV_64F);

   // Calculate alf
   cv::Mat alf = posWeakOrder-negWeakOrder;

   // Calculate percentage of pixels using weak order
   double percent = 100.0*cv::countNonZero(alf)/alf.total();
   std::cout<<round(percent)<<"% pixel use weak order"<<std::endl;

   return alf;
}

cv::Mat TMOLu14::calculatePolyGrad(const cv::Mat& img, int order) {
   // Check for empty input
   if (img.empty()) {
      throw std::runtime_error("Input image is empty");
   }

   // Resize image to roughly constant pixel count
   double factor = 64.0/std::sqrt(img.rows*img.cols);

   cv::Mat resized;
   cv::resize(img, resized, cv::Size(), factor, factor, cv::INTER_NEAREST);

   // Split color channels
   std::vector<cv::Mat> channels(3);
   cv::split(resized, channels);

   cv::Mat ImR = channels[0];  // Note: OpenCV loads images in BGR format
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

   // Transform featuresMatrices to a single rows (reshape(0,1))
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


   cv::split(img, channels);

   // Initialize P1 to be 
   cv::Mat P1;
   // Compute polynomial gradients for each combination
   for (int r = 0; r<=order; ++r) {
      for (int g = 0; g<=order; ++g) {
         for (int b = 0; b<=order; ++b) {
            if ((r+g+b)<=order&&(r+g+b)>0) {
               cv::Mat curIm = cv::Mat::ones(img.size(), CV_64F);

               cv::pow(channels[0], r, curIm); // Apply power to Red channel

               if (g>0) {
                  cv::Mat temp;
                  cv::pow(channels[1], g, temp); // Apply power to Green channel
                  curIm = curIm.mul(temp);
               }
               if (b>0) {
                  cv::Mat temp;
                  cv::pow(channels[2], b, temp); // Apply power to Blue channel
                  curIm = curIm.mul(temp);
               }

               // Compute gradients
               cv::Mat gradX = computeHorizontalGradient(curIm);
               cv::Mat gradY = computeVerticalGradient(curIm);

               // Flatten rx to single column

               gradX = gradX.reshape(0, 1).t(); // Flatten Rx to single column
               gradY = gradY.reshape(0, 1).t(); // Flatten Ry to single column


               cv::Mat tmp;
               // Append to P1 merged gradX and gradY
               cv::vconcat(gradX, gradY, tmp);

               // Add as new column to  P1
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

   // Transpose poly to get B
   cv::Mat B = poly.t();

   // Compute (poly' * poly) which is the same as (poly.t() * poly)
   cv::Mat polyTransposed = poly.t();
   cv::Mat A = polyTransposed*poly;

   // Solve the linear system A * Mt = B
   cv::Mat Mt;
   cv::solve(A, B, Mt, cv::DECOMP_SVD); // Using SVD decomposition to solve the system

   return Mt;
}

double energyCalcu(const cv::Mat& Cg, const cv::Mat& polyGrad, const cv::Mat& wei, double sigma) {
   // Matrix multiplication
   cv::Mat polyGradWei = polyGrad*wei;

   // Element-wise operations: Calculate polyGrad*wei - Cg and polyGrad*wei + Cg
   cv::Mat negDiff = polyGradWei-Cg;
   cv::Mat posDiff = polyGradWei+Cg;

   // Squaring, scaling, and exponential
   negDiff = negDiff.mul(negDiff)/sigma;  // (polyGrad*wei - Cg).^2 / sigma
   posDiff = posDiff.mul(posDiff)/sigma;  // (polyGrad*wei + Cg).^2 / sigma

   cv::exp(-negDiff, negDiff);  // exp(-(polyGrad*wei - Cg).^2 / sigma)
   cv::exp(-posDiff, posDiff);  // exp(-(polyGrad*wei + Cg).^2 / sigma)

   // Sum of exponentials
   cv::Mat expSum = negDiff+posDiff;

   // Logarithmic transformation
   cv::Mat P;
   cv::log(expSum, P);
   P = -P;

   // Mean of P
   cv::Scalar meanValue = cv::mean(P);

   return meanValue[0];  // Return the mean as double
}

cv::Mat grayImContruct(const std::vector<double>& wei, const cv::Mat& Im, int order) {
   // Assuming Im is a CV_8UC3 image (standard colored image)
   cv::Mat channels[3];
   cv::split(Im, channels);  // Split into R, G, B

   cv::Mat ImR = channels[0];  // Note: OpenCV loads images in BGR format
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

               cv::pow(ImR, r, curIm); // Apply power to Red channel

               if (g>0) {
                  cv::Mat temp;
                  cv::pow(ImG, g, temp); // Apply power to Green channel
                  curIm = curIm.mul(temp);
               }
               if (b>0) {
                  cv::Mat temp;
                  cv::pow(ImB, b, temp); // Apply power to Blue channel
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


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOLu14::Transform()
{
   // Parameters and initial setup
   double pre_E = std::numeric_limits<double>::infinity();
   double E = 0;
   double tol = std::pow(10, -4);  // Tolerance
   int maxIter = 15;

   // Source image is stored in local parameter pSrc
   // Destination image is in pDst

   pSrc->Convert(TMO_RGB, false);
   pDst->Convert(TMO_RGB, false);

   double* pSourceData = pSrc->GetData();
   double* pDestinationData = pDst->GetData();

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();

   cv::Mat sourceImg(height, width, CV_64FC3, pSourceData);

   cv::Mat polyGrad = calculatePolyGrad(sourceImg, 2);

   cv::Mat Cg = calculateCg(polyGrad);
   auto thr = 0.05;

   cv::Mat filteredPolyGrad(0, polyGrad.cols, polyGrad.type());
   cv::Mat filteredCg(0, Cg.cols, Cg.type());
   for (int i = 0; i<Cg.rows; i++) {
      if (Cg.at<double>(i)>=thr) {
         filteredPolyGrad.push_back(polyGrad.row(i));
         filteredCg.push_back(Cg.row(i));
      }
   }

   cv::Mat alf = weak_order(filteredPolyGrad);


   cv::Mat Mt = wei_update_matrix(filteredPolyGrad);

   cv::Mat wei = cv::Mat::zeros(9, 1, CV_64F);  // Creates a 9x1 matrix of doubles, all initialized to zero
   std::vector<int> wei_index = { 0, 1, 0, 1, 1, 0, 1, 1, 1 };
   std::vector<double> values = { 0.1140, 0.5870, 0.2989 };

   // Reverse the vector 'values'
   std::reverse(values.begin(), values.end());

   // Assign reversed values to 'wei' based on 'wei_index'
   for (size_t i = 0, j = 0; i<wei_index.size(); ++i) {
      if (wei_index[i]==0) {
         if (j<values.size()) {
            wei.at<double>(i, 0) = values[j++];
         }
      }
   }

   // Selecting columns 1, 3, 6 from polyGrad (0-based indexing in C++)
   std::vector<int> indices = { 0, 2, 5 };  // Corresponds to MATLAB's 1, 3, 6
   cv::Mat selectedColumns = cv::Mat(filteredPolyGrad.rows, indices.size(), filteredPolyGrad.type());
   for (size_t i = 0; i<indices.size(); ++i) {
      filteredPolyGrad.col(indices[i]).copyTo(selectedColumns.col(i));
   }

   // Select weights
   cv::Mat selectedWeights = cv::Mat(indices.size(), 1, wei.type());
   for (size_t i = 0; i<indices.size(); ++i) {
      selectedWeights.at<double>(i, 0) = wei.at<double>(indices[i], 0);
   }

   // Perform the weighted sum operation
   cv::Mat weightedSum = selectedColumns*selectedWeights;

   // Compute Cg1 and Cg2
   cv::Mat Cg1 = filteredCg-weightedSum;
   cv::Mat Cg2 = filteredCg+weightedSum;

   // While TODO: Implement while loop

   // Main loop
   int iterCount = 0;
   while (std::abs(E-pre_E)>tol) {
      iterCount++;
      pre_E = E;

      // Calculate polyGrad * wei
      weightedSum = filteredPolyGrad*wei;

      // Compute differences and sums
      cv::Mat diff = weightedSum-filteredCg;
      cv::Mat summation = weightedSum+filteredCg;

      // Square, scale by -0.5 / sigma^2
      diff = diff.mul(diff); // Element-wise square
      summation = summation.mul(summation); // Element-wise square

      double scale = -0.5/(sigmaParameter*sigmaParameter);
      diff *= scale;
      summation *= scale;

      // Exponential
      cv::exp(diff, diff); // Element-wise exp
      cv::exp(summation, summation); // Element-wise exp

      // Combine with alf
      cv::Mat G_pos = (1.0+alf).mul(0.5*diff); // Element-wise multiplication
      cv::Mat G_neg = (1.0-alf).mul(0.5*summation); // Element-wise multiplication

      cv::Mat ExpSum = G_pos+G_neg;
      cv::Mat ExpTerm;

      cv::Mat numerator = G_pos.mul(Cg1)-G_neg.mul(Cg2);  // Element-wise multiplication and subtraction

      cv::Mat zeros = ExpSum==0;
      zeros.convertTo(zeros, CV_64F);

      cv::Mat denominator = ExpSum+zeros;  // Avoid division by zero
      cv::divide(numerator, denominator, ExpTerm);   // Element-wise division

      cv::Mat weiTemp = Mt*ExpTerm;

      // Conditional assignment based on wei_index
      int j = 0; // Index for accessing elements in wei_temp
      for (size_t i = 0; i<wei_index.size(); ++i) {
         if (wei_index[i]==1) {
            if (j<weiTemp.rows) {
               wei.at<double>(i, 0) = weiTemp.at<double>(j++, 0);
            }
         }
      }

      E = energyCalcu(filteredCg, filteredPolyGrad, wei, sigmaParameter);

      // Print the current energy
      std::cout<<"Energy: "<<E<<std::endl;

      if (iterCount>maxIter) {
         break;
      }
   }

   // Construct the grayscale image
   cv::Mat grayIm = grayImContruct(wei, sourceImg, 2);

   // Copy the grayscale image to the destination image
   std::memcpy(pDestinationData, grayIm.data, width*height*3*sizeof(double));

   int N = width * height;
    for (int i = 0; i < N; i++) {
        pDestinationData[i * 3] = grayIm.at<double>(i / width, i % width);
        pDestinationData[i * 3 + 1] = grayIm.at<double>(i / width, i % width);
        pDestinationData[i * 3 + 2] = grayIm.at<double>(i / width, i % width);
    }

   return 0;
}
