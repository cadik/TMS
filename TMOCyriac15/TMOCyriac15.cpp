/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   A Tone Mapping Operator Based on Neural and                *
 *                    Psychophysical Models of Visual Perception                *
 * 																			                    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMOCyriac15.h"

TMOCyriac15::TMOCyriac15()
{
   SetName(L"Cyriac15");                             // TODO - Insert operator name
   SetDescription(L"Add your TMO description here"); // TODO - Insert description

   dParameter.SetName(L"ParameterName");               // TODO - Insert parameters names
   dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
   dParameter.SetDefault(1);                           // TODO - Add default values
   dParameter = 1.;
   dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
   this->Register(dParameter);
}

TMOCyriac15::~TMOCyriac15()
{
}

cv::Mat TMOCyriac15::convertToMat()
{
   cv::Mat mat(pDst->GetWidth(), pDst->GetHeight(), CV_64FC3);

   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         double *pixel = pDst->GetPixel(x, y);
         mat.at<cv::Vec3d>(x, y) = cv::Vec3d(pixel[0], pixel[1], pixel[2]);
      }
   }

   return mat;
}

void TMOCyriac15::convertToPDst(cv::Mat &img)
{
   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         cv::Vec3d pixel = img.at<cv::Vec3d>(x, y);
         pDst->GetPixel(x, y)[0] = pixel[0];
         pDst->GetPixel(x, y)[1] = pixel[1];
         pDst->GetPixel(x, y)[2] = pixel[2];
      }
   }
}

cv::Mat TMOCyriac15::normalizeAndLuminance(std::vector<double> *RGBmax)
{
   cv::Mat luminanceMat(pDst->GetWidth(), pDst->GetHeight(), CV_64FC1);

   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         double *pixel = pDst->GetPixel(x, y);
         double R = pixel[0];
         double G = pixel[1];
         double B = pixel[2];
         pixel[0] = std::max(R, 0.0) / RGBmax->at(0);
         pixel[1] = std::max(G, 0.0) / RGBmax->at(1);
         pixel[2] = std::max(B, 0.0) / RGBmax->at(2);
         luminanceMat.at<double>(x, y) = 0.2126 * R + 0.7152 * G + 0.0722 * B;
      }
   }

   return luminanceMat;
}

void TMOCyriac15::clip(cv::Mat *luminanceMat)
{
   double maxLuminance = 0.0;
   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         double L = luminanceMat->at<double>(x, y);
         if (L > maxLuminance)
         {
            maxLuminance = L;
         }
      }
   }

   double clipValue = 0.99 * maxLuminance;

   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         luminanceMat->at<double>(x, y) = std::min(luminanceMat->at<double>(x, y) / clipValue, 1.0);
         double *pixel = pDst->GetPixel(x, y);
         pixel[0] = std::min(pixel[0] / clipValue, 1.0);
         pixel[1] = std::min(pixel[1] / clipValue, 1.0);
         pixel[2] = std::min(pixel[2] / clipValue, 1.0);
      }
   }
}

void TMOCyriac15::addToCumulativeHistogram(std::vector<double> *cumulativeHistogram, double value)
{
   for (int i = std::round(value * (histogramBins - 1)); i < cumulativeHistogram->size(); i++)
   {
      cumulativeHistogram->at(i) += 1.0;
   }
}

double TMOCyriac15::findGammaEnc(cv::Mat *luminanceMat)
{
   double gammaSys = 0.4;
   double gammaPsy = 0.4;
   double gammaDec = 2.2;
   double previousF = 0.0;
   double difference = 1.0;
   cv::Mat tmpMat = luminanceMat->clone();

   while ((difference > 0.0) && (gammaSys < 2.2))
   {
      gammaSys += 0.1;

      // compute cumulative histogram
      for (int y = 0; y < pDst->GetHeight(); y++)
      {
         for (int x = 0; x < pDst->GetWidth(); x++)
         {
            double L = luminanceMat->at<double>(x, y);
            tmpMat.at<double>(x, y) = std::pow(L, gammaSys * gammaPsy);
         }
      }

      tmpMat.convertTo(tmpMat, CV_32F);
      float range[] = {0, 1};
      const float *histogramRange = {range};
      cv::Mat histogram;
      cv::calcHist(&tmpMat, 1, 0, cv::Mat(), histogram, 1, &histogramBins, &histogramRange, true, false);

      // compute sum
      double sum = 0.0;
      for (int i = 0; i < histogramBins; i++)
      {
         histogram.at<float>(i) += histogram.at<float>(i - 1); // cumulative histogram
         sum += std::pow(histogram.at<float>(i) - i, 2);
      }

      double F = 1 - std::sqrt(sum) / histogramBins;
      difference = F - previousF;
      previousF = F;
   }

   return gammaSys / gammaDec;
}

cv::Mat TMOCyriac15::createGaussianKernel(int size, double sigma)
{
   size = size | 1;
   cv::Mat kernel = cv::getGaussianKernel(size, sigma);
   cv::Mat kernel2D = kernel * kernel.t();
   return kernel2D;
}

cv::Mat TMOCyriac15::computeLocalMean(cv::Mat &img)
{
   cv::Mat G1 = createGaussianKernel(31, 10.0);
   cv::Mat G2 = createGaussianKernel(31, 250.0);
   double n1 = 1.0;
   double n2 = 0.5;
   cv::Mat K = n1 * G1 + n2 * G2;
   cv::Mat result;
   cv::filter2D(img, result, CV_64F, K);

   return result;
}

double TMOCyriac15::sgn(double x)
{
   if (x > 0)
   {
      return 1.0;
   }
   else if (x < 0)
   {
      return -1.0;
   }
   else
   {
      return 0.0;
   }
}

cv::Mat TMOCyriac15::computeGradient(cv::Mat &I, cv::Mat &I0, cv::Mat &localMean)
{
   cv::Mat gradient = cv::Mat::zeros(I.size(), CV_64F);
   cv::Mat paddedI;
   cv::copyMakeBorder(I, paddedI, 1, 1, 1, 1, cv::BORDER_REFLECT);

   for (int y = 0; y < I.rows; y++)
   {
      for (int x = 0; x < I.cols; x++)
      {
         double sum = 0.0;

         for (int dy = -1; dy <= 1; dy++)
         {
            for (int dx = -1; dx <= 1; dx++)
            {
               if (dx == 0 && dy == 0)
                  continue;

               double w = exp(-(dx * dx + dy * dy) / (2 * sigma_w * sigma_w));
               double diff = paddedI.at<double>(y + 1, x + 1) - paddedI.at<double>(y + 1 + dy, x + 1 + dx);
               sum += w * sgn(diff) * diff;
            }
         }

         gradient.at<double>(y, x) = -alpha * (I.at<double>(y, x) - localMean.at<double>(y, x)) + gamma * sum - beta * (I.at<double>(y, x) - I0.at<double>(y, x));
      }
   }

   return gradient;
}

int TMOCyriac15::Transform()
{
   // copy source image to destination image and find min and max RGB values
   std::vector<double> RGBmax = {0.0, 0.0, 0.0};

   for (int y = 0; y < pSrc->GetHeight(); y++)
   {
      for (int x = 0; x < pSrc->GetWidth(); x++)
      {
         double *pixel = pSrc->GetPixel(x, y);
         pDst->GetPixel(x, y)[0] = pixel[0];
         pDst->GetPixel(x, y)[1] = pixel[1];
         pDst->GetPixel(x, y)[2] = pixel[2];

         RGBmax[0] = std::max(RGBmax[0], pixel[0]);
         RGBmax[1] = std::max(RGBmax[1], pixel[1]);
         RGBmax[2] = std::max(RGBmax[2], pixel[2]);
      }
   }
   cv::Mat luminanceMat = normalizeAndLuminance(&RGBmax);
   clip(&luminanceMat);
   double gammaEnc = findGammaEnc(&luminanceMat);

   // apply gamma encoding
   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         double *pixel = pDst->GetPixel(x, y);
         pixel[0] = std::pow(pixel[0], gammaEnc);
         pixel[1] = std::pow(pixel[1], gammaEnc);
         pixel[2] = std::pow(pixel[2], gammaEnc);
         luminanceMat.at<double>(x, y) = std::pow(luminanceMat.at<double>(x, y), gammaEnc);
      }
   }

   pDst->Convert(TMO_Yxy);

   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         pDst->GetPixel(x, y)[0] = luminanceMat.at<double>(x, y);
      }
   }

   pDst->Convert(TMO_RGB);

   cv::Mat img = convertToMat();
   cv::Mat I;
   img.convertTo(I, CV_64F, 1.0 / 255.0);
   cv::Mat I0 = I.clone();
   cv::Mat previousI;
   int iteration = 0;
   double maxDiff;

   do
   {
      previousI = I.clone();

      cv::Mat localMean = computeLocalMean(I);

      cv::Mat gradient = computeGradient(I, I0, localMean);

      cv::Mat gradientStep;
      cv::multiply(gradient, deltaT, gradientStep);
      cv::add(I, gradientStep, I);

      // cv::min(cv::max(I, 0.0), 1.0, I);

      // cv::Mat diff;
      // cv::absdiff(I, previousI, diff);
      // cv::minMaxLoc(diff, nullptr, &maxDiff);

      iteration++;
   } while (maxDiff > convergenceThreshold && iteration < 1000);

   // convertToPDst(I);

   return 0;
}
