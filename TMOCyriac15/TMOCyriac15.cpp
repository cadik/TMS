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
   std::vector<double> luminancesVector;

   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         luminancesVector.push_back(luminanceMat->at<double>(x, y));
      }
   }

   std::sort(luminancesVector.begin(), luminancesVector.end());

   double clipValue = luminancesVector.at(std::round(0.99 * luminancesVector.size()));

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

cv::Mat TMOCyriac15::createMat(int colorChannel)
{
   cv::Mat mat(pDst->GetWidth(), pDst->GetHeight(), CV_64FC1);
   for (int y = 0; y < pDst->GetHeight(); y++)
   {
      for (int x = 0; x < pDst->GetWidth(); x++)
      {
         mat.at<double>(x, y) = pDst->GetPixel(x, y)[colorChannel];
      }
   }

   return mat;
}

cv::Mat TMOCyriac15::get2DGaussianKernel(int size, double sigma)
{
   cv::Mat kernel(size, size, CV_64FC1);
   double mean = size / 2;
   double sum = 0.0;
   for (int x = 0; x < size; x++)
   {
      for (int y = 0; y < size; y++)
      {
         kernel.at<double>(x, y) = std::exp(-0.5 * (std::pow(x - mean, 2) + std::pow(y - mean, 2)) / std::pow(sigma, 2));
         sum += kernel.at<double>(x, y);
      }
   }

   return kernel / sum;
}

cv::Mat TMOCyriac15::computeLocalMean(cv::Mat *originalIntensity)
{
   cv::Mat G1 = get2DGaussianKernel(std::max(pDst->GetWidth(), pDst->GetHeight()), 10);
   cv::Mat G2 = get2DGaussianKernel(std::max(pDst->GetWidth(), pDst->GetHeight()), 250);
   double n1 = 1.0;
   double n2 = 0.5;

   cv::Mat K = n1 * G1 + n2 * G2;

   cv::filter2D(*originalIntensity, K, -1, K);

   return K;
}

double TMOCyriac15::sgn(double value)
{
   if (value > 0)
   {
      return 1.0;
   }
   else if (value < 0)
   {
      return -1.0;
   }
   else
   {
      return 0.0;
   }
}

double TMOCyriac15::gradientDescent(int x, int y, cv::Mat *intensity, cv::Mat *localMean, cv::Mat *kernel, cv::Mat *originalIntensity)
{
   double I0 = originalIntensity->at<double>(x, y);
   double mu = localMean->at<double>(x, y);
   double Ix = intensity->at<double>(x, y);
   double Iy;
   double sum = 0.0;

   for (int i = 0; i < pDst->GetWidth(); i++)
   {
      for (int j = 0; j < pDst->GetHeight(); j++)
      {
         Iy = intensity->at<double>(i, j);
         sum += kernel->at<double>(i, j) * sgn(Ix - Iy);
      }
   }

   Ix = -alpha * (Ix - mu) + gamma * sum - beta * (Ix - I0);
   intensity->at<double>(x, y) = Ix;
   return Ix;
}

double TMOCyriac15::findMinEnergy(int x, int y, cv::Mat *intensity, cv::Mat *localMean, cv::Mat *kernel, cv::Mat *originalIntensity)
{
   double Iprev = intensity->at<double>(x, y);
   double Inew = Iprev;

   while (Iprev - Inew >= convergenceThreshold)
   {
      Iprev = Inew;
      Inew = Iprev + deltaT * gradientDescent(x, y, intensity, localMean, kernel, originalIntensity);
   }

   return Inew;
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

   // second stage

   cv::Mat intensityR = createMat(0);
   cv::Mat intensityG = createMat(1);
   cv::Mat intensityB = createMat(2);
   cv::Mat originalIntensityR = intensityR.clone();
   cv::Mat originalIntensityG = intensityG.clone();
   cv::Mat originalIntensityB = intensityB.clone();
   cv::Mat localMeanR = computeLocalMean(&originalIntensityR);
   cv::Mat localMeanG = computeLocalMean(&originalIntensityG);
   cv::Mat localMeanB = computeLocalMean(&originalIntensityB);
   cv::Mat kernel = get2DGaussianKernel(std::max(pDst->GetWidth(), pDst->GetHeight()), sigma_w);

   for (int x = 0; x < pDst->GetWidth(); x++)
   {
      for (int y = 0; y < pDst->GetHeight(); y++)
      {
         pDst->GetPixel(x, y)[0] = findMinEnergy(x, y, &intensityR, &localMeanR, &kernel, &originalIntensityR);
         pDst->GetPixel(x, y)[1] = findMinEnergy(x, y, &intensityG, &localMeanG, &kernel, &originalIntensityG);
         pDst->GetPixel(x, y)[2] = findMinEnergy(x, y, &intensityB, &localMeanB, &kernel, &originalIntensityB);
      }
   }

   return 0;
}
