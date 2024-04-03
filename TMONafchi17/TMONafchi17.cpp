/* --------------------------------------------------------------------------- *
 * TMONafchi17.cpp: implementation of the TMONafchi17 class.   *
 * --------------------------------------------------------------------------- */

#include "TMONafchi17.h"
#include <algorithm>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMONafchi17::TMONafchi17()
{
	SetName(L"TMONafchi17");
	SetDescription(L"Color to Gray Conversion by Correlation.");

	dParameter.SetName(L"inv_cor");
	dParameter.SetDescription(L"Controls the level of contribution of the inverse correlations");
	dParameter.SetDefault(0.5);
	dParameter = 0.5;
	dParameter.SetRange(0.0, 1.0);
	this->Register(dParameter);

   iParameter.SetName(L"r");
	iParameter.SetDescription(L"Constant r to use to calculate downsampling factor (f = r / min(height, width))");
	iParameter.SetDefault(0);
	iParameter = 0;
	iParameter.SetRange(0, 512);
	this->Register(iParameter);

   bParameter.SetName(L"com");
	bParameter.SetDescription(L"Standard deviation image will be subsituted for it`s complement");
	bParameter.SetDefault(false);
	bParameter = false;
	this->Register(bParameter);
}

TMONafchi17::~TMONafchi17()
{
}

double getCorr(std::vector<double> *X, double meanX, std::vector<double> *Y, double meanY)
{
   double sumOfDiffs = 0.0, sumOfDiffX = 0.0, sumOfDiffY = 0.0;
   double downPart;

   for (int i = 0; i < X->size(); i++) {
      sumOfDiffs += (X->at(i) - meanX) * (Y->at(i) - meanY);
      sumOfDiffX += std::pow(X->at(i) - meanX, 2);
      sumOfDiffY += std::pow(Y->at(i) - meanY, 2);
   }

   downPart = std::sqrt(sumOfDiffX * sumOfDiffY);

   return sumOfDiffs / downPart;
}

cv::Mat getRGBImage(TMOImage *image)
{
   double *data = image->GetData();
   cv::Mat RGB(image->GetHeight(), image->GetWidth(), CV_64FC3);

   for (int j = 0; j < image->GetHeight(); j++)
	{
		for (int i = 0; i < image->GetWidth(); i++)
		{
         RGB.at<cv::Vec3d>(j, i)[0] = *data++;
			RGB.at<cv::Vec3d>(j, i)[1] = *data++;
			RGB.at<cv::Vec3d>(j, i)[2] = *data++;
      }
   }

   return RGB;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMONafchi17::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format
   pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

   cv::Mat srcRGB = getRGBImage(pSrc);
   cv::Mat srcRGBResized;
   double *pDestinationData = pDst->GetData();
   double minD = (double) std::min(srcRGB.rows, srcRGB.cols);

   if (iParameter > 0 && iParameter <= 512 && minD > iParameter) {
      double f = ((double) iParameter) / minD;

      cv::resize(srcRGB, srcRGBResized, cv::Size((int) (srcRGB.rows * f), (int) (srcRGB.cols * f)), cv::INTER_AREA);
   } else {
      srcRGBResized = srcRGB;
   }

   std::vector<double> cSums(3, 0.0);
   std::vector<double> cMean(3, 0.0);
   std::vector<double> cBeta(3, 0.0);
   std::vector<double> cWeight(3, 0.0);
   std::vector<double> cLambda(3, 0.0);
   std::vector<double> pearsCQ(3, 0.0);

	double u, d, diff, g;
   double meanQ = 0.0, sumQ = 0.0;
   double pearsMin = 0.0, pearsMax = 0.0, pearsSum = 0.0, pearsDiff = 0.0;

   double lambdaSum = 0.0;

   int totalPixelsOriginal = srcRGB.rows * srcRGB.cols;
   int totalPixels = srcRGBResized.rows * srcRGBResized.cols;

   std::vector<double> qVect(totalPixels, 0.0);
   std::vector<double> rVect(totalPixels, 0.0);
   std::vector<double> gVect(totalPixels, 0.0);
   std::vector<double> bVect(totalPixels, 0.0);

   const double dThrs = (147.2243d / 255.0d);

	int j = 0, index = 0;

   for (j = 0; j < srcRGBResized.rows; ++j) {
      pSrc->ProgressBar(j, srcRGBResized.rows);

      for (int i = 0; i < srcRGBResized.cols; ++i) {
         index = (j * srcRGBResized.cols) + i;

         cv::Vec3d pixel = srcRGBResized.at<cv::Vec3d>(j, i);

         cSums[0] += rVect[index] = pixel[0];
         cSums[1] += gVect[index] = pixel[1];
         cSums[2] += bVect[index] = pixel[2];

         u = (pixel[0] + pixel[1] + pixel[2]) / 3.0d;
         d = std::sqrt((std::pow(pixel[0] -  u, 2) + std::pow(pixel[1] -  u, 2) + std::pow(pixel[2] -  u, 2)) / 2.0d) / dThrs;

         if (bParameter) {
            d = 1 - d;
         }

         sumQ += qVect[index] = u * d;
      }
   }

   for (int i = 0; i < 3; ++i) {
      cMean[i] = cSums[i] / totalPixels;
   }

   meanQ = sumQ / totalPixels;

   pearsCQ[0] = getCorr(&qVect, meanQ, &rVect, cMean[0]);
   pearsCQ[1] = getCorr(&qVect, meanQ, &gVect, cMean[1]);
   pearsCQ[2] = getCorr(&qVect, meanQ, &bVect, cMean[2]);

   pearsMax = *std::max_element(pearsCQ.begin(), pearsCQ.end());
   pearsMin = *std::min_element(pearsCQ.begin(), pearsCQ.end());
   pearsSum = std::abs(pearsCQ[0]) + std::abs(pearsCQ[1])  + std::abs(pearsCQ[2]);

   pearsDiff = pearsMax - pearsMin;

   for (int i = 0; i < 3; ++i) {
      cBeta[i] = std::abs(pearsCQ[i]) / pearsSum;
      cWeight[i] = ((pearsCQ[i] - pearsMin) / pearsDiff) - dParameter;
      cLambda[i] = std::abs(cBeta[i] + std::min(cBeta[i], cWeight[i]));
   }

   lambdaSum = cLambda[0] + cLambda[1] + cLambda[2];

   for (int i = 0; i < 3; ++i) {
      cLambda[i] /= lambdaSum;
   }

   for (j = 0; j < srcRGB.rows; ++j) {
		pSrc->ProgressBar(j, srcRGB.rows);

      for (int i = 0; i < srcRGB.cols; ++i) {
         cv::Vec3d pixel = srcRGB.at<cv::Vec3d>(j, i);

         g = ((cLambda[0] * pixel[0]) + (cLambda[1] * pixel[1]) + (cLambda[2] * pixel[2]));

			*pDestinationData++ = g;
         *pDestinationData++ = g;
         *pDestinationData++ = g;
		}
	}

	pDst->Convert(TMO_RGB);
	return 0;
}
