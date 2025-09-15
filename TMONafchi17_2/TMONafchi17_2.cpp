/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio	                                       *
*                                                                                   *
*                       Project to subject Computational photography (VYF)          *
*                       Author: Boris Strbak [xstrba05 AT stud.fit.vutbr.cz]        *
*                       Brno 2024                                                   *
*                                                                                   *
*                      Converting color images to grayscale using correlation       *
*                      method by Nafchi, H. Z., Shahkolaei, A.,                     *
*                      Hedjam, R., and Cheriet, M.                                  *
*                                                                                   *
************************************************************************************/
/**
 * @file TMONafchi17_2.cpp
 * @brief Converting color images to grayscale using correlation method by Nafchi, H. Z., Shahkolaei, A., Hedjam, R., and Cheriet, M.
 * @author Boris Strbak
 * @class TMONafchi17_2.cpp
 */

#include "TMONafchi17_2.h"
#include <algorithm>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMONafchi17_2::TMONafchi17_2()
{
	SetName(L"Nafchi17_2");
	SetDescription(L"Color to Gray Conversion by Correlation.");

	invCorrParam.SetName(L"inv_cor");
	invCorrParam.SetDescription(L"Controls the level of contribution of the inverse correlations");
	invCorrParam.SetDefault(0.5);
	invCorrParam = 0.5;
	invCorrParam.SetRange(0.0, 1.0);
	this->Register(invCorrParam);

   downSampleParam.SetName(L"r");
	downSampleParam.SetDescription(L"Constant r to use to calculate downsampling factor (f = r / min(height, width))");
	downSampleParam.SetDefault(0);
	downSampleParam = 0;
	downSampleParam.SetRange(0, 512);
	this->Register(downSampleParam);

   comParam.SetName(L"com");
	comParam.SetDescription(L"Standard deviation image will be subsituted for it`s complement");
	comParam.SetDefault(false);
	comParam = false;
	this->Register(comParam);
}

TMONafchi17_2::~TMONafchi17_2()
{
}

/**
  *  @brief Return pearson`s correlation of two datasets
  *
  *  @param  X dataset X
  *  @param  meanX precomputed mean of dataset X
  *  @param  Y dataset Y
  *  @param  meanY precomputed mean of dataset Y
  *
  *  @return pearson`s correlation
  */
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

/**
  *  @brief Return image converted to cv::Mat structure
  *   where each item in 2d array has 3 values for each RGB
  *   channel
  *
  *  @param  image original image loaded by TMO library
  *
  *  @return cv::MAT represantion of image
  */
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
 * This overloaded function is an implementation of tone mapping operator      *
 * --------------------------------------------------------------------------- */
int TMONafchi17_2::Transform()
{
   // make sure images are in RGB
   pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

   cv::Mat srcRGB = getRGBImage(pSrc);
   cv::Mat srcRGBResized;

   double *pDestinationData = pDst->GetData();

   // create resized image for computing alpha parameters
   double minD = (double) std::min(srcRGB.rows, srcRGB.cols);
   if (downSampleParam > 0 && downSampleParam <= 512 && minD > downSampleParam) {
      double f = ((double) downSampleParam) / minD;

      cv::resize(srcRGB, srcRGBResized, cv::Size((int) (srcRGB.rows * f), (int) (srcRGB.cols * f)), cv::INTER_AREA);
   } else {
      srcRGBResized = srcRGB;
   }

   int totalPixelsOriginal = srcRGB.rows * srcRGB.cols;
   int totalPixels = srcRGBResized.rows * srcRGBResized.cols;

   // vectors for storing some data from computations that will be used for
   // other computations

   // here three values in each vector represent value for each RGB channel
   std::vector<double> cSums(3, 0.0);
   std::vector<double> cMean(3, 0.0);
   std::vector<double> cBeta(3, 0.0);
   std::vector<double> cWeight(3, 0.0);
   std::vector<double> cLambda(3, 0.0);
   std::vector<double> pearsCQ(3, 0.0);

   // here each value represent value for one pixel from src image
   std::vector<double> qVect(totalPixels, 0.0);
   std::vector<double> rVect(totalPixels, 0.0);
   std::vector<double> gVect(totalPixels, 0.0);
   std::vector<double> bVect(totalPixels, 0.0);

	double u, d, diff, g;
   double meanQ = 0.0, sumQ = 0.0;
   double pearsMin = 0.0, pearsMax = 0.0, pearsSum = 0.0, pearsDiff = 0.0;

   double lambdaSum = 0.0;

   // this should be max value of standard deviation in pixel
   // values for standard deviation are normalized with this value
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

         if (comParam) {
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
      cWeight[i] = ((pearsCQ[i] - pearsMin) / pearsDiff) - invCorrParam;
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
