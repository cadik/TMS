/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2022                                              *
*                                                                              *
*                       Implementation of the TMOShan10 class                  *
*                                                                              *
*******************************************************************************/

#include "TMOShan10.h"

#define CHANNELSCNT 3
#define PADARR_REPLICATE 0
#define PADARR_CONSTANT 1

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOShan10::TMOShan10()
{
   SetName(L"Shan10");
   SetDescription(L"Globally Optimized Linear Windowed Tone-Mapping");

   winSize.SetName(L"winSize");
   winSize.SetDescription(L"The local window size about each pixel.");
   winSize.SetDefault(3);
   winSize = 3;
   winSize.SetRange(3, 15);
   this->Register(winSize);

   sSat.SetName(L"sSat");
   sSat.SetDescription(L"The saturation factor in Eq. (5).");
   sSat.SetDefault(0.6);
   sSat = 0.6;
   sSat.SetRange(0.4, 0.6);
   this->Register(sSat);

   beta2.SetName(L"beta2");
   beta2.SetDescription(L"Beta2 parameter. Article Eq. (6).");
   //beta2.SetDefault(0.1);
   //beta2 = 0.1;
   //beta2.SetRange(0.1, 0.4);
   beta2.SetDefault(0.0);
   beta2 = 0.0;
   beta2.SetRange(0.0, 0.4);
   this->Register(beta2);

   beta1.SetName(L"beta1");
   beta1.SetDescription(L"Beta1 parameter. Article Eq. (6).");
   beta1.SetDefault(0.4);
   beta1 = 0.4;
   beta1.SetRange(0.4, 0.9);
   this->Register(beta1);
}

TMOShan10::~TMOShan10()
{
}

/*
 * Own implementation of MATLAB padarray function using OpenCV copyMakeBorder().
 *
 * Original MATLAB function: https://www.mathworks.com/help/images/ref/padarray.html
 */
Mat padArray (Mat srcMat, int rowPad, int colPad, uint8_t borderType)
{
   int height, width;

   assert(borderType == PADARR_REPLICATE || borderType == PADARR_CONSTANT);

   height = srcMat.rows;
   width = srcMat.cols;

   Mat tmp = Mat::zeros(height + rowPad*2, width + colPad*2, srcMat.type());

   if (borderType == PADARR_REPLICATE)
   {
      copyMakeBorder(srcMat, tmp, rowPad, rowPad, colPad, colPad, BORDER_REPLICATE);
   }
   else if (borderType == PADARR_CONSTANT)
   {
      copyMakeBorder(srcMat, tmp, rowPad, rowPad, colPad, colPad, BORDER_CONSTANT, 0.0);
   }

   return tmp;
}

/*
 * Own implementation of MATLAB colfilt function.
 *
 * colfilt function settings:
 * - Sliding Neighborhood Operations
 * - @std function (Standard deviation)
 *
 * Original MATLAB function: https://www.mathworks.com/help/images/ref/colfilt.html
 */
Mat colFiltSlidingStd(Mat srcMat, int blockRowsCnt, int blockColsCnt)
{
   int srcHeight, srcWidth;

   srcHeight = srcMat.rows;
   srcWidth = srcMat.cols;

   assert(srcMat.type() == CV_64FC1);

   Mat result(srcHeight, srcWidth, srcMat.type());

   int rowsStepsCnt = srcHeight - blockRowsCnt + 1;
   int colsStepsCnt = srcWidth - blockColsCnt + 1;

   Mat tmp(blockRowsCnt * blockColsCnt, rowsStepsCnt * colsStepsCnt, CV_64FC1);

   // im2col function
   // https://www.mathworks.com/help/images/ref/im2col.html
   for (int x = 0; x < colsStepsCnt; x++)
   {
      for (int y = 0; y < rowsStepsCnt; y++)
      {
         for (int i = x; i < x + blockColsCnt; i++)
         {
            for (int j = y; j < y + blockRowsCnt; j++)
            {
                tmp.at<double>((i-x) * blockRowsCnt + (j-y), x * rowsStepsCnt + y) = srcMat.at<double>(j, i);
            }
         }
      }
   }

   // std function - standard deviation
   // https://www.mathworks.com/help/matlab/ref/std.html

   double *colData = new double[tmp.rows];
   double *stdDev = new double[tmp.cols];

   for (int x = 0; x < tmp.cols; x++)
   {
      // Get column data
      for (int y = 0; y < tmp.rows; y++)
      {
         colData[y] = tmp.at<double>(y, x);
      }

      // Calculate mean of column data
      double mean = 0;
      for (int i = 0; i < tmp.rows; i++)
      {
         mean += colData[i];
      }

      mean /= (double)tmp.rows;

      double sum = 0;

      for (int i = 0; i < tmp.rows; i++)
      {
         sum += pow(colData[i]-mean, 2);
      }

      // Calculate result standard deviation
      stdDev[x] = sqrt(sum / (double)tmp.rows-1);
   }

   // reshape function
   // https://www.mathworks.com/help/matlab/ref/reshape.html

   for (int x = 0; x < srcWidth; x++)
   {
      for (int y = 0; y < srcHeight; y++)
      {
         result.at<double>(y, x) = stdDev[x*srcHeight + y];
      }
   }

   delete[] colData;
   delete[] stdDev;

   return result;
}

/*
 * Convert TMOImage to cv::Mat
 */
Mat TMOImage2Mat(TMOImage *pSrc)
{
   double* pSourceData;
   int rowsCnt, colsCnt;

   pSourceData = pSrc->GetData();
   rowsCnt = pSrc->GetHeight();
   colsCnt = pSrc->GetWidth();

   Mat srcConvMat(rowsCnt, colsCnt, CV_64FC3);

   if (srcConvMat.isContinuous())
   {
      colsCnt *= rowsCnt;
      rowsCnt = 1;
   }

   for (int y = 0; y < rowsCnt; y++)
   {
      for (int x = 0; x < colsCnt; x++)
      {
         for (int c = 0; c < CHANNELSCNT; c++)
         {
            srcConvMat.at<double>(y, x*CHANNELSCNT+c) = pSourceData[CHANNELSCNT-c-1];
         }

         /* Add count of channels (RGB) to pointer */
         pSourceData += CHANNELSCNT;
      }
   }

   return srcConvMat;
}

/*
 * Set value by HSV color model.
 *
 * This function does not implement whole algorithm, but only Value part.
 * Article: Smith, A. R. “Color Gamut Transform Pairs”. SIGGRAPH 78 Conference Proceedings. 1978, pp. 12–19.
 */
void setValueHSV(Mat &srcMat, Mat &lum)
{
   double R, G, B, V;
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   int y = 0;
   for (y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         R = srcMat.at<double>(y, x*CHANNELSCNT);
         G = srcMat.at<double>(y, x*CHANNELSCNT+1);
         B = srcMat.at<double>(y, x*CHANNELSCNT+2);

         V = max(R, G);
         V = max(V, B);

         lum.at<double>(y, x) = V;
      }
   }
}

/*
 * Check luminance value.
 *
 * If the value is lower than 1e-3, to actual luminance value is add 1e-3.
 */
void checkLumo(Mat &lumo)
{
   int height, width;

   height = lumo.rows;
   width = lumo.cols;

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         if (lumo.at<double>(y, x) < 1e-3)
         {
            lumo.at<double>(y, x) = lumo.at<double>(y, x) + 1e-3;
         }
      }
   }
}

/*
 * Create guidance map.
 */
Mat TMOShan10::generateGuidanceMap(Mat lumo,
                                   int winSize,
                                   double beta1,
                                   double beta2,
                                   const double beta3,
                                   const double kappa)
{
   int height, width;
   int halfWinSize;

   height = lumo.rows;
   width = lumo.cols;

   halfWinSize = (winSize - 1) / 2;

   Mat meanFilter = Mat::ones(winSize, winSize, lumo.type());

   /*
    * Source: https://stackoverflow.com/questions/21874774/sum-of-elements-in-a-matrix-in-opencv
    * Author: Shai (https://stackoverflow.com/users/1714410/shai)
    */
   meanFilter = meanFilter / cv::sum(meanFilter)[0];

   Mat meanLum(height, width, lumo.type());
   Mat tLum(height, width, lumo.type());

   filter2D(lumo, meanLum, -1, meanFilter, Point(-1, -1), 0, BORDER_REPLICATE); // TODO jestli je to korelace

   tLum = padArray(lumo, halfWinSize, halfWinSize, PADARR_REPLICATE);

   Mat tmpStdLum;

   tmpStdLum = colFiltSlidingStd(tLum, winSize, winSize);

   Mat stdLum(height, width, lumo.type());

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         stdLum.at<double>(y, x) = tmpStdLum.at<double>(y + halfWinSize, x + halfWinSize);
      }
   }

   Mat ci(height, width, lumo.type());

   pow(meanLum, beta1, meanLum);
   pow(stdLum, beta2, stdLum);
   pow(lumo, beta3, lumo);

   ci = meanLum.mul(stdLum);
   ci = ci.mul(lumo);

   ci = ci + kappa;
   ci = ci * 0.2;

   ci = 1.0 / ci;

   return ci;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOShan10::Transform()
{
   int height, width;
   int pWinSize;
   double pBeta1, pBeta2;
   double pSSat;

   const uint8_t levelNum = 3;
   const double kappa = 0.5;
   const double beta3 = 0.1;
   const double epsilon = 0.1;
   const double lambda = 0.0;

   /***********************/
   /* Get user parameters */
   /***********************/

   pBeta1 = beta1.GetDouble();
   pBeta2 = beta2.GetDouble();
   pSSat = sSat.GetDouble();
   pWinSize = winSize.GetInt();

   pSrc->convertToGrayscale();

   Mat srcMat = TMOImage2Mat(pSrc);

   height = srcMat.rows;
   width = srcMat.cols;

   Mat tLum(height, width, CV_64FC1);
   Mat lumo(height, width, CV_64FC1);

   setValueHSV(srcMat, tLum);

   Scalar tLumMean = mean(mean(mean(tLum)));

   // Because tLum has only one channel, the tLumMean scalar
   // has only value at first element position -> val[0]
   srcMat = srcMat * (200.0 / tLumMean.val[0]);

   setValueHSV(srcMat, lumo);
   checkLumo(lumo);

   Mat map = generateGuidanceMap(lumo, pWinSize, pBeta1, pBeta2, beta3, kappa);

   /*-----------------------------------------------------*/
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);

	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			//*pDestinationData++ = *pSourceData++;
         //*pDestinationData++ = *pSourceData++;
         //*pDestinationData++ = *pSourceData++;
         *pDestinationData++ = map.at<double>(j, i);
         *pDestinationData++ = map.at<double>(j, i);
         *pDestinationData++ = map.at<double>(j, i);
		}
	}

	pSrc->ProgressBar(j, pSrc->GetHeight());

	return 0;
}
