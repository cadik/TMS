/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2022                                              *
*                                                                              *
*                       Implementation of the MatUtil class                    *
*                                                                              *
*******************************************************************************/

#include "MatUtil.h"

/*
 * Own implementation of MATLAB padarray function using OpenCV copyMakeBorder().
 *
 * Original MATLAB function: https://www.mathworks.com/help/images/ref/padarray.html
 */
Mat MatUtil::padArray (Mat &srcMat, int rowPad, int colPad, uint8_t borderType)
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
Mat MatUtil::colFiltSlidingStd(Mat &srcMat, int blockRowsCnt, int blockColsCnt)
{
   int srcHeight, srcWidth;

   srcHeight = srcMat.rows;
   srcWidth = srcMat.cols;

   assert(srcMat.type() == CV_64FC1);

   Mat result(srcHeight, srcWidth, srcMat.type());

   int verticalPadSize = blockRowsCnt / 2;
   int horizontalPadSize = blockColsCnt / 2;   

   Mat zeroPad(srcHeight + 2*verticalPadSize, srcWidth + 2*horizontalPadSize, CV_64FC1);
   Mat tmp(blockRowsCnt * blockColsCnt, srcHeight * srcWidth, CV_64FC1);

   copyMakeBorder(srcMat, zeroPad, verticalPadSize, verticalPadSize, horizontalPadSize, horizontalPadSize, BORDER_CONSTANT, 0.0);

   // im2col function - sliding
   // http://matrix.etseq.urv.es/manuals/matlab/toolbox/images/block9.html   
   for (int x = 0; x < srcWidth; x++)
   {
      for (int y = 0; y < srcHeight; y++)
      {
         for (int i = x; i < x + blockColsCnt; i++)
         {
            for (int j = y; j < y + blockRowsCnt; j++)
            {                
                tmp.at<double>((i-x) * blockRowsCnt + (j-y), x * srcHeight + y) = zeroPad.at<double>(j, i);
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
      double mean = 0.0;
      for (int i = 0; i < tmp.rows; i++)
      {
         mean += colData[i];
      }

      mean /= (double)tmp.rows;

      double sum = 0.0;

      for (int i = 0; i < tmp.rows; i++)
      {
         sum += pow(colData[i]-mean, 2.0);
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

   tmp.deallocate();
   zeroPad.deallocate();

   return result;
}


/*
 * Set value by HSV color model.
 *
 * This function does not implement whole algorithm, but only Value part.
 * Article: Smith, A. R. “Color Gamut Transform Pairs”. SIGGRAPH 78 
 * Conference Proceedings. 1978, pp. 12–19.
 */
void MatUtil::setValueHSV(Mat &srcMat, Mat &lum)
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
 * If the value is lower than 1e-3, to actual luminance 
 * value is add 1e-3.
 */
void MatUtil::checkLumo(Mat &lumo)
{
   int height, width;

   height = lumo.rows;
   width = lumo.cols;

   double min;
   double max;
   Point minLoc;
   Point maxLoc;

   minMaxLoc(lumo, &min, &max, &minLoc, &maxLoc);

   if (min < 1e-3)
   {
      lumo += 1e-3;
   }
}

/*
 * Calculates the assigned percentile of matrix data.
 */
double MatUtil::percentile(Mat &data, double percentile)
{
   int index = (int) ceil(percentile / 100.0 * data.rows);

   return data.at<double>(index, 0);
}

/*
 * Cuts off unreasonable peaks and Valleys. 
 */
Mat MatUtil::cutPeakValley(Mat &rgb,
                           double lowEndCutPercentage,
                           double highEndCutPercentage)
{
   double cmpVal;
   double prcLowEnd = lowEndCutPercentage;
   double prcHighEnd = 100.0 - highEndCutPercentage;      

   double height = rgb.rows;
   double width = rgb.cols;   

   Mat linearRgb(height*width*CHANNELSCNT, 1, CV_64FC1);

   for (int y = 0; y < rgb.rows; y++)
   {
      for (int x = 0; x < rgb.cols; x++)
      {
         for (int c = 0; c < CHANNELSCNT; c++)
         {
            linearRgb.at<double>(y * rgb.cols * CHANNELSCNT + x * CHANNELSCNT + c, 0) =
               rgb.at<double>(y, x*CHANNELSCNT+c);
         }                  
      }      
   }
   
   Mat sorted(height*width*CHANNELSCNT, 1, CV_64FC1);
   cv::sort(linearRgb, sorted, SORT_EVERY_COLUMN + SORT_ASCENDING);

   double lowEnd = percentile(sorted, prcLowEnd);
   double highEnd = percentile(sorted, prcHighEnd);

   rgb -= lowEnd;
   rgb /= highEnd-lowEnd; 

   for (int y = 0; y < rgb.rows; y++)
   {
      for (int x = 0; x < rgb.cols; x++)
      {
         cmpVal = rgb.at<double>(y, x*CHANNELSCNT+2);

         if (cmpVal > 1.0) {
            rgb.at<double>(y, x*CHANNELSCNT+2) = 1.0;
         }
         else if (cmpVal < 0.0)
         {
            rgb.at<double>(y, x*CHANNELSCNT+2) = 0.0;
         }

         cmpVal = rgb.at<double>(y, x*CHANNELSCNT+1);

         if (cmpVal > 1.0) {
            rgb.at<double>(y, x*CHANNELSCNT+1) = 1.0;
         }
         else if (cmpVal < 0.0)
         {
            rgb.at<double>(y, x*CHANNELSCNT+1) = 0.0;
         }

         cmpVal = rgb.at<double>(y, x*CHANNELSCNT);

         if (cmpVal > 1.0) {
            rgb.at<double>(y, x*CHANNELSCNT) = 1.0;
         }
         else if (cmpVal < 0.0)
         {
            rgb.at<double>(y, x*CHANNELSCNT) = 0.0;
         }
      }      
   }
   
   return rgb;
}

/*
 * Normalize result image.
 */
void MatUtil::normImage(Mat &img)
{
  double min;
  double max;
  Point minLoc;
  Point maxLoc;

  minMaxLoc(img, &min, &max, &minLoc, &maxLoc);

  img = (img-min)/(max-min);
}
