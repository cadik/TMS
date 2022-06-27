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
 * Function for calculate differences.
 */
void MatUtil::diff(Mat &srcMat, Mat &dstMat, uint8_t dim)
{
   int height, width, channelsCnt;

   assert(dim == 2 || dim == 1);

   height      = srcMat.rows;
   width       = srcMat.cols;
   channelsCnt = srcMat.channels();

   if (dim == 2)
   {
      for (int y = 0; y < height; y++)
      {
         for (int x = 0; x < width - 1; x++)
         {
            for (int c = 0; c < channelsCnt; c++)
            {
               dstMat.at<double>(y, x * channelsCnt + c) = srcMat.at<double>(y, (x+1)*channelsCnt+c) - srcMat.at<double>(y, x*channelsCnt+c);
            }
         }
      }
   }
   else if (dim == 1)
   {
      for (int y = 0; y < height - 1; y++)
      {
         for (int x = 0; x < width; x++)
         {
            for (int c = 0; c < channelsCnt; c++)
            {
               dstMat.at<double>(y, x*channelsCnt + c) = srcMat.at<double>(y+1, x*channelsCnt + c) - srcMat.at<double>(y, x*channelsCnt + c);
            }
         }
      }
   }
}

/*
 * Function for convert subscripts to linear indices (for one dimension).
 */
Mat MatUtil::sub2indOneDim(Size indSize, Size resSize, Mat rows, Mat cols)
{
   int height, width;
   int indHeight, indWidth;
   Mat tmp;

   height = resSize.height;
   width = resSize.width;

   indHeight = indSize.height;
   indWidth = indSize.width;

   assert(rows.rows == height &&
          rows.cols == width &&
          cols.rows == height &&
          cols.cols == width);

   Mat resMat(height, width, CV_32SC1);

   multiply(cols, indHeight, tmp);

   resMat = tmp + rows;

   return resMat;
}

/*
 * Function for find first greater value than reference value.
 * If greater element was found, returns his index. Otherwise returns 0.
 */
int MatUtil::findFirstGtr(Mat srcMat, double cmpVal, int startIdx)
{
   int i;
   int width;

   width = srcMat.size().width;

   i = startIdx;

   while ((i < width) && (srcMat.at<double>(0, i) <= cmpVal))
   {
      i++;
   }

   return (i < width) ? (i - startIdx) : 0;
}

/*
 * Own implementation of MATLAB padarray function using OpenCV copyMakeBorder().
 *
 * Original MATLAB function: https://www.mathworks.com/help/images/ref/padarray.html
 */
Mat MatUtil::padArray (Mat srcMat, int rowPad, int colPad, uint8_t borderType)
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
