/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2022                                              *
*                                                                              *
*                       Implementation of the TMOGastal11 class                *
*                                                                              *
*******************************************************************************/

/*
 * Following code is inspired of original public MATLAB implementation
 * of article:
 *
 *    Domain Transform for Edge-Aware Image and Video Processing
 *    Eduardo S. L. Gastal  and  Manuel M. Oliveira
 *    ACM Transactions on Graphics. Volume 30 (2011), Number 4.
 *    Proceedings of SIGGRAPH 2011, Article 69.
 *
 * Source: http://inf.ufrgs.br/~eslgastal/DomainTransform/
 * Authors: Eduardo S. L. Gastal, Manuel M. Oliveira
 *
 * All rights of original implementation reserved to the authors
 * and any possible similarities are taken as citation of mentioned source.
 */

#include <limits>
#include <float.h>

#include "TMOGastal11.h"

#define CHANNELSCNT 3
#define NC 1
#define IC 2
#define RF 3
#define PADARR_REPLICATE 0
#define PADARR_CONSTANT 1

using namespace std;
using namespace cv;

TMOGastal11::TMOGastal11()
{
	SetName(L"Gastal11");
	SetDescription(L"Domain Transform for Edge-Aware Image and Video Processing");

   /* filterType parameter */
   filterType.SetName(L"filterType");
   filterType.SetDescription(L"Type of filter. Possible values are: 1 (NC), 2 (IC), 3 (RF)");
   filterType.SetDefault(1);
   filterType = 1;
   filterType.SetRange(1, 3);
   this->Register(filterType);

   /* sigma_r parameter */
   sigma_r.SetName(L"sigma_r");
   sigma_r.SetDescription(L"Sigma_r parameter: <0.0, 3.0>");
   sigma_r.SetDefault(0.4);
   sigma_r = 0.4;
   sigma_r.SetRange(0.0, 3.0);
   this->Register(sigma_r);

   /* sigma_s parameter */
   sigma_s.SetName(L"sigma_s");
   sigma_s.SetDescription(L"Sigma_s parameter: <0.0, 300.0>");
   sigma_s.SetDefault(60.0);
   sigma_s = 60.0;
   sigma_s.SetRange(0.0, 300.0);
   this->Register(sigma_s);

   /* numIter parameter */
   numIter.SetName(L"numIter");
   numIter.SetDescription(L"Number of iterations: <1, 10>");
   numIter.SetDefault(3);
   numIter = 3;
   numIter.SetRange(1, 10);
   this->Register(numIter);
}

TMOGastal11::~TMOGastal11()
{
}

/*
 * Convert TMOImage to cv::Mat
 */
Mat TMOImage2Mat(TMOImage* pSrc)
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
 * Function for calculate differences.
 */
void diff(Mat &srcMat, Mat &dstMat, uint8_t dim)
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
Mat sub2indOneDim(Size indSize, Size resSize, Mat rows, Mat cols)
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
int findFirstGtr(Mat srcMat, double cmpVal, int startIdx = 0)
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
 * Compute the domain transform (Equation 11).
 */
void domainTransform (Mat srcMat,
                      Mat &dHdx,
                      Mat &dVdy,
                      Mat &ct_H,
                      Mat &ct_V,
                      double sigma_s,
                      double sigma_r,
                      uint8_t filterOper)
{
   int height, width, channelsCnt;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();

   /* Horizontal and vertical partial derivatives using finite differences */
   Mat dIcdx(height, width-1, CV_64FC3);
   Mat dIcdy(height-1, width, CV_64FC3);

   diff(srcMat, dIcdx, 2);
   diff(srcMat, dIcdy, 1);

   Mat dIdx = Mat::zeros(height, width, CV_64FC1);
   Mat dIdy = Mat::zeros(height, width, CV_64FC1);

   /* Compute the l1-norm distance of neighbor pixels */
   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width - 1; x++)
      {
         for (int c = 0; c < channelsCnt; c++)
         {
            dIdx.at<double>(y, x+1) = dIdx.at<double>(y, x+1) + abs(dIcdx.at<double>(y, x*channelsCnt + c));
         }
      }
   }

   for (int y = 0; y < height-1; y++)
   {
      for (int x = 0; x < width; x++)
      {
         for (int c = 0; c < channelsCnt; c++)
         {
            dIdy.at<double>(y+1, x) = dIdy.at<double>(y+1, x) + abs(dIcdy.at<double>(y, x*channelsCnt + c));
         }
      }
   }

   /* Compute the derivatives of the horizontal and vertical domain transforms */
   multiply(dIdx, sigma_s/sigma_r, dHdx);
   multiply(dIdy, sigma_s/sigma_r, dVdy);

   dHdx = dHdx + 1.0;
   dVdy = dVdy + 1.0;

   if (filterOper == RF)
   {
       dVdy = dVdy.t();
   }
   else
   {
       /* NC or IC */

       /* Integrate the domain transforms */

       /* Cumulative sum dHdx, dim 2 */
       for (int y = 0; y < height; y++)
       {
          ct_H.at<double>(y, 0) = dHdx.at<double>(y, 0);

          for (int x = 1; x < width; x++)
          {
             ct_H.at<double>(y, x) = ct_H.at<double>(y, x-1) + dHdx.at<double>(y, x);
          }
       }

       /* Cumulative sum dVdy, dim 1 */
       for (int y = 0; y < width; y++)
       {
          ct_V.at<double>(0, y) = dVdy.at<double>(0, y);

          for (int x = 1; x < height; x++)
          {
             ct_V.at<double>(x, y) = ct_V.at<double>(x-1, y) + dVdy.at<double>(x, y);
          }
       }

       /* The vertical pass is performed using a transposed image */
       ct_V = ct_V.t();
   }
}

/*
 * Function for computer box filter.
 * Used for computing normalized or interpolated convolution.
 */
void boxFilter (Mat &srcMat,
                Mat &domainPosition,
                Mat &lowerPos,
                Mat &upperPos,
                Mat &lowerIdx,
                Mat &upperIdx,
                double boxRadius)
{
   int height, width, channelsCnt;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();

   lowerPos = domainPosition - boxRadius;
   upperPos = domainPosition + boxRadius;

   Mat domainPosRow(1, width+1, CV_64FC1);

   for (int y = 0; y < height; y++)
   {
      for (int i = 0; i < width; i++)
      {
         domainPosRow.at<double>(0, i) = domainPosition.at<double>(y, i);
      }

      domainPosRow.at<double>(0, width) = numeric_limits<double>::infinity();

      Mat lowerPosRow(1, width, CV_64FC1);
      Mat upperPosRow(1, width, CV_64FC1);

      for (int i = 0; i < width; i++)
      {
         lowerPosRow.at<double>(0, i) = lowerPos.at<double>(y, i);
         upperPosRow.at<double>(0, i) = upperPos.at<double>(y, i);
      }

      Mat localLowerIdx = Mat::zeros(1, width, CV_32SC1);
      Mat localUpperIdx = Mat::zeros(1, width, CV_32SC1);

      localLowerIdx.at<int>(0, 0) = findFirstGtr(domainPosRow, lowerPosRow.at<double>(0, 0));
      localUpperIdx.at<int>(0, 0) = findFirstGtr(domainPosRow, upperPosRow.at<double>(0, 0));

      for (int x = 1; x < width; x++)
      {
         localLowerIdx.at<int>(0, x) = localLowerIdx.at<int>(0, x-1) +
                                       findFirstGtr(domainPosRow,
                                                    lowerPosRow.at<double>(0, x),
                                                    localLowerIdx.at<int>(0, x-1));
         localUpperIdx.at<int>(0, x) = localUpperIdx.at<int>(0, x-1) +
                                       findFirstGtr(domainPosRow,
                                                    upperPosRow.at<double>(0, x),
                                                    localUpperIdx.at<int>(0, x-1));
      }

      for (int i = 0; i < width; i++)
      {
         lowerIdx.at<int>(y, i) = localLowerIdx.at<int>(0, i);
         upperIdx.at<int>(y, i) = localUpperIdx.at<int>(0, i);
      }
   }
}

/*
 * Function for compute normalized convolution.
 */
void NCfilter (Mat srcMat, Mat &outMat, Mat &domainPosition, double radius)
{
   int height, width, channelsCnt;
   int bVal, aVal;
   int bCol, bRow;
   int aCol, aRow;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();

   Mat resMat(height, width, CV_64FC3);

   Mat lowerPos(height, width, CV_64FC1);
   Mat upperPos(height, width, CV_64FC1);

   Mat lowerIdx = Mat::zeros(height, width, CV_32SC1);
   Mat upperIdx = Mat::zeros(height, width, CV_32SC1);

   boxFilter(srcMat, domainPosition, lowerPos, upperPos, lowerIdx, upperIdx, radius);

   // Compute the box filter using a summed area table
   Mat SAT = Mat::zeros(height, width+1, CV_64FC3);

   /* cumsum */
   for (int i = 0; i < height; i++)
   {
      SAT.at<double>(i, channelsCnt)   = srcMat.at<double>(i, 0);
      SAT.at<double>(i, channelsCnt+1) = srcMat.at<double>(i, 1);
      SAT.at<double>(i, channelsCnt+2) = srcMat.at<double>(i, 2);

      for (int j = 2; j < SAT.size().width; j++)
      {
         for (int c = 0; c < channelsCnt; c++)
         {
            SAT.at<double>(i, j*channelsCnt+c) = srcMat.at<double>(i, (j-1)*channelsCnt+c) + SAT.at<double>(i, (j-1)*channelsCnt+c);

         }
      }
   }

   Mat rowIndices(height, width, CV_32SC1);

   /* repmat */
   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         rowIndices.at<int>(y, x) = y;
      }
   }

   Mat a(height, width, CV_32SC1);
   Mat b(height, width, CV_32SC1);

   Size SATSize = SAT.size();
   Size resSize = srcMat.size();

   a = sub2indOneDim(SATSize, resSize, rowIndices, lowerIdx);
   b = sub2indOneDim(SATSize, resSize, rowIndices, upperIdx);

   for (int c = 0; c < channelsCnt; c++)
   {
      for (int y = 0; y < height; y++)
      {
         for (int x = 0; x < width; x++)
         {
            bVal = b.at<int>(y,x);

            bCol = bVal / height;
            bRow = bVal % height;

            aVal = a.at<int>(y, x);

            aCol = aVal / height;
            aRow = aVal % height;

            resMat.at<double>(y, x*channelsCnt+c) = (SAT.at<double>(bRow, bCol*channelsCnt+c) - SAT.at<double>(aRow, aCol*channelsCnt+c)) /
                    (double)(upperIdx.at<int>(y, x) - lowerIdx.at<int>(y, x));
         }
      }
   }

   resMat.copyTo(outMat);
}

/*
 * Function for compute interpolated convolution.
 */
void ICfilter (Mat srcMat, Mat &outMat, Mat &domainPosition, double radius)
{
   double C, L, R;
   double alpha, yi;
   int l1_val, l1_col, l1_row;
   int u0_val, u0_col, u0_row;
   int l0_val, l0_col, l0_row;
   int u1_val, u1_col, u1_row;
   int height, width, channelsCnt;
   int domainPosHeight, domainPosWidth;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();

   Mat resMat(height, width, CV_64FC3);

   Mat lowerPos(height, width, CV_64FC1);
   Mat upperPos(height, width, CV_64FC1);

   Mat lowerIdx = Mat::zeros(height, width, CV_32SC1);
   Mat upperIdx = Mat::zeros(height, width, CV_32SC1);

   lowerPos = domainPosition - radius;
   upperPos = domainPosition + radius;

   boxFilter(srcMat, domainPosition, lowerPos, upperPos, lowerIdx, upperIdx, radius);

   // Compute the box filter using a summed area table. This SAT is built using
   // the area under the graph (in the transformed domain) of the interpolated
   // signal. We use linear interpolation and compute the area using the
   // trapezoidal rule.

   Mat img(height, width-1, CV_64FC3);
   Mat domain(height, width-1, CV_64FC3);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width-1; x++)
      {
         double domainVal = domainPosition.at<double>(y,x+1) - domainPosition.at<double>(y, x);

         for (int c = 0; c < channelsCnt; c++)
         {
            img.at<double>(y, x*channelsCnt+c) = 0.5 * (srcMat.at<double>(y, (x+1)*channelsCnt+c) + srcMat.at<double>(y, x*channelsCnt+c));
            domain.at<double>(y, x*channelsCnt+c) = domainVal;
         }
      }
   }

   Mat areas(height, width-1, CV_64FC3);

   multiply(img, domain, areas);

   Mat SAT = Mat::zeros(height, width, CV_64FC3);

   /* cumsum */
   for (int i = 0; i < height; i++)
   {
      SAT.at<double>(i, channelsCnt)   = areas.at<double>(i, 0);
      SAT.at<double>(i, channelsCnt+1) = areas.at<double>(i, 1);
      SAT.at<double>(i, channelsCnt+2) = areas.at<double>(i, 2);

      for (int j = 2; j < width; j++)
      {
         for (int c = 0; c < channelsCnt; c++)
         {
            SAT.at<double>(i, j*channelsCnt+c) = areas.at<double>(i, (j-1)*channelsCnt+c) + SAT.at<double>(i, (j-1)*channelsCnt+c);
         }
      }
   }

   Mat rowIndices(height, width, CV_32SC1);

   /* repmat */
   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         rowIndices.at<int>(y, x) = y;
      }
   }

   srcMat = padArray(srcMat, 0, 1, PADARR_REPLICATE);
   SAT = padArray(SAT, 0, 1, PADARR_CONSTANT);
   domainPosition = padArray(domainPosition, 0, 1, PADARR_REPLICATE);

   // Pixel values outside the bounds of the image are assumed to equal the
   // nearest pixel border value.

   domainPosHeight = domainPosition.rows;
   domainPosWidth = domainPosition.cols;

   for (int y = 0; y < height; y++) {
      domainPosition.at<double>(y, 0) = domainPosition.at<double>(y, 0) - 1.2 * radius;
      domainPosition.at<double>(y, domainPosWidth-1) = domainPosition.at<double>(y, domainPosWidth-1) + 1.2*radius;
   }

   lowerIdx = lowerIdx + 1;

   Mat l1 = sub2indOneDim(SAT.size(), resMat.size(), rowIndices, lowerIdx);
   Mat u0 = sub2indOneDim(SAT.size(), resMat.size(), rowIndices, upperIdx);

   Mat l0 = sub2indOneDim(SAT.size(), resMat.size(), rowIndices, lowerIdx-1);
   Mat u1 = sub2indOneDim(SAT.size(), resMat.size(), rowIndices, upperIdx+1);

   for (int c = 0; c < channelsCnt; c++)
   {
      for (int y = 0; y < height; y++)
      {
         for (int x = 0; x < width; x++)
         {
            l1_val = l1.at<int>(y, x);
            l1_col = l1_val / height;
            l1_row = l1_val % height;

            u0_val = u0.at<int>(y, x);
            u0_col = u0_val / height;
            u0_row = u0_val % height;

            l0_val = l0.at<int>(y, x);
            l0_col = l0_val / height;
            l0_row = l0_val % height;

            u1_val = u1.at<int>(y, x);
            u1_col = u1_val / height;
            u1_row = u1_val % height;

            C = SAT.at<double>(u0_row, u0_col*channelsCnt+c) - SAT.at<double>(l1_row, l1_col*channelsCnt+c);

            // Left fractional areas
            alpha = (lowerPos.at<double>(y, x) - domainPosition.at<double>(l0_row, l0_col)) /
               (domainPosition.at<double>(l1_row, l1_col) - domainPosition.at<double>(l0_row, l0_col));

            yi = srcMat.at<double>(l0_row, l0_col*channelsCnt+c) + alpha *
               (srcMat.at<double>(l1_row, l1_col*channelsCnt+c) - srcMat.at<double>(l0_row, l0_col*channelsCnt+c));

            L = 0.5 * (yi + srcMat.at<double>(l1_row, l1_col*channelsCnt+c)) * (1.0 - alpha) *
               (domainPosition.at<double>(l1_row, l1_col) - domainPosition.at<double>(l0_row, l0_col));

            // Right fractional areas
            alpha = (upperPos.at<double>(y, x) - domainPosition.at<double>(u0_row, u0_col)) /
               (domainPosition.at<double>(u1_row, u1_col) - domainPosition.at<double>(u0_row, u0_col));

            yi = srcMat.at<double>(u0_row, u0_col*channelsCnt+c) + alpha *
               (srcMat.at<double>(u1_row, u1_col*channelsCnt+c) - srcMat.at<double>(u0_row, u0_col*channelsCnt+c));

            R = 0.5 * (yi + srcMat.at<double>(u0_row, u0_col*channelsCnt+c)) * alpha *
               (domainPosition.at<double>(u1_row, u1_col) - domainPosition.at<double>(u0_row, u0_col));

            resMat.at<double>(y, x*channelsCnt+c) = (L + C + R) / (2 * radius);
         }
      }
   }

   resMat.copyTo(outMat);
}

/*
 * Function for compute recursive filter.
 */
void Rfilter(Mat srcMat, Mat &outMat, Mat &domainPosition, double sigma_H)
{
   int height, width, channelsCnt;

   height = srcMat.size().height;
   width = srcMat.size().width;
   channelsCnt = srcMat.channels();

   double a = exp((-1.0 * sqrt(2.0)) / sigma_H);

   Mat resMat(height, width, CV_64FC3);
   Mat V = Mat(height, width, CV_64FC1);

   srcMat.copyTo(resMat);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
          V.at<double>(y, x) = pow(a, domainPosition.at<double>(y, x));
      }
   }

   for(int y = 0; y < height; y++)
   {
      for(int x = 1; x < width; x++)
      {
         for(int c = 0; c < channelsCnt; c++)
         {
            resMat.at<double>(y, x*channelsCnt+c) = resMat.at<double>(y, x*channelsCnt+c) +
               V.at<double>(y,x) * (resMat.at<double>(y, (x-1)*channelsCnt+c) - resMat.at<double>(y, x*channelsCnt+c));
         }
      }
   }

   for(int y = 0; y < height; y++)
   {
      for(int x = width-2; x >= 0; x--)
      {
         for(int c = 0; c < channelsCnt; c++)
         {
            resMat.at<double>(y, x*channelsCnt+c) = resMat.at<double>(y, x*channelsCnt+c) +
               V.at<double>(y,x+1) * (resMat.at<double>(y, (x+1)*channelsCnt+c) - resMat.at<double>(y, x*channelsCnt+c));
         }
      }
   }

   resMat.copyTo(outMat);
}

/*
 * Function for handle NC filter type.
 */
void filterOperationNC(Mat srcMat, Mat &resMat, double sigma_s, double sigma_r, uint8_t numIter)
{
    int height, width;
    double multiplier, divider;
    double sigma_H_i;
    double boxRadius;
    double sqrtOf3 = sqrt(3.0);

    height = srcMat.size().height;
    width  = srcMat.size().width;

    Mat matT(height, width, CV_64FC3);
    Mat outMat(height, width, CV_64FC3);

    Mat dHdx(height, width, CV_64FC1);
    Mat dVdy(height, width, CV_64FC1);

    Mat ct_H(height, width, CV_64FC1);
    Mat ct_V(height, width, CV_64FC1);

    domainTransform (srcMat, dHdx, dVdy, ct_H, ct_V, sigma_s, sigma_r, NC);

    multiplier = sigma_s * sqrt(3.0); /* sigma_H = sigma_s */
    divider = sqrt(pow(4.0, numIter) - 1.0);

    matT = srcMat;

    for (int i = 0; i < numIter; i++)
    {
       sigma_H_i = multiplier * pow(2.0, numIter - (i+1)) / divider;
       boxRadius = sqrtOf3 * sigma_H_i;

       NCfilter(matT, outMat, ct_H, boxRadius);
       matT = outMat.t();

       NCfilter(matT, outMat, ct_V, boxRadius);
       matT = outMat.t();
    }

    resMat = matT.clone();
}

/*
 * Function for handle IC filter type.
 */
void filterOperationIC(Mat srcMat,
                       Mat &resMat,
                       double sigma_s,
                       double sigma_r,
                       uint8_t numIter)
{
   int height, width;
   double multiplier, divider;
   double sigma_H_i;
   double boxRadius;
   double sqrtOf3 = sqrt(3.0);

   height = srcMat.size().height;
   width  = srcMat.size().width;

   Mat matT(height, width, CV_64FC3);
   Mat outMat(height, width, CV_64FC3);

   Mat dHdx(height, width, CV_64FC1);
   Mat dVdy(height, width, CV_64FC1);

   Mat ct_H(height, width, CV_64FC1);
   Mat ct_V(height, width, CV_64FC1);

   domainTransform (srcMat, dHdx, dVdy, ct_H, ct_V, sigma_s, sigma_r, NC);

   multiplier = sigma_s * sqrt(3.0); // sigma_H = sigma_s
   divider = sqrt(pow(4.0, numIter) - 1.0);

   matT = srcMat;

   for (int i = 0; i < numIter; i++)
   {
      sigma_H_i = multiplier * pow(2.0, numIter - (i+1)) / divider;
      boxRadius = sqrtOf3 * sigma_H_i;

      ICfilter(matT, outMat, ct_H, boxRadius);
      matT = outMat.t();

      ICfilter(matT, outMat, ct_V, boxRadius);
      matT = outMat.t();
   }

   resMat = matT.clone();
}

/*
 * Function for handle RF filter type.
 */
void filterOperationRF(Mat srcMat,
                       Mat &resMat,
                       double sigma_s,
                       double sigma_r,
                       uint8_t numIter)
{
    int height, width;
    double multiplier, divider;
    double sigma_H_i;

    height = srcMat.size().height;
    width  = srcMat.size().width;

    Mat matT(height, width, CV_64FC3);
    Mat outMat(height, width, CV_64FC3);

    Mat dHdx(height, width, CV_64FC1);
    Mat dVdy(height, width, CV_64FC1);

    Mat ct_H, ct_V;

    domainTransform (srcMat, dHdx, dVdy, ct_H, ct_V, sigma_s, sigma_r, RF);

    multiplier = sigma_s * sqrt(3.0); /* sigma_H = sigma_s */
    divider = sqrt (pow(4.0, numIter) - 1.0);
    matT = srcMat;

    for (int i = 0; i < numIter; i++)
    {
        sigma_H_i = multiplier * pow(2.0, numIter - (i+1)) / divider;

        Rfilter(matT, outMat, dHdx, sigma_H_i);
        matT = outMat.t();

        Rfilter(matT, outMat, dVdy, sigma_H_i);
        matT = outMat.t();
    }

    resMat = matT.clone();
}

/*
 * Function for handle type of filter (NC, IC or RF).
 */
void filterOperation(Mat srcMat, Mat &resMat, double sigma_s, double sigma_r, uint8_t filterType, uint8_t numIter)
{
    switch (filterType) {
        case NC:
            filterOperationNC(srcMat, resMat, sigma_s, sigma_r, numIter);
            break;
        case IC:
            filterOperationIC(srcMat, resMat, sigma_s, sigma_r, numIter);
            break;
        case RF:
            filterOperationRF(srcMat, resMat, sigma_s, sigma_r, numIter);
            break;
        default:
            break;
    }
}

/*
 * Main method for Domain Transform for Edge-Aware Image and Video Processing operator
 */
int TMOGastal11::Transform()
{
   int height, width, channelsCnt;
   uint8_t pNumIter, pFilterType;
   double pSigmaS, pSigmaR;
   double* pDestData;

   /***********************/
   /* Get user parameters */
   /***********************/

   pNumIter = numIter.GetInt();
   pFilterType = filterType.GetInt();
   pSigmaS = sigma_s.GetDouble();
   pSigmaR = sigma_r.GetDouble();

   Mat srcMat = TMOImage2Mat(pSrc);

   Mat result(srcMat.size().height, srcMat.size().width, CV_64FC3);

   filterOperation(srcMat, result, pSigmaS, pSigmaR, pFilterType, pNumIter);

   width = result.size().width;
   height = result.size().height;
   channelsCnt = result.channels();

   pDestData = pDst->GetData();

   /*
	 * Save result to the destination image
	 */
	int y = 0;
	for (; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			*pDestData++ = result.at<double>(y, x*channelsCnt+2);
			*pDestData++ = result.at<double>(y, x*channelsCnt+1);
			*pDestData++ = result.at<double>(y, x*channelsCnt);
		}
	}

	return 0;
}
