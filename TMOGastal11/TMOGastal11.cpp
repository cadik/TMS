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

using namespace std;
using namespace cv;

TMOGastal11::TMOGastal11()
{
	SetName(L"Gastal11");
	SetDescription(L"Domain Transform for Edge-Aware Image and Video Processing");

   /* pencilColor parameter */
   pencilColor.SetName(L"pencilColor");
   pencilColor.SetDescription(L"Color version of pencil drawing.");
   pencilColor.SetDefault(false);
   pencilColor = false;
   this->Register(pencilColor);

   /* shadeFactor parameter */
   shadeFactor.SetName(L"shadeFactor");
   shadeFactor.SetDescription(L"Used only for pencil drawing.\n"
                              "\t\tShade factor scale the output image intensity.\n"
                              "\t\tThe higher the value, the brighter is the result (0.0 - 0.1).");
   shadeFactor.SetDefault(0.05);
   shadeFactor = 0.05;
   shadeFactor.SetRange(0.0, 0.1);
   this->Register(shadeFactor);

   /* filter application */
   /*
    * 0 - edge aware smoothing (usage of basic filters)
    * 1 - detail enhancement
    * 2 - stylization
    * 3 - pencil drawing
    * 4 - tone mapping
    */
   filterAppl.SetName(L"filterAppl");
   filterAppl.SetDescription(L"Application of filter. "
                             "Possible values are:\n"
                             "\t\t0 (Edge aware smoothing - basic filtering),\n"
                             "\t\t1 (Detail enhancement),\n"
                             "\t\t2 (Stylization),\n"
                             "\t\t3 (Pencil drawing - filterType is unused, "
                             "recommended sigma_s = 10.0, sigma_r = 0.1),\n"
                             "\t\t4 (Tone Mapping) - only numIter parameter is used");
   filterAppl.SetDefault(4);
   filterAppl = 4;
   filterAppl.SetRange(0, 4);
   this->Register(filterAppl);

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
Mat TMOGastal11::TMOImage2Mat(TMOImage* pSrc)
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
            srcConvMat.at<double>(y, x*CHANNELSCNT+c) =
               pSourceData[CHANNELSCNT-c-1];
         }

         /* Add count of channels (RGB) to pointer */
         pSourceData += CHANNELSCNT;
      }
   }

   return srcConvMat;
}

/*
 * Compute the domain transform (Equation 11).
 */
void TMOGastal11::domainTransform (Mat srcMat,
                                   Mat &dHdx,
                                   Mat &dVdy,
                                   Mat &ct_H,
                                   Mat &ct_V,
                                   double sigma_s,
                                   double sigma_r,
                                   uint8_t filterOper)
{
   int height, width, channelsCnt, type;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();
   type = srcMat.type();

   /* Horizontal and vertical partial derivatives using finite differences */
   Mat dIcdx(height, width-1, type);
   Mat dIcdy(height-1, width, type);

   matUtil.diff(srcMat, dIcdx, 2);
   matUtil.diff(srcMat, dIcdy, 1);

   Mat dIdx = Mat::zeros(height, width, CV_64FC1);
   Mat dIdy = Mat::zeros(height, width, CV_64FC1);

   /* Compute the l1-norm distance of neighbor pixels */
   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width - 1; x++)
      {
         for (int c = 0; c < channelsCnt; c++)
         {
            dIdx.at<double>(y, x+1) = dIdx.at<double>(y, x+1) +
               abs(dIcdx.at<double>(y, x*channelsCnt + c));
         }
      }
   }

   for (int y = 0; y < height-1; y++)
   {
      for (int x = 0; x < width; x++)
      {
         for (int c = 0; c < channelsCnt; c++)
         {
            dIdy.at<double>(y+1, x) = dIdy.at<double>(y+1, x) +
               abs(dIcdy.at<double>(y, x*channelsCnt + c));
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
             ct_H.at<double>(y, x) = ct_H.at<double>(y, x-1) +
               dHdx.at<double>(y, x);
          }
       }

       /* Cumulative sum dVdy, dim 1 */
       for (int y = 0; y < width; y++)
       {
          ct_V.at<double>(0, y) = dVdy.at<double>(0, y);

          for (int x = 1; x < height; x++)
          {
             ct_V.at<double>(x, y) = ct_V.at<double>(x-1, y) +
               dVdy.at<double>(x, y);
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
void TMOGastal11::boxFilter (Mat &srcMat,
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

      localLowerIdx.at<int>(0, 0) =
         matUtil.findFirstGtr(domainPosRow, lowerPosRow.at<double>(0, 0));
      localUpperIdx.at<int>(0, 0) =
         matUtil.findFirstGtr(domainPosRow, upperPosRow.at<double>(0, 0));

      for (int x = 1; x < width; x++)
      {
         localLowerIdx.at<int>(0, x) = localLowerIdx.at<int>(0, x-1) +
                                       matUtil.findFirstGtr(domainPosRow,
                                                    lowerPosRow.at<double>(0, x),
                                                    localLowerIdx.at<int>(0, x-1));
         localUpperIdx.at<int>(0, x) = localUpperIdx.at<int>(0, x-1) +
                                       matUtil.findFirstGtr(domainPosRow,
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
void TMOGastal11::NCfilter (Mat srcMat, Mat &outMat, Mat &domainPosition, double radius)
{
   int height, width, channelsCnt, type;
   int bVal, aVal;
   int bCol, bRow;
   int aCol, aRow;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();
   type = srcMat.type();

   Mat resMat(height, width, type);

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
            SAT.at<double>(i, j*channelsCnt+c) =srcMat.at<double>(i, (j-1)*channelsCnt+c) +
               SAT.at<double>(i, (j-1)*channelsCnt+c);

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

   a = matUtil.sub2indOneDim(SATSize, resSize, rowIndices, lowerIdx);
   b = matUtil.sub2indOneDim(SATSize, resSize, rowIndices, upperIdx);

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

            resMat.at<double>(y, x*channelsCnt+c) =
               (SAT.at<double>(bRow, bCol*channelsCnt+c) - SAT.at<double>(aRow, aCol*channelsCnt+c)) /
                  (double)(upperIdx.at<int>(y, x) - lowerIdx.at<int>(y, x));
         }
      }
   }

   resMat.copyTo(outMat);
}

/*
 * Function for compute interpolated convolution.
 */
void TMOGastal11::ICfilter (Mat srcMat, Mat &outMat, Mat &domainPosition, double radius)
{
   double C, L, R;
   double alpha, yi;
   int l1_val, l1_col, l1_row;
   int u0_val, u0_col, u0_row;
   int l0_val, l0_col, l0_row;
   int u1_val, u1_col, u1_row;
   int height, width, channelsCnt, type;
   int domainPosHeight, domainPosWidth;

   width  = srcMat.size().width;
   height = srcMat.size().height;
   channelsCnt = srcMat.channels();
   type = srcMat.type();

   Mat resMat(height, width, type);

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

   Mat img(height, width-1, type);
   Mat domain(height, width-1, type);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width-1; x++)
      {
         double domainVal = domainPosition.at<double>(y,x+1) - domainPosition.at<double>(y, x);

         for (int c = 0; c < channelsCnt; c++)
         {
            img.at<double>(y, x*channelsCnt+c) =
               0.5 * (srcMat.at<double>(y, (x+1)*channelsCnt+c) + srcMat.at<double>(y, x*channelsCnt+c));
            domain.at<double>(y, x*channelsCnt+c) = domainVal;
         }
      }
   }

   Mat areas(height, width-1, type);

   multiply(img, domain, areas);

   Mat SAT = Mat::zeros(height, width, type);

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
            SAT.at<double>(i, j*channelsCnt+c) = areas.at<double>(i, (j-1)*channelsCnt+c) +
               SAT.at<double>(i, (j-1)*channelsCnt+c);
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

   srcMat = matUtil.padArray(srcMat, 0, 1, PADARR_REPLICATE);
   SAT = matUtil.padArray(SAT, 0, 1, PADARR_CONSTANT);
   domainPosition = matUtil.padArray(domainPosition, 0, 1, PADARR_REPLICATE);

   // Pixel values outside the bounds of the image are assumed to equal the
   // nearest pixel border value.

   domainPosHeight = domainPosition.rows;
   domainPosWidth = domainPosition.cols;

   for (int y = 0; y < height; y++) {
      domainPosition.at<double>(y, 0) =
         domainPosition.at<double>(y, 0) - 1.2 * radius;
      domainPosition.at<double>(y, domainPosWidth-1) =
         domainPosition.at<double>(y, domainPosWidth-1) + 1.2 * radius;
   }

   lowerIdx = lowerIdx + 1;

   Mat l1 = matUtil.sub2indOneDim(SAT.size(), resMat.size(), rowIndices, lowerIdx);
   Mat u0 = matUtil.sub2indOneDim(SAT.size(), resMat.size(), rowIndices, upperIdx);

   Mat l0 = matUtil.sub2indOneDim(SAT.size(), resMat.size(), rowIndices, lowerIdx-1);
   Mat u1 = matUtil.sub2indOneDim(SAT.size(), resMat.size(), rowIndices, upperIdx+1);

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
void TMOGastal11::Rfilter(Mat srcMat, Mat &outMat, Mat &domainPosition, double sigma_H)
{
   int height, width, channelsCnt, type;

   height = srcMat.size().height;
   width = srcMat.size().width;
   channelsCnt = srcMat.channels();
   type = srcMat.type();

   double a = exp((-1.0 * sqrt(2.0)) / sigma_H);

   Mat resMat(height, width, type);
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
void TMOGastal11::filterOperationNC(Mat srcMat,
                                    Mat &resMat,
                                    double sigma_s,
                                    double sigma_r,
                                    uint8_t numIter)
{
    int height, width, type;
    double multiplier, divider;
    double sigma_H_i;
    double boxRadius;
    double sqrtOf3 = sqrt(3.0);

    height = srcMat.size().height;
    width  = srcMat.size().width;
    type   = srcMat.type();

    Mat matT(height, width, type);
    Mat outMat(height, width, type);

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
void TMOGastal11::filterOperationIC(Mat srcMat,
                                    Mat &resMat,
                                    double sigma_s,
                                    double sigma_r,
                                    uint8_t numIter)
{
   int height, width, type;
   double multiplier, divider;
   double sigma_H_i;
   double boxRadius;
   double sqrtOf3 = sqrt(3.0);

   height = srcMat.size().height;
   width  = srcMat.size().width;
   type   = srcMat.type();

   Mat matT(height, width, type);
   Mat outMat(height, width, type);

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
void TMOGastal11::filterOperationRF(Mat srcMat,
                                    Mat &resMat,
                                    double sigma_s,
                                    double sigma_r,
                                    uint8_t numIter)
{
    int height, width, type;
    double multiplier, divider;
    double sigma_H_i;

    height = srcMat.size().height;
    width  = srcMat.size().width;
    type   = srcMat.type();

    Mat matT(height, width, type);
    Mat outMat(height, width, type);

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
void TMOGastal11::filterOperation(Mat srcMat,
                                  Mat &resMat,
                                  double sigma_s,
                                  double sigma_r,
                                  uint8_t filterType,
                                  uint8_t numIter)
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
 * Function for compute edge aware smoothing.
 *
 * Function computes basic filtering with NC, IC or RF.
 */
void TMOGastal11::edgeAwareSmoothing(Mat srcMat,
                                     Mat &resMat,
                                     double sigma_s,
                                     double sigma_r,
                                     uint8_t filterType,
                                     uint8_t numIter)
{
   filterOperation(srcMat, resMat, sigma_s, sigma_r, filterType, numIter);
}

/*
 * Function for compute detail enhancement as a concrete usage of filters.
 *
 * Function is based on an OpenCV implementation of detail enhancement
 * Source: https://github.com/opencv/opencv/blob/17234f82d025e3bbfbf611089637e5aa2038e7b8/modules/photo/src/npr.cpp
 * Author: Siddharthk (https://github.com/Siddharthk)
 * Original implementation page: https://www.inf.ufrgs.br/~eslgastal/DomainTransform/Detail_Manipulation/index.html
 */
void TMOGastal11::detailEnhancement(Mat srcMat,
                                    Mat &resMat,
                                    double sigma_s,
                                    double sigma_r,
                                    uint8_t filterType,
                                    uint8_t numIter)
{
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   Mat srcMatLab(height, width, CV_64FC3);
   vector<Mat> sepChannel;

   converter.bgr2Lab64FType(srcMat, srcMatLab);

   split(srcMatLab, sepChannel);

   Mat L = sepChannel[0];

   Mat outMat;

   filterOperation(L, outMat, sigma_s, sigma_r, filterType, numIter);

   Mat D = L - outMat;

   multiply(D, 3.0, D);
   L = outMat + D;

   sepChannel[0] = L;

   merge(sepChannel, srcMatLab);

   converter.lab2Bgr64FType(srcMatLab, resMat);
}

/*
 * Function for compute magnitute.
 *
 * Function is based on an OpenCV implementation of stylization
 * Source: https://github.com/opencv/opencv/blob/17234f82d025e3bbfbf611089637e5aa2038e7b8/modules/photo/src/npr.hpp
 * Author: Siddharthk (https://github.com/Siddharthk)
 */
void TMOGastal11::findMagnitude(Mat srcMat, Mat &mag)
{
   int height = srcMat.rows;
   int width = srcMat.cols;

   vector<Mat> planes;
   split(srcMat, planes);

   Mat magXR = Mat(height, width, CV_64FC1);
   Mat magYR = Mat(height, width, CV_64FC1);

   Mat magXG = Mat(height, width, CV_64FC1);
   Mat magYG = Mat(height, width, CV_64FC1);

   Mat magXB = Mat(height, width, CV_64FC1);
   Mat magYB = Mat(height, width, CV_64FC1);

   Sobel(planes[0], magXR, CV_64FC1, 1, 0, 3);
   Sobel(planes[0], magYR, CV_64FC1, 0, 1, 3);

   Sobel(planes[1], magXG, CV_64FC1, 1, 0, 3);
   Sobel(planes[1], magYG, CV_64FC1, 0, 1, 3);

   Sobel(planes[2], magXB, CV_64FC1, 1, 0, 3);
   Sobel(planes[2], magYB, CV_64FC1, 0, 1, 3);

   Mat mag1 = Mat(height,width,CV_64FC1);
   Mat mag2 = Mat(height,width,CV_64FC1);
   Mat mag3 = Mat(height,width,CV_64FC1);

   magnitude(magXR,magYR,mag1);
   magnitude(magXG,magYG,mag2);
   magnitude(magXB,magYB,mag3);

   mag = mag1 + mag2 + mag3;
   mag = 1.0f - mag;
}

/*
 * Function for compute stylization as a concrete usage of filters.
 *
 * Function is based on an OpenCV implementation of stylization
 * Source: https://github.com/opencv/opencv/blob/17234f82d025e3bbfbf611089637e5aa2038e7b8/modules/photo/src/npr.hpp
 * Author: Siddharthk (https://github.com/Siddharthk)
 * Original implementation page: https://www.inf.ufrgs.br/~eslgastal/DomainTransform/Stylization/index.html
 */
void TMOGastal11::stylization(Mat srcMat,
                              Mat &resMat,
                              double sigma_s,
                              double sigma_r,
                              uint8_t filterType,
                              uint8_t numIter)
{
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   Mat out(height, width, CV_64FC3);

   filterOperation(srcMat, out, sigma_s, sigma_r, filterType, numIter);

   Mat mag(height, width, CV_64FC1);

   findMagnitude(out, mag);

   Mat stylized;

   vector <Mat> tmp;
   split(out, tmp);

   multiply(tmp[0],mag,tmp[0]);
   multiply(tmp[1],mag,tmp[1]);
   multiply(tmp[2],mag,tmp[2]);
   merge(tmp,stylized);

   stylized.copyTo(resMat);
}

/*
 * Function for compute pencil sketch as a concrete usage of filters.
 *
 * Function is based on an OpenCV implementation of pencil sketch
 * Source: https://github.com/opencv/opencv/blob/17234f82d025e3bbfbf611089637e5aa2038e7b8/modules/photo/src/npr.hpp
 * Author: Siddharthk (https://github.com/Siddharthk)
 * Original implementation page: https://www.inf.ufrgs.br/~eslgastal/DomainTransform/Pencil_Drawing/index.html
 */
void TMOGastal11::pencilSketch(Mat srcMat,
                               Mat &sketch,
                               Mat &colorSketch,
                               double sigma_s,
                               double sigma_r,
                               double shadeFactor,
                               bool color,
                               uint8_t numIter)
{
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   Mat colorSketchTmp(height, width, CV_64FC3);

   Mat dHdx(height, width, CV_64FC1);
   Mat dVdy(height, width, CV_64FC1);

   Mat ct_H(height, width, CV_64FC1);
   Mat ct_V(height, width, CV_64FC1);

   domainTransform(srcMat, dHdx, dVdy, ct_H, ct_V, sigma_s, sigma_r, NC);

   converter.bgr2YCrCb64FType(srcMat, colorSketchTmp);

   vector <Mat> YUV_channel;

   Mat penx = Mat(height, width, CV_32SC1);
   Mat peny = Mat(width, height, CV_32SC1);
   Mat penRes = Mat::zeros(height, width, CV_64FC1);
   Mat penyT;

   double boxRadius;

   double sqrtOf3 = sqrt(3.0);

   Mat lowerPos(height, width, CV_64FC1);
   Mat upperPos(height, width, CV_64FC1);

   Mat lowerIdx = Mat::zeros(height, width, CV_32SC1);
   Mat upperIdx = Mat::zeros(height, width, CV_32SC1);
   Mat lowerIdxV = Mat::zeros(width, height, CV_32SC1);
   Mat upperIdxV = Mat::zeros(width, height, CV_32SC1);

   double multiplier = sigma_s * sqrtOf3; /* sigma_H = sigma_s */

   double divider = sqrt (pow(4.0, numIter) - 1.0);
   double sigma_H_i;

   Mat matT;
   matT = srcMat;

   sigma_H_i = (multiplier * pow(2.0, numIter - 1)) / divider;
   boxRadius = sqrtOf3 * sigma_H_i;

   boxFilter(matT, ct_H, lowerPos, upperPos, lowerIdx, upperIdx, boxRadius);
   matT = matT.t();
   penx = upperIdx - lowerIdx;

   boxFilter(matT, ct_V, lowerPos, upperPos, lowerIdxV, upperIdxV, boxRadius);
   matT = matT.t();
   peny = upperIdxV - lowerIdxV;
   penyT = peny.t();

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         penRes.at<double>(y, x) = shadeFactor *
            (double)(penx.at<int>(y, x) + penyT.at<int>(y, x));
      }
   }

   sketch = penRes.clone();
   split(colorSketchTmp, YUV_channel);
   penRes.copyTo(YUV_channel[0]);
   merge(YUV_channel, colorSketchTmp);
   converter.yCrCb2Bgr64FType(colorSketchTmp, colorSketch);
}

/*
 * Function for compute pencil drawing.
 *
 * The main task of this function is choose version of pencil sketch (colored or grayscale).
 */
void TMOGastal11::pencilDrawing(Mat srcMat,
                                Mat &resMat,
                                double sigma_s,
                                double sigma_r,
                                double shadeFactor,
                                bool color,
                                uint8_t numIter)
{
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   Mat sketch(height, width, CV_64FC1);
   Mat colorSketch(height, width, CV_64FC3);

   pencilSketch(srcMat, sketch, colorSketch, sigma_s, sigma_r, shadeFactor, color, numIter);

   if (color)
   {
      colorSketch.copyTo(resMat);
   }
   else
   {
      sketch.copyTo(resMat);
   }
}

/*
 * According to the article's Tone Mapping description
 * in the part on Real-Time Applications.
 */
void TMOGastal11::toneMappingFilter(Mat srcMat,
                                    Mat &resMat,
                                    uint8_t numIter)
{
   int height, width;
   double sigma_s, sigma_r;
   double min, max;
   Point minLoc, maxLoc;
   Mat L, j1, j2, j3;

   height = srcMat.rows;
   width = srcMat.cols; 
   
   L = srcMat;

   sigma_s = 20.0;
   sigma_r = 0.33;

   // J1
   filterOperation(L, j1, sigma_s, sigma_r, RF, numIter);

   sigma_s = 50.0;
   sigma_r = 0.67;

   // J2
   filterOperation(j1, j2, sigma_s, sigma_r, RF, numIter);

   sigma_s = 100.0;
   sigma_r = 1.34;

   // J3
   filterOperation(j2, j3, sigma_s, sigma_r, RF, numIter);
   
   minMaxLoc(j3, &min, &max, &minLoc, &maxLoc);

   Mat B = (j3 - min) / (max - min);
   Scalar tmpMean = mean(mean(mean(B)));
   double mean = tmpMean.val[0];

   Mat d0 = L - j1;
   Mat d1 = j1 - j2;
   Mat d2 = j2 - j3;

   Mat lc = 0.12 + mean + 0.9 * (B - mean) +
      0.3 * d0 + 0.2 * d1 + 0.2 * d2;
      
   resMat = lc;
}

/*
 * Real-Time Application Tone mapping.
 * This method provides tone mapping
 * of the input HDR image with the use
 * of the implemented RF filter.
 */
void TMOGastal11::toneMapping(Mat srcMat,
                              Mat &resMat,
                              uint8_t numIter)
{
   const double eps = pow(2,-52);

   vector<Mat> sepChannelSource;
   split(srcMat, sepChannelSource);

   // Luminance channel.
   Mat lum = 0.299 * sepChannelSource[2] +
             0.587 * sepChannelSource[1] +
             0.114 * sepChannelSource[0];

   merge(sepChannelSource, srcMat);     

   // Using of log-luminance.
   Mat logLum;
   log(lum + eps, logLum);     

   // Filter luminance channel.
   Mat res;
   toneMappingFilter(logLum, res, numIter);

   // Fill the all result math channels
   // with the result luminance channel.
   vector<Mat> sepChannelResult;
   split(resMat, sepChannelResult);          

   sepChannelResult[0] = res;
   sepChannelResult[1] = res;
   sepChannelResult[2] = res;

   merge(sepChannelResult, resMat);

   Mat divider(srcMat.rows, srcMat.cols, CV_64FC3);

   vector<Mat> sepChannelDivider;
   split(divider, sepChannelDivider);

   lum += eps;

   sepChannelDivider[0] = lum;
   sepChannelDivider[1] = lum;
   sepChannelDivider[2] = lum;

   merge(sepChannelDivider, divider);

   srcMat /= divider;
   resMat = resMat.mul(srcMat);

   // Correct Gamma
   pow(resMat, 1.0/1.8, resMat);
}

/*
 * Main method for Domain Transform for Edge-Aware Image and Video Processing operator
 */
int TMOGastal11::Transform()
{
   int height, width, channelsCnt;
   uint8_t pNumIter, pFilterType, pFilterAppl;
   double pSigmaS, pSigmaR, pShadeFactor;
   double* pDestData;
   bool pPencilColor;

   /***********************/
   /* Get user parameters */
   /***********************/

   pNumIter = numIter.GetInt();
   pFilterType = filterType.GetInt();
   pFilterAppl = filterAppl.GetInt();
   pSigmaS = sigma_s.GetDouble();
   pSigmaR = sigma_r.GetDouble();
   pShadeFactor = shadeFactor.GetDouble();
   pPencilColor = pencilColor.GetBool();

   Mat srcMat = TMOImage2Mat(pSrc);

   Mat result(srcMat.rows, srcMat.cols, CV_64FC3);        

   switch (pFilterAppl) {
      case EDGEAWARESMOOTH:
         edgeAwareSmoothing(srcMat, result, pSigmaS, pSigmaR, pFilterType, pNumIter);
         break;
      case DETAILENHC:
         detailEnhancement(srcMat, result, pSigmaS, pSigmaR, pFilterType, pNumIter);
         break;
      case STYLIZATION:
         stylization(srcMat, result, pSigmaS, pSigmaR, pFilterType, pNumIter);
         break;
      case PENCILSKETCH:
         pencilDrawing(srcMat, result, pSigmaS, pSigmaR, pShadeFactor, pPencilColor, pNumIter);
         break;
      case TONEMAPPING:
         toneMapping(srcMat, result, pNumIter);
         break;
      default:
         break;
   }

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
