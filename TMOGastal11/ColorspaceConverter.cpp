/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2022                                              *
*                                                                              *
*                       Implementation of the ColorspaceConverter class        *
*                                                                              *
*******************************************************************************/

#include "ColorspaceConverter.h"

/*
 * sRGB Companding.
 */
double ColorspaceConverter::srgbCompanding(double color)
{
   return (color <= 0.0031308) ? (12.92 * color) : (1.055 * pow(color, 1.0 / 2.4) - 0.055);
}

/*
 * Inverse sRGB Companding.
 */
double ColorspaceConverter::inverseSrgbCompanding(double color)
{
   return (color <= 0.04045) ? (color / 12.92) : pow((color + 0.055) / 1.055, 2.4);
}

/*
 * Mat conversion function from BGR colorspace to XYZ.
 */
void ColorspaceConverter::bgr2Xyz64FType(Mat srcMat, Mat &outMat)
{
   int height, width, channelsCnt;

   height = srcMat.rows;
   width  = srcMat.cols;
   channelsCnt = srcMat.channels();

   assert(channelsCnt == CHANNELSCNT);
   assert(srcMat.type() == CV_64FC3);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         outMat.at<double>(y, x*channelsCnt) = 0.0;
         outMat.at<double>(y, x*channelsCnt+1) = 0.0;
         outMat.at<double>(y, x*channelsCnt+2) = 0.0;

         for (int c = 0; c < channelsCnt; c++)
         {
            /* opencv mat is in order BGR and convert matrix is in RGB */
            double b = inverseSrgbCompanding(srcMat.at<double>(y, x*channelsCnt));
            double g = inverseSrgbCompanding(srcMat.at<double>(y, x*channelsCnt+1));
            double r = inverseSrgbCompanding(srcMat.at<double>(y, x*channelsCnt+2));

            outMat.at<double>(y, x*channelsCnt+c) += RGB2XYZ[c][2] * b;
            outMat.at<double>(y, x*channelsCnt+c) += RGB2XYZ[c][1] * g;
            outMat.at<double>(y, x*channelsCnt+c) += RGB2XYZ[c][0] * r;
         }
      }
   }

   multiply(outMat, 100.0, outMat);
}

/*
 * Mat conversion function from XYZ colorspace to LAB.
 */
void ColorspaceConverter::xyz2Lab(Mat srcMat, Mat &outMat)
{
   int height, width, channelsCnt;

   height = srcMat.rows;
   width = srcMat.cols;
   channelsCnt = srcMat.channels();

   assert(channelsCnt == CHANNELSCNT);
   assert(srcMat.type() == CV_64FC3);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         double _x = srcMat.at<double>(y, x*channelsCnt);
         double _y = srcMat.at<double>(y, x*channelsCnt+1);
         double _z = srcMat.at<double>(y, x*channelsCnt+2);

         double xr = _x / XYZ_WHITE_X;
         double yr = _y / XYZ_WHITE_Y;
         double zr = _z / XYZ_WHITE_Z;

         double fx = (xr > EPS) ? pow(xr, 1.0 / 3.0) : ((KAPPA * xr + 16.0) / 116.0);
         double fy = (yr > EPS) ? pow(yr, 1.0 / 3.0) : ((KAPPA * yr + 16.0) / 116.0);
         double fz = (zr > EPS) ? pow(zr, 1.0 / 3.0) : ((KAPPA * zr + 16.0) / 116.0);

         double L = (116.0 * fy - 16.0) > 0.0 ? 116.0 * fy - 16.0 : 0.0;

         double a = 500 * (fx - fy);
         double b = 200 * (fy - fz);

         outMat.at<double>(y, x*channelsCnt) = L;
         outMat.at<double>(y, x*channelsCnt+1) = a;
         outMat.at<double>(y, x*channelsCnt+2) = b;
      }
   }
}

/*
 * Mat conversion function from LAB colorspace to XYZ.
 */
void ColorspaceConverter::lab2Xyz(Mat srcMat, Mat &outMat)
{
   int height, width, channelsCnt;

   height = srcMat.rows;
   width = srcMat.cols;
   channelsCnt = srcMat.channels();

   assert(channelsCnt == CHANNELSCNT);
   assert(srcMat.type() == CV_64FC3);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         double L = srcMat.at<double>(y, x*channelsCnt);
         double a = srcMat.at<double>(y, x*channelsCnt+1);
         double b = srcMat.at<double>(y, x*channelsCnt+2);

         double fy = (L + 16.0) / 116.0;
         double fx = a / 500.0 + fy;
         double fz = fy - b / 200.0;

         double xr = (pow(fx, 3) > EPS) ? pow(fx, 3) : (116.0 * fx - 16.0) / KAPPA;
         double yr = (L > KAPPA * EPS) ? pow((L + 16.0) / 116.0, 3) : L / KAPPA;
         double zr = (pow(fz, 3) > EPS) ? pow(fz, 3) : (116.0 * fz - 16.0) / KAPPA;

         outMat.at<double>(y, x*channelsCnt) = xr * XYZ_WHITE_X;
         outMat.at<double>(y, x*channelsCnt+1) = yr * XYZ_WHITE_Y;
         outMat.at<double>(y, x*channelsCnt+2) = zr * XYZ_WHITE_Z;
      }
   }
}

/*
 * Mat conversion function from XYZ colorspace to BGR.
 */
void ColorspaceConverter::xyz2Bgr(Mat srcMat, Mat &outMat)
{
   int height, width, channelsCnt;

   height = srcMat.rows;
   width = srcMat.cols;
   channelsCnt = srcMat.channels();

   assert(channelsCnt == CHANNELSCNT);
   assert(srcMat.type() == CV_64FC3);

   double* pix = new double[3];

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         pix[0] = 0.0;
         pix[1] = 0.0;
         pix[2] = 0.0;

         for (int c = 0; c < channelsCnt; c++)
         {
            pix[c] += XYZ2RGB[c][0] * (srcMat.at<double>(y,x*channelsCnt)/100.0);
            pix[c] += XYZ2RGB[c][1] * (srcMat.at<double>(y,x*channelsCnt+1)/100.0);
            pix[c] += XYZ2RGB[c][2] * (srcMat.at<double>(y,x*channelsCnt+2)/100.0);
         }

         outMat.at<double>(y, x*channelsCnt) = srgbCompanding(pix[2]);
         outMat.at<double>(y, x*channelsCnt+1) = srgbCompanding(pix[1]);
         outMat.at<double>(y, x*channelsCnt+2) = srgbCompanding(pix[0]);
      }
   }

   delete [] pix;
}

/*
 * Mat conversion function from BGR colorspace to LAB for CV_64FC3 type.
 */
void ColorspaceConverter::bgr2Lab64FType(Mat srcMat, Mat &outMat)
{
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   Mat resMat(height, width, CV_64FC3);

   bgr2Xyz64FType(srcMat, resMat);
   xyz2Lab(resMat, outMat);
}

/*
 * Mat conversion function from LAB colorspace to BGR for CV_64FC3 type.
 */
void ColorspaceConverter::lab2Bgr64FType(Mat srcMat, Mat &outMat)
{
   int height, width;

   height = srcMat.rows;
   width = srcMat.cols;

   Mat resMat(height, width, CV_64FC3);

   lab2Xyz(srcMat, resMat);
   xyz2Bgr(resMat, outMat);
}

/*
 * Mat conversion function from BGR colorspace to YCrCb for CV_64FC3 type.
 *
 * Based on OpenCV conversion calculation: https://docs.opencv.org/3.4/de/d25/imgproc_color_conversions.html
 */
void ColorspaceConverter::bgr2YCrCb64FType(Mat srcMat, Mat &outMat)
{
   int height, width, channelsCnt;
   const double delta = 0.5;

   height = srcMat.rows;
   width  = srcMat.cols;
   channelsCnt = srcMat.channels();

   assert(channelsCnt == CHANNELSCNT);
   assert(srcMat.type() == CV_64FC3);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         double r = srcMat.at<double>(y, x*channelsCnt+2);
         double g = srcMat.at<double>(y, x*channelsCnt+1);
         double b = srcMat.at<double>(y, x*channelsCnt);
         // BGR
         double Y = 0.299 * r + 0.587 * g + 0.114 * b;
         double Cr = (r - Y) * 0.713 + delta;
         double Cb = (b - Y) * 0.564 + delta;

         outMat.at<double>(y, x*channelsCnt)   = Y;
         outMat.at<double>(y, x*channelsCnt+1) = Cr;
         outMat.at<double>(y, x*channelsCnt+2) = Cb;
      }
   }
}

/*
 * Mat conversion function from YCrCb colorspace to BGR for CV_64FC3 type.
 *
 * Based on OpenCV conversion calculation: https://docs.opencv.org/3.4/de/d25/imgproc_color_conversions.html
 */
void ColorspaceConverter::yCrCb2Bgr64FType(Mat srcMat, Mat &outMat)
{
   int height, width, channelsCnt;
   const double delta = 0.5;

   height = srcMat.rows;
   width = srcMat.cols;
   channelsCnt = srcMat.channels();

   assert(channelsCnt == CHANNELSCNT);
   assert(srcMat.type() == CV_64FC3);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         double Y = srcMat.at<double>(y, x*channelsCnt);
         double Cr = srcMat.at<double>(y, x*channelsCnt+1);
         double Cb = srcMat.at<double>(y, x*channelsCnt+2);

         // BGR
         double r = Y + 1.403*(Cr - delta);
         double g = Y - 0.714 * (Cr - delta) - 0.344 * (Cb - delta);
         double b = Y + 1.773 * (Cb - delta);

         outMat.at<double>(y, x*channelsCnt) = b;
         outMat.at<double>(y, x*channelsCnt+1) = g;
         outMat.at<double>(y, x*channelsCnt+2) = r;
      }
   }
}
