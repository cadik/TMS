#ifndef TMS_COLORSPACECONVERTER_H
#define TMS_COLORSPACECONVERTER_H

#include "opencv2/opencv.hpp"

#define EPS (216.0 / 24389.0)
#define KAPPA (24389.0 / 27.0)

// XYZ reference white
#define XYZ_WHITE_X 95.047
#define XYZ_WHITE_Y 100.0
#define XYZ_WHITE_Z 108.8830

#define CHANNELSCNT 3

using namespace cv;

class ColorspaceConverter
{
public:
    void bgr2Lab64FType(Mat srcMat, Mat &outMat);

    void lab2Bgr64FType(Mat srcMat, Mat &outMat);

    void bgr2YCrCb64FType(Mat srcMat, Mat &outMat);

    void yCrCb2Bgr64FType(Mat srcMat, Mat &outMat);

private:
   double XYZ2RGB[3][3] =
   {
           {3.2404542, -1.5371385, -0.4985314},
           {-0.9692660, 1.8760108, 0.0415560},
           {0.0556434, -0.2040259, 1.0572252}
   };

   double RGB2XYZ[3][3] =
   {
           {0.4124564, 0.3575761, 0.1804375},
           {0.2126729, 0.7151522, 0.0721750},
           {0.0193339, 0.1191920, 0.9503041}
   };

    double srgbCompanding(double color);

    double inverseSrgbCompanding(double color);

    void bgr2Xyz64FType(Mat srcMat, Mat &outMat);

    void xyz2Lab(Mat srcMat, Mat &outMat);

    void lab2Xyz(Mat srcMat, Mat &outMat);

    void xyz2Bgr(Mat srcMat, Mat &outMat);
};


#endif //TMS_COLORSPACECONVERTER_H
