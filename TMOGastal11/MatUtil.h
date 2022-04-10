#ifndef TMS_MATUTIL_H
#define TMS_MATUTIL_H

#include "opencv2/opencv.hpp"

#define PADARR_REPLICATE 0
#define PADARR_CONSTANT 1

using namespace cv;

class MatUtil
{
public:
    void diff(Mat &srcMat, Mat &dstMat, uint8_t dim);

    Mat sub2indOneDim(Size indSize, Size resSize, Mat rows, Mat cols);

    int findFirstGtr(Mat srcMat, double cmpVal, int startIdx = 0);

    Mat padArray (Mat srcMat, int rowPad, int colPad, uint8_t borderType);
};


#endif //TMS_MATUTIL_H
