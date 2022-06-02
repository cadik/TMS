#ifndef TMS_MATUTIL_H
#define TMS_MATUTIL_H

#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

#define CHANNELSCNT 3
#define PADARR_REPLICATE 0
#define PADARR_CONSTANT 1

class MatUtil
{
public:
    /*
     * Own implementation of MATLAB padarray function using OpenCV copyMakeBorder().
     *
     * Original MATLAB function: https://www.mathworks.com/help/images/ref/padarray.html
     */
    Mat padArray (Mat &srcMat, int rowPad, int colPad, uint8_t borderType);

    /*
     * Own implementation of MATLAB colfilt function.
     *
     * colfilt function settings:
     * - Sliding Neighborhood Operations
     * - @std function (Standard deviation)
     *
     * Original MATLAB function: https://www.mathworks.com/help/images/ref/colfilt.html
     */
    Mat colFiltSlidingStd(Mat &srcMat, int blockRowsCnt, int blockColsCnt);

    /*
     * Set value by HSV color model.
     *
     * This function does not implement whole algorithm, but only Value part.
     * Article: Smith, A. R. “Color Gamut Transform Pairs”. SIGGRAPH 78 Conference Proceedings. 1978, pp. 12–19.
     */
    void setValueHSV(Mat &srcMat, Mat &lum);

    /*
     * Check luminance value.
     *
     * If the value is lower than 1e-3, to actual luminance value is add 1e-3.
     */
    void checkLumo(Mat &lumo);    

    /*
     * Cuts off unreasonable peaks and Valleys. 
     */
    Mat cutPeakValley(Mat &rgb, double lowEndCutPercentage, double highEndCutPercentage);

    /*
     * Normalize result image.
     */
    void normImage(Mat &img);

private:
    /*
     * Calculates the assigned percentile of matrix data.
     */
    double percentile(Mat &data, double percentile);

};


#endif //TMS_MATUTIL_H
