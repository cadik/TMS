/*
* Author of the code: Mária Nováková (xnovak2w@stud.fit.vutbr.cz)
* Author of the paper: Zhang L., Wan Y.
* Date of implementation: 28.04.2024
* Description: Implementation of color-to-gray image conversion using salient colors and radial basis functions.
*/

#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

using namespace std;
using namespace cv;

class TMOZhangWan24 : public TMO
{
public:
	TMOZhangWan24();
	virtual ~TMOZhangWan24();
	virtual int Transform();

	void quantizeColors(Mat imageLab, Mat image, Mat grayImage, int &k, double theta_0, double theta_1);
	void ordering();
	Mat createGrayScale(Mat imageLab, Mat image);
	double getGreyValue(Vec3d img_color, vector<double> a);


	vector<Vec3d> centers; 
    vector<vector<pair<Point, Vec3d>>> clusters;
    vector<double> grayvalues;

protected:
	TMODouble max_k;
	TMODouble sigma;
};
