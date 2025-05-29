#include "TMO.h"
#include <cmath>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Sparse>

class TMOHsin11 : public TMO
{
public:
	TMOHsin11();
	virtual ~TMOHsin11();
	virtual int Transform();
	cv::Mat globalMapping(cv::Mat& input);
	cv::Mat localLightness(cv::Mat& grayscale, cv::Mat& colorInput);
	cv::Mat conjugateGrad(cv::Mat& input);

protected:
	TMODouble dParameter;
};
