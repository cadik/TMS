#include "osqp.h"
#include <cmath>
#include <opencv2/opencv.hpp>

extern std::vector<cv::Mat> optimizationWithOsqp(cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels);
