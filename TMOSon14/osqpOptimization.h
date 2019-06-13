/* -------------------------------------------------------------------------- *
 * osqpOptimization.h: solve quadratic programming problem with OSQP library  *
 * Authors: Tomas Hudziec (2019)                                              *
 *          - based on Pavel Sedlar's code with qpOASES                       *
 *          - replaced QP library qpOASES with OSQP                           *
 * -------------------------------------------------------------------------- */
#include "osqp.h"
#include <cmath>
#include <opencv2/opencv.hpp>

extern std::vector<cv::Mat> optimizationWithOsqp(cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels);
