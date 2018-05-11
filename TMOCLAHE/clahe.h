#include <opencv2/opencv.hpp>
#include <stdlib.h>     /* srand, rand */
#include <cmath>

#include <cfloat>

#include <stdio.h>

cv::Mat histogramEqualization(	cv::Mat matrix, 
								int height, 
								int width, 
								int gridRegions, 
								double cl);
															