#include <opencv2/opencv.hpp>
#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <fstream>
#include "glog/logging.h"
#include "pallas/simulated_annealing.h"

#include <cfloat>

#include <stdio.h>
#include <time.h>

extern cv::Mat getFirstTerm(   cv::Mat image, 
                                            const double* x, 
                                            int rows, 
                                            int cols);

extern cv::Mat getSecondTerm (  const double* x, 
                                int rows, 
                                int cols, 
                                float ypsilon);    

extern cv::Mat optimizeForSigma(int height, 
                                int width, 
                                cv::Mat originalImage, 
                                cv::Mat smoothenedImage);     

extern std::vector<cv::Mat> optimizeForGettingSAndTparameters(  int height, 
                                                                int width, 
                                                                cv::Mat detailImage, 
                                                                cv::Mat weight1, 
                                                                cv::Mat weight2, 
                                                                std::vector<cv::Mat> baseChannels, 
                                                                std::vector<cv::Mat> detailChannels);                               