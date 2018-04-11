#include <opencv2/opencv.hpp>
#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <fstream>
#include "glog/logging.h"
#include "pallas/simulated_annealing.h"

#include <cfloat>

#include <stdio.h>
#include <time.h>

extern cv::Mat myOwn2DFilterOptimization(   cv::Mat image, 
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