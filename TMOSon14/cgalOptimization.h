// outputDetailMaximalization: construct a nonnegative quadratic program from data
// the QP below is the first nonnegative quadratic program outputDetailMaximalization 
// in the user manual

#include <iostream>
#include <regex>

#include <sstream>
#include <string>
#include <iostream>
#include <cassert>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include<gmp.h>
// choose exact integral type

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif
// #include <cstdio>
#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <fstream>

#include <opencv2/opencv.hpp>
#include <cmath>
// program and solution types
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;
using namespace cv;
using namespace std;

extern void tryIt();

extern std::vector<cv::Mat> optimizeForGettingSAndTparametersWithCgal(  int height, 
                                                                        int width, 
                                                                        cv::Mat detailImage, 
                                                                        cv::Mat weight1, 
                                                                        cv::Mat weight2, 
                                                                        std::vector<cv::Mat> baseChannels, 
                                                                        std::vector<cv::Mat> detailChannels);
