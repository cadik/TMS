/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*             Fast Local Laplacian Filters: Theory and Applications (2014)     *
* by Mathieu Aubry, Sylvain Paris, Samuel W. Hasinoff, Jan Kautz, Fredo Durand *
*                         ACM Transactions on Graphics                         *
*                                                                              *
*             Author: Tomas Hudziec [xhudzi01 AT stud.fit.vutbr.cz]            *
*         Term project for Computational Photography course - 2018             *
*       Part of master thesis (HDR support, code reorganization) - 2019        *
*                                                                              *
*******************************************************************************/

#ifndef FAST_LOCAL_LAPL_FILT
#define FAST_LOCAL_LAPL_FILT

#include "TMO.h"
#include <vector>
#include "opencv2/opencv.hpp"

cv::Mat FastLocalLaplFilt(cv::Mat I, double sigma, double fact, int N, TMOImage *pSrc);

static std::vector<double> linspace(double min, double max, int n);

#endif /* end of include guard: FAST_LOCAL_LAPL_FILT */
