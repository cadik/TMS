/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*          Log-Euclidean Metrics for Contrast Preserving Decolorization        *
*                                                                              *
*             Author: Peter Zdravecký [xzdrav00 AT stud.fit.vutbr.cz]          *
*                                    Brno 2024                                 *
*                                                                              *
*******************************************************************************/
#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <algorithm>
#include <random>
#include <cmath>

using namespace cv;

class TMOLiu17 : public TMO
{
public:
	TMOLiu17();
	virtual ~TMOLiu17();
	virtual int Transform();
	Mat wei();

protected:
	TMODouble gamma;
};
