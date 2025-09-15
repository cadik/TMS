/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*          						   GcsDecolor2s                                *
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

class TMOLiu15 : public TMO
{
public:
	TMOLiu15();
	virtual ~TMOLiu15();
	virtual int Transform();
	Mat wei();

protected:
	TMODouble Lpp;
};
