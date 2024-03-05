/***********************************************************************************
* Author: Bc. Tomas Krivanek (xkriva29)   xkriva29@stud.fit.vutbr.cz               *
************************************************************************************/
#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cmath>
#include <fstream>

class TMONafchi17 : public TMO
{
public:
	TMONafchi17();
	virtual ~TMONafchi17();
	virtual int Transform();

protected:
	TMODouble r;
	TMOBool imageType;

private:
	cv::Mat calculateMeanImage(cv::Mat &in01, int width, int height);
	cv::Mat calculateStdDevImage(cv::Mat &in01, cv::Mat &meanImage, int width, int height);
	std::tuple<double,double,double> calculatePearsonCoeff(cv::Mat &in01, cv::Mat &contrastMap, int width, int height);
	std::tuple<double,double,double> calculateLambda(cv::Mat &in01, std::tuple<double,double,double> rhos, int width, int height);
};
