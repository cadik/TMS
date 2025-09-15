/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*             Author: Tomas Krivanek [xkriva29 AT stud.fit.vutbr.cz]           *
*                                    Brno 2024                                 *
*                                                                              *
*******************************************************************************/

#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

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

	bool haveImageGreaterValuesThen(cv::Mat &in01, double value);
	cv::Mat convertToDouble(cv::Mat &in01);
};
