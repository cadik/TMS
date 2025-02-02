/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   A Tone Mapping Operator Based on Neural and                *
 *                    Psychophysical Models of Visual Perception                *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMO.h"

class TMOCyriac15 : public TMO
{
public:
	TMOCyriac15();
	virtual ~TMOCyriac15();
	cv::Mat convertToMat();
	cv::Mat normalizeAndLuminance(std::vector<double> *RGBmax);
	void clip(cv::Mat *luminanceMat);
	void addToCumulativeHistogram(std::vector<double> *cumulativeHistogram, double value);
	double findGammaEnc(cv::Mat *luminanceMat);
	cv::Mat createGaussianKernel(int size, double sigma);
	cv::Mat computeLocalMean(cv::Mat &img);
	double sgn(double x);
	cv::Mat computeGradient(cv::Mat &I, cv::Mat &I0, cv::Mat &localMean);
	virtual int Transform();

protected:
	int histogramBins = std::pow(2, 16);
	double alpha = 1.0;
	double beta = 1.0;
	double gamma = 1.0;
	double deltaT = 0.15;
	double convergenceThreshold = 0.005;
	double sigma_w = 200.0;
	TMODouble dParameter;
};
