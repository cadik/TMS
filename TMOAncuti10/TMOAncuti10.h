#include "TMO.h"
#include <cmath>
#include <iostream>
#include <opencv2/opencv.hpp>

class TMOAncuti10 : public TMO
{
public:
	TMOAncuti10();
	virtual ~TMOAncuti10();
	virtual int Transform();
	void convertRGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z);
	void convertXYZtoCIELAB(double X, double Y, double Z, double &L, double &a, double &b);
	void convertCIELABtoCIELCh(double L, double a, double b, double &C, double &h);
	void convertRGBtoHSL(double R, double G, double B, double &H, double &S, double &L);
	double computeL_HK(double L, double C, double H);
	double calculateAPV(double* channel, int width, int height);
	void applySeparableBinomialKernel(double* input, double* output, int width, int height);
	void computeSaliencyMap(double* channel, double* saliencyMap, int width, int height);
	void computeExposednessMap(double* channel, double* exposednessMap, int width, int height);
	void computeChromaticMap(double* channel, double* chromaticMap, double* saturation, int width, int height);
	void computeLaplacianPyramid(const cv::Mat& input, std::vector<cv::Mat>& pyramid, int levels);
	void computeGaussianPyramid(const cv::Mat& input, std::vector<cv::Mat>& pyramid, int levels);
	void normalizeWeightMaps(std::vector<cv::Mat>& weightMaps);
	void fusePyramids(const std::vector<std::vector<cv::Mat>>& laplacianPyramids, const std::vector<std::vector<cv::Mat>>& gaussianPyramids, std::vector<cv::Mat>& fusedPyramid);
	cv::Mat reconstructFromPyramid(const std::vector<cv::Mat>& pyramid);
protected:
	TMODouble dParameter;
};
