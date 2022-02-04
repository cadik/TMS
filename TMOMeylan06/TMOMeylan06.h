/*
author: Jan Kohút (xkohut08)
date: 9. 4. 2019
*/
#include "TMO.h"
#include <math.h>
#include <algorithm>
#include <string>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

/*
High Dynamic Range Image Rendering with a Retinex-Based Adaptive Filter

This is implementation of MATLAB code that can be found at: https://infoscience.epfl.ch/record/50199
*/

class TMOMeylan06 : public TMO
{
public:
	TMOMeylan06();
	virtual ~TMOMeylan06();
	virtual int Transform();

protected:
	TMODouble kernelRadiusParameter;
	TMODouble sigmaOrigParameter;
	TMODouble sigmaEdgeParameter;
	TMODouble saturationParameter;

	int numberOfPixels;
	int numberOfPixelsRGB;
	cv::PCA pca;
	double kernelRadius;
	double sigmaOrig;
	double sigmaEdge;

	/* PCA ANALYSIS */
	/****************************************************************************/
	// Implementation of "myRGB2PCA.m".
	cv::Mat RGBToPCA(cv::Mat &data);
	// Implementation of "myPCA2RGB.m".
	cv::Mat PCAToRGB(cv::Mat &PCAProjection);
	cv::Mat GetLuminance(cv::Mat &PCAProjection);
	void ScaleSaturation(cv::Mat &data, double saturationEnhancement);
	/****************************************************************************/

	/* GLOBAL TONE MAPPING */
	/****************************************************************************/
	// Implementation of "globalOpPow.m".
	void GlobalToneMap(cv::Mat &data, std::string type);
	void ScaleRGB(cv::Mat &data, double RScale, double GScale, double BScale);
	double ComputeAL(cv::Mat &data, double scale);
	/****************************************************************************/

	/* EDGE DETECTION */
	/****************************************************************************/
	// Implementation of "edgeImage.m".
	cv::Mat GetEdges(cv::Mat &luminance, double upperThresholdRatio);
	cv::Mat DilatateEdges(cv::Mat &edges);
	cv::Mat ResizeLuminance(cv::Mat &luminance, int maskMaxSize);
	/****************************************************************************/

	/* LOCAL TONE MAPPING */
	/****************************************************************************/
	// Refactorization of "runIterRetinex.c".
	cv::Mat GetMask(cv::Mat &luminance, cv::Mat &edges);
	double GetMaskVal(cv::Mat &luminance, cv::Mat &edges, cv::Mat &crossCounter, int x, int y, int counter);
	bool IsAnEdge(cv::Mat &edges, cv::Mat &crossCounter, int x, int y, int counter);
	double GaussDist(double d, double s);
	/****************************************************************************/

	/* BETA FACTOR */
	/****************************************************************************/
	// Implementation of "sigmFac(max(min(H_Ylog01,1),0)*255,10)".
	cv::Mat GetBetaFactor(cv::Mat &luminance, double c, double a);
	/****************************************************************************/

	/* HELPER METHODS */
	/****************************************************************************/
	// Implementation of "H_lin01/max(H_lin01(:))".
	void Normalize(cv::Mat &data);

	// Implementation of "histoClip.m".
	void HistoClip(cv::Mat &data, int numberOfBuckets, double minThreshold, double maxThreshold);

	// Implementation of "log(max(1,rgb1*100))/log(100)".
	void LogMaxScale(cv::Mat &data, double max, double scale);

	void Normalize(cv::Mat &data, double lowerBound, double upperBound);
	cv::Mat ElementWiseMul(cv::Mat &first, cv::Mat &second);
	cv::Mat ElementWiseSub(cv::Mat &first, cv::Mat &second);
	void Pow(cv::Mat &data, double exponent);
	void Max(cv::Mat &data, double max);
	void Min(cv::Mat &data, double min);
	double GetMax(double *data, int dataLength);
	double GetMin(double *data, int dataLength);
	void SaveImg(std::string name, cv::Mat &data);
	/****************************************************************************/
};
