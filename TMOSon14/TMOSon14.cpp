/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio	                                        *
*                                                                                   *
*                       Authors: Pavel Sedler (initial code),                       *
*                       Tomas Hudziec (added debug code,                            *
*                       replaced QP library qpOASES with OSQP, code cleanup)        *
*                       Brno 2018                                                   *
*                                                                                   *
*                       Implementation of "Art-Photographic Detail Enhancement"     *
*                       Minjung Son, Yunjin Lee, Henry Kang, Seungyong Lee          *
*                       Computer Graphics Forum 2014                                *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOSon14.cpp
 * @brief Implementation of "Art-Photographic Detail Enhancement" Minjung Son, Yunjin Lee, Henry Kang, Seungyong Lee, Computer Graphics Forum 2014
 * @author Pavel Sedler
 * @author Tomas Hudziec
 * @class TMOSon14.cpp
 */ 

#include "TMOSon14.h"
// #include <fftw3.h>
#include <iostream>
#include <fstream>
#include <cfloat>

// #include <stdlib.h>     /* srand, rand */
// #include <time.h>

using namespace std;
using namespace cv;
using namespace Eigen;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */

TMOSon14::TMOSon14()
{
	SetName(L"Son14");
	SetDescription(L"Art-Photography detail enhancement");

	debugFlag.SetName(L"debug");
	debugFlag.SetDescription(L"enable to store intermediate results");
	debugFlag.SetDefault(false);
	debugFlag = false;
	this->Register(debugFlag);

	/**
	 * Mu - Parameter
	 **/
	mu.SetName(L"Mu");
	mu.SetDescription(L"Rate Mu for detail maximization");
	mu.SetDefault(0.5);
	mu=0.5;
	mu.SetRange(0.0,1.0);
	this->Register(mu);
	/**
	 * Iteration of optimizing sigma control - Parameter
	 **/
	optim1Iteration.SetName(L"Sigma");
	optim1Iteration.SetDescription(L"Number of iterations to repeat for getting sigma map");
	optim1Iteration.SetDefault(50);
	optim1Iteration=10;
	optim1Iteration.SetRange(1, 1000);
	this->Register(optim1Iteration);
}

TMOSon14::~TMOSon14()
{
}

// https://stackoverflow.com/a/6417908
std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

/* --------------------------------------------------------------------------- *
 * This overloaded function
 * -------------------------------------------------------- */
int TMOSon14::Transform()
{
	ofstream myfile;
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();

	/*
	 * Base matrix
	 **/
	cv::Mat r;
	cv::Mat g;
	cv::Mat b;

	r = cv::Mat::zeros(height, width, CV_32F);
	g = cv::Mat::zeros(height, width, CV_32F);
	b = cv::Mat::zeros(height, width, CV_32F);

	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_RGB);								// This is format of Y as luminance
	// pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array
											// of three doubles representing

	// get filename
	const char *dstFilename = pDst->GetFilename();
	std::string fnWoExt = remove_extension(std::string(dstFilename));

	/*
	 * Fill base matrix
	 * */
	 int j;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			r.at<float>(j,i) = float(*pSourceData++);
			g.at<float>(j,i) = float(*pSourceData++);  //getting separate RGB channels
			b.at<float>(j,i) = float(*pSourceData++);
		}
	}

	pSrc->ProgressBar(j, pSrc->GetHeight());
    /*
	 * For L0 smoothing
	 * */

	std::vector<cv::Mat> array_to_merge;

    array_to_merge.push_back(b);
    array_to_merge.push_back(g);
    array_to_merge.push_back(r);

    cv::Mat originalImage;

    cv::merge(array_to_merge, originalImage);

	/*
	 * Base Decomposition
	 **/

	/*
	 * Phase 1 - Original L0 Smoothing
	 **/
	std::cout << "Base Phase1" << std::endl;
	cv::Mat basePhase1 = minimizeL0Gradient1(originalImage);

	/*
	 * Phase 2 - L0 smooting with adaptive lambda matrix
	 **/
	std::cout << "Base Phase2" << std::endl;
	cv::Mat gradientFrom1stSmoothing = getGradientMagnitude(basePhase1);
	basePhase1.release();
	cv::Mat adaptiveLambdaMatrix1 = getAdaptiveLambdaMatrix(gradientFrom1stSmoothing, height, width);
	// cv::Mat::zeros(height, width, CV_32F); //
	/*for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			adaptiveLambdaMatrix1.at<float>(j,i) = 0.00001;
		}
	}*/

	cv::Mat basePhase2 = minimizeL0GradientSecondFaze(originalImage, adaptiveLambdaMatrix1, height, width);
	adaptiveLambdaMatrix1.release();
	gradientFrom1stSmoothing.release();
	/*
     * split basePhase2;
     **/

	cv::Mat basePhase2Chan[3];
	cv::split(basePhase2, basePhase2Chan);

    (basePhase2Chan[0]).convertTo(basePhase2Chan[0], CV_32F);
    (basePhase2Chan[1]).convertTo(basePhase2Chan[1], CV_32F);
    (basePhase2Chan[2]).convertTo(basePhase2Chan[2], CV_32F);

    /*
     * Phase 3 -- getting final base layer
     **/

    /*
	cv::Mat sumOfCostsBase = basePhase1Chan[2] + basePhase2Chan[1] + basePhase2Chan[0];

	cv::Mat sumOfCostsOriginal = r + g + b;
	std::cout << "Base Phase3" << std::endl;
	cv::Mat sigmaMap = stochasticOptimizationForGetSigma(sumOfCostsBase/256.0, sumOfCostsOriginal, height, width, 50000);

	cv::Mat basePhase3R = myOwn2DFilter(r, sigmaMap, height, width);
	cv::Mat basePhase3G = myOwn2DFilter(g, sigmaMap, height, width);
	cv::Mat basePhase3B = myOwn2DFilter(b, sigmaMap, height, width);
	*/
	std::cout << "Base phase -- COMPLETED" << std::endl;

	cv::Mat &basePhase2R = basePhase2Chan[2];
	cv::Mat &basePhase2G = basePhase2Chan[1];
	cv::Mat &basePhase2B = basePhase2Chan[0];
	cv::Mat detailLayerR = r - basePhase2R;
	cv::Mat detailLayerG = g - basePhase2G;
	cv::Mat detailLayerB = b - basePhase2B;

	cv::Mat sumOfDetail = detailLayerR + detailLayerG + detailLayerB;
	cv::Mat sumOfBase = basePhase2R + basePhase2G + basePhase2B;

    std::vector<cv::Mat> base_channels;

    base_channels.push_back(basePhase2R);
    base_channels.push_back(basePhase2G);
    base_channels.push_back(basePhase2B);

    cv::Mat baseImage;

    cv::merge(base_channels, baseImage);
	/*
		Getting weights
	*/

    cv::Mat gradientOfBaseLayer = getGradientMagnitude(baseImage);

	// FIXME call only once and then adjust with r1 and r2
	cv::Mat r1Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 200);
	cv::Mat r2Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 500);

	if(debugFlag) {
		cv::Mat showWeight;
		cv::normalize(r1Layer, showWeight, 0, 255, cv::NORM_MINMAX, r1Layer.type());
		cv::imwrite(fnWoExt + "_Weight.png", showWeight);
	}

	std::vector<cv::Mat> detail;
	detail.push_back((detailLayerR.clone()));
	detail.push_back((detailLayerG.clone()));
	detail.push_back((detailLayerB.clone()));

	cv::normalize(sumOfDetail, sumOfDetail, 0, 1, cv::NORM_MINMAX, sumOfDetail.type());
	std::vector<cv::Mat> ST = optimizationWithOsqp(sumOfDetail, r1Layer, r2Layer, base_channels, detail);

	// some error occured
	if(ST.empty()) {
		std::cerr << "some error occured during QP computation, its output is empty" << '\n';
		return 1;
	}

	detail.clear();
	base_channels.clear();
	sumOfDetail.release();
	sumOfBase.release();
	r1Layer.release();
	r2Layer.release();
	std::cout << "Detail maximalization -- COMPLETED" << std::endl;

	cv::Mat &s = ST[0];
	cv::Mat &t = ST[1];

	if(debugFlag) {
		cv::Mat showScale, showshifT;
		cv::normalize(s, showScale, 0, 255, cv::NORM_MINMAX, s.type());
		cv::imwrite(fnWoExt + "_Scale.png", showScale);
		cv::normalize(t, showshifT, 0, 255, cv::NORM_MINMAX, t.type());
		cv::imwrite(fnWoExt + "_shifT.png", showshifT);
	}

	/*
	 * Function for control details enhancement of picture
	 **/
	// detailMaximizedLayer = (mu*s + (1-mu))*D + B + mu*t
	cv::Mat detailMaximizedLayerR = getDetailControl(basePhase2R, detailLayerR, s, t, mu, height, width);
	cv::Mat detailMaximizedLayerG = getDetailControl(basePhase2G, detailLayerG, s, t, mu, height, width);
	cv::Mat detailMaximizedLayerB = getDetailControl(basePhase2B, detailLayerB, s, t, mu, height, width);

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			*pDestinationData++ = ((detailMaximizedLayerR).at<float>(j,i));
			*pDestinationData++ = ((detailMaximizedLayerG).at<float>(j,i));
			*pDestinationData++ = ((detailMaximizedLayerB).at<float>(j,i));
		}
	}
	pDst->Convert(TMO_RGB);
	return 0;
}
