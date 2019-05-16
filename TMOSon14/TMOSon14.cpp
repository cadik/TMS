/* --------------------------------------------------------------------------- *
 * TMOSon14.cpp: implementation of the TMOSon14 class.   *
 * --------------------------------------------------------------------------- */
#include "TMOSon14.h"
// #include <fftw3.h>
//#include <opencv2/core/core.hpp>
// #include "L0minimization.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;
using namespace Eigen;
/* --------------------------------------------------------------------------- *
 * Constructor se<ves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
/*
	pokus pallas
*/
#include <opencv2/opencv.hpp>
#include <stdlib.h>     /* srand, rand */

#include <cfloat>

#include <stdio.h>
#include <time.h>

void showImg(std::string name, const cv::Mat &img)
{
	cv::Mat show_img = img.clone();
	cv::normalize(show_img, show_img, 0, 1, cv::NORM_MINMAX, show_img.type());

	// cv::cvtColor(show_img, show_img, cv::COLOR_RGB2BGR);
	// cv::cvtColor(show_img, show_img, cv::COLOR_Lab2BGR);
	cv::imshow(name, show_img);
	cv::waitKey(0);
	cv::destroyAllWindows();
}

void printMatRange(std::string name, const cv::Mat &mat)
{
	double min, max;
	cv::minMaxLoc(mat, &min, &max);
	std::cout << "name: " << name ;
	std::cout << ", min: " << min ;
	std::cout << ", max: " << max << "\n\n";
}

/*
	pokus pallas
*/
TMOSon14::TMOSon14()
{
	SetName(L"TMOSon14");						// TODO - Insert operator name
	SetDescription(L"Art-Photography detail enhancement");	// TODO - Insert description
	/**
	 * Mu - Parameter
	 **/
	mu.SetName(L"Mu");				// TODO - Insert parameters names
	mu.SetDescription(L"Represents rate Mu for detail maximization");	// TODO - Insert parameter descriptions
	mu.SetDefault(1.0);							// TODO - Add default values
	mu=1.0;
	mu.SetRange(0.0,1.0);				// TODO - Add acceptable range if needed
	this->Register(mu);
	/**
	 * Iteration of optimizing sigma control - Parameter
	 **/
	optim1Iteration.SetName(L"Sigma Iteration Control");				// TODO - Insert parameters names
	optim1Iteration.SetDescription(L"Represents number of iteration to repeat for getting sigma map");	// TODO - Insert parameter descriptions
	optim1Iteration.SetDefault(50);							// TODO - Add default values
	optim1Iteration=10;
	optim1Iteration.SetRange(1, 1000);				// TODO - Add acceptable range if needed
	this->Register(optim1Iteration);	
}

TMOSon14::~TMOSon14()
{
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
	/*
	 * Fill base matrix
	 * */
	 int j;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			r.at<float>(j,i) = *pSourceData++; 
			g.at<float>(j,i) = *pSourceData++;  //getting separate RGB channels
			b.at<float>(j,i) = *pSourceData++;
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
	
    /*cv::Mat sumOfCostsBase = basePhase1Chan[2] + basePhase2Chan[1] + basePhase2Chan[0];
	
	cv::Mat sumOfCostsOriginal = r + g + b;
	std::cout << "Base Phase3" << std::endl;*/
	//cv::Mat sigmaMap = optimizeForSigma(height, width, sumOfCostsOriginal/255.0, sumOfCostsBase/255.0, optim1Iteration);
	/////cv::Mat sigmaMap = stochasticOptimizationForGetSigma(sumOfCostsBase/256.0, sumOfCostsOriginal, height, width, 50000);
	
	/*cv::Mat basePhase3R = myOwn2DFilter(r, sigmaMap, height, width);
	cv::Mat basePhase3G = myOwn2DFilter(g, sigmaMap, height, width);
	cv::Mat basePhase3B = myOwn2DFilter(b, sigmaMap, height, width);*/
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
	
	// printMatRange("detail", sumOfDetail);
	// printMatRange("w1", r1Layer);
	// printMatRange("w2", r2Layer);
	
	std::vector<cv::Mat> detail;
	detail.push_back((detailLayerR.clone()));
	detail.push_back((detailLayerG.clone()));
	detail.push_back((detailLayerB.clone()));
	
	//std::vector<cv::Mat> SdsaT = detailMaximalization(sumOfBase, sumOfDetail, r1Layer, r2Layer, height, width, 50, detail); 
	
	
	cv::normalize(sumOfDetail, sumOfDetail, 0, 1, cv::NORM_MINMAX, sumOfDetail.type());
	// sumOfDetail /= 255.0;
	printMatRange("detail normalized", sumOfDetail);
	printMatRange("r1Layer", r1Layer);
	printMatRange("r2Layer", r2Layer);
	std::vector<cv::Mat> ST = optimizationWithOsqp(sumOfDetail, r1Layer, r2Layer, base_channels, detail);

	// some error occured
	if(ST.empty())
		return 1;

	cv::Mat showW1, showW2;
	cv::Size size(256, 256);
	cv::resize(r1Layer, showW1, size, 0, 0, cv::INTER_NEAREST);
	cv::normalize(showW1, showW1, 0, 255, cv::NORM_MINMAX, showW1.type());
	cv::imwrite("./img/showW1.png", showW1);
	cv::resize(r2Layer, showW2, size, 0, 0, cv::INTER_NEAREST);
	cv::normalize(showW2, showW2, 0, 255, cv::NORM_MINMAX, showW2.type());
	cv::imwrite("./img/showW2.png", showW2);

	// base_channels.release();
	sumOfDetail.release();
	sumOfBase.release();
	r1Layer.release();
	r2Layer.release();
	std::cout << "Detail maximalization -- COMPLETED" << std::endl;
	// std::cout << detailMaximizedLayerY*255 << std::endl;

	cv::Mat &s = ST[0];
	cv::Mat &t = ST[1];
	// s += 1;

	std::ofstream sfile("s-sparse.txt", std::ios::out);
	if (sfile.is_open()) {
		sfile << s;
		sfile.close();
	}
	
	std::ofstream tfile("t-sparse.txt", std::ios::out);
	if (tfile.is_open()) {
		tfile << t;
		tfile.close();
	}

	printMatRange("s", s);
	printMatRange("t", t);

	cv::Mat showS, showT;
	cv::resize(s, showS, size, 0, 0, cv::INTER_NEAREST);
	cv::resize(t, showT, size, 0, 0, cv::INTER_NEAREST);
	cv::resize(originalImage, originalImage, size, 0, 0, cv::INTER_NEAREST);
	// showImg("s", showS);
	// showImg("t", showT);
	cv::normalize(showS, showS, 0, 255, cv::NORM_MINMAX, showS.type());
	cv::normalize(showT, showT, 0, 255, cv::NORM_MINMAX, showT.type());
	cv::normalize(originalImage, originalImage, 0, 255, cv::NORM_MINMAX, originalImage.type());
	cv::imwrite("./img/showS.png", showS);
	cv::imwrite("./img/showT.png", showT);
	cv::imwrite("./img/showI.png", originalImage);
	// (s).convertTo(s, CV_32F);
	// (t).convertTo(t, CV_32F);
	
	// (mu*s + (1-mu))*D + B + mu*t
	// cv::Mat detailMaximizedLayerR = (mu*s + (1-mu))*detailLayerR/255.0 + basePhase2R/255.0 + mu*t;
    // cv::Mat detailMaximizedLayerG = (mu*s + (1-mu))*detailLayerG/255.0 + basePhase2G/255.0 + mu*t;
    // cv::Mat detailMaximizedLayerB = (mu*s + (1-mu))*detailLayerB/255.0 + basePhase2B/255.0 + mu*t;

	cv::Mat detailMaximizedLayerR = getDetailControl(basePhase2R, detailLayerR, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerG = getDetailControl(basePhase2G, detailLayerG, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerB = getDetailControl(basePhase2B, detailLayerB, ST[0], ST[1], mu, height, width);
	
	std::vector<cv::Mat> detailMaximizedLayers;
	detailMaximizedLayers.push_back(detailMaximizedLayerB);
	detailMaximizedLayers.push_back(detailMaximizedLayerG);
	detailMaximizedLayers.push_back(detailMaximizedLayerR);
	cv::Mat showE;
	cv::merge(detailMaximizedLayers, showE);
	cv::resize(showE, showE, size, 0, 0, cv::INTER_NEAREST);
	cv::normalize(showE, showE, 0, 255, cv::NORM_MINMAX, showE.type());
	cv::imwrite("./img/showE.png", showE);
	
	/*
	 * Function for control details enhancement of picture 
	 **/
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{				
			// *pDestinationData++ = basePhase1.at<cv::Vec3b>(j,i)[2];
			// *pDestinationData++ = basePhase1.at<cv::Vec3b>(j,i)[1];
			// *pDestinationData++ = basePhase1.at<cv::Vec3b>(j,i)[0];
			*pDestinationData++ = ((detailMaximizedLayerR).at<float>(j,i));
			*pDestinationData++ = ((detailMaximizedLayerG).at<float>(j,i));
			*pDestinationData++ = ((detailMaximizedLayerB).at<float>(j,i));
		}
	}
	pDst->Convert(TMO_RGB);
	return 0;
}
