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
	mu.SetDefault(0.5);							// TODO - Add default values
	mu=0.5;
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
	// tryIt();
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
     * split basePhase1
     **/
	cv::Mat basePhase1Chan[3];
	cv::split(basePhase1, basePhase1Chan);
	
    (basePhase1Chan[0]).convertTo(basePhase1Chan[0], CV_32F);
    (basePhase1Chan[1]).convertTo(basePhase1Chan[1], CV_32F);
    (basePhase1Chan[2]).convertTo(basePhase1Chan[2], CV_32F);  
	/*
	 * Phase 2 - L0 smooting with adaptive lambda matrix
	 **/
	std::cout << "Base Phase2" << std::endl;
	cv::Mat gradientFrom1stSmoothing = getGradientMagnitude(basePhase1);
	cv::Mat adaptiveLambdaMatrix1 = getAdaptiveLambdaMatrix(gradientFrom1stSmoothing, height, width);
	// cv::Mat::zeros(height, width, CV_32F); //
	/*for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			adaptiveLambdaMatrix1.at<float>(j,i) = 0.00001;
		}
	}*/
	
	cv::Mat basePhase2 = minimizeL0GradientSecondFaze(originalImage, adaptiveLambdaMatrix1, height, width);
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
     
    /*cv::Mat sumOfCostsBase = getSumOfCostsForSigmaOptimization(basePhase1Chan[2], basePhase2Chan[1], basePhase2Chan[0], height, width);
	
	cv::Mat sumOfCostsOriginal = getSumOfCostsForSigmaOptimization(r, g, b, height, width);
	std::cout << "Base Phase3" << std::endl;*/
	//cv::Mat sigmaMap = optimizeForSigma(height, width, sumOfCostsOriginal/255.0, sumOfCostsBase/255.0, optim1Iteration);
	/////cv::Mat sigmaMap = stochasticOptimizationForGetSigma(sumOfCostsBase/256.0, sumOfCostsOriginal, height, width, 50000);
	
	/*cv::Mat basePhase3R = myOwn2DFilter(r, sigmaMap, height, width);
	cv::Mat basePhase3G = myOwn2DFilter(g, sigmaMap, height, width);
	cv::Mat basePhase3B = myOwn2DFilter(b, sigmaMap, height, width);*/
	std::cout << "Base phase -- COMPLETED" << std::endl;

     
    cv::Mat detailLayerR = getDetailLayer(r, basePhase2Chan[2], height, width);
    cv::Mat detailLayerG = getDetailLayer(g, basePhase2Chan[1], height, width);
    cv::Mat detailLayerB = getDetailLayer(b, basePhase2Chan[0], height, width);
    cv::Mat sumOfDetail = getSumOfCosts(detailLayerR, detailLayerG, detailLayerB, height, width);
    cv::Mat sumOfBase = getSumOfCosts(basePhase2Chan[0], basePhase2Chan[1], basePhase2Chan[2], height, width);

    std::vector<cv::Mat> array_to_merge1;

    array_to_merge1.push_back(basePhase2Chan[2]);
    array_to_merge1.push_back(basePhase2Chan[1]);
    array_to_merge1.push_back(basePhase2Chan[0]);

    cv::Mat baseImage;
    
    cv::merge(array_to_merge1, baseImage);
	/*
		Getting weights
	*/

    cv::Mat gradientOfBaseLayer = getGradientMagnitude(baseImage);
    
	cv::Mat r1Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 200);
    cv::Mat r2Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 500);

	std::vector<cv::Mat> detail;
	detail.push_back((detailLayerR.clone()));
	detail.push_back((detailLayerG.clone()));
	detail.push_back((detailLayerB.clone()));
	std::vector<cv::Mat> ST = detailMaximalization(sumOfBase, sumOfDetail, r1Layer, r2Layer, height, width, 50, detail); 
	// std::vector<cv::Mat> ST = optimizeForGettingSAndTparametersWithCgal(height, width, sumOfDetail, r1Layer, r2Layer, array_to_merge1, detail);
	std::cout << "Detail maximalization -- COMPLETED" << std::endl;

	cv::Mat detailMaximizedLayerR = getDetailControl(basePhase2Chan[2], detailLayerR, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerG = getDetailControl(basePhase2Chan[1], detailLayerG, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerB = getDetailControl(basePhase2Chan[0], detailLayerB, ST[0], ST[1], mu, height, width);
	
	/*
	 * Function for control details enhancement of picture 
	 **/
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{				
			*pDestinationData++ = (detailMaximizedLayerR).at<float>(j,i);// + (detailChan[2]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = (detailMaximizedLayerG).at<float>(j,i);// + (detailChan[1]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = (detailMaximizedLayerB).at<float>(j,i);//
		}
	}
	pDst->Convert(TMO_RGB);
	return 0;
}

