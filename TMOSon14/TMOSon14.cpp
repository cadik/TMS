/* --------------------------------------------------------------------------- *
 * TMOSon14.cpp: implementation of the TMOSon14 class.   *
 * --------------------------------------------------------------------------- */
#include "TMOSon14.h"
#include <fftw3.h>
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
TMOSon14::TMOSon14()
{
	SetName(L"Son14");						// TODO - Insert operator name
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

    array_to_merge.push_back(b*256);
    array_to_merge.push_back(g*256);
    array_to_merge.push_back(r*256);

    cv::Mat originalImage;
    
    cv::merge(array_to_merge, originalImage);
   
	/*
	 * Base Decomposition
	 **/
	 
	/*
	 * Phase 1 - Original L0 Smoothing
	 **/
	  
	// std::vector<cv::Mat> result = minimizeL0Gradient(color);
	cv::Mat basePhase1 = minimizeL0Gradient1(originalImage);
	 /*
     * split basePhase1, experiment; 
     **/
     
	cv::Mat basePhase1Chan[3];
	cv::split(basePhase1, basePhase1Chan);
	
    (basePhase1Chan[0]).convertTo(basePhase1Chan[0], CV_32F);
    (basePhase1Chan[1]).convertTo(basePhase1Chan[1], CV_32F);
    (basePhase1Chan[2]).convertTo(basePhase1Chan[2], CV_32F);  
	/*
	 * Phase 2 - L0 smooting with adaptive lambda matrix
	 **/	
	cv::Mat gradientFrom1stSmoothing = getGradientMagnitude(basePhase1);
	cv::Mat adaptiveLambdaMatrix1 = getAdaptiveLambdaMatrix(gradientFrom1stSmoothing, height, width);
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
     
    cv::Mat sumOfCostsBase = getSumOfCosts(basePhase2Chan[0], basePhase2Chan[1], basePhase2Chan[2], height, width);
    /*myfile.open ("./TMOSon14/txt/sumOfCostsBase.txt");
	myfile << sumOfCostsBase;
	myfile.close();*/
	
	cv::Mat sumOfCostsOriginal = getSumOfCosts(r, g, b, height, width);
	
	cv::Mat sigmaMap = stochasticOptimizationForGetSigma(sumOfCostsBase/256.0, sumOfCostsOriginal, height, width, 50000);
	
	cv::Mat basePhase3R = myOwn2DFilter(r*256, sigmaMap, height, width);
	cv::Mat basePhase3G = myOwn2DFilter(g*256, sigmaMap, height, width);
	cv::Mat basePhase3B = myOwn2DFilter(b*256, sigmaMap, height, width);
	
	/*
	 * DETAIL MAXIMALIZATION
	 **/
	/*
     * Getting Detail Layer
     **/
     
    cv::Mat detailLayerR = getDetailLayer(r*256, basePhase3R, height, width);
    cv::Mat detailLayerG = getDetailLayer(g*256, basePhase3G, height, width);
    cv::Mat detailLayerB = getDetailLayer(b*256, basePhase3B, height, width);
    
    /*
     * Getting Sum of Cost of detailLayer
     * */
    cv::Mat sumOfDetail = getSumOfCosts(detailLayerR, detailLayerG, detailLayerB, height, width);
    cv::Mat sumOfBase = getSumOfCosts(basePhase3R, basePhase3G, basePhase3B, height, width);
    
    
    /*
     * Gettig gradient from base image
     * */
    std::vector<cv::Mat> array_to_merge1;

    array_to_merge1.push_back(basePhase3R);
    array_to_merge1.push_back(basePhase3G);
    array_to_merge1.push_back(basePhase3B);

    cv::Mat baseImage;
    
    cv::merge(array_to_merge1, baseImage);
    
    cv::Mat gradientOfBaseLayer = getGradientMagnitude(baseImage);
    
    /*
     * Getting weigth r1 and r2 layers 
     **/
     
    cv::Mat r1Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 200);
    cv::Mat r2Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 500);
	
	std::vector<cv::Mat> detail;
	detail.push_back((detailLayerR.clone()) / 256.0);
	detail.push_back((detailLayerG.clone()) / 256.0);
	detail.push_back((detailLayerB.clone()) / 256.0);
	
	std::vector<cv::Mat> ST = detailMaximalization(sumOfBase/256.0, sumOfDetail/256.0, r1Layer, r2Layer, height, width, 30000, detail);	
	/*
	Replace function above, with quadtratic solution
	*/
	cv::Mat detailMaximizedLayerR = getDetailControl(basePhase3R, detailLayerR, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerG = getDetailControl(basePhase3G, detailLayerG, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerB = getDetailControl(basePhase3B, detailLayerB, ST[0], ST[1], mu, height, width);
myfile.open ("./TMOSon14/txt/sumOfCostsBase.txt");
	myfile << sumOfCostsBase;
	myfile.close();
	/*
	 * Function for control details enhancement of picture 
	 **/

	/*
	 * Showing picture (shows blurred picture)
	 **/
	/* for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{													// simple variables
			// and store results to the destination image
			*pDestinationData++ = ((basePhase3R).at<float>(j,i) + (detailLayerR).at<float>(j,i)) / 256.0;
			*pDestinationData++ = ((basePhase3G).at<float>(j,i) + (detailLayerG).at<float>(j,i)) / 256.0;
			*pDestinationData++ = ((basePhase3B).at<float>(j,i) + (detailLayerB).at<float>(j,i)) / 256.0;
		}
	}*/
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{													// simple variables
			// and store results to the destination image
			*pDestinationData++ = (detailMaximizedLayerR).at<float>(j,i) / 256.0;// + (detailChan[2]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = (detailMaximizedLayerG).at<float>(j,i) / 256.0;// + (detailChan[1]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = (detailMaximizedLayerB).at<float>(j,i) / 256.0;// + (detailChan[0]).at<float>(j,i)) / 256.0;
		}
	}
	pDst->Convert(TMO_RGB);
return 0;
}

