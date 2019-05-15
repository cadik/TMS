/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                Content Adaptive Image Detail Enhancement (2015)              *
*              by Fei Kou, Weihai Chen, Zhengguo Li, Changyun Wen              *
*                                                                              *
*                            part of diploma thesis                            *
*             Author: Tomas Hudziec [xhudzi01 AT stud.fit.vutbr.cz]            *
*                                    Brno 2019                                 *
*                                                                              *
*******************************************************************************/

/* --------------------------------------------------------------------------- *
 * TMOKou15.cpp: implementation of the TMOKou15 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOKou15.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOKou15::TMOKou15()
{
	SetName(L"Kou15");
	SetDescription(L"Content Adaptive Image Detail Enhancement");

	etaParameter.SetName(L"eta");
	etaParameter.SetDescription(L"Parameter controlling how many times of detail adding to the input image.");
	etaParameter.SetDefault(4.0);
	etaParameter=4.0;
	etaParameter.SetRange(0.1, 10.0);

	lambdaParameter.SetName(L"lambda");
	lambdaParameter.SetDescription(L"Lagrangian factor - smoothing parameter controlling the degree of enhance. Powers of 2 are mostly used.");
	lambdaParameter.SetDefault(0.16);
	lambdaParameter=0.16;
	lambdaParameter.SetRange(0.016, 0.64);

	kappaParameter.SetName(L"kappa");
	kappaParameter.SetDescription(L"Rate control parameter. Small kappa results in more iterations and with sharper edges.");
	kappaParameter.SetDefault(2.0);
	kappaParameter=2.0;
	kappaParameter.SetRange(1.1, 3.0);

	this->Register(kappaParameter);
	this->Register(lambdaParameter);
	this->Register(etaParameter);
}

TMOKou15::~TMOKou15()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOKou15::Transform()
{
	// get parameters
	float eta = etaParameter;
	// recomputation of lambda from original L0 smoothing alg. to modified L0 enhancing alg.
	float lambda = lambdaParameter/(eta*eta);
	float kappa = kappaParameter;
	std::cout << "original lambda from L0 smoothing algorithm = " << lambdaParameter << std::endl;
	std::cout << "recomputed lambda for modified L0 enhancing algorithm = " << lambda << " (= lambda/(eta^2))" << std::endl;

	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_RGB);								// This is format of Y as luminance
	// pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components

	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();

	cv::Mat I(height, width, CV_32FC3); // INPUT IMAGE IN RGB 
	double r, g, b;

	int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			I.at<cv::Vec3f>(j,i)[0] = r = *pSourceData++;
			I.at<cv::Vec3f>(j,i)[1] = g = *pSourceData++;
			I.at<cv::Vec3f>(j,i)[2] = b = *pSourceData++;
		}
	}
	
	// L0 smoothing
	// std::vector<cv::Mat> B_layers = minimizeL0Gradient(I);
	// cv::Mat B;
	// cv::merge(B_layers, B);
	// cv::normalize(I, I, 0, 1, cv::NORM_MINMAX, I.type());
	// L0 enhancing
	cv::Mat outputImage = minimizeL0Gradient(inputImage, eta, lambda, kappa);

	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			// store results to the destination image
			*pDestinationData++ = B.at<cv::Vec3b>(j,i)[0];
			*pDestinationData++ = B.at<cv::Vec3b>(j,i)[1];
			*pDestinationData++ = B.at<cv::Vec3b>(j,i)[2];
		}
	}

	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}
