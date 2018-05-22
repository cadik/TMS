/* --------------------------------------------------------------------------- *
 * TMOCLAHE.cpp: implementation of the TMOCLAHE class.   *
 * --------------------------------------------------------------------------- */
#include "TMOCLAHE.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;
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

TMOCLAHE::TMOCLAHE()
{
	SetName(L"TMOCLAHE");						// TODO - Insert operator name
	SetDescription(L"Contrast Limited Adaptive Histogram Equalization");	// TODO - Insert description
	/**
	 * clip limit
	 **/
	cl.SetName(L"clipLimit");				// TODO - Insert parameters names
	cl.SetDescription(L"Represents clip limit for histogram contrast limiting");	// TODO - Insert parameter descriptions
	cl.SetDefault(0.5);							// TODO - Add default values
	cl=0.01;
	cl.SetRange(0.0, 15.0);				// TODO - Add acceptable range if needed
	this->Register(cl);
	/**
	 * region Size
	 **/
	gridRegions.SetName(L"regionSize");				// TODO - Insert parameters names
	gridRegions.SetDescription(L"Represents size of cotextual region in the picture.");	// TODO - Insert parameter descriptions
	gridRegions.SetDefault(9);							// TODO - Add default values
	gridRegions=2;
	gridRegions.SetRange(9, 999);				// TODO - Add acceptable range if needed
	this->Register(gridRegions);
}

TMOCLAHE::~TMOCLAHE()
{
}


/* --------------------------------------------------------------------------- *
 * This overloaded function
 * -------------------------------------------------------- */
int TMOCLAHE::Transform()
{ 
	ofstream myfile;
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	/*
	 * Base matrix 
	 **/
    cv::Mat Y;
	cv::Mat x;
	cv::Mat y;
	
	Y = cv::Mat::zeros(height, width, CV_32F);
	x = cv::Mat::zeros(height, width, CV_32F);
	y = cv::Mat::zeros(height, width, CV_32F);

	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
											// of three doubles representing


	/*
	 * Fill matrix
	 * */
	 int j;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			Y.at<float>(j,i) = *pSourceData++; 
			x.at<float>(j,i) = *pSourceData++;  //getting separate RGB channels
			y.at<float>(j,i) = *pSourceData++;
		}
	}
	cv::Mat newImage;
	newImage = cv::Mat::zeros(height, width, CV_32F);	
    newImage = histogramEqualization(Y, height, width, gridRegions, cl);


	pSrc->ProgressBar(j, pSrc->GetHeight());
	/*
	 * Function for control details enhancement of picture 
	 **/
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{											// simple variables		
			*pDestinationData++ = newImage.at<float>(j,i);// + (detailChan[2]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = x.at<float>(j,i);// + (detailChan[1]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = y.at<float>(j,i);
		}
	}
	pDst->Convert(TMO_RGB);

	return 0;
}


