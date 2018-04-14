/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                 Bachelor thesis                              *
*             Author: Filip Brezna [xbrezn00 AT stud.fit.vutbr.cz]             *
*                                    Brno 2018                                 *
*                                                                              *
*******************************************************************************/


/* --------------------------------------------------------------------------- *
 * TMOJin17.cpp: implementation of the TMOJin17 class.   		  			   *
 * --------------------------------------------------------------------------- */

#include "TMOJin17.h"
#include <stdlib.h>
#include <math.h>

//tmolib contains definition for EPS,
//so does openCV so we undefine it and then put it back. After link
#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#undef EPS

using namespace cv;

#define BETA 0.5
#define CIELAB_NUM_CHANNELS 3
#define SHIFT_TO_L 0
#define SHIFT_TO_A 1
#define SHIFT_TO_B 2			

unsigned int IMAGE_WIDTH; //global variable for width of loaded image
unsigned int IMAGE_HEIGHT; //global variable for height of loaded image


#include <fstream>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOJin17::TMOJin17()
{
	SetName(L"Jin17");						
	SetDescription(L"Preserving perceptual contrast in decolorization with optimized color orders");	

	sigma.SetName(L"Sigma");
	sigma.SetDescription(L"Sigma for image blur, removing noise ; <0,10>");	
	sigma.SetDefault(0);
	sigma=0;		
	sigma.SetRange(0,10);
	this->Register(sigma);
}

TMOJin17::~TMOJin17()
{
}

unsigned int SignFunction(Mat centersRGB, int x, int y)
{

	//if r_x >= r_y AND g_x >= g_y AND b_x >= b_y
	if  (
			//SHIFT_TO_L == SHIFT_TO_R
			( centersRGB.at <float> (x, SHIFT_TO_L) >= centersRGB.at <float> (y, SHIFT_TO_L) ) &&
			//SHIFT_TO_L == SHIFT_TO_G
			( centersRGB.at <float> (x, SHIFT_TO_L) >= centersRGB.at <float> (y, SHIFT_TO_A) ) &&
		    //SHIFT_TO_L == SHIFT_TO_B
			( centersRGB.at <float> (x, SHIFT_TO_L) >= centersRGB.at <float> (y, SHIFT_TO_B) ) 
		)
		return 1;

	//if r_x < r_y AND g_x < g_y AND b_x < b_y
	else if
		(
			//SHIFT_TO_L == SHIFT_TO_R
			( centersRGB.at <float> (x, SHIFT_TO_L) < centersRGB.at <float> (y, SHIFT_TO_L) ) &&
			//SHIFT_TO_L == SHIFT_TO_G
			( centersRGB.at <float> (x, SHIFT_TO_L) < centersRGB.at <float> (y, SHIFT_TO_A) ) &&
		    //SHIFT_TO_L == SHIFT_TO_B
			( centersRGB.at <float> (x, SHIFT_TO_L) < centersRGB.at <float> (y, SHIFT_TO_B) ) 
		)
		return 2;

	//else
	return 3;
}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOJin17::Transform()
{
	std::ofstream myfile;
	myfile.open ("test.txt");


 
	// Source image is stored in local parameter pSrc
	pSrc->Convert(TMO_LAB);								// This is format of CIELab
	// Destination image is in pDst
	pDst->Convert(TMO_LAB);								// This is format of CIELab

	double *sourceImage = pSrc->GetData();				// You can work at low level data
	double *destinationImage = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components

	//saving values of image width and height to global variables
	IMAGE_WIDTH = pSrc->GetWidth();
	IMAGE_HEIGHT = pSrc->GetHeight();

	unsigned int row, col, CIELab_channel;


	Mat sourceImage_2D (IMAGE_HEIGHT * IMAGE_WIDTH, CIELAB_NUM_CHANNELS, CV_32F);


	//go through all image pixels
	for (row = 0; row < IMAGE_HEIGHT; ++row)
	{
		for (col = 0; col < IMAGE_WIDTH; ++col)
		{
			//we copy and transform source image into 2D picture [row and col][color channel]
			for (CIELab_channel = 0; CIELab_channel < CIELAB_NUM_CHANNELS; ++CIELab_channel)
				sourceImage_2D.at <float> (row + col * IMAGE_HEIGHT, CIELab_channel) = *(sourceImage + ((col + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + CIELab_channel);	
		
		}
	}		

	//Number of clusters to split the set by.
	int clusterCount = 30;

	// Flag to specify the number of times the algorithm is executed using different initial labellings.
	//The algorithm returns the labels that yield the best compactness (see the last function parameter).
	int attempts = 5;
	
	//Input/output integer array that stores the cluster indices for every sample.
	Mat labels;

	//Output matrix of the cluster centers, one row per each cluster center.
	Mat centers;

	//function for kmeans clustering. Image color reduction. For faster computations and memory save up
	kmeans(sourceImage_2D, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centers);

	//myfile << "labels:" << labels << std::endl;
	//myfile << "centers:" << centers << std::endl;
	unsigned int cluserID;

	//go through all image pixels
	for (row = 0; row < IMAGE_HEIGHT; ++row)
	{	
		for (col = 0; col < IMAGE_WIDTH; ++col)
	    { 
	    	//get id of cluster on row and col
	    	//zero because it is col with zero ID. We have only one anyway
			cluserID = labels.at <int> (row + col * IMAGE_HEIGHT, 0);

			//return cluster from center on clusterID and save it to destination image
			*(destinationImage + ((col + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_TO_L) = centers.at <float> (cluserID, SHIFT_TO_L);
			*(destinationImage + ((col + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_TO_A) = centers.at <float> (cluserID, SHIFT_TO_A);
			*(destinationImage + ((col + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_TO_B) = centers.at <float> (cluserID, SHIFT_TO_B);
		}
	}
		
	
	//contrast loss computation

	int labelsSize = labels.size().height;
	int centersSize = centers.size().height;

	unsigned int *PixNumberInClusters = new unsigned int [centersSize];
	//fill array PixNumberInClusters with zeroes - starting value
	std::fill_n(PixNumberInClusters, centersSize, 0);

	unsigned int indexToArray;
	
	//computation for number of pixels in specific clusters
	for (unsigned int i = 0; i < labelsSize; ++i)
	{
		//find value in labels, take it as index
		indexToArray = labels.at <int> (i, 0);
		
		//on that index increment number of occurrence
		PixNumberInClusters[indexToArray] += 1;
	}

	double contrastLossSum = 0.0;
	double weightClusterPair;

	Mat centersXYZ;	

	//LAB -> XYZ; XYZ -> RGB
	//cvtColor(centers, centersXYZ, cv::COLOR_RGB2XYZ);

	Mat centersGray;
		myfile << "tady anno";
	//cvtColor(centersRGB, centersGray, cv::COLOR_RGB2GRAY);
	//myfile << centersGray;
	unsigned int signFunReturned;
	double temporary;
/*
	//main cycle for Contrast Loss computation
	for (unsigned int Ci = 0; Ci < centersSize; ++Ci)
	{
		for (unsigned int Cj = 0; Cj < centersSize; ++Cj)
		{
					
			temporary = 
			(  ( ( 1 / log(1 + BETA) ) * ( BETA / ( BETA * (centersGray.at <float> (Cj, 0)) + 1) ) * ( (centersGray.at <float> (Ci, 0)) - (centersGray.at <float> (Cj, 0)) ) )  - 
				//second square root for this three channels
				sqrt(
					//(L_Ci - L_Cj)^2 +
					pow((centers.at <float> (Ci, SHIFT_TO_L) - (centers.at <float> (Cj, SHIFT_TO_L))), 2) +
					//(A_Ci - A_Cj)^2 +
					pow((centers.at <float> (Ci, SHIFT_TO_L) - (centers.at <float> (Cj, SHIFT_TO_L))), 2) +
					//(B_Ci - B_Cj)^2
					pow((centers.at <float> (Ci, SHIFT_TO_L) - (centers.at <float> (Cj, SHIFT_TO_L))), 2)
				 	)
	
			);

			signFunReturned = SignFunction(centersRGB, Ci, Cj);

			if (signFunReturned == 2)
				temporary *= -1;

			else if (signFunReturned == 3)
				temporary = std::abs(temporary);

			//^2 for the rest
			pow(temporary, 2);

			weightClusterPair = (PixNumberInClusters[Ci] * PixNumberInClusters[Cj]) / pow((0.01 * IMAGE_WIDTH * IMAGE_HEIGHT), 2);

			contrastLossSum += weightClusterPair * temporary;
		}
	}

	myfile << "sum:" << contrastLossSum << std::endl;
*/
 	//dealocation of array for storage
	delete[] PixNumberInClusters;

  myfile.close();
 
	return 0;
}