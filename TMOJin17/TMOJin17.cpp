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
#include <cstring>

#include <mlpack/core.hpp>
#include <mlpack/core/optimizers/adam/adam.hpp>
#include <mlpack/core/optimizers/function.hpp>

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
#define OMEGA_SIZE 9

unsigned int IMAGE_WIDTH; //global variable for width of loaded image
unsigned int IMAGE_HEIGHT; //global variable for height of loaded image

//prototype of function used in class under this
double contrastLossComputation(unsigned int *PixNumberInClusters, Mat centersLAB, Mat centersRGB, Mat centersGray, int centersSize, arma::mat Omega, arma::umat ShuffledIndex);


class AdamClassType 
{
public:
	unsigned int *PixNumberInClusters;
 	Mat CentersLAB, CentersRGB, CentersGray;
 	int CentersSize;
 	arma::umat ShuffledIndex;
 	arma::mat oldCoordinates;
 	int MaxIter;
 	double oldContrastLoss;
 	double contrastLoss;

 	void setPixNumberInCluster(unsigned int * MyPixNumberInClusters);
 	void setCentersLAB(Mat MyCentersLAB);
 	void setCentersRGB(Mat MyCentersRGB);
 	void setCentersGray(Mat MyCentersGray);
 	void setCentersSize(int MyCentersSize);
 	void setOldCoordinates(arma::mat MyOldCoordinates);
 	void setMaxIter(int MyMaxIter);
 	void setOldContrastLoss(double MyOldContrastLoss);
 	void setContrastLoss(double MyContrastLoss);

 	//9 indexes of array for shuffling.
	void Shuffle(void) {ShuffledIndex = arma::shuffle(arma::uvec("0 1 2 3 4 5 6 7 8"));}

	//function that computes new Omega values in variable substitute coordinates
	double Evaluate (const arma::mat &coordinates)
	{
		contrastLoss = contrastLossComputation(PixNumberInClusters, CentersLAB, CentersRGB, CentersGray, CentersSize, coordinates, ShuffledIndex);
		return contrastLoss;
	}

	void Gradient (const arma::mat &coordinates, arma::mat &gradient)
	{
		//mlpack:: already deleted Gradient in this moment
		if (MaxIter >= 1000)
			return;
	
		//check changes 
		for (unsigned int i = 0; i < OMEGA_SIZE; ++i)
			gradient[ShuffledIndex[i]] = (oldContrastLoss - contrastLoss) / (oldCoordinates[ShuffledIndex[i]] - coordinates[ShuffledIndex[i]]);

		MaxIter++;

	}

	//function that calls other functions from this class AdamClassType
	double EvaluateWithGradient (const arma::mat &coordinates, const size_t begin, arma::mat &gradient, const size_t batchSize)
	{
		Shuffle();
		Evaluate(coordinates);
		Gradient(coordinates, gradient);

		oldCoordinates = coordinates;
		oldContrastLoss = contrastLoss;
	}

	size_t NumFunctions(void) {return 1;}

};

//functions for setting starting values of variables for AdamClassType
void AdamClassType::setPixNumberInCluster(unsigned int * MyPixNumberInClusters)
{
	PixNumberInClusters = MyPixNumberInClusters;
}

void AdamClassType::setCentersLAB(Mat MyCentersLAB)
{
	CentersLAB = MyCentersLAB;
}

void AdamClassType::setCentersRGB(Mat MyCentersRGB)
{
	CentersRGB = MyCentersRGB;
}

void AdamClassType::setCentersGray(Mat MyCentersGray)
{
	CentersGray = MyCentersGray;
}

void AdamClassType::setCentersSize(int MyCentersSize)
{
	CentersSize = MyCentersSize;
}


void AdamClassType::setOldCoordinates(arma::mat MyOldCoordinates)
{
	oldCoordinates = MyOldCoordinates;
}

void AdamClassType::setMaxIter(int MyMaxIter)
{
	MaxIter = MyMaxIter;
}


void AdamClassType::setOldContrastLoss(double MyOldContrastLoss)
{
	oldContrastLoss = MyOldContrastLoss;
}

void AdamClassType::setContrastLoss(double MyContrastLoss)
{
	contrastLoss = MyContrastLoss;
}

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOJin17::TMOJin17()
{
	SetName(L"Jin17");
	SetDescription(L"Preserving perceptual contrast in decolorization with optimized color orders");	

	K_clusters.SetName(L"K");
	K_clusters.SetDescription(L"Value of K for K-means clustering ; <5,50>");
	K_clusters.SetDefault(30);
	K_clusters=30;		
	K_clusters.SetRange(5,50);
	this->Register(K_clusters);

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


//color conversion from CIELab to RGB
//Always must go through XYZ color space
//i am not author of this function: inpired by
//https://stackoverflow.com/questions/9372626/converting-lab-values-to-rgb-values-in-opencv?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
void convertLAB2RGB(Mat srcLab, Mat &dstRGB)
{
	float L, a, b;
	float X, Y, Z;
	float R, G, B;

	// Lab -> normalized XYZ (X,Y,Z are all in 0...1)
	for (unsigned int i = 0; i < srcLab.size().height; ++i)
	{
		L = srcLab.at <float> (i, SHIFT_TO_L);
		a = srcLab.at <float> (i, SHIFT_TO_A);
		b = srcLab.at <float> (i, SHIFT_TO_B);


		Y = L * (1.0/116.0) + 16.0/116.0;
		X = a * (1.0/500.0) + Y;
		Z = b * (-1.0/200.0) + Y;

		X = X > 6.0/29.0 ? X * X * X : X * (108.0/841.0) - 432.0/24389.0;
		Y = L > 8.0 ? Y * Y * Y : L * (27.0/24389.0);
		Z = Z > 6.0/29.0 ? Z * Z * Z : Z * (108.0/841.0) - 432.0/24389.0;

		// normalized XYZ -> linear sRGB (in 0...1)

		R = X * (1219569.0/395920.0)     + Y * (-608687.0/395920.0)    + Z * (-107481.0/197960.0);
		G = X * (-80960619.0/87888100.0) + Y * (82435961.0/43944050.0) + Z * (3976797.0/87888100.0);
		B = X * (93813.0/1774030.0)      + Y * (-180961.0/887015.0)    + Z * (107481.0/93370.0);

		// linear sRGB -> gamma-compressed sRGB (in 0...1)

		R = R > 0.0031308 ? pow(R, 1.0 / 2.4) * 1.055 - 0.055 : R * 12.92;
		G = G > 0.0031308 ? pow(G, 1.0 / 2.4) * 1.055 - 0.055 : G * 12.92;
		B = B > 0.0031308 ? pow(B, 1.0 / 2.4) * 1.055 - 0.055 : B * 12.92;

		dstRGB.at <float> (i, SHIFT_TO_L) = R * 255; 
		dstRGB.at <float> (i, SHIFT_TO_A) = G * 255; 
		dstRGB.at <float> (i, SHIFT_TO_B) = B * 255; 

	}

}


//decolorization method via order set
void colorToGrayscaleConversion(Mat sourceImageRGB, Mat &destImageGray, int Mat_height, int Mat_width, arma::mat Omega, arma::umat ShuffledIndex)
{

	unsigned int Omega_array_index = 0;
	double sourceImageRGB_r, sourceImageRGB_g, sourceImageRGB_b;
	double result;

	for (unsigned int r = 0; r <= 2; ++r)
	{
		for (unsigned int g = 0; g <= 2; ++g)
		{
			for (unsigned int b = 0; b <= 2; ++b)
			{
				
				if ( ((r + g + b) <= 2) &&  ((r + g + b) > 0) )
				{
					
					for (unsigned int row = 0; row < Mat_height; ++row)
					{
				
						for (unsigned int col = 0; col < Mat_width; ++col)
						{
							//taking all RGB values for one pixel and computing new ones to destination image
							sourceImageRGB_r = sourceImageRGB.at <float> (col + row * Mat_width, SHIFT_TO_L);
							sourceImageRGB_g = sourceImageRGB.at <float> (col + row * Mat_width, SHIFT_TO_A);
							sourceImageRGB_b = sourceImageRGB.at <float> (col + row * Mat_width, SHIFT_TO_B);

							//multiplying values with computed weight coefficients
							result = Omega[ShuffledIndex[Omega_array_index]] * pow(sourceImageRGB_r, r) * pow(sourceImageRGB_g, g) * pow(sourceImageRGB_b, b);

							//actualization of destination image for current pixel
							destImageGray.at <float> (row, col) += result;

						}
					}
					
					Omega_array_index++;
				}
			}	
		}
	}

}

//function that computes Contrast Loss between lab and grayscale. We are looking for minimum returned double with Adam Gradient Descent solver
double contrastLossComputation(unsigned int *PixNumberInClusters, Mat centersLAB, Mat centersRGB, Mat centersGray, int centersSize, arma::mat Omega, arma::umat ShuffledIndex)
{
	double weightClusterPair;

	unsigned int signFunReturned;
	double midStepCalculation;

	double contrastLossSum = 0.0;

	colorToGrayscaleConversion(centersLAB, centersGray, centersLAB.size().height, 1, Omega, ShuffledIndex);

	//main cycle for Contrast Loss computation
	for (unsigned int Ci = 0; Ci < centersSize; ++Ci)
	{
		for (unsigned int Cj = 0; Cj < centersSize; ++Cj)
		{

			midStepCalculation = ( 1 / log(1 + BETA) ) * ( BETA / ( BETA * (centersGray.at <float> (Cj, 0)) + 1) ) * ( (centersGray.at <float> (Ci, 0)) - (centersGray.at <float> (Cj, 0)) ); 

			signFunReturned = SignFunction(centersRGB, Ci, Cj);

			if (signFunReturned == 2)
				midStepCalculation *= -1;

			else if (signFunReturned == 3)
				midStepCalculation = std::abs(midStepCalculation);

			//second square root for this three channels
			midStepCalculation -= sqrt(
				//(L_Ci - L_Cj)^2 +
				pow( (centersLAB.at <float> (Ci, SHIFT_TO_L) - (centersLAB.at <float> (Cj, SHIFT_TO_L)) ), 2) +
				//(A_Ci - A_Cj)^2 +
				pow( (centersLAB.at <float> (Ci, SHIFT_TO_A) - (centersLAB.at <float> (Cj, SHIFT_TO_A)) ), 2) +
				//(B_Ci - B_Cj)^2
				pow( (centersLAB.at <float> (Ci, SHIFT_TO_B) - (centersLAB.at <float> (Cj, SHIFT_TO_B)) ), 2)
		 	);

			//^2 for the rest
			midStepCalculation = pow(midStepCalculation, 2);

			weightClusterPair = (PixNumberInClusters[Ci] * PixNumberInClusters[Cj]) / pow((0.01 * IMAGE_WIDTH * IMAGE_HEIGHT), 2);

			contrastLossSum += weightClusterPair * midStepCalculation;
		}
	}

	return contrastLossSum;	
}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOJin17::Transform()
{

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

	// Flag to specify the number of times the algorithm is executed using different initial labellings.
	//The algorithm returns the labels that yield the best compactness (see the last function parameter).
	int attempts = 5;
	
	//Input/output integer array that stores the cluster indices for every sample.
	Mat labelsLAB, labelsGray;

	//Output matrix of the cluster centers, one row per each cluster center.
	Mat centersLAB;

	//function for kmeans clustering. Image color reduction. For faster computations and memory save up
	kmeans(sourceImage_2D, K_clusters, labelsLAB, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centersLAB);


	unsigned int cluserID;

	//###CONTRAST LOSS COMPUTATION###

	int labelsSize = labelsLAB.size().height;
	int centersSize = centersLAB.size().height;

	unsigned int *PixNumberInClusters = new unsigned int [centersSize];
	//fill array PixNumberInClusters with zeroes - starting value
	std::fill_n(PixNumberInClusters, centersSize, 0);

	unsigned int indexToArray;
	
	//computation for number of pixels in specific clusters
	for (unsigned int i = 0; i < labelsSize; ++i)
	{
		//find value in labels, take it as index
		indexToArray = labelsLAB.at <int> (i, 0);
		
		//on that index increment number of occurrence
		PixNumberInClusters[indexToArray] += 1;
	}

	Mat centersRGB = Mat::zeros(centersLAB.size().height, centersLAB.size().width, CV_32F);	
	convertLAB2RGB(centersLAB, centersRGB);
	
	Mat centersGray = Mat::zeros(centersLAB.size().height, 1, CV_32F);


	//					w_b  w_b^2  w_g  w_gb w_g^2 w_r w_rb w_rg w_r^2
	arma::mat Omega = {0.33, 0.01, 0.33, 0.03, 0.02, 0.33, 0.05, 0.04, 0.06};

	//Initialization of Adam Gradient Descent with parameters StepSize and MaxIterations
	mlpack::optimization::AdamType<mlpack::optimization::AdamUpdate> AdamGradientDescent;

	double LearnRate = 0.0005;


	double &stepSizeRef = AdamGradientDescent.StepSize();
	stepSizeRef = LearnRate;

	size_t &MaxIterationsRef = AdamGradientDescent.MaxIterations();
	MaxIterationsRef = 1000;


	double &Tolerance = AdamGradientDescent.Tolerance();
	Tolerance = 0.0;	

	//Instanciation of class AdamClassType
	AdamClassType AdamClassTypeInit;
	AdamClassTypeInit.setPixNumberInCluster(PixNumberInClusters);
	AdamClassTypeInit.setCentersLAB(centersLAB);
	AdamClassTypeInit.setCentersRGB(centersRGB);
	AdamClassTypeInit.setCentersGray(centersGray);
	AdamClassTypeInit.setCentersSize(centersSize);
	AdamClassTypeInit.setOldCoordinates({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
	AdamClassTypeInit.setMaxIter(0);
	AdamClassTypeInit.setOldContrastLoss(0.0);
	AdamClassTypeInit.setContrastLoss(0.0);

	AdamGradientDescent.Optimize(AdamClassTypeInit, Omega);


 	//dealocation of array for storage
	delete[] PixNumberInClusters;


	//###Final decolorization###

	//convert color spaces back to RGB
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	//save backup pointer to image arrays
	double *destinationImage_P_backup;
	double *sourceImage_P_backup;
	destinationImage_P_backup = destinationImage;
	sourceImage_P_backup = sourceImage;

	double sourceImage_r, sourceImage_g, sourceImage_b;
	double result;
	unsigned int Omega_array_index = 0;

	for (unsigned int r = 0; r <= 2; ++r)
	{
		for (unsigned int g = 0; g <= 2; ++g)
		{
			for (unsigned int b = 0; b <= 2; ++b)
			{
				//it takes exactly what we expect in Omega.
				//If exponent for color is 2, other colors have zero etc.
				if ( ((r + g + b) <= 2) &&  ((r + g + b) > 0) )
				{
					//setting pointer to source and destination image again on starting position.
					destinationImage = destinationImage_P_backup;
					sourceImage = sourceImage_P_backup;

					//going through all pixels, first x - cols then y - rows
					for (unsigned int y = 0; y < IMAGE_HEIGHT; ++y)
					{
						pSrc->ProgressBar(y, pSrc->GetHeight());
						for (unsigned int x = 0; x < IMAGE_WIDTH; ++x)
						{
							//taking all RGB values for one pixel and computing new ones to destination image
							sourceImage_r = *(sourceImage++);
							sourceImage_g = *(sourceImage++);
							sourceImage_b = *(sourceImage++);

							//multiplying values with computed weight coefficients
							result = Omega[AdamClassTypeInit.ShuffledIndex[Omega_array_index]] * pow(sourceImage_r, r) * pow(sourceImage_g, g) * pow(sourceImage_b, b);

							//actualization of destination image for current pixel
							*(destinationImage++) += result;
							*(destinationImage++) += result;
							*(destinationImage++) += result;

						}
					}

					//color recomputed for all pixels, lets take another color combination from Omega
					Omega_array_index++;
						
				}
			}	
		}
	}


	return 0;
}