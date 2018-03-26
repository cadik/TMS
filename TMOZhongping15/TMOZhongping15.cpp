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
 * TMOZhongping15.cpp: implementation of the TMOZhongping15 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOZhongping15.h"
#include <math.h>
#include <cmath> 

 #include <iostream>
 #include <fstream>

#define CIELAB_NUM_CHANNELS 3
#define GAUSSIAN_FILTER_SIZE 5

const double EULER = 2.71828182845904523536;

unsigned int IMAGE_WIDTH; //global variable for width of loaded image
unsigned int IMAGE_HEIGHT; //global variable for height of loaded image


std::ofstream myfile;


/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOZhongping15::TMOZhongping15()
{
	SetName(L"Zhongping15");						
	SetDescription(L"Efficient decolorization preserving dominant distinctions");	

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOZhongping15::~TMOZhongping15()
{
}


//compitation of gradients, we won't consider X and Y axis as separated
//so we are counting Gradient magnitude (size of summed vector) for every CIELab channel - L, A, B
void computeGradientPlane (double  *image, double **gradientPlane)
{
	double gradient_cur, gradient_col_neighbor, gradient_row_neighbor;
	double gradient_col, gradient_row;
	//go through all pixels
	for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
	{
		for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
		{

			//SHIFT_CIELAB = 0 = Luminance channel
			//SHIFT_CIELAB = 1 = A channel 
			//SHIFT_CIELAB = 2 = B channel 
			for (unsigned int SHIFT_CIELAB = 0; SHIFT_CIELAB < CIELAB_NUM_CHANNELS; ++SHIFT_CIELAB)
			{

				//gradient for current cell
				gradient_cur = *(image + ((col + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_CIELAB);


				//x axis

				//if this is last column or last row
				if ((col == IMAGE_WIDTH-1) || (row == IMAGE_HEIGHT-1))
					gradient_col = 0;

				//it is not:
				else
				{
					//neigbor cell value on x axis
					gradient_col_neighbor = *(image + ((col + 1 + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_CIELAB);

					//value X difference for [x,y] - [x,y+1]
					gradient_col = gradient_cur - gradient_col_neighbor;

				}

				//y axis
				//if we are on last row
				if (row == IMAGE_HEIGHT-1)
					gradient_row = gradient_cur;

			
				//it is not last row, last col does not matter for y axis
				else
				{
					//neigbor cell value on y axis
					gradient_row_neighbor = *(image + ((col + IMAGE_WIDTH + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_CIELAB);

					//value y difference for [x,y] - [x+1,y]
					gradient_row = gradient_cur - gradient_row_neighbor;
				}	
				
				gradientPlane[row][col + IMAGE_WIDTH*SHIFT_CIELAB] = sqrt(pow(gradient_row, 2) + pow(gradient_col, 2));
			
			}	

		}
		
	}

	return;
}


void computeChromaticGradient(double **chromaticGradient, double **gradientPlane)
{
	double sum_A_vector, sum_B_vector, gradient_diff, weight_function;
	short int sign = 1;

	for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
	{		
		for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
		{

			//SHIFT_CIELAB = 0 = Luminance channel
			//SHIFT_CIELAB = 1 = A channel 
			//SHIFT_CIELAB = 2 = B channel 
			for (unsigned int SHIFT_CIELAB = 0; SHIFT_CIELAB < CIELAB_NUM_CHANNELS; ++SHIFT_CIELAB)
			{

				//border cell. Skip it
				if ((row == 0) || (row == IMAGE_HEIGHT-1) || (col == 0) || (col == IMAGE_WIDTH-1))
					continue;

				sum_A_vector = sum_B_vector = 0.0;

				//we have here Np which neighbor pixel pair set.
				//Combinations shown as ascii images below 
				/*
				|   |   |   |
				| x |   | x |
				|   |   |   |
				*/

				gradient_diff = gradientPlane[row][(col - 1) + IMAGE_WIDTH * SHIFT_CIELAB] - gradientPlane[row][(col + 1) + IMAGE_WIDTH * SHIFT_CIELAB];
				
				//if we are working with luminance channel currently
				if (SHIFT_CIELAB == 0)
				{
					//1 + x^2; for us 1 + luminance_difference^2
					weight_function = 1 + pow(gradient_diff, 2);

					//if gradient difference for luminance channel is smaller than zero	
					if (gradient_diff < 0)
							sign = -1;
						
				}

				//A channel or B channel
				(SHIFT_CIELAB == 1) ? sum_A_vector += sign * weight_function * gradient_diff : sum_B_vector += sign * weight_function * gradient_diff;
					
				sign = 1;
			}


			for (unsigned int SHIFT_CIELAB = 0; SHIFT_CIELAB < CIELAB_NUM_CHANNELS; ++SHIFT_CIELAB)
			{

				//border cell. Skip it
				if ((row == 0) || (row == IMAGE_HEIGHT-1) || (col == 0) || (col == IMAGE_WIDTH-1))
					continue;

				sum_A_vector = sum_B_vector = 0.0;

				//we have here Np which neighbor pixel pair set.
				//Combinations shown as ascii images below 
				/*
				|   | x |   |
				|   |   |   |
				|   | x |   |
				*/

				gradient_diff = gradientPlane[row - 1][col + IMAGE_WIDTH * SHIFT_CIELAB] - gradientPlane[row + 1][col + IMAGE_WIDTH * SHIFT_CIELAB];
				
				//if we are working with luminance channel currently
				if (SHIFT_CIELAB == 0)
				{
					//1 + x^2; for us 1 + luminance_difference^2
					weight_function = 1 + pow(gradient_diff, 2);

					//if gradient difference for luminance channel is smaller than zero	
					if (gradient_diff < 0)
							sign = -1;
						
				}

				//A channel or B channel
				(SHIFT_CIELAB == 1) ? sum_A_vector += sign * weight_function * gradient_diff : sum_B_vector += sign * weight_function * gradient_diff;
					
				sign = 1;
			}


			for (unsigned int SHIFT_CIELAB = 0; SHIFT_CIELAB < CIELAB_NUM_CHANNELS; ++SHIFT_CIELAB)
			{

				//border cell. Skip it
				if ((row == 0) || (row == IMAGE_HEIGHT-1) || (col == 0) || (col == IMAGE_WIDTH-1))
					continue;

				sum_A_vector = sum_B_vector = 0.0;

				//we have here Np which neighbor pixel pair set.
				//Combinations shown as ascii images below 
				/*
				| x |   |   |
				|   |   |   |
				|   |   | x |
				*/

				gradient_diff = gradientPlane[row - 1][(col - 1) + IMAGE_WIDTH * SHIFT_CIELAB] - gradientPlane[row + 1][(col + 1) + IMAGE_WIDTH * SHIFT_CIELAB];
				
				//if we are working with luminance channel currently
				if (SHIFT_CIELAB == 0)
				{
					//1 + x^2; for us 1 + luminance_difference^2
					weight_function = 1 + pow(gradient_diff, 2);

					//if gradient difference for luminance channel is smaller than zero	
					if (gradient_diff < 0)
							sign = -1;
						
				}

				//A channel or B channel
				(SHIFT_CIELAB == 1) ? sum_A_vector += sign * weight_function * gradient_diff : sum_B_vector += sign * weight_function * gradient_diff;
					
				sign = 1;
			}

			for (unsigned int SHIFT_CIELAB = 0; SHIFT_CIELAB < CIELAB_NUM_CHANNELS; ++SHIFT_CIELAB)
			{

				//border cell. Skip it
				if ((row == 0) || (row == IMAGE_HEIGHT-1) || (col == 0) || (col == IMAGE_WIDTH-1))
					continue;

				sum_A_vector = sum_B_vector = 0.0;

				//we have here Np which neighbor pixel pair set.
				//Combinations shown as ascii images below 
				/*
				|   |   | x |
				|   |   |   |
				| x |   |   |
				*/

				gradient_diff = gradientPlane[row - 1][(col + 1) + IMAGE_WIDTH * SHIFT_CIELAB] - gradientPlane[row + 1][(col - 1) + IMAGE_WIDTH * SHIFT_CIELAB];
				
				//if we are working with luminance channel currently
				if (SHIFT_CIELAB == 0)
				{
					//1 + x^2; for us 1 + luminance_difference^2
					weight_function = 1 + pow(gradient_diff, 2);

					//if gradient difference for luminance channel is smaller than zero	
					if (gradient_diff < 0)
							sign = -1;
						
				}

				//A channel or B channel
				(SHIFT_CIELAB == 1) ? sum_A_vector += sign * weight_function * gradient_diff : sum_B_vector += sign * weight_function * gradient_diff;
					
				sign = 1;
			}

			
			chromaticGradient[row][col] = sum_A_vector;
			chromaticGradient[row][col + IMAGE_WIDTH] = sum_B_vector;
				
		}
		
	}

	return;
}


void chromaticGlobalOrientation(double **chromaticGradient, double *O_a, double *O_b)
{
	double sum_A_vector, sum_B_vector = 0.0;
	//going through all pixels
	for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
	{
		for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
		{
			sum_A_vector += chromaticGradient[row][col];
			sum_B_vector += chromaticGradient[row][col + IMAGE_WIDTH];
		}
	}

	*O_a = sum_A_vector / (std::abs(sum_A_vector));
	*O_b = sum_B_vector / (std::abs(sum_B_vector));
}

void computeGaussianFilter(double **gaussianFilter, unsigned int sigma)
{
	//counts distance from border to center pixel, we know that GAUSSIAN_FILTER_SIZE is odd number
	unsigned int distanceFromCenter = (GAUSSIAN_FILTER_SIZE - 1) / 2;
	double soucet = 0.0;
	int shifted_row, shifted_col = 0;
	
	for (unsigned int row = 0; row < GAUSSIAN_FILTER_SIZE; ++row)
	{
		for (unsigned int col = 0; col < GAUSSIAN_FILTER_SIZE; ++col)
		{
			//because gaussian mean (peak) is equal to (0,0)
			//indexing: -2 -1 0 1 2 for 5x5
			shifted_row = row - distanceFromCenter;
			shifted_col = col - distanceFromCenter;
 
			//Gaussian blur equation:
			// (1 / (2 * pi * sigma^2)) * e^-((x^2 + y^2) / (2 * sigma^2))
			gaussianFilter[row][col] = (1 / (2 * M_PI * pow(sigma, 2)) ) * pow(EULER, -1*( (pow(shifted_row, 2) + pow(shifted_col, 2)) / (2 * pow(sigma, 2)) ));
			myfile << "[" << row << " | " << col << "] "<< gaussianFilter[row][col] << std::endl;

			soucet += gaussianFilter[row][col];
		}
	}
	myfile << "soucet:" << soucet << std::endl;
}

void differenceOfGaussians(double **firstGaussianFilter, double **secondGaussianFilter, double **resultGaussianFilter)
{
	for (unsigned int row = 0; row < GAUSSIAN_FILTER_SIZE; ++row)
	{
		for (unsigned int col = 0; col < GAUSSIAN_FILTER_SIZE; ++col)
			resultGaussianFilter[row][col] = firstGaussianFilter[row][col] - secondGaussianFilter[row][col];
	}
}


void convolutionWithGaussianFilter(double *sourceImage, double *destinationImage, double **gaussianFilter)
{
	//counts distance from border to center pixel, we know that GAUSSIAN_FILTER_SIZE is odd number
	unsigned int distanceFromCenter = (GAUSSIAN_FILTER_SIZE - 1) / 2;
	int shifted_gaussianFilter_row, shifted_gaussianFilter_col;

	//double *destinationImage_P_backup;
	//destinationImage_P_backup = destinationImage;

	//going through all pixels
	for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
	{
		for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
		{
			myfile << "distanceFromCenter" << distanceFromCenter << std::endl;
			//we are on border, skip it, no mercy. It makes no difference on big images
			if ( (row < distanceFromCenter) || (row >= (IMAGE_HEIGHT - distanceFromCenter))
				|| (col < distanceFromCenter) || (col >= (IMAGE_WIDTH - distanceFromCenter)) )
				continue;

			//concatenation will be for all CIELAB channels
			for (unsigned int SHIFT_CIELAB = 0; SHIFT_CIELAB < CIELAB_NUM_CHANNELS; ++SHIFT_CIELAB)
			{
				//go through whole Gaussasian kernel and concatenate it with sourceImage values
				for (unsigned int gaussianFilter_row = 0; gaussianFilter_row < GAUSSIAN_FILTER_SIZE; ++gaussianFilter_row)
				{
					for (unsigned int gaussianFilter_col = 0; gaussianFilter_col < GAUSSIAN_FILTER_SIZE; ++gaussianFilter_col)
					{
						//shift indexes even to negative values
						shifted_gaussianFilter_row = gaussianFilter_row - distanceFromCenter;
						shifted_gaussianFilter_col = gaussianFilter_col - distanceFromCenter;

						//destination image pixel on row, col
						*(destinationImage + ((col + IMAGE_WIDTH  * row ) * CIELAB_NUM_CHANNELS) + SHIFT_CIELAB) +=
						//source image pixel on row, col or shifted by filter
						*(sourceImage + ((col + shifted_gaussianFilter_col + IMAGE_WIDTH  * (row + shifted_gaussianFilter_row))
						* CIELAB_NUM_CHANNELS) + SHIFT_CIELAB) 
						//multiplied by gaussianFilter
						* gaussianFilter[gaussianFilter_row][gaussianFilter_col];
					
					}

				}

			}
		}
	}

//destinationImage = destinationImage_P_backup;

}

void computeDirectionalDistance(double *destinationImage, double O_a, double O_b)
{
	unsigned short int SHIFT_TO_L = 0;
	unsigned short int SHIFT_TO_A = 1;
	unsigned short int SHIFT_TO_B = 2;
	double result;
	//going through all pixels
	for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
	{
		for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
		{
			//actual pixel on row and col
			result = (*(destinationImage + SHIFT_TO_L)) + (*(destinationImage + SHIFT_TO_A)) * O_a + (*(destinationImage + SHIFT_TO_B)) * O_b;

			*(destinationImage++) = result;
			*(destinationImage++) = result;
			*(destinationImage++) = result;
		}
	}

}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOZhongping15::Transform()
{

	myfile.open ("example.txt");
  


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

	unsigned int i;

	//represents double laplacianPlane[IMAGE_HEIGHT][IMAGE_WIDTH*CIELAB_NUM_CHANNELS];
	double **gradientPlane = new double*[IMAGE_HEIGHT];
	for(i = 0; i < IMAGE_HEIGHT; ++i)
		//*3 because we have 3 channels, L; A; B
		gradientPlane[i] = new double[IMAGE_WIDTH*CIELAB_NUM_CHANNELS];

	computeGradientPlane(sourceImage, gradientPlane);

	//represents double chromaticGradient[IMAGE_HEIGHT][IMAGE_WIDTH*2];
	double **chromaticGradient = new double*[IMAGE_HEIGHT];
	for(i = 0; i < IMAGE_HEIGHT; ++i)
		//*2 because we have [a b] vector 
		chromaticGradient[i] = new double[IMAGE_WIDTH*2];

	computeChromaticGradient(chromaticGradient, gradientPlane);


	double O_a, O_b;
	chromaticGlobalOrientation(chromaticGradient, &O_a, &O_b);

	//represents double gaussianFilter[GAUSSIAN_FILTER_SIZE][GAUSSIAN_FILTER_SIZE];
	//this filter size should be odd, it will usually be with size 5x5
	//if (GAUSSIAN_FILTER_SIZE % 2 == 0)
		//error

	double **gaussianFilter = new double*[GAUSSIAN_FILTER_SIZE];
	for(i = 0; i < GAUSSIAN_FILTER_SIZE; ++i)
		gaussianFilter[i] = new double[GAUSSIAN_FILTER_SIZE];

	unsigned int sigma = 1;
	computeGaussianFilter(gaussianFilter, sigma);


	convolutionWithGaussianFilter(sourceImage, destinationImage, gaussianFilter);


	//computeDirectionalDistance(destinationImage, O_a, O_b);

	double sourceImage_L, sourceImage_A, sourceImage_B;
	//going through all pixels,
	/*
	for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
	{
		pSrc->ProgressBar(row, pSrc->GetHeight());
		for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
		{
			//taking all LAB values in CIELab for one pixel and computing new ones to destination image
			sourceImage_L = *(sourceImage++);
			sourceImage_A = *(sourceImage++);
			sourceImage_B = *(sourceImage++);

			
			//actualization of destination image for current pixel
			*(destinationImage++) = 0;
			*(destinationImage++) = chromaticGradient[row][col];
			*(destinationImage++) = chromaticGradient[row][col + IMAGE_WIDTH];
		}
	}

*/
  myfile.close();


	return 0;
}