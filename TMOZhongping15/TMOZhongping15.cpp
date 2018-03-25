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

#define CIELAB_NUM_CHANNELS 3

unsigned int IMAGE_WIDTH; //global variable for width of loaded image
unsigned int IMAGE_HEIGHT; //global variable for height of loaded image


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


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOZhongping15::Transform()
{

	// Source image is stored in local parameter pSrc
	pSrc->Convert(TMO_LAB);								// This is format of CIELab
	// Destination image is in pDst
	pDst->Convert(TMO_LAB);								// This is format of CIELab

	double* sourceImage = pSrc->GetData();				// You can work at low level data
	double* destinationImage = pDst->GetData();			// Data are stored in form of array 
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


	double sourceImage_L, sourceImage_A, sourceImage_B;
	//going through all pixels,
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
			*(destinationImage++) = 1;
			*(destinationImage++) = chromaticGradient[row][col];
			*(destinationImage++) = chromaticGradient[row][col + IMAGE_WIDTH];
		}
	}


	return 0;
}