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
 * TMOXiong17.cpp: implementation of the TMOXiong17 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOXiong17.h"
#include <stdlib.h> 
#include <math.h> //power
#include <iostream>
#include <string>

#include <Eigen/Dense>

using namespace Eigen;

#define RED_CHOSEN 0
#define GREEN_CHOSEN 1
#define BLUE_CHOSEN 2
#define COLOR_NOT_CHOSEN 3

unsigned int IMAGE_WIDTH; //global variable for width of loaded image
unsigned int IMAGE_HEIGHT; //global variable for height of loaded image
unsigned int SIZE_OF_Z2 = 1; //global variable for size of Z2 polynomial space, dependent on order

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOXiong17::TMOXiong17()
{
	SetName(L"Xiong17");						
	SetDescription(L"Parametric ratio-based method for efficient contrast-preserving decolorization");	

	//input optional parameters 

	//sets how many time we will recompute weight coeficients; later set from command line
	maximum_of_iterations.SetName(L"MaxIter");
	maximum_of_iterations.SetDescription(L"Maximum of iterations, iterative method ; <1,200>");	
	maximum_of_iterations.SetDefault(50);
	maximum_of_iterations=50;		
	maximum_of_iterations.SetRange(1,200);
	this->Register(maximum_of_iterations);

	order.SetName(L"Order");
	order.SetDescription(L"Order of method, Z2 size is dependent on it ; <1,3>");	
	order.SetDefault(2);
	order=2;		
	order.SetRange(1,3);	
	this->Register(order);

}

TMOXiong17::~TMOXiong17()
{
}


//function that computes gradient for our equations
//it takes color intensity from neighbor pixels and compares it with given X Y pixel
void singleChannelGradient(double *image, double *gradient, int GRADIENT_COLOR_OPTION)
{
	unsigned int image_cell_shift, image_color_shift;

	//color was not chosen, we have image with only one channel
	if (GRADIENT_COLOR_OPTION == COLOR_NOT_CHOSEN)
	{
		image_cell_shift = 1;
		image_color_shift = 0; 
	}
	
	//multichannel image, RGB expected
	else
	{
		image_cell_shift = 3;
		image_color_shift = GRADIENT_COLOR_OPTION;
	}

	double gradient_cur, gradient_col_neighbor, gradient_row_neighbor;
	double gradient_col, gradient_row;

	//go through all pixels
	for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
	{
		for(unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
		{
		
			//gradient for current cell
			gradient_cur = *(image + ((col + IMAGE_WIDTH  * row ) * image_cell_shift) + image_color_shift);

			//x axis

			//if this is last column or last row
			if ((col == IMAGE_WIDTH-1) || (row == IMAGE_HEIGHT-1))
				gradient[col * IMAGE_HEIGHT + row] = 0;
		
			//it is not:
			else
			{
				//neigbor cell value on x axis
				gradient_col_neighbor = *(image + ((col + 1 + IMAGE_WIDTH  * row ) * image_cell_shift) + image_color_shift);

				//value X difference for [x,y] - [x,y+1]
				gradient_col = gradient_cur - gradient_col_neighbor;


				gradient[col * IMAGE_HEIGHT + row] = gradient_col;

			}	
		
			//y axis
			//if we are on last row
			if (row == IMAGE_HEIGHT-1)
				gradient[IMAGE_HEIGHT*IMAGE_WIDTH + col * IMAGE_HEIGHT + row] = gradient_cur;

		
			//it is not last row, last col does not matter for y axis
			else
			{
				//neigbor cell value on y axis
				gradient_row_neighbor = *(image + ((col + IMAGE_WIDTH + IMAGE_WIDTH  * row ) * image_cell_shift) + image_color_shift);

				//value y difference for [x,y] - [x+1,y]
				gradient_row = gradient_cur - gradient_row_neighbor;
				
				//vector with row gradients will be added behind all col gradients
				gradient[IMAGE_HEIGHT*IMAGE_WIDTH + col * IMAGE_HEIGHT + row] = gradient_row;
			}

		}
	}
				
}

//function that computes polygrad = array of gradients for all Z2 channel combinations
//						 combination = stores channel combinations, represents Z2 
void gradSystem(double **polygrad, double *image, int order, int combination[][3])
{
	//allocation of arrays with new operation, needed for big resolution images
	double *curIm = new double [IMAGE_HEIGHT*IMAGE_WIDTH];
	double *curGrad = new double [IMAGE_HEIGHT*IMAGE_WIDTH*order];

	int iteration_Number = 0;

	for (int r = 0; r <= order; ++r)
	{
		for (int g = 0; g <= order; ++g)
		{
			for (int b = 0; b <= order; ++b)
			{
				//it takes exactly what we expect in Z2.
				//If exponent for color is 2, other colors have zero etc.
				if ( ((r + g + b) <= order) &&  ((r + g + b) > 0) )
				{
					//we will later need it. Saves useful combinations represented by Z2 linear combination of colors
					combination[iteration_Number][0] = r;
					combination[iteration_Number][1] = g;
					combination[iteration_Number][2] = b;
					
					//going through all pixels
					for (unsigned int row = 0; row < IMAGE_HEIGHT; ++row)
					{
						for (unsigned int col = 0; col < IMAGE_WIDTH; ++col)
						{			
							//going through image with Z2 editation; Ib, Ib*Ig, Ir*
							curIm[col + IMAGE_WIDTH  * row] = pow ( (*(image + (col + IMAGE_WIDTH  * row) * 3 + RED_CHOSEN) ) , r ) 
							* pow ( (*(image + (col + IMAGE_WIDTH  * row) * 3 + GREEN_CHOSEN) ), g )
							* pow ( (*(image + (col + IMAGE_WIDTH  * row) * 3 + BLUE_CHOSEN) ), b );
						}
					}

					//compute gradient for image 
					singleChannelGradient(curIm, curGrad, COLOR_NOT_CHOSEN);

					//save all gradients to polygrad
					for (unsigned int i = 0; i < IMAGE_HEIGHT*IMAGE_WIDTH*order; ++i)
						//copy array of gradient values to array that stores it for all iterations
						polygrad[i][iteration_Number] = curGrad[i];
						
					
					iteration_Number++;
				}
			}	
		}
	}

	//dealocation of storage for arrays
	delete [] curIm;
	delete [] curGrad;
}


void colorGradient(double *image, double *colorGradient, double *max_val)
{
	//allocation of arrays with new operation, needed for big resolution images
	double *ImL = new double [IMAGE_WIDTH*IMAGE_HEIGHT*2];
	double *ImA = new double [IMAGE_WIDTH*IMAGE_HEIGHT*2];
	double *ImB = new double [IMAGE_WIDTH*IMAGE_HEIGHT*2];

	double result;

	singleChannelGradient(image, ImL, RED_CHOSEN);
	singleChannelGradient(image, ImA, GREEN_CHOSEN);
	singleChannelGradient(image, ImB, BLUE_CHOSEN);


	for (unsigned int i = 0; i < IMAGE_HEIGHT*IMAGE_WIDTH*2; ++i)
	{
		//
		result = (sqrt(pow(ImL[i], 2) + pow(ImA[i], 2) + pow(ImB[i], 2))) / 100.0;

		//looking for biggest value in Color gradient delta_xy, will be used for constant K for S_xy computation
		if (*max_val < result)
			*max_val = result;
		
		//saving result of operation for current cell
		colorGradient[i] = result;
	}

	//dealocation of storage for arrays
	delete [] ImL;
	delete [] ImA;
	delete [] ImB;
}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOXiong17::Transform()
{
	//source image is stored in local parameter pSrc
	//destination image is in pDst
	//initialy images are in RGB format, but you can convert it into other format
	pSrc->Convert(TMO_RGB);								// This is format of RGB
	pDst->Convert(TMO_RGB);								// Ir (x), Ig (y), Ib (z) as color information


	//saving values of image width and height to global variables
	IMAGE_WIDTH = pSrc->GetWidth();
	IMAGE_HEIGHT = pSrc->GetHeight();

	//size of Z2 polynomial space
	//2^order
	if (order == 1)
		SIZE_OF_Z2 = 3;

	else if (order == 2)
		SIZE_OF_Z2 = 9;

	else if (order == 3)
		SIZE_OF_Z2 = 27;

	double* sourceImage = pSrc->GetData();				// You can work at low level data
	double* destinationImage = pDst->GetData();			// Data are stored in form of array 
		
	//avoids numerical instability
	//Gamma and Epsilon values based on dataset research. by Doc. Ing. Martin Čadík, Ph.D. Looking for best CCPR
	//where CCPR stands for Color Contrast Preserving Ratio
	double epsilon = 0.015;
	unsigned int gamma = 50;

	unsigned int j;
	unsigned int i;


	//###### POLYNOMIAL INITIALIZATION #######
	
	//represents double polygrad[IMAGE_WIDTH*IMAGE_HEIGHT*order][SIZE_OF_Z2];
	double **polygrad =	new double*[IMAGE_WIDTH*IMAGE_HEIGHT*order];
	for(i = 0; i < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++i)
		polygrad[i] = new double[SIZE_OF_Z2];


	int combination[SIZE_OF_Z2][3];

	gradSystem (polygrad, sourceImage, order, combination);

	//weight coefficients; they will concatenate to what we are searching for
	
	double *W_l = new double[SIZE_OF_Z2];
	//starting values, later on it will change
	for (i = 0; i < SIZE_OF_Z2; ++i)
	{
		
	//RESULT for common order 2:	
	//			w_b  w_b^2  w_g  w_gb w_g^2 w_r w_rb w_rg w_r^2
	//W_l[9] = {0.33, 0.0, 0.33, 0.0, 0.0, 0.33, 0.0, 0.0, 0.0};

		//if we have basic channel R or G or B, not their combination
		if ((combination[i][0] + combination[i][1] + combination[i][2]) == 1)
			W_l[i] = 0.33;
		else
			W_l[i] = 0;
	}


	double *sourceImageBW = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order];
	double *sourceImage_P_backup = sourceImage; 
	double *sourceImageBW_P_backup = sourceImageBW;

	double pixel_r, pixel_g, pixel_b;

	//go through all pixels and count grayscale image will be used for starting values for S_xy
	//then values are updated for bettter results
	for (j = 0; j < IMAGE_HEIGHT; ++j)
	{
		for (i = 0; i < IMAGE_WIDTH; ++i)
		{
			pixel_r = *(sourceImage++);
			pixel_g = *(sourceImage++);
			pixel_b = *(sourceImage++);

			//making image gray with common easy method
			sourceImageBW[i + j*IMAGE_WIDTH] = 0.299*pixel_r+0.587*pixel_g+0.114*pixel_b;
			
		}
	}

	//returning starting pointer value
	sourceImage = sourceImage_P_backup;
	sourceImageBW = sourceImageBW_P_backup;

	//gradient for gray image
	double *S_xy = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order];

	singleChannelGradient(sourceImageBW, S_xy, COLOR_NOT_CHOSEN);

	//Langrange multiplier
	double *Langrange_multiplier = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order];

	//fill array Langrange_multiplier with zeroes - starting value
	std::fill_n(Langrange_multiplier, IMAGE_WIDTH*IMAGE_HEIGHT*order, 0);

	//color gradient for image
	double *delta_xy = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order];

	//transfer color model to CIELab
	pSrc->Convert(TMO_LAB);
	//LAB <> L for lightness and a and b for the color opponents green–red and blue–yellow

	double *sourceImageCIELab = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order*3];
	//save starting pointer
	double *sourceImageCIELab_P_backup = sourceImageCIELab;

	//we have to copy content, cannot use just pSrc->GetData()
	// because pSrc is going to be converted to RGB again, and it is just pointer
	for (j = 0; j < IMAGE_HEIGHT; ++j)
	{
		for (i = 0; i < IMAGE_WIDTH; ++i)
		{
			//L
			sourceImageCIELab[(i + j*IMAGE_WIDTH)*3] = *(sourceImage++);
			//A
			sourceImageCIELab[(i + j*IMAGE_WIDTH)*3 + 1] = *(sourceImage++);
			//B
			sourceImageCIELab[(i + j*IMAGE_WIDTH)*3 + 2] = *(sourceImage++);
		}
	}
	
	sourceImageCIELab = sourceImageCIELab_P_backup;
	sourceImage = sourceImage_P_backup;
	//convert it back
	pSrc->Convert(TMO_RGB);


	double max_val;
	colorGradient(sourceImageCIELab, delta_xy, &max_val);

	// (1 / max) value from color gradient delta_xy
	double k = 1 / max_val;

	double row_sum = 0.0;

	double *top_part_fraction = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order];
	double *t = new double [IMAGE_WIDTH*IMAGE_HEIGHT*order];

	//B is for transformed polygrad matrix
	//represents double B[SIZE_OF_Z2][IMAGE_WIDTH*IMAGE_HEIGHT*order];
	double **B = new double*[SIZE_OF_Z2];
	for(i = 0; i < SIZE_OF_Z2; ++i)
		B[i] = new double[IMAGE_WIDTH*IMAGE_HEIGHT*order];

	VectorXf sum_B(SIZE_OF_Z2);
	VectorXf Mt(SIZE_OF_Z2);
	MatrixXf polygrad_polygradT(SIZE_OF_Z2, SIZE_OF_Z2);
	
	double col_sum = 0.0;

	//Below is main cycle of events; all iterations:
	/*	Initialization: w0 , s0 0, 0 0;
		3: For k 0 to K 1 do
		4: update wk 1 according to Eq. (8)
		5: update
		1
		,
		k
		x y s according to Eq. (12)
		6: update k 1
		according to Eq. (7)
		7: End (For)
	*/

	for (unsigned int iterationCount = 0; iterationCount < maximum_of_iterations; ++iterationCount)
	{	
		
		//progress bar shows how far we are in iterations
		pSrc->ProgressBar(iterationCount, maximum_of_iterations);

	
		//computation of help variable S_xy
		//equation 12 
		//this part is top part of fraction:
		for (unsigned int i = 0; i < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++i)
		{
			row_sum = 0.0;
			//Matrix multiplied by vector. Ax9 * 9x1 = Ax1 size. Multiplied colums by weight are summed to one.
			for (unsigned int w_array_index = 0; w_array_index < SIZE_OF_Z2; ++w_array_index)
				row_sum += polygrad[i][w_array_index] * W_l[w_array_index];


			top_part_fraction[i] = gamma * row_sum + Langrange_multiplier[i];
		}


		//why four times? Algorithm authors figured out it has better solution than (2-5 times)
		for (unsigned int s_counter = 0; s_counter < 4; ++s_counter)
		{
			for (unsigned int i = 0; i < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++i)
			{
				//bottom part of fraction:
				
				//epsilon is small parameter to prevent numerical instability
				t[i] = 1 / (std::abs(S_xy[i]) + epsilon);

				//new computation for our S_xy for this iteration. Needed for every weight recompute
				S_xy[i] = top_part_fraction[i] / (2 + gamma - 2*k*t[i]*delta_xy[i]);

			} 

		}
	

		//computation of weight vector
		//equation 8

		for (unsigned int column = 0; column < SIZE_OF_Z2; ++column)
		{
			//radeji zkonrolovat
			for (unsigned int row = 0; row < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++row)	
				//B has size of transformed polygrad
				B[column][row] = polygrad[row][column]*(S_xy[row] - Langrange_multiplier[row] / gamma);
			
		}


		//Sum B contains sums of rows. So result has only one row and same number of columns
		for (unsigned int column = 0; column < SIZE_OF_Z2; ++column)
		{
			col_sum = 0.0;

			for (unsigned int row = 0; row < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++row)
				col_sum += B[column][row];

			sum_B [column] = col_sum;

		}

		//poly'*poly not doing Transposition, Faking it, so this is the way to go
		//multiplying cols with cols
		for (unsigned int col_polygrad_trans = 0; col_polygrad_trans < SIZE_OF_Z2; ++col_polygrad_trans)
		{
			for (unsigned int col_polygrad = 0; col_polygrad < SIZE_OF_Z2; ++col_polygrad)
			{
				col_sum = 0.0;

				for (unsigned int row_polygrad = 0; row_polygrad < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++row_polygrad)
					col_sum += polygrad[row_polygrad][col_polygrad_trans]*polygrad[row_polygrad][col_polygrad];

				polygrad_polygradT (col_polygrad_trans, col_polygrad) = col_sum;

			}
		}

		//Basic linear solving
		//You have a system of equations, that you have written as a single matrix equation:
		//polygrad_polygradT * Mt = sum_B
		//Where polygrad_polygradT and sum_B are matrices (sum_B could be a vector, as a special case and it is).
		//You want to find a solution Mt.
		Mt = polygrad_polygradT.colPivHouseholderQr().solve(sum_B);

		//computation of Langrange multiplier
		//equation 7

		for (unsigned int i = 0; i < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++i)
		{
			row_sum = 0.0;
			//Matrix multiplied by vector. Ax9 * 9x1 = Ax1 size. Multiplied colums by weight are summed to one.
			for (unsigned int w_array_index = 0; w_array_index < SIZE_OF_Z2; ++w_array_index)
				row_sum += polygrad[i][w_array_index] * W_l[w_array_index];
			

			Langrange_multiplier[i] += gamma*(row_sum - S_xy[i]);
		}

	}

	//delete rrays now, we won't need it anymore
	for(int counter = 0; counter < IMAGE_WIDTH*IMAGE_HEIGHT*order; ++counter)
	{
	    delete [] polygrad[counter];
	}
	delete [] polygrad;

	//array delete
	for(int counter = 0; counter < SIZE_OF_Z2; ++counter)
	{
	    delete [] B[counter];
	}
	delete [] B;

	delete [] top_part_fraction;
	delete [] t;
	delete [] S_xy;
	delete [] Langrange_multiplier;
	delete [] delta_xy;

	//######## color weights computed - final version after all iterations ########
	unsigned int w_array_index = 0;
	double sourceImage_r, sourceImage_g, sourceImage_b;
	double result;
	double *destinationImage_P_backup;

	//saving pointer to starting position in destination image, for later easier manipulation with shifting
	destinationImage_P_backup = destinationImage;

	for (unsigned int r = 0; r <= order; ++r)
	{
		for (unsigned int g = 0; g <= order; ++g)
		{
			for (unsigned int b = 0; b <= order; ++b)
			{
				//it takes exactly what we expect in Z2.
				//If exponent for color is 2, other colors have zero etc.
				if ( ((r + g + b) <= order) &&  ((r + g + b) > 0) )
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
							result = W_l[w_array_index] * pow(sourceImage_r, r) * pow(sourceImage_g, g) * pow(sourceImage_b, b);

							//actualization of destination image for current pixel
							*(destinationImage++) += result;
							*(destinationImage++) += result;
							*(destinationImage++) += result;

						}
					}

					//color recomputed for all pixels, lets take another color combination from Z2
					w_array_index++;
						
				}
			}	
		}
	}

	//delete final weight array
	delete [] W_l;


	destinationImage = destinationImage_P_backup;
	//stores max and min value for any canal in that image
	double minValDestImage = 255.0;
	double maxValDestImage = 0.0;
	//we will search for max and min value, it will be needed for histogram cut
	for (unsigned int y = 0; y < IMAGE_HEIGHT; ++y)
	{
		pSrc->ProgressBar(y, pSrc->GetHeight());
		for (unsigned int x = 0; x < IMAGE_WIDTH; ++x)
		{
							
			if ((*(destinationImage)) < minValDestImage)
				minValDestImage = (*(destinationImage));

			if ((*(destinationImage)) > maxValDestImage)
				maxValDestImage = (*(destinationImage));

			destinationImage++;
		}
	}



	//change of image scale histogram for better contrast
	//(Gray - minValue) / (maxValue - minValue)
	//cutting useless histogram edges
	//going through all pixels, first x - cols then y - rows
	destinationImage = destinationImage_P_backup;
	for (unsigned int y = 0; y < IMAGE_HEIGHT; ++y)
	{
		pSrc->ProgressBar(y, pSrc->GetHeight());
		for (unsigned int x = 0; x < IMAGE_WIDTH; ++x)
		{
			result = (*(destinationImage + 1) - minValDestImage) / (maxValDestImage - minValDestImage);
			*(destinationImage++) = result;
			*(destinationImage++) = result;
			*(destinationImage++) = result;
		}
	}


return 0;
}
