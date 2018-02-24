/* --------------------------------------------------------------------------- *
 * TMOXiong17.cpp: implementation of the TMOXiong17 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOXiong17.h"
#include <stdlib.h> 
#include <math.h>

#include <opencv2/imgproc/imgproc.hpp> //cvtColor
#include "optim.hpp" //ConjGradSolver


//avoids numerical instability
//Gamma and Epsilon values based on dataset research. by Doc. Ing. Martin Čadík, Ph.D. Looking for best CCPR
//where CCPR stands for Color Contrast Preserving Ratio
#define IDEAL_CCPR_EPSILON 4
#define IDEAL_CCPR_GAMMA 3
#define NOT_SET -1.0

#define RED_CHOSEN 0
#define GREEN_CHOSEN 1
#define BLUE_CHOSEN 2

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOXiong17::TMOXiong17()
{
	SetName(L"Xiong17");						
	SetDescription(L"Parametric ratio-based method for efficient contrast-preserving decolorization");	

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOXiong17::~TMOXiong17()
{
}


//function that computes gradient for our equations
//it takes color intensity from neighbor pixels and compares it with given X Y pixel
double gradientCompute(int x, int y, double I_c int GRADIENT_COLOR_OPTION)
{

	double gradient_x_part1, gradient_x_part2, gradient_y_part1, gradient_y_part2 = NOT_SET;
	double gradient_x, gradient_y, gradient;


	if (y > 0) //if this is not first row
		//funtion GetPixel returns pointer to area starting on  [x, y] 
		//GRADIENT_COLOR_OPTION is color offset for moving in array of colors
		gradient_y_part1 = abs(I_c - (*pSrc->GetPixel(x, y-1) + GRADIENT_COLOR_OPTION));

	if (y < pSrc->GetHeight()-1) //if this is not last row
		gradient_y_part2 = abs(I_c - (*pSrc->GetPixel(x, y+1) + GRADIENT_COLOR_OPTION));


	if (gradient_y_part1 == NOT_SET)
		gradient_y = gradient_y_part2;

	gradient_y_part2 != NOT_SET ? gradient_y = ((gradient_y_part1 + gradient_y_part1) / 2) : gradient_y = gradient_y_part1;


	//x axis
	if (x > 0) //if this is not first column
		//funtion GetPixel returns pointer to area starting on  [x, y] 
		gradient_x_part1 = abs(I_c - (*pSrc->GetPixel(x-1, y) + GRADIENT_COLOR_OPTION));

	if (x < pSrc->GetWidth()-1) //if this is not last column
		gradient_x_part2 = abs(I_c - (*pSrc->GetPixel(x+1, y) + GRADIENT_COLOR_OPTION));


	if (gradient_x_part1 == NOT_SET)
		gradient_x = gradient_x_part2;

	gradient_x_part2 != NOT_SET ? gradient_yx = ((gradient_x_part1 + gradient_x_part1) / 2) : gradient_x = gradient_x_part1;

	//gradient magnitude
	gradient = sqrt(pow(gradient_x,2.0) + pow(gradient_y,2.0));	
	
	gradient_x_part1 = gradient_x_part2 = gradient_y_part1 = gradient_y_part2 = NOT_SET;

	return gradient;
				
}

//function for 	computation of distance between Ix (first_pixel) and Iy (second_pixel) in CIELab model
//color contrast absolute variable; signed value
double deltaXY(int first_pixel_x, int first_pixel_y, int second_pixel_x, int second_y, double *pSourceDataCIELab)
{
	double first_pixel_L = (*pSourceDataCIELab->GetPixel(first_pixel_x, first_pixel_y)); //Lx
	double first_pixel_a = (*(pSourceDataCIELab->GetPixel(first_pixel_x, first_pixel_y)+1)); //ax
	double first_pixel_b = (*(pSourceDataCIELab->GetPixel(first_pixel_x, first_pixel_y)+2)); //bx

	double second_pixel_L = (*pSourceDataCIELab->GetPixel(second_pixel_x, second_y)); //Ly
	double second_pixel_a = (*(pSourceDataCIELab->GetPixel(second_pixel_x, second_y)+1)); //ay
	double second_pixel_b = (*(pSourceDataCIELab->GetPixel(second_pixel_x, second_y)+2)); //by

	//deltaXY = (sqrt(Lx-Ly))^2 + (ax - ay)^2 + (bx - by)^2 
	return ( pow(sqrt(first_pixel_L - second_pixel_L), 2.0) + pow(first_pixel_a - second_pixel_a, 2.0) + pow(first_pixel_b - second_pixel_b, 2.0) );
}



class findMinWeight : public cv::MinProblemSolver::Function
{
	double calc(nejake promenne!!)const
	{
		gradient = gradientCompute(x, y, I_c, proccessed_color);
		sum_of_gradient = W_l[proccessed_color] * gradient;

		//equation number 8; new compute of weight
		//W_k_next = W_k+1
		w_k_next = gamma[IDEAL_CCPR_GAMMA] / 2 ( sum_of_gradient - (Sxy_k - (lambda_k / gamma[IDEAL_CCPR_GAMMA]) ) );
	}
};

//Function for optimalization of weight compute
static void optimalizationFunction(cv::Ptr<cv::ConjGradSolver> optimalizator,cv::Ptr<cv::MinProblemSolver::Function> ptr_F, double Sxy_k, double lambda_k)
{
	//Getter for the optimized function
	optimalizator->setFunction(ptr_F);
	
	//actually runs the algorithm and performs the minimization.
	double res = optimalizator->minimize(x);

}



/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOXiong17::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_RGB);								// This is format of RGB
	pDst->Convert(TMO_RGB);								// Ir (x), Ig (y), Ib (z) as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components
	//value of red, green, blue, and general value
	double I_r, I_g, I_b, I_c;

	//avoids numerical instability
	//Gamma and Epsilon values based on dataset research. by Doc. Ing. Martin Čadík, Ph.D. Looking for best CCPR
	//where CCPR stands for Color Contrast Preserving Ratio
	double epsilon[8] = {0.00015, 0.0012, 0.0017, 0.001, 0.0015, 0.002, 0.015, 0.15};
	unsigned int gamma[8] = {20, 80, 85, 90, 100, 120, 150, 200};

	//weight coefficients in iterations of two initializations
	//starting values, later on it will change
				//w_r	w_g	  w_b	w_rg w_gb w_rb w_r^2 w_g^2 w_b^2
	double W_l = {0.33, 0.33, 0.33, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,}

	//for image weights 
	double w_k, w_k_next = 0;

	//Minimum of weight variable
	double argmin_w = w_k;

	//Auxiliary variable. Starting value at first then value of last iteration
	double Sxy_k = 0.0; 


	//Langrange multiplier constant. Starting value at first then value of last iteration
	double lambda_k = 0.0;
	
	int k = 1;


	double sum_of_gradient, gradient;



	double* pSourceDataCIELab = pSrc->GetData();

	//change values from 0..255 to 0..1
	pSourceDataCIELab *= 1./255;
	//transfer color model to CIELab
	cvtColor(pSourceDataCIELab, pSourceDataCIELab, COLOR_RGB2Lab);
	//LAB <> L for lightness and a and b for the color opponents green–red and blue–yellow


	//we set pointer to funtion for optimalization
	cv::Ptr<cv::MinProblemSolver::Function> ptr_F(new findMinWeight());

    unsigned int y=0;
	for (y = 0; y < pSrc->GetHeight(); ++y)
	{
		pSrc->ProgressBar(y, pSrc->GetHeight());	// You can provide progress bar
		for (unsigned int x = 0; x < pSrc->GetWidth(); ++x)
		{
			I_r = *pSourceData++;
			I_g = *pSourceData++;
			I_b = *pSourceData++;

			// Here you can use your transform 
			// expressions and techniques...
			//pY *= dParameter;							// Parameters can be used like
														// simple variables


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

			//where proccessed_color = 0 representes RED; = 1 is GREEN; = 2 is BLUE
			for (unsigned int proccessed_color = 0; proccessed_color < 2; ++proccessed_color)
			{
				for (unsigned int counter_for_K = 0; counter_for_K < K-1; ++counter_for_K)
				{
					//sum_of_gradient of w_l_k and gradient of I_l where I2 belongs to Z2 which is span family of monomials
					sum_of_gradient = 0.0;

					//proccessed_color == 0
					if (proccessed_color == RED_CHOSEN)
						I_c = I_r;
					else if (proccessed_color == GREEN_CHOSEN)
						I_c = I_g;
					//processed color is BLUE
					else
						I_c = I_b;


					//w_k_next = findMinWeight(x, y, Sxy_k, lambda_k);
					optimalizationFunction(optimalizator, ptr_F);

					//updating last iteration value for weight variable 
					w_k = w_k_next;

					//we change value of min founded weight if new one is less
					if (argmin_w > w_k_next)
						argmin_w = w_k_next;


					//equation number 12; new compute of auxiliary variable S_xy
					//updating last iteration value for auxiliary variable  S_xy

					//t = 1/ + epsilon[IDEAL_CCPR_EPSILON]; 
					//delta_xy = deltaXY(x1, y1, x2, y2, pSourceDataCIELab)
					S_xy_next = (gamma[IDEAL_CCPR_GAMMA] * sum_of_gradient + lambda_k) / (2 + gamma[IDEAL_CCPR_GAMMA] - 2 * K * t * delta_xy);


					Sxy_k = S_xy_next;


					//equation number 7; new compute of Langrange multiplier
					lambda_k_next = lambda_k + gamma[IDEAL_CCPR_GAMMA] * (sum_of_gradient - S_xy_next);

					//updating last iteration value for Langrange multiplier variable 
					lambda_k = lambda_k_next;

				}

			// and store results for a color to the destination image
			*pDestinationData++ *= W_k_next;
				
			}

		}
	}

	//pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0;

}

