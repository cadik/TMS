/* --------------------------------------------------------------------------- *
 * TMOEfficientC2GConversionForDIInGD.cpp: implementation of the TMOEfficientC2GConversionForDIInGD class.   *
 * 
 * Dependencies: FFTW library and boost::multi_array.
 * --------------------------------------------------------------------------- */

#include "TMOEfficientC2GConversionForDIInGD.h"

// matrix library
#include "../matrix.h"

// libraries for PES
#include <boost/multi_array.hpp>
#include "../poisson_pde/laplace.h"

// constant c is used to ensure that the largest chromatic difference will not be completely scaled down
#define CONSTANT_C 2.0

using namespace pde;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOEfficientC2GConversionForDIInGD::TMOEfficientC2GConversionForDIInGD()
{
	SetName(L"EfficientC2GConversionForDIInGD");
	SetDescription(L"Efficient Color-to-Gray Conversion for Digital Images in Gradinet Domain");

	alpha.SetName(L"alpha");
	alpha.SetDescription(L"Strength of the chromatic difference affecting the sign of the gradient");
	alpha.SetDefault(0.2);
	alpha=0.2;
	alpha.SetRange(0.0, 1.0);
	this->Register(alpha);
	
	beta.SetName(L"beta");
	beta.SetDescription(L"Key parameter to reduce grayscale distortion");
	beta.SetDefault(0.2);
	beta=0.2;
	beta.SetRange(0.0, std::numeric_limits<double>::max());
	this->Register(beta);

	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Key parameter to reduce grayscale distortion");
	gamma.SetDefault(0.5);
	gamma=0.5;
	gamma.SetRange(0.0, std::numeric_limits<double>::max());
	this->Register(gamma);

	theta.SetName(L"theta");
	theta.SetDescription(L"Direction in a-b plane");
	theta.SetDefault(0.0);
	theta=0.0;
	theta.SetRange(0.0, 360.0);
	this->Register(theta);	
	
	x_max = 0.0;
}

TMOEfficientC2GConversionForDIInGD::~TMOEfficientC2GConversionForDIInGD()
{
}

/**
 * attenuation function scals down the input signal
 * correspondents to part of Eq.(1) in papper
 * 
 * @param x - input value
 * @return scaled value
 */
double TMOEfficientC2GConversionForDIInGD::AttenuationFunction(double x){
	x_max = (x > x_max) ? x : x_max;		
	//std::cerr << "A in: " << x << "A out: " << x * (beta * (1 - pow(x / (CONSTANT_C * x_max), gamma))) << std::endl;	
	return x * (beta * (1 - pow(x / (CONSTANT_C * x_max), gamma)));	
}

/**
 * gets modulated difference of 2 colors
 * 
 * @param delta_l - L2 - L1
 * @param delta_a - a2 - a1
 * @param delta_b - b2 - b1
 * 
 * @return delta_E - modulated difference of 2 colors
 */
double TMOEfficientC2GConversionForDIInGD::ModulatedColorDifference(double delta_l, double delta_a, double delta_b){
	return sqrt(pow(delta_l, 2) + pow(AttenuationFunction(sqrt(pow(delta_a, 2) + pow(delta_b, 2))), 2));
}

/**
 * sign function for sing of gradinet for correct color ordering
 * correspondents to part of Eq.(4) in papper
 * 
 * @param delta_l - L2 - L1
 * @param a2 - component a of color 2 in LAB color space
 * @param a1 - component a of color 1 in LAB color space
 * @param b2 - component b of color 2 in LAB color space
 * @param b1 - component b of color 1 in LAB color space
 * 
 * @return sign of gradient
 */
double TMOEfficientC2GConversionForDIInGD::Sign(double delta_l, double a2, double a1, double b2, double b1){	
	// v_theta and delta_c are vectors (2 x 1) and (1 x 2)
	theta = TMOImage::DegreesToRadians(theta);		
	
	double value = delta_l + alpha * (cos(theta) * (a2 - a1) + sin(theta) * (b2 - b1));			
	
	if (value >= 0.0){		
		return 1.0;
	}else{		
		return -1.0;
	}
}

/**
 * get Cx or Cy - chromatic gradient component of gradinet field C
 * correspondents to part of Eq.(3) in papper
 * 
 * @param a2 - component a of color 2 in LAB color space
 * @param a1 - component a of color 1 in LAB color space
 * @param b2 - component b of color 2 in LAB color space
 * @param b1 - component b of color 1 in LAB color space
 * 
 * @return Cx or Cy
 */
double TMOEfficientC2GConversionForDIInGD::Chromatic_gradient_component(double a2, double a1, double b2, double b1){
	return sqrt(pow(a2 - a1, 2) + pow(b2 - b1, 2));
}

/**
 * compute one component (x-direction or y-direction) of gradient_field_C
 * correspondents to first part of Eq.(4) in papper
 * 
 * @param delta_l - Lx or Ly, depeneds on direction
 * @param a2 - component a of color 2 in LAB color space
 * @param a1 - component a of color 1 in LAB color space
 * @param b2 - component b of color 2 in LAB color space
 * @param b1 - component b of color 1 in LAB color space
 * @param c_component - chromatic gradient component in direction x or y (Cx or Cy)
 * 
 * @return component x or of gradient_field_C
 */
double TMOEfficientC2GConversionForDIInGD::Gradient_filed_component(
	double delta_l, double a2, double a1, double b2, double b1, double c_component){
	
//	return Sign(delta_l, a2, a1, b2, b1) * sqrt(pow(delta_l, 2) + pow(AttenuationFunction(c_component), 2));
	return Sign(delta_l, a2, a1, b2, b1) * sqrt(pow(delta_l, 2) + AttenuationFunction(pow(c_component, 2)));
//	return Sign(delta_l, a2, a1, b2, b1) * sqrt(pow(0.0, 2) + AttenuationFunction(pow(c_component, 2)));
}

/*void TMOEfficientC2GConversionForDIInGD::Get_grayscale_field_simple(Matrix m1, Matrix m2){
	
}*/

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOEfficientC2GConversionForDIInGD::Transform(){
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_LAB);
	pDst->Convert(TMO_LAB);

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double L, a, b, delta_l, a2, a1, b2, b1, Cx, Cy;
	double min_L = std::numeric_limits<double>::max();
	double max_L = std::numeric_limits<double>::min();
	
	// gradient field C for reconstruct grayscale image using PES
	// Matrix gradient_field_C(pSrc->GetWidth(), pSrc->GetHeight());
	// struct for x and y direction, because matrix.h doesnt support 3D matrices
	struct gradient_field_C_type
	{
		mtx::Matrix x;
		mtx::Matrix y;
	} gradient_field_C;
			
	// set dimensions of gradient_field_C
	gradient_field_C.x = mtx::Matrix(pSrc->GetWidth(), pSrc->GetHeight());
	gradient_field_C.y = mtx::Matrix(pSrc->GetWidth(), pSrc->GetHeight());
	
	// create resulting grayscale field
	mtx::Matrix grayscale_field(pSrc->GetWidth(), pSrc->GetHeight(), 0.0);
        
	// step 1: compute gradient field C
	for (int j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++){
			L = *pSourceData++;
			a = *pSourceData++;
			b = *pSourceData++;
			
			// 1a: x-direction
			// get delta_l, a2, a1, b2 and b1
			//delta_l = pSrc->GetPixel(i + 1, j)[0] - pSrc->GetPixel(i - 1, j)[0];
			delta_l = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[0] : pSrc->GetPixel(i + 1, j)[0];
			delta_l -= pSrc->GetPixel(i, j)[0];
			a2 = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i + 1, j)[1];
			a1 = pSrc->GetPixel(i, j)[1];
			b2 = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i + 1, j)[2];
			b1 = pSrc->GetPixel(i, j)[2];
			
			// get Cx
			Cx = Chromatic_gradient_component(a2, a1, b2, b1);
			//Cx *= 10000000.0;			
			
			// fill gradinet fild C (components x-direction and y-direction)
			gradient_field_C.x(i, j) = Gradient_filed_component(delta_l, a2, a1, b2, b1, Cx);
			
			// 1b: y-direction
			// get delta_l, a2, a1, b2 and b1
			//delta_l = pSrc->GetPixel(i, j + 1)[0] - pSrc->GetPixel(i, j - 1)[0];
			delta_l = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[0] : pSrc->GetPixel(i, j + 1)[0];
			delta_l -= pSrc->GetPixel(i, j)[0];
			a2 = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i, j + 1)[1];
			a1 = pSrc->GetPixel(i, j)[1]; 		// redundant
			b2 = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i, j + 1)[2];
			b1 = pSrc->GetPixel(i, j)[2];			// redundant
			
			// get Cy
			Cy = Chromatic_gradient_component(a2, a1, b2, b1);			
			//Cy *= 10000000.00;
			
			// fill gradinet fild C (components x-direction and y-direction)
			gradient_field_C.y(i, j) = Gradient_filed_component(delta_l, a2, a1, b2, b1, Cy);
		}
	}
	
	// step 2 SOLUION I (without PES)
	//compute grayscale image from gradient_field_C
	//Get_grayscale_field_simple();	
	/*
	// convert data from gradient_field_C.x|y to grayscale_field
	// step 2a: x axe	
	double previous = 0.0;
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) = previous + gradient_field_C.x(i, j);
			previous = grayscale_field(i, j);
			
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);
		}
	}
	for (int j = 0; j < pSrc->GetWidth(); j++){
		for (int i = 0; i < pSrc->GetHeight(); i++){
			grayscale_field(i, j) -= gradient_field_C.y(j, i);
			previous = grayscale_field(i, j);
			
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);
		}
	}	
	
	// step 2c: normalize values in grayscale_field
	double difference = abs(min_L) + abs(max_L);
	double divider = difference / 100.0;
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) += abs(min_L); // move to range [0; max_l]
			grayscale_field(i, j) /= divider; // normalize
			
			// show resulted lightness for debug
			//std::cerr << "L: " << grayscale_field(i, j) << std::endl;
		}
	}*/
	
	// step 2 SOLUION II (use PES solver)
	
	
	// step 2a: prepare variables	
	unsigned int n1 = pSrc->GetHeight(), n2 = pSrc->GetWidth();
	double h1=1.0, h2=1.0, a1p=1.0, a2p=1.0;		
	
	boost::multi_array<double,2> F(boost::extents[n2][n1]);	
	boost::multi_array<double,2> U(boost::extents[n2][n1]);		
	boost::multi_array<double,2> U2(boost::extents[n2][n1]);
	boost::multi_array<double,2> F2(boost::extents[n2][n1]);
	
	// copy data from gradient_field_C to F
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){			
			F[i][j] = gradient_field_C.x(i, j) + gradient_field_C.y(i, j);			
		}
	}	
	
	//pde::types::boundary bdtype=pde::types::Neumann;
	pde::types::boundary bdtype=pde::types::Dirichlet;	
	double bdvalue=pde::neumann_compat(F,a1p,a2p,h1,h2);								
	//double bdvalue = 0.0;
	std::cerr << "bdvalue: " << bdvalue << std::endl;	
	
	// set number of threads
	pde::fftw_threads(2);		
	
	// step 2b: run PES
	pde::poisolve(U,F,a1p,a2p,h1,h2,bdvalue,bdtype,false);	
	//pde::laplace(U,F,a1,a2,h1,h2,bdvalue,bdtype);  
	
	// step 2c: copy result to grayscale_field variable
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) = U[i][j];
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);			
		}
	}	
	
	// setp 2d: normalize grayscale_field
	double difference = abs(max_L - min_L);
	double divider = difference / 100.0;				
	std::cerr << "DIVIDER " << divider << std::endl;	
	
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) += (min_L > 0) ? -abs(min_L) : abs(min_L); // move to range [0; max_l]
			grayscale_field(i, j) = grayscale_field(i, j) / divider; // normalize
			
			// show resulted lightness for debug
			std::cerr << "L: " << grayscale_field(i, j) << std::endl;
		}
	}			
	
	// step 3: copy grayscale map to pDestinationData
	for (int j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){

			//std::cerr << gradient_field_C.x(i, j) << " " << gradient_field_C.y(i, j) << std::endl;						
			
			*pDestinationData++ = grayscale_field(i, j);						
			*pDestinationData++ = 0.0;
			*pDestinationData++ = 0.0;
		}
	}	
	
	//std::cerr << "min_L: " << min_L << ", max_L: " << max_L << std::endl;
	
	pDst->Convert(TMO_RGB);
	return 0;
}

