/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                 diploma thesis                               *
*             Author: Petr Pospisil [xpospi68 AT stud.fit.vutbr.cz]            *
*                                    Brno 2016                                 *
*                                                                              *
*******************************************************************************/

/* --------------------------------------------------------------------------- *
 * TMOZheng15.cpp: implementation of the TMOZheng15 class.                     *
 *                 Efficient Color-to-Gray Conversion for Digital Images in    *
 *                 Gradinet Domain                                             *
 * Dependencies: FFTW library and boost::multi_array.                          *
 * Method number: 2                                                            *
 * --------------------------------------------------------------------------- */

#include "TMOZheng15.h"

// matrix library
#include "../tmolib/matrix.h"

// libraries for PES
#include <boost/multi_array.hpp>
#include "../tmolib/poisson_pde/laplace.h"

// removed
// #include "../pfstools/pde.h"

using namespace pde;

/**
 * constructor, prepare parameters
 */
TMOZheng15::TMOZheng15()
{
	SetName(L"Zheng15");
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
	
	// show/hide debug info
	verbose.SetName(L"v");
	verbose.SetDescription(L"Verbose output");
	verbose.SetDefault(false);
	verbose=false;	
	this->Register(verbose);	
	
	x_max = 0.0;
}

TMOZheng15::~TMOZheng15()
{
}

/**
 * attenuation function scales down the input signal
 * correspondents to part of Eq.(1) in papper
 * 
 * @param x - input value
 * @return scaled value
 */
double TMOZheng15::AttenuationFunction(double x){
	x_max = (x > x_max) ? x : x_max;		
	//std::cerr << "A in: " << x << "A out: " << x * (beta * (1 - pow(x / (CONSTANT_C * x_max), gamma))) << std::endl;	
	return x * (beta * (1 - pow(x / (CONSTANT_C * x_max), gamma)));	
}

/**
 * gets modulated difference of 2 colors
 * 
 * not used, integrated in GradientFiledComponent
 * 
 * @param delta_l - L2 - L1
 * @param delta_a - a2 - a1
 * @param delta_b - b2 - b1
 * 
 * @return delta_E - modulated difference of 2 colors
 */
/*double TMOZheng15::ModulatedColorDifference(double delta_l, double delta_a, double delta_b){
	return sqrt(pow(delta_l, 2) + pow(AttenuationFunction(sqrt(pow(delta_a, 2) + pow(delta_b, 2))), 2));
}*/

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
double TMOZheng15::Sign(double delta_l, double a2, double a1, double b2, double b1){	
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
double TMOZheng15::ChromaticGradientComponent(double a2, double a1, double b2, double b1){
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
double TMOZheng15::GradientFiledComponent(
	double delta_l, double a2, double a1, double b2, double b1, double c_component){		
	double result = Sign(delta_l, a2, a1, b2, b1) * sqrt(pow(delta_l, 2) + pow(AttenuationFunction(c_component), 2));
	
	// check if result is nan
	return (result != result) ? 0.0 : result;
}

/**
 * transformation function
 * @return exit code
 */
int TMOZheng15::Transform(){
	pSrc->Convert(TMO_LAB);
	pDst->Convert(TMO_LAB);

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double L, a, b, Lx, Ly, a2, a1, b2, b1, Cx, Cy;
	double min_L = std::numeric_limits<double>::max();
	double max_L = std::numeric_limits<double>::min();
	
	// gradient field C for reconstruct grayscale image using PES
	// struct for x and y direction, because matrix.h doesnt support 3D matrices
	// possibly could be replaced by boost.multi_array in next version
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
        
	// step 1: compute gradient field C SOLUTION I
	for (int j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){
			L = *pSourceData++;
			a = *pSourceData++;
			b = *pSourceData++;
			
			// 1a: x-direction
			// get Lx, a2, a1, b2, b1 and Cx
			//Lx = (i == 0) ? 0.0 : pSrc->GetPixel(i, j)[0] - pSrc->GetPixel(i - 1, j)[0];					// backward differenc
			Lx = (i == pSrc->GetWidth() - 1) ? 0.0 : pSrc->GetPixel(i + 1, j)[0] - pSrc->GetPixel(i, j)[0];			// forward difference
			//Lx = (i == pSrc->GetWidth() - 1 || i == 0) ? 0.0 : pSrc->GetPixel(i + 1, j)[0] - pSrc->GetPixel(i - 1, j)[0];	// central difference
			
			a2 = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i + 1, j)[1];
			a1 = pSrc->GetPixel(i, j)[1];
			b2 = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i + 1, j)[2];
			b1 = pSrc->GetPixel(i, j)[2];
			Cx = ChromaticGradientComponent(a2, a1, b2, b1);
			
			// fill gradinet fild C (components x-direction and y-direction)
			gradient_field_C.x(i, j) = GradientFiledComponent(Lx, a2, a1, b2, b1, Cx);
			
			// basic gradient
			//gradient_field_C.x(i, j) = Lx;
			
			// 1b: y-direction
			// get Ly, a2, a1, b2, b1 and Cy			
			//Ly = (j == 0) ? 0.0 : pSrc->GetPixel(i, j)[0] - pSrc->GetPixel(i, j - 1)[0];		 			 // backward difference
			Ly = (j == pSrc->GetHeight() - 1) ? 0.0 : pSrc->GetPixel(i, j + 1)[0] - pSrc->GetPixel(i, j)[0];		 // forward difference
			//Ly = (j == pSrc->GetHeight() - 1 || j == 0) ? 0.0 : pSrc->GetPixel(i, j + 1)[0] - pSrc->GetPixel(i, j - 1)[0]; // central difference
						
			a2 = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i, j + 1)[1];
			a1 = pSrc->GetPixel(i, j)[1]; 		// redundant
			b2 = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i, j + 1)[2];
			b1 = pSrc->GetPixel(i, j)[2];			// redundant
			Cy = ChromaticGradientComponent(a2, a1, b2, b1);
			
			// fill gradinet fild C (components x-direction and y-direction)
			gradient_field_C.y(i, j) = GradientFiledComponent(Ly, a2, a1, b2, b1, Cy);
			
			// basic gradient
			//gradient_field_C.y(i, j) = Ly;
		}
	}
	
	// if verbose, print gradient_field_C.x, gradient_field_C.y
	if (verbose){
		for (int j = 0; j < pSrc->GetHeight(); j++){		
			for (int i = 0; i < pSrc->GetWidth(); i++){	
				std::cerr << gradient_field_C.x(i, j) << ", " << gradient_field_C.y(i, j) << std::endl;
			}
		}
	}
	
	// step 1: compute gradient field C SOLUTION II
	/*for (int j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){
			gradient_field_C.x(i, j)
		}
	}*/
	
	// step 2 SOLUION I (without PES - simple integration)
	//compute grayscale image from gradient_field_C
		
	// convert data from gradient_field_C.x|y to grayscale_field
	// step 2a: x axe	
	/*double previous = 0.0;	
	for (int j = 0; j < pSrc->GetHeight(); j++){
		previous = pSrc->GetPixel(0, j)[0];		
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) = previous + gradient_field_C.x(i, j);
			previous = grayscale_field(i, j);
			
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);
		}
	}*/
	
	// step 2b: y axe			
 	/*for (int i = 0; i < pSrc->GetWidth(); i++){
		previous = pSrc->GetPixel(i, 0)[0];
		for (int j = 0; j < pSrc->GetHeight(); j++){			
			grayscale_field(i, j) += previous + gradient_field_C.y(i, j);						
			previous = gradient_field_C.y(i, j);
			
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);
		}
	}*/
	
	// step 2ab x+y combinated
	/*double previousX = 0.0;
	double previousY = 0.0;	
	for (int j = 0; j < pSrc->GetHeight(); j++){
		previousX = 0.0;
		previousY = 0.0;
		for (int i = 0; i < pSrc->GetWidth(); i++){			
			grayscale_field(i, j) = previousX + gradient_field_C.x(i, j) + previousY + gradient_field_C.y(i, j);			
			previousX += gradient_field_C.x(i, j);
			previousY += gradient_field_C.y(i, j);
			
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);
		}
	}*/
	
	// step 2c: normalize values in grayscale_field
	/*double difference = abs(max_L - min_L);
	double divider = difference / 100.0;
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) += (min_L > 0) ? -abs(min_L) : abs(min_L); // move to range [0; max_l]					
			grayscale_field(i, j) /= divider; // normalize
			
			// show resulted lightness for debug
			//std::cerr << "L: " << grayscale_field(i, j) << std::endl;
		}
	}*/
	// step 2 SOLUION II (use PES solver from http://kluge.in-chemnitz.de/opensource/poisson_pde/)
	
	// step 2a: prepare variables	
	unsigned int n1 = pSrc->GetHeight(), n2 = pSrc->GetWidth();	
	double h1=1.0, h2=1.0, a1p=1.0, a2p=1.0;		
	
	boost::multi_array<double,2> F(boost::extents[n2][n1]);
	boost::multi_array<double,2> U(boost::extents[n2][n1]);
	std::vector<double> bd1a (pSrc->GetHeight(), 0.0);
	std::vector<double> bd1b (pSrc->GetHeight(), 0.0);
	std::vector<double> bd2a (pSrc->GetWidth(), 0.0);
	std::vector<double> bd2b (pSrc->GetWidth(), 0.0);
	
	// copy data from gradient_field_C to F
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			// backward difference
			F[i][j] = (i == 0) ? 0.0 : gradient_field_C.x(i, j) - gradient_field_C.x(i - 1, j);
			F[i][j] += (j == 0) ? 0.0 : gradient_field_C.y(i, j) - gradient_field_C.y(i, j - 1);
			
			// central difference
			/*F[i][j] = (i == 0 || i == pSrc->GetWidth() - 1) ? 0.0 : gradient_field_C.x(i + 1, j) - gradient_field_C.x(i - 1, j);
			F[i][j] += (j == 0 || j == pSrc->GetHeight() - 1 ) ? 0.0 : gradient_field_C.y(i, j + 1) - gradient_field_C.y(i, j - 1);*/
			
			// forward difference
			/*F[i][j] = (i == pSrc->GetWidth() - 1) ? 0.0 : gradient_field_C.x(i + 1, j) - gradient_field_C.x(i, j);
			F[i][j] += (j == pSrc->GetHeight() - 1) ? 0.0 : gradient_field_C.y(i, j + 1) - gradient_field_C.y(i, j);*/
		}
	}		
	
	// set boundary
	for (int i = 0; i < pSrc->GetWidth(); i++){
		bd2a[i] = pSrc->GetPixel(i, 0)[0];		
		bd2b[i] = pSrc->GetPixel(i, pSrc->GetHeight() - 1)[0];		
	}
	for (int j = 0; j < pSrc->GetHeight(); j++){
		bd1a[j] = pSrc->GetPixel(0, j)[0];		
		bd1b[j] = pSrc->GetPixel(pSrc->GetWidth() - 1, j)[0];		
	}
	
	//pde::types::boundary bdtype=pde::types::Neumann;
	pde::types::boundary bdtype=pde::types::Dirichlet;		
	double bdvalue = 0.0;
	if (bdtype==pde::types::Neumann){      
		bdvalue=pde::neumann_compat(F,a1p,a2p,h1,h2);      
	}	
	if (verbose) std::cerr << "bdvalue: " << bdvalue << std::endl;	
	
	// set number of threads
	pde::fftw_threads(1);		
	
	// step 2b: run PES
	//double trunc = pde::poisolve(U, F, a1p, a2p, h1, h2, bdvalue, bdtype, false);		
	double trunc = pde::poisolve(U, F, a1p, a2p, h1, h2, bd1a, bd1b, bd2a, bd2b, bdtype, false);	
			
	// step 2c: copy result to grayscale_field variable
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) = U[i][j];
			if (U[i][j] > max_L) max_L = U[i][j];
			if (U[i][j] < min_L) min_L = U[i][j];			
		}
	}	
	
	// setp 2d: normalize grayscale_field
	double difference = abs(max_L - min_L);
	double divider = difference / 100;	
	if (verbose) std::cerr << "min: " << min_L << ", max: " << max_L << std::endl;
	if (verbose) std::cerr << "DIVIDER " << divider << std::endl;	
	
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) -= min_L; 			// move to range [0; max_l]					
			grayscale_field(i, j) /= divider; 			// normalize			
		}
	}
	
	// step 2 SOLUION III (use pfstools)
	// step 2a: prepare data
	/*pfstmo::Array2D in = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	pfstmo::Array2D out = pfstmo::Array2D(pSrc->GetWidth(),pSrc->GetHeight());
	
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){			
			in(i, j) = (i == 0) ? 0.0 : gradient_field_C.x(i, j) - gradient_field_C.x(i - 1, j);
			in(i, j) += (j == 0) ? 0.0 : gradient_field_C.y(i, j) - gradient_field_C.y(i, j - 1);
			
			// ~ same results
			//in(i, j) = (i == pSrc->GetWidth() - 1) ? 0.0 : gradient_field_C.x(i + 1, j) - gradient_field_C.x(i, j);
			//in(i, j) += (j == pSrc->GetHeight() - 1) ? 0.0 : gradient_field_C.y(i, j + 1) - gradient_field_C.y(i, j);
			
			// worse results
			//in(i, j) = (i == pSrc->GetWidth() - 1 || i == 0) ? 0.0 : gradient_field_C.x(i + 1, j) - gradient_field_C.x(i - 1, j);
			//in(i, j) += (j == pSrc->GetHeight() - 1 || j == 0) ? 0.0 : gradient_field_C.y(i, j + 1) - gradient_field_C.y(i, j - 1);
		}
	}	
	
	// step 2b: call solver	
	solve_pde_multigrid( &in, &out);
	//solve_pde_sor( &in, &out, 50);
	//exact_sollution( &in, &out);	
	//solve_pde_fft(&in, &out, true);
	
	// step 2c: copy output data
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			grayscale_field(i, j) = out(i, j);			
			
			if (grayscale_field(i, j) > max_L) max_L = grayscale_field(i, j);
			if (grayscale_field(i, j) < min_L) min_L = grayscale_field(i, j);			
		}
	}
	
	// setp 2d: normalize grayscale_field
	double difference = (max_L - min_L);
	double divider = difference / 100.0;		
	std::cerr << "min: " << min_L << ", max: " << max_L << std::endl;
	std::cerr << "DIVIDER " << divider << std::endl;	
	
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){			
			grayscale_field(i, j) -= min_L; // move to range [0; max_l]
			grayscale_field(i, j) /= divider; // normalize			
		}
	}*/	
	
	double comulatedError = 0.0;
	
	// step 3: copy grayscale map to pDestinationData
	for (int j = 0; j < pSrc->GetHeight(); j++){
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){			
			if (verbose) std::cerr << j << ", " << i << ": Original: " << pSrc->GetPixel(i, j)[0] << "  Restored: " << grayscale_field(i, j) << std::endl;
			comulatedError += abs(pSrc->GetPixel(i, j)[0] - grayscale_field(i, j));
						
			*pDestinationData++ = grayscale_field(i, j);		
			*pDestinationData++ = 0.0;
			*pDestinationData++ = 0.0;
		}
	}	
	
	comulatedError /= pSrc->GetHeight() * pSrc->GetWidth();
	
	if (verbose) std::cerr << "trunc: " << trunc << std::endl;
	if (verbose) std::cerr << "AVERAGE ERROR: " << comulatedError << std::endl;
	pDst->Convert(TMO_RGB);	
	return 0;
}
