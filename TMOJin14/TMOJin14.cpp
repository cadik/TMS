/* --------------------------------------------------------------------------- *
 * TMOJin14.cpp: implementation of the TMOJin14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOJin14.h"
#include <cmath>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;

using namespace std;
#define PI 3.14159265
#define R_I 0
#define G_I 1
#define B_I 2
#define TAU_NOT_SET -1.
#define MU_NOT_SET -1.

typedef CGAL::Quadratic_program_from_iterators
<double**,                                             // for A
 double*,                                              // for b
 CGAL::Comparison_result*,                             // for r
 bool*,                                                // for fl
 double*,                                              // for l
 bool*,                                                // for fu
 double*,                                              // for u
 double**,                                             // for D
 double*>                                              // for c 
Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOJin14::TMOJin14()
{
	SetName(L"Jin14");
	SetDescription(L"A contrast maximization method for color-to-grayscale conversion"); 
	tau.SetName(L"tau");
	tau.SetDescription(L"...");
	tau.SetDefault(TAU_NOT_SET);
	tau=TAU_NOT_SET;
	this->Register(tau);
	mu.SetName(L"mu");
	mu.SetDescription(L"...");
	mu.SetDefault(MU_NOT_SET);
	mu=MU_NOT_SET;
	this->Register(mu);
	rescale.SetName(L"rescale");
	rescale.SetDescription(L"If image color should be rescaled to interval [0,1]");
	rescale.SetDefault(true);
	rescale=true;
	this->Register(rescale);
}

TMOJin14::~TMOJin14()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOJin14::Transform()
{
	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();	
	
	//compute mean value of R, G, B channels
	double sum_rgb[3] = {0,0,0};
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){
			sum_rgb[R_I] += *pSourceData++;
			sum_rgb[G_I] += *pSourceData++;
			sum_rgb[B_I] += *pSourceData++;
		}
	}
	double mean_rgb[3] = {sum_rgb[R_I]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                      sum_rgb[G_I]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                      sum_rgb[B_I]/(pSrc->GetHeight()*pSrc->GetWidth())};
	
	if(tau == TAU_NOT_SET){
		//compute tau based on average color
		tau=0.025;
		if((mean_rgb[R_I]+mean_rgb[G_I]+mean_rgb[B_I]) < 1) tau=0.5;
	}
	//compute covariance between color channels
	//https://en.wikipedia.org/wiki/Covariance
	pSourceData = pSrc->GetData();
	double sum_6[6] = {0,0,0,0,0,0};
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){
			double r = *pSourceData++ - mean_rgb[R_I];
			double g = *pSourceData++ - mean_rgb[G_I];
			double b = *pSourceData++ - mean_rgb[B_I];
			sum_6[0] += r*r;
			sum_6[1] += r*g;
			sum_6[2] += g*g;
			sum_6[3] += r*b;
			sum_6[4] += g*b;
			sum_6[5] += b*b;
		}
	}
	double cov_6[6] = {sum_6[0]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                 sum_6[1]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                 sum_6[2]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                 sum_6[3]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                 sum_6[4]/(pSrc->GetHeight()*pSrc->GetWidth()),
	                 sum_6[5]/(pSrc->GetHeight()*pSrc->GetWidth())};
	
	double mean_6[6] = {mean_rgb[R_I]*mean_rgb[R_I],
	                   mean_rgb[R_I]*mean_rgb[G_I],
	                   mean_rgb[G_I]*mean_rgb[G_I],
	                   mean_rgb[R_I]*mean_rgb[B_I],
	                   mean_rgb[G_I]*mean_rgb[B_I],
	                   mean_rgb[B_I]*mean_rgb[B_I]};
	
	double h_6[6] = {-2*cov_6[0] + tau*mean_6[0],
	                -2*cov_6[1] + tau*mean_6[1],
	                -2*cov_6[2] + tau*mean_6[2],
	                -2*cov_6[3] + tau*mean_6[3],
	                -2*cov_6[4] + tau*mean_6[4],
	                -2*cov_6[5] + tau*mean_6[5]};
	
	if(mu == MU_NOT_SET){
		//compute mu based on eigen value of hessian    
	
		//compute minimal eigen value
		//https://en.wikipedia.org/wiki/Eigenvalue_algorithm
		double eigen;
		double p1 = h_6[1]*h_6[1] + h_6[3]*h_6[3] + h_6[4]*h_6[4];
		if (p1 == 0){ 
			eigen = min(min(h_6[0],h_6[2]),h_6[5]);
		}else{
			double q = (h_6[0]+h_6[2]+h_6[5])/3;
			double p2 = (h_6[0] - q)*(h_6[0] - q) + (h_6[2] - q)*(h_6[2] - q) + (h_6[5] - q)*(h_6[5] - q) + 2 * p1;
			double p = sqrt(p2 / 6);
			double b_6[6] = {(1 / p) * (h_6[0] - q),
		                    (1 / p) * (h_6[1]),
		                    (1 / p) * (h_6[2] - q),
		                    (1 / p) * (h_6[3]),
		                    (1 / p) * (h_6[4]),
		                    (1 / p) * (h_6[5] - q)};
			double r = (b_6[0]*b_6[2]*b_6[5]+2*b_6[1]*b_6[3]*b_6[4]-b_6[0]*b_6[4]*b_6[4]-b_6[1]*b_6[1]*b_6[5]-b_6[2]*b_6[3]*b_6[3]) / 2;
			// In exact arithmetic for a symmetric matrix  -1 <= r <= 1
			// but computation error can leave it slightly outside this range.
			double phi;
			if (r <= -1) 
				phi = PI / 3;
			else if (r >= 1)
				phi = 0;
			else
				phi = acos(r) / 3;
		
			eigen = q + 2 * p * cos(phi + (2*PI/3));
		}
		if(eigen<0)
			mu=-eigen+0.0001;
		else 
			mu=0;
	}
	
	h_6[0]+=mu;
	h_6[2]+=mu;
	h_6[5]+=mu;
	
	double A1[] = {1};
	double A2[] = {1};
	double A3[] = {1};
	double* A[] = {A1, A2, A3}; // A comes columnwise
	double   b[] = {1};         // right-hand side
	CGAL::Comparison_result r[] = {CGAL::EQUAL};
	bool fl[] = {true, true, true}; // all x, y, z are lower-bounded
	double   l[] = {0, 0, 0};
	bool fu[] = {true, true, true}; // all x, y, z are upper-bounded
	double   u[] = {1, 1, 1};
	double  D1[] = {h_6[0]};                       // 2D_{0,0}
	double  D2[] = {h_6[1],h_6[2]};                 // 2D_{1,0}, 2D_{1,1}
	double  D3[] = {h_6[3],h_6[4],h_6[5]};           // 2D_{2,0}, 2D_{2,1}, 2D_{2,2}
	double*  D[] = {D1, D2, D3};                  // D-entries on/below diagonal 
	double c[] = {-tau*mean_rgb[R_I],-tau*mean_rgb[G_I],-tau*mean_rgb[B_I]};
	double c0 = 0;
	
	// now construct the quadratic program; the first two parameters are
	// the number of variables and the number of constraints (rows of A)
	Program qp (3, 1, A, b, r, fl, l, fu, u, D, c, c0);
	// solve the program, using ET as the exact type
	Solution s = CGAL::solve_quadratic_program(qp, ET());
	// output solution
	
	double alpha = CGAL::to_double(*(s.variable_values_begin()));
	double beta = CGAL::to_double(*(s.variable_values_begin()+1));
	double gama = CGAL::to_double(*(s.variable_values_begin()+2));
	
	pSourceData = pSrc->GetData();
	double min=1000000;
	double max=-1000000;
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++){
			double r = *pSourceData++;
			double g = *pSourceData++;
			double b = *pSourceData++;
			double out = r*alpha+g*beta+b*gama;
			if(out<min) min = out;
			if(out>max) max = out;
			*pDestinationData++ = out;
			*pDestinationData++ = out;
			*pDestinationData++ = out;
		}
	}
	
	if(rescale){
		pDestinationData = pDst->GetData();	
		for (int j = 0; j < pSrc->GetHeight(); j++)
		{
			pSrc->ProgressBar(j, pSrc->GetHeight());
			for (int i = 0; i < pSrc->GetWidth(); i++){
				double out = (*pDestinationData-min)/(max-min);
				*pDestinationData++ = out;
				*pDestinationData++ = out;
				*pDestinationData++ = out;
			}
		}
	}   
	return 0;
}

