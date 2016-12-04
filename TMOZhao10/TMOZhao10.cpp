/* --------------------------------------------------------------------------- *
 * TMOZhao10.cpp: implementation of the TMOZhao10 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOZhao10.h"
#include <fftw3.h>
#include <math.h>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOZhao10::TMOZhao10()
{
	SetName(L"Zhao10");						// TODO - Insert operator name
	SetDescription(L"Add your TMO description here");	// TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOZhao10::~TMOZhao10()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOZhao10::Transform()
{
	double* data;
	int size = pSrc->GetHeight()*pSrc->GetWidth();
	int spec_size = pSrc->GetHeight()*(pSrc->GetWidth()/2+1);
	double *rgb[] = {fftw_alloc_real(size),fftw_alloc_real(size),fftw_alloc_real(size)};
	double *lab[] = {fftw_alloc_real(size),fftw_alloc_real(size),fftw_alloc_real(size)};
	fftw_complex *spec_rgb[] = {fftw_alloc_complex(spec_size),fftw_alloc_complex(spec_size),fftw_alloc_complex(spec_size)};
	fftw_complex *spec_lab[] = {fftw_alloc_complex(spec_size),fftw_alloc_complex(spec_size),fftw_alloc_complex(spec_size)};
	
	double *theta = fftw_alloc_real(spec_size);
	double *phi = fftw_alloc_real(spec_size);
	
	fftw_complex *spec_gray = fftw_alloc_complex(spec_size);
	double *gray = fftw_alloc_real(size);
	
	fftw_plan p = fftw_plan_dft_r2c_2d(pSrc->GetHeight(), pSrc->GetWidth(),rgb[0], spec_rgb[0], FFTW_ESTIMATE);
	
	//copy data channels to r,g,b arrays
	data=pSrc->GetData();
	for(int i=0;i<size;++i){
		rgb[0][i] = *data++;
		rgb[1][i] = *data++;
		rgb[2][i] = *data++;
	}
	
	//transform to Lab space
	pSrc->Convert(TMO_LAB);
	pDst->Convert(TMO_LAB);
	
	
	//copy data channels to l,a,b array
	data=pSrc->GetData();
	for(int i=0;i<size;++i){
		lab[0][i] = *data++/100;
		lab[1][i] = *data++/100;
		lab[2][i] = *data++/100;
		
		//fprintf(stderr,"%f %f %f\n",lab[0][i],lab[1][i],lab[2][i]);
	}
	
	//compute fft of all channels
	fftw_execute_dft_r2c(p,rgb[0], spec_rgb[0]);
	fftw_execute_dft_r2c(p,rgb[1], spec_rgb[1]);
	fftw_execute_dft_r2c(p,rgb[2], spec_rgb[2]);
	fftw_execute_dft_r2c(p,lab[0], spec_lab[0]);
	fftw_execute_dft_r2c(p,lab[1], spec_lab[1]);
	fftw_execute_dft_r2c(p,lab[2], spec_lab[2]);
	
	fftw_destroy_plan(p);
	p = fftw_plan_dft_c2r_2d(pSrc->GetHeight(), pSrc->GetWidth(),spec_gray, gray, FFTW_ESTIMATE);
	
	//compute phi and theta coefficient
	double thetasum=0;
	double phisum=0;
	for(int i=0;i<spec_size;++i){
		double a2 = spec_lab[1][i][0]*spec_lab[1][i][0]+spec_lab[1][i][1]*spec_lab[1][i][1];
		double b2 = spec_lab[2][i][0]*spec_lab[2][i][0]+spec_lab[2][i][1]*spec_lab[2][i][1];
		phi[i] = a2/(a2+b2);
		phisum += phi[i];
		
		double rr2 = spec_rgb[0][i][0]*spec_rgb[0][i][0]+spec_rgb[0][i][1]*spec_rgb[0][i][1];
		double gg2 = spec_rgb[1][i][0]*spec_rgb[1][i][0]+spec_rgb[1][i][1]*spec_rgb[1][i][1];
		double bb2 = spec_rgb[2][i][0]*spec_rgb[2][i][0]+spec_rgb[2][i][1]*spec_rgb[2][i][1];
		double l2 = spec_lab[0][i][0]*spec_lab[0][i][0]+spec_lab[0][i][1]*spec_lab[0][i][1];
		double rgb2 = rr2+gg2+bb2;
		theta[i] = (rgb2-l2)/rgb2;
		thetasum += theta[i];
	}
	
	thetasum=0.6*spec_size;
	phisum=0.9*spec_size;
	for(int i=0;i<spec_size;++i){
		spec_gray[i][0] = ((1-thetasum/spec_size)*spec_lab[0][i][0]+thetasum/spec_size*(phisum/spec_size*spec_lab[1][i][0]+(1-phisum/spec_size)*spec_lab[2][i][0]))/size;
		//spec_gray[i][0] = ((1-theta[i])*spec_lab[0][i][0]+theta[i]*(phi[i]*spec_lab[1][i][0]+(1-phi[i])*spec_lab[2][i][0]))/size;
		spec_gray[i][1] = ((1-thetasum/spec_size)*spec_lab[0][i][1]+thetasum/spec_size*(phisum/spec_size*spec_lab[1][i][1]+(1-phisum/spec_size)*spec_lab[2][i][1]))/size; 
		//spec_gray[i][1] = ((1-theta[i])*spec_lab[0][i][1]+theta[i]*(phi[i]*spec_lab[1][i][1]+(1-phi[i])*spec_lab[2][i][1]))/size; 
	}
	
	fftw_execute(p);
	
	double minimum = 99999999999999999;
	double maximum = -99999999999999999;
	data=pDst->GetData();
	for(int i=0;i<size;++i){
		if(gray[i]>maximum) maximum=gray[i];
		if(gray[i]<minimum) minimum=gray[i];
		*data++ = gray[i]*100;
		*data++ = 0;//gray[i];
		*data++ = 0;//gray[i];
	}
	
	fprintf(stderr,"%f %f %d %d\n",minimum,maximum,size,spec_size);
	data=pDst->GetData();
	for(int i=0;i<size;++i){
		*data = 100*(gray[i]-minimum)/(maximum-minimum);
		data += 3;
	}
	
	fftw_destroy_plan(p);
	fftw_free(gray); fftw_free(spec_gray);
	fftw_free(phi); fftw_free(theta);
	fftw_free(rgb[0]);fftw_free(rgb[1]);fftw_free(rgb[2]);
	fftw_free(lab[0]);fftw_free(lab[1]);fftw_free(lab[2]);
	fftw_free(spec_rgb[0]);fftw_free(spec_rgb[1]);fftw_free(spec_rgb[2]);
	fftw_free(spec_lab[0]);fftw_free(spec_lab[1]);fftw_free(spec_lab[2]);
	
	pDst->Convert(TMO_RGB);
	return 0;
}

