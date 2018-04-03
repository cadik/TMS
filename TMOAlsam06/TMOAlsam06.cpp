/* --------------------------------------------------------------------------- *
 * TMOBiswas05.cpp: implementation of the Biswas K. K., Pattanaik S.		   *
 * A Simple Spatial Tone Mapping Operator for High Dynamic Range Images.   	   *
 * --------------------------------------------------------------------------- */

#include "TMOAlsam06.h"

#include <algorithm>
#include <math.h> 

#include <iostream>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAlsam06::TMOAlsam06()
{
	SetName(L"Alsam06");								
	SetDescription(L"Generates greyscale image with color separation and texture enhancement.");	
	
	alpha.SetName(L"alpha");				
	alpha.SetDescription(L"Weight of the RED color in final image.");	
	alpha.SetDefault(0.30);
	alpha=0.30;
	alpha.SetRange(0.0,1.0);
	this->Register(alpha);

	beta.SetName(L"beta");				
	beta.SetDescription(L"Weight of the GREEN color in final image.");	
	beta.SetDefault(0.59);
	beta=0.59;
	beta.SetRange(0.0,1.0);				
	this->Register(beta);

	gamma.SetName(L"gamma");				
	gamma.SetDescription(L"Weight of the BLUE color in final image.");	
	gamma.SetDefault(0.11);
	gamma=0.11;
	gamma.SetRange(0.0,1.0);				
	this->Register(gamma);
}

TMOAlsam06::~TMOAlsam06()
{
}

void TMOAlsam06::RecalculateWeights() {
	double sum = alpha + beta + gamma;

	if (sum != 1) {
		alpha /= sum;
		beta /= sum;
		gamma /= sum;
	}
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAlsam06::Transform()
{
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			

	RecalculateWeights();

	int x, y;
	double r, g, b;
	double out;

	for (y = 0; y < pSrc->GetHeight(); y++) {
		pSrc->ProgressBar(y, pSrc->GetHeight());

		for (x = 0; x < pSrc->GetWidth(); x++) {
			r = *pSourceData++;
			g = *pSourceData++;
			b = *pSourceData++;

			out = alpha * r + beta * g + gamma * b;
			std::cout << out << std::endl;

			*pDestinationData++ = out;
			*pDestinationData++ = out;
			*pDestinationData++ = out;
		}
	}

	pSrc->ProgressBar(y, pSrc->GetHeight());
	return 0;
}

