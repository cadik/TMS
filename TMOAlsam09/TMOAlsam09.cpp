/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                bachelor thesis                               *
*             Author: Martin Molek [xmolek00 AT stud.fit.vutbr.cz]             *
*                                    Brno 2017                                 *
*                                                                              *
*******************************************************************************/

/*---------------------------------------------------------------------------- *
 * TMOAlsam09.cpp: Fast Multispectral2Gray				       *
 *                 Ali Alsam, Mark S. Drew				       *
 * Method number: 2MM                                                          *
 * --------------------------------------------------------------------------- */


#include "TMOAlsam09.h"
#include <fftw3.h>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAlsam09::TMOAlsam09()
{
	SetName(L"Alsam09");						
	SetDescription(L"Add your TMO description here");	

	dParameter.SetName(L"ParameterName");				
	dParameter.SetDescription(L"ParameterDescription");	
	dParameter.SetDefault(1);							
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				
	this->Register(dParameter);
}

TMOAlsam09::~TMOAlsam09()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int maxChannel (int x, int y, int maxx, int maxy, double* pStartSourceData) 
{
	if ((pStartSourceData[3 * (y * maxy + x)] > pStartSourceData[3 * (y * maxy + x)+1]) && (pStartSourceData[3 * (y * maxy + x)] > pStartSourceData[3 * (y * maxy + x)+2])) //R>G && R>B
		return 0;
	if ((pStartSourceData[3 * (y * maxy + x)+1] > pStartSourceData[3 * (y * maxy + x)]) && (pStartSourceData[3 * (y * maxy + x)+1] > pStartSourceData[3 * (y * maxy + x)+2])) //G>R && G>R
		return 1;
	return 2;
}

double pixelMax (int x, int y, int maxx, int maxy, double* pStartSourceData) 
{
	return std::max(std::max(pStartSourceData[3 * (y * maxy + x)], pStartSourceData[3 * (y * maxy + x)+1]), pStartSourceData[3 * (y * maxy + x)+2]);
}

double pixelMaxSum (int x, int y, int maxx, int maxy, double* pStartSourceData) 
{
	double sum = 0.0;
	int outOfPicture = 0;
	
	if (x-1 > 0)
		sum += std::max(std::max(pStartSourceData[3 * (y * maxy + x-1)], pStartSourceData[3 * (y * maxy + x-1)+1]), pStartSourceData[3 * (y * maxy + x-1)+2]);
	else
		outOfPicture++;
	if (y-1 > 0)
		sum += std::max(std::max(pStartSourceData[3 * ((y-1) * maxy + x)], pStartSourceData[3 * ((y-1) * maxy + x)+1]), pStartSourceData[3 * ((y-1) * maxy + x)+2]);
	else
		outOfPicture++;
	if (x+1 < maxx)
		sum += std::max(std::max(pStartSourceData[3 * (y * maxy + x+1)], pStartSourceData[3 * (y * maxy + x+1)+1]), pStartSourceData[3 * (y * maxy + x+1)+2]);
	else
		outOfPicture++;
	if (y+1 < maxy)
		sum += std::max(std::max(pStartSourceData[3 * ((y+1) * maxy + x)], pStartSourceData[3 * ((y+1) * maxy + x)+1]), pStartSourceData[3 * ((y+1) * maxy + x)+2]);
	else
		outOfPicture++;
	sum += outOfPicture * std::max(std::max(pStartSourceData[3 * (y * maxy + x)], pStartSourceData[3 * (y * maxy + x)+1]), pStartSourceData[3 * (y * maxy + x)+2]);
	return sum;
}

double gradSum (int x, int y, int maxx, int maxy, int channel, double* pStartSourceData) 
{
	double sum = 0.0;
	double me = pStartSourceData[3 * (y * maxy + x)+channel];
	int outOfPicture = 0;
	if (x-1 > 0)
		sum += pStartSourceData[3 * (y * maxy + x-1)+channel] - me;
	if (y-1 > 0)
		sum += pStartSourceData[3 * ((y-1) * maxy + x)+channel] - me;
	if (x+1 < maxx)
		sum += pStartSourceData[3 * (y * maxy + x+1)+channel] - me;
	if (y+1 < maxy)
		sum += pStartSourceData[3 * ((y+1) * maxy + x)+channel] - me;
	return sum;
}

int TMOAlsam09::Transform()
{
	
	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();
	double r, g, b, shade;
	int mC, i;
	int size = pSrc->GetHeight()*pSrc->GetWidth()*3;
	double* tmp1 = fftw_alloc_real(size);
	double* tmp2 = fftw_alloc_real(size);
	for (i = 0; i<size; i++)
		tmp2[i] = pSourceData[i];

	for (int omega = 0; omega < 5; omega++)	
	{
		for (i = 0; i<size; i++)
			tmp1[i] = tmp2[i];
		i = 0;
		for (int y = 0; y < pSrc->GetHeight(); y++){
			for (int x = 0; x < pSrc->GetWidth(); x++){
				
				mC = maxChannel(x,y,pSrc->GetHeight(),pSrc->GetWidth(), pSrc->GetData());
				tmp2[i+mC] = tmp1[i+mC] + gradSum(x, y, pSrc->GetHeight(), pSrc->GetWidth(), mC, tmp1)/4.0;

				// if (x >	303 && y == 303)
				//	std::cerr << "GS: " << gradSum(x, y, pSrc->GetHeight(), pSrc->GetWidth(), mC, tmp1) << std::endl; 
				// shade = r * 0.299 + g * 0.587 + b * 0.114;

				if (tmp2[i+mC] > 1.0) {
					tmp2[i+mC] = 1.0;
				} else if (tmp2[i+mC] < 0.0) {
					tmp2[i+mC] = 0.0;
				}
				i += 3;
			}	
			pSrc->ProgressBar(y, pSrc->GetHeight());
		}
	}

	for (i = 0; i<size; i += 3) 
	{
		shade = tmp2[i] * 0.299 + tmp2[i+1] * 0.587 + tmp2[i+2] * 0.114;
		pDestinationData[i] = shade;
		pDestinationData[i+1] = shade;
		pDestinationData[i+2] = shade;
	}

	//pDst->Convert(TMO_RGB);
	fftw_free(tmp1);fftw_free(tmp2);
	return 0;
}

