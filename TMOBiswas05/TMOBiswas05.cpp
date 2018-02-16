/* --------------------------------------------------------------------------- *
 * TMOBiswas05.cpp: implementation of the TMOBiswas05 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOBiswas05.h"

#include <algorithm>
#include <math.h> 

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOBiswas05::TMOBiswas05()
{
	SetName(L"Biswas05");								
	SetDescription(L"This simple operator takes into account local and global average of the luminance");	
	
	dParameter.SetName(L"Average luminance multiplier");				
	dParameter.SetDescription(L"Used to adjust the value of the average luminance of the whole image");	
	dParameter.SetDefault(0.15);							
	dParameter=1.;
	dParameter.SetRange(0.1,0.3);				
	this->Register(dParameter);
}

TMOBiswas05::~TMOBiswas05()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOBiswas05::Transform()
{
	pSrc->Convert(TMO_Yxy);								
	pDst->Convert(TMO_Yxy);								

	double* pSourceData = pSrc->GetData();	
	double* pSourceDataArray = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			

	double pY, px, py;
	
	//average luminance
	double ya = 0;
	int pixelCount = pSrc->GetHeight()*pSrc->GetWidth(); 
	for (int i = 0; i < pixelCount; i++)
		ya += pSourceDataArray[i*3];
	ya /= pixelCount;

    int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			//global contrast
			double gc = dParameter * ya;
			
			//local luminance computed with median filter
			double window[MEDIAN_DIMENSION*MEDIAN_DIMENSION];
			int halfDim = MEDIAN_DIMENSION/2;
			for(int y=-halfDim; y<=halfDim; y++)
				for(int x=-halfDim; x<halfDim; x++)
				{
					int ix = x+i;
					int iy = y+j;
					
					//boundaries
					if(ix<0)
						ix=0;
					else if(ix>=pSrc->GetWidth())
						ix = pSrc->GetWidth()-1;
						
					if(iy<0)
						iy=0;
					else if(iy>=pSrc->GetHeight())
						iy = pSrc->GetHeight()-1;
				
					window[y*MEDIAN_DIMENSION + x] = pSourceDataArray[iy*pSrc->GetWidth()*3 + ix*3];
				}
				
			std::sort(window, window + MEDIAN_DIMENSION*MEDIAN_DIMENSION);
			double yl = window[halfDim+1];

			
			//offset to avoid singularity when computing log
			double offset = 0.00000001;
			double cl = yl*log(offset + yl/pY) + gc; 
			pY = pY/(pY + cl);						
			
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}

