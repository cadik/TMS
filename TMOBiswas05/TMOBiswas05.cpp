/* --------------------------------------------------------------------------- *
 * TMOBiswas05.cpp: implementation of the Biswas K. K., Pattanaik S.		   *
 * A Simple Spatial Tone Mapping Operator for High Dynamic Range Images.   	   *
 * --------------------------------------------------------------------------- */

#include "TMOBiswas05.h"

#include <algorithm>
#include <math.h> 

#include <iostream>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOBiswas05::TMOBiswas05()
{
	SetName(L"Biswas05");								
	SetDescription(L"This simple operator takes into account local and global average of the luminance");	
	
	dParameter.SetName(L"lum");				
	dParameter.SetDescription(L"Used to adjust the value of the average luminance of the whole image");	
	dParameter.SetDefault(0.15);							
	dParameter=0.2;
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
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			

	double pY, px, py;
	
	//average luminance and copy colors
	double ya = 0; 
	for (int x = 0; x < pSrc->GetWidth(); x++)
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			ya += pSrc->GetLuminance(x,y);
			*pDestinationData++ = *pSourceData++;
			*pDestinationData++ = *pSourceData++;
			*pDestinationData++ = *pSourceData++;
		}
	ya /= pSrc->GetHeight()*pSrc->GetWidth();
	
	std::cerr << "Global average luminance: " << ya << std::endl;

	//global contrast
	double gc = dParameter * ya;

    int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{		
			pY = pSrc->GetLuminance(i,j);
			
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
				
					window[y*MEDIAN_DIMENSION + x] = pSrc->GetLuminance(ix,iy);
				}
				
			std::sort(window, window + MEDIAN_DIMENSION*MEDIAN_DIMENSION);
			double yl = window[halfDim+1];

			
			//offset to avoid singularity when computing log
			double offset = 0.00001;
			double cl = yl*log(offset + yl/pY) + gc; 
			pY = pY/(pY + cl);		
			if (pY > 1.0)
				pY = 1.0;
				
			pDst->SetLuminance(i,j, pY);
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	return 0;
}

