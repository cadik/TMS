/* --------------------------------------------------------------------------- *
 * TMONguy15.cpp: implementation of the TMONguy15 class.   *
 * --------------------------------------------------------------------------- */

#include "TMONguy15.h"
#include "mainprepare.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMONguy15::TMONguy15()
{
	SetName(L"Nguy15");						// TODO - Insert operator name
	SetDescription(L"Color to grayscale");	// TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	//this->Register(dParameter);
}

TMONguy15::~TMONguy15()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMONguy15::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	//pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	//pDst->Convert(TMO_Yxy);								// x, y as color information

    double* pSourceData = pSrc->GetData();				// You can work at low level data

    double* pDestinationData = pDst->GetData();			// Data are stored in form of array
    double* pom = mainprepare(pSourceData,pSrc->GetWidth(),pSrc->GetHeight());		// entering main calculations
    // three colour components
	double pY, px, py;

	int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		//pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pom++;
			px = *pom++;
			py = *pom++;

			// Here you can use your transform 
			// expressions and techniques...
			//pY *= dParameter;							// Parameters can be used like
														// simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	//pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0;
}

