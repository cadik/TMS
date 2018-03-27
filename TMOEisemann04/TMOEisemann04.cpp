/* --------------------------------------------------------------------------- *
 * TMOEisemann04.cpp: implementation of the TMOEisemann04 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOEisemann04.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOEisemann04::TMOEisemann04()
{
	SetName(L"Eisemann04");						// TODO - Insert operator name
	SetDescription(L"Add your TMO description here");	// TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOEisemann04::~TMOEisemann04()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
TMOEisemann04::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components
	double pY, px, py;

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform 
			// expressions and techniques...
			pY *= dParameter;							// Parameters can be used like
														// simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}