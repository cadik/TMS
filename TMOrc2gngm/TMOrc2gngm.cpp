/* --------------------------------------------------------------------------- *
 * TMOrc2gngm.cpp: implementation of the TMOrc2gngm class.   *
 * 	rc2gngm = Robust Color-to-gray via Nonlinear Global Mapping
 * --------------------------------------------------------------------------- */

#include "TMOrc2gngm.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOrc2gngm::TMOrc2gngm()
{
	SetName(L"rc2gngm");
	SetDescription(L"Robust Color-to-gray via Nonlinear Global Mapping");

	//dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	//dParameter.SetDescription(L"ParameterDescription");		// TODO - Insert parameter descriptions
	//dParameter.SetDefault(1);					// TODO - Add default values
	//dParameter=1.;
	//dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	//this->Register(dParameter);
}

TMOrc2gngm::~TMOrc2gngm()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOrc2gngm::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can convert it into other format
	pSrc->Convert(TMO_LCH);
	pDst->Convert(TMO_LCH);

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
									// of three doubles representing
									// three colour components
	double l, c, h, g, f;

        int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());		// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			l = *pSourceData++;
			c = *pSourceData++;
			h = *pSourceData++;

			f = 1.0; // TODO!
			g = l + f*c;					// global mapping

			// store results to the destination image
			*pDestinationData++ = g;
			*pDestinationData++ = 1.0/3.0;
			*pDestinationData++ = 1.0/3.0;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}

