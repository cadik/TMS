/* --------------------------------------------------------------------------- *
 * TMOMeylan06.cpp: implementation of the TMOMeylan06 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMeylan06.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMeylan06::TMOMeylan06()
{
	SetName(L"Meylan06");						// TODO - Insert operator name
	SetDescription(L"Add your TMO description here");	// TODO - Insert description
	/*
	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
	*/
}

TMOMeylan06::~TMOMeylan06()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMeylan06::Transform()
{

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array
														// of three doubles representing
														// three colour components
	double pY, px, py;

        int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
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

	pDst->Convert(TMO_RGB);
	return 0;
}
