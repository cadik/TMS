/* --------------------------------------------------------------------------- *
 * TMOMikamo14.cpp: implementation of the TMOMikamo14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMikamo14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMikamo14::TMOMikamo14()
{
	SetName(L"Mikamo14");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOMikamo14::~TMOMikamo14()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMikamo14::Transform()
{
	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();

	double luminanceAverage = 0;

	for (int x = 0; x < pSrc->GetWidth(); x++)
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			luminanceAverage += pSrc->GetLuminance(x, y);
		}

	luminanceAverage /= pSrc->GetHeight() * pSrc->GetWidth();

	double lambda_l, lambda_m, lambda_s, k1, k2, k3, k4, k5, k6;

	lambda_l = -18.3 / (1 + 7.2 * pow(I, -0.7)) - 0.9;
	lambda_m = -44.6 / (1 + 35.4 * pow(I, -1.2)) + 22.0;
	lambda_s = 43.0 / (1 + 9.0 * pow(I, -1.5)) + 28.0;
	k1 = 6.69 / (1 + 2500 * pow(I, -2.65)) + 0.80;
	k2 = -6.24 / (1 + 2500 * pow(I, -2.5)) - 0.77;
	k3 = 0.36 / (1 + 50.02 * pow(I, -1.5)) + 0.04;
	k4 = 0.24 / (1 + 50.04 * pow(I, -1.7)) + 0.03;
	k5 = 0.42 / (1 + 1.76 * pow(I, -0.02)) + 0.14;
	k6 = 0.15 / (1 + 2.8 * pow(I, -0.46)) - 0.27;

	

	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py;

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			pY *= dParameter; // Parameters can be used like
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
