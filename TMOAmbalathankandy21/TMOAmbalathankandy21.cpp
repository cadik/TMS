/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAmbalathankandy21.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAmbalathankandy21::TMOAmbalathankandy21()
{
	SetName(L"Ambalathankandy21");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

//	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
//	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
//	dParameter.SetDefault(1);							// TODO - Add default values
//	dParameter = 1.;
//	 dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
//	this->Register(dParameter);
}

TMOAmbalathankandy21::~TMOAmbalathankandy21()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAmbalathankandy21::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	//double pY, px, py;
   double R, G, B, LWhite, LRG, LWarm, LCool, LB, L; 
   double parameter = 0.75;

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			R = *pSourceData++;
			G = *pSourceData++;
			B = *pSourceData++;

         LWhite = sqrt((R*R + G*G + B*B)/3);
         LB = B;
         LRG = sqrt(parameter*R*R + (1-parameter)*G*G);
         LWarm = R/(R+G+B)*LRG + (1-R/(R+G+B))*LWhite;
         LCool = (1-B/(R+G+B))*LB + B/(R+G+B)*LWhite;
         L = sqrt(parameter*LWarm*LWarm + (1-parameter)*LCool*LCool);
			// Here you can use your transform
			// expressions and techniques...
			//pY *= dParameter; // Parameters can be used like
							  // simple variables

			// and store results to the destination image
			*pDestinationData++ = L;
			*pDestinationData++ = L;
			*pDestinationData++ = L;
		}
	}
//	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0; 
}
