/* --------------------------------------------------------------------------- *
 * TMOSlomp12.cpp: implementation of the TMOSlomp12 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOSlomp12.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOSlomp12::TMOSlomp12()
{
	SetName(L"YourOperatorName");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOSlomp12::~TMOSlomp12()
{
}

cv::Mat TMOSlomp12::TMOImageToLogLuminanceMat(TMOImage *pSrc)
{
	pSrc->Convert(TMO_XYZ);
	double *pSourceData = pSrc->GetData();
	cv::Mat logLuminanceMat(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC1);

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			*pSourceData++;
			double Y = *pSourceData++;
			*pSourceData++;

			logLuminanceMat.at<double>(y, x) = log10(Y);
		}
	}
}

int TMOSlomp12::Transform()
{

	return 0;
}
