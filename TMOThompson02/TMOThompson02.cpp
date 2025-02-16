/* --------------------------------------------------------------------------- *
 * TMOThompson02.cpp: implementation of the TMOThompson02 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOThompson02.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOThompson02::TMOThompson02()
{
	SetName(L"Thompson02");							  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOThompson02::~TMOThompson02()
{
}

TMOThompson02::matAndDouble TMOThompson02::getLuminanceMat()
{
	matAndDouble result;
	double max = 0.0;
	cv::Mat luminanceMat(pSrc->GetWidth(), pSrc->GetHeight(), CV_64F);

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double *pixel = pSrc->GetPixel(x, y);
			double L = 0.27 * pixel[0] + 0.67 * pixel[1] + 0.06 * pixel[2];
			luminanceMat.at<double>(x, y) = L;
			if (L > max)
			{
				max = L;
			}
		}
	}

	result.mat = luminanceMat;
	result.d = max;
	return result;
}

double TMOThompson02::getDisplayableLuminance(double L, double maxLuminance)
{
	return (L * (1 + L / (maxLuminance * maxLuminance))) / (1 + L);
}

void TMOThompson02::mapLuminance(cv::Mat &luminanceMat, double maxLuminance)
{
	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double L = luminanceMat.at<double>(x, y);
			double Ld = getDisplayableLuminance(L, maxLuminance);
			luminanceMat.at<double>(x, y) = Ld;
		}
	}
}

int TMOThompson02::Transform()
{
	TMOThompson02::matAndDouble result = getLuminanceMat();
	cv::Mat luminanceMat = result.mat;
	double maxLuminance = result.d;

	mapLuminance(luminanceMat, maxLuminance);

	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);

	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *srcPixel = pSrc->GetPixel(x, y);
			double *dstPixel = pDst->GetPixel(x, y);

			dstPixel[0] = luminanceMat.at<double>(x, y);
			dstPixel[1] = srcPixel[1];
			dstPixel[2] = srcPixel[2];
		}
	}

	pDst->Convert(TMO_RGB);

	return 0;
}
