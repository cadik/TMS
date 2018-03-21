/* --------------------------------------------------------------------------- *
 * TMOAubry14.cpp: implementation of the TMOAubry14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAubry14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAubry14::TMOAubry14()
{
	SetName(L"Aubry14");						// TODO - Insert operator name
	SetDescription(L"Tone mapping using fast local Laplacian filters");	// TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(0.0,1.0);				// TODO - Add acceptable range if needed
	// this->Register(dParameter);
}

TMOAubry14::~TMOAubry14()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAubry14::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	// pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	// pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	int height = pSrc->GetHeight();
	int width  = pSrc->GetWidth();

	double red, green, blue, gray;

	cv::Mat rgbMat, grayMat;
	rgbMat = cv::Mat::zeros (height, width, CV_32FC3);
	grayMat = cv::Mat::zeros (height, width, CV_32FC1);


	// convert to grayscale
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// need to store rgb in mat to calculate ratio later
			rgbMat.at<cv::Vec3f>(j,i)[0] = *pSourceData++;
			rgbMat.at<cv::Vec3f>(j,i)[1] = *pSourceData++;
			rgbMat.at<cv::Vec3f>(j,i)[2] = *pSourceData++;
		}
	}

	// cvtColor handles max 32b floats
	cv::cvtColor(rgbMat, grayMat, CV_RGB2GRAY);

	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// store results to the destination image
			*pDestinationData++ = grayMat.at<float>(j,i);
			*pDestinationData++ = grayMat.at<float>(j,i);
			*pDestinationData++ = grayMat.at<float>(j,i);
		}
	}

	// calculate ratio for converting to rgb at the end
	// ...

	// calculate LLF
	// ...

	// multiply result with ratio
	// ...

	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}
