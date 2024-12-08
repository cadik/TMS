/* --------------------------------------------------------------------------- *
 * TMOSlomp12.cpp: implementation of the TMOSlomp12 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOSlomp12.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOSlomp12::TMOSlomp12()
{
	SetName(L"Slomp12");							  // TODO - Insert operator name
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
	double delta = 0.00001;

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			*pSourceData++;
			double Y = *pSourceData++;
			*pSourceData++;

			logLuminanceMat.at<double>(y, x) = log10(Y) + delta;
		}
	}

	return logLuminanceMat;
}

void TMOSlomp12::logLuminanceImage(cv::Mat srcMat)
{
	double *pDestinationData = pDst->GetData();
	for (int i = 0; i < pDst->GetHeight(); i++)
	{
		for (int j = 0; j < pDst->GetWidth(); j++)
		{
			*pDestinationData++ = srcMat.at<double>(i, j);
			*pDestinationData++ = srcMat.at<double>(i, j);
			*pDestinationData++ = srcMat.at<double>(i, j);
		}
	}
	pDst->Convert(TMO_RGB);
}

cv::Mat TMOSlomp12::mipmap(cv::Mat srcMat, int levels)
{
	if (levels == -1)
	{
		levels = (int)log2(std::min(srcMat.cols, srcMat.rows));
	}
	cv::Mat dstMat = srcMat.clone();
	for (int i = 0; i < levels; i++)
	{
		cv::pyrDown(dstMat, dstMat);
	}
	return dstMat;
}

int TMOSlomp12::Transform()
{
	cv::Mat srcMat = TMOImageToLogLuminanceMat(pSrc);
	cv::Mat mipmapMat = mipmap(srcMat, -1);
	std::cerr << "Mipmap matrix size: " << mipmapMat.size() << std::endl;
	std::cerr << "Final value: " << mipmapMat.at<double>(0, 0) << std::endl;

	// uncomment the following line to see the logarithmic luminance image
	logLuminanceImage(srcMat);
	return 0;

	return 0;
}
