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

cv::Mat TMOSlomp12::TMOImageToLogLuminanceMat()
{
	cv::Mat logLuminanceMat(pSrc->GetWidth(), pSrc->GetHeight(), CV_64FC1);
	double delta = 0.00001;

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double Y = pSrc->GetPixel(x, y)[1];
			logLuminanceMat.at<double>(x, y) = log10(Y + delta);
		}
	}

	return logLuminanceMat;
}

void TMOSlomp12::logLuminanceImage(cv::Mat srcMat)
{
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pSrcPixel = pSrc->GetPixel(x, y);
			double *pDstPixel = pDst->GetPixel(x, y);
			double value = srcMat.at<double>(x, y);
			pDstPixel[0] = value;
			pDstPixel[1] = pSrcPixel[1];
			pDstPixel[2] = pSrcPixel[2];
		}
	}
	pDst->Convert(TMO_RGB);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pPixel = pDst->GetPixel(x, y);
			double value = 0.299 * pPixel[0] + 0.587 * pPixel[1] + 0.114 * pPixel[2];
			pPixel[0] = value;
			pPixel[1] = value;
			pPixel[2] = value;
		}
	}
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

void TMOSlomp12::scaleLuminance(cv::Mat *srcMat, double keyValue, double alpha)
{
	for (int i = 0; i < srcMat->rows; i++)
	{
		for (int j = 0; j < srcMat->cols; j++)
		{
			srcMat->at<double>(i, j) = pSrc->GetPixel(i, j)[1] * (alpha / keyValue);
		}
	}
}

void TMOSlomp12::scaledLuminanceImage(cv::Mat srcMat)
{
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pSrcPixel = pSrc->GetPixel(x, y);
			double *pDstPixel = pDst->GetPixel(x, y);
			double value = srcMat.at<double>(x, y) / (1 + srcMat.at<double>(x, y));
			pDstPixel[0] = value;
			pDstPixel[1] = pSrcPixel[1];
			pDstPixel[2] = pSrcPixel[2];
		}
	}
	pDst->Convert(TMO_RGB);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pPixel = pDst->GetPixel(x, y);
			double value = 0.299 * pPixel[0] + 0.587 * pPixel[1] + 0.114 * pPixel[2];
			pPixel[0] = value;
			pPixel[1] = value;
			pPixel[2] = value;
		}
	}
}

double TMOSlomp12::redResponseValue(double illuminance)
{
	return 70 / (1 + pow(10 / illuminance, 0.383)) + 22;
}

double TMOSlomp12::arithLuminanceAverage()
{
	double sum = .0;
	for (int i = 0; i < pSrc->GetHeight(); i++)
	{
		for (int j = 0; j < pSrc->GetWidth(); j++)
		{
			sum += pSrc->GetPixel(i, j)[1];
		}
	}
	return sum / (pSrc->GetHeight() * pSrc->GetWidth());
}

int TMOSlomp12::Transform()
{
	pSrc->Convert(TMO_XYZ);
	cv::Mat srcMat = TMOImageToLogLuminanceMat();

	// uncomment the following lines to see the logarithmic luminance image
	// logLuminanceImage(srcMat);
	// return 0;

	cv::Mat mipmapMat = mipmap(srcMat, -1);

	double keyValue;
	int numberOfCells = mipmapMat.rows * mipmapMat.cols;
	if (numberOfCells == 1)
	{
		keyValue = mipmapMat.at<double>(0, 0);
	}
	else
	{
		double sum = .0;
		for (int i = 0; i < mipmapMat.rows; i++)
		{
			for (int j = 0; j < mipmapMat.cols; j++)
			{
				sum += mipmapMat.at<double>(i, j);
			}
		}
		keyValue = sum / numberOfCells;
	}

	keyValue = pow(10, keyValue);

	double alpha = 0.18;
	// scaled luminances are stored in srcMat
	scaleLuminance(&srcMat, keyValue, alpha);

	cv::Mat globalLuminanceOperatorMat = srcMat.clone();

	// uncomment the following lines to see the scaled luminance image
	// scaledLuminanceImage(srcMat);
	// return 0;

	mipmapMat = mipmap(srcMat, -1);
	double averageValue;
	numberOfCells = mipmapMat.rows * mipmapMat.cols;
	if (numberOfCells == 1)
	{
		averageValue = mipmapMat.at<double>(0, 0);
	}
	else
	{
		double sum = .0;
		for (int i = 0; i < mipmapMat.rows; i++)
		{
			for (int j = 0; j < mipmapMat.cols; j++)
			{
				sum += mipmapMat.at<double>(i, j);
			}
		}
		averageValue = sum / numberOfCells;
	}

	for (int x = 0; x < srcMat.rows; x++)
	{
		for (int y = 0; y < srcMat.cols; y++)
		{
			srcMat.at<double>(x, y) = srcMat.at<double>(x, y) - averageValue;
		}
	}

	cv::Mat SAT;
	cv::integral(srcMat, SAT, CV_64F);

	double mesopicLightness = redResponseValue(10);

	if (true) // Spatially-uniform mesopic vision reproduction operator
	{
		double arithLuminanceAvg = arithLuminanceAverage();
		double coefficientRo = redResponseValue(arithLuminanceAvg) / mesopicLightness;

		pSrc->Convert(TMO_RGB);
		pSrc->Convert(TMO_LAB);
		pDst->Convert(TMO_LAB);
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			for (int x = 0; x < pSrc->GetWidth(); x++)
			{
				double *srcPixel = pSrc->GetPixel(x, y);
				double *dstPixel = pDst->GetPixel(x, y);
				if (srcPixel[1] > 0)
				{
					dstPixel[1] = srcPixel[1] * coefficientRo;
				}
				else
				{
					dstPixel[1] = srcPixel[1];
				}
				dstPixel[0] = srcPixel[0];
				dstPixel[2] = srcPixel[2];
			}
		}

		pDst->Convert(TMO_RGB);
		pDst->Convert(TMO_Yxy);
		for (int y = 0; y < pDst->GetHeight(); y++)
		{
			for (int x = 0; x < pDst->GetWidth(); x++)
			{
				pDst->GetPixel(x, y)[0] = globalLuminanceOperatorMat.at<double>(x, y) * (1 + globalLuminanceOperatorMat.at<double>(x, y));
			}
		}
		pDst->Convert(TMO_RGB);
	}

	return 0;
}
