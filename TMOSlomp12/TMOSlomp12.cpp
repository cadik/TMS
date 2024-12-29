/* --------------------------------------------------------------------------- *
 * TMOSlomp12.cpp: implementation of the TMOSlomp12 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOSlomp12.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOSlomp12::TMOSlomp12()
{
	SetName(L"Slomp12");
	SetDescription(L"Fast Local Tone Mapping, Summed-Area Tables and Mesopic Vision Simulation");

	local.SetName(L"local");
	local.SetDescription(L"Turn on or off local or global luminance compression (tone mapping).");
	local.SetDefault(false);
	local = false;
	this->Register(local);

	mesopic.SetName(L"mesopic");
	mesopic.SetDescription(L"Turn on or off the filter which simulates mesopic vision color change.");
	mesopic.SetDefault(true);
	mesopic = true;
	this->Register(mesopic);

	varying.SetName(L"varying");
	varying.SetDescription(L"Choose between spatially-varying and spatially-uniform mesopic vision reproduction operator.");
	varying.SetDefault(false);
	varying = false;
	this->Register(varying);
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

void TMOSlomp12::logLuminanceImage(cv::Mat luminanceMat)
{
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pSrcPixel = pSrc->GetPixel(x, y);
			double *pDstPixel = pDst->GetPixel(x, y);
			double value = luminanceMat.at<double>(x, y);
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

cv::Mat TMOSlomp12::mipmap(cv::Mat mat, int levels)
{
	if (levels == -1)
	{
		levels = (int)log2(std::min(mat.cols, mat.rows));
	}
	cv::Mat mipmapMat = mat.clone();
	for (int i = 0; i < levels; i++)
	{
		cv::pyrDown(mipmapMat, mipmapMat);
	}
	return mipmapMat;
}

void TMOSlomp12::scaleLuminance(cv::Mat *luminanceMat, double keyValue)
{
	for (int i = 0; i < luminanceMat->rows; i++)
	{
		for (int j = 0; j < luminanceMat->cols; j++)
		{
			luminanceMat->at<double>(i, j) = pSrc->GetPixel(i, j)[1] * (alpha / keyValue);
		}
	}
}

void TMOSlomp12::scaledLuminanceImage(cv::Mat luminanceMat)
{
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pSrcPixel = pSrc->GetPixel(x, y);
			double *pDstPixel = pDst->GetPixel(x, y);
			double value = luminanceMat.at<double>(x, y) / (1 + luminanceMat.at<double>(x, y));
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

double TMOSlomp12::boxFilter(cv::Mat *SAT, int x, int y, int s)
{
	int x0 = std::max(0, x - s);
	int y0 = std::max(0, y - s);
	int x1 = std::min(SAT->cols - 1, x + s);
	int y1 = std::min(SAT->rows - 1, y + s);

	double sum = SAT->at<double>(x1, y1) - SAT->at<double>(x0, y1) - SAT->at<double>(x1, y0) + SAT->at<double>(x0, y0);
	return sum / std::pow(2 * s + 1, 2);
}

double TMOSlomp12::getNormalizedDifference(double conv0, double conv1, int s)
{
	return (conv0 - conv1) / (std::pow(phi, 2) * (alpha / std::pow(s, 2)) + conv0);
}

int TMOSlomp12::getMaxScale(cv::Mat *SAT, int x, int y)
{
	int s = 0;
	double convolution0 = SAT->at<double>(x, y);
	double normalizedDifference = 0.;
	int maxScale = 0;

	while (normalizedDifference < epsilon)
	{
		maxScale = s;
		s++;
		double convolution1 = boxFilter(SAT, x, y, s);
		normalizedDifference = getNormalizedDifference(convolution0, convolution1, s);
		convolution0 = convolution1;
	}

	return maxScale;
}

double TMOSlomp12::redResponseValue(double illuminance)
{
	return 70 / (1 + pow(10 / illuminance, 0.383)) + 22;
}

double TMOSlomp12::arithLuminanceAverage()
{
	double sum = .0;
	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			sum += pSrc->GetPixel(x, y)[1];
		}
	}
	return sum / (pSrc->GetHeight() * pSrc->GetWidth());
}

int TMOSlomp12::Transform()
{
	pSrc->Convert(TMO_XYZ);
	cv::Mat luminanceMat = TMOImageToLogLuminanceMat();

	// uncomment the following lines to see the logarithmic luminance image
	// logLuminanceImage(luminanceMat);
	// return 0;

	cv::Mat mipmapMat = mipmap(luminanceMat, -1);

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

	// scaled luminances are stored in luminanceMat
	scaleLuminance(&luminanceMat, keyValue);

	// uncomment the following lines to see the scaled luminance image
	// scaledLuminanceImage(luminanceMat);
	// return 0;

	if (!local)
	{
		for (int y = 0; y < luminanceMat.cols; y++)
		{
			for (int x = 0; x < luminanceMat.rows; x++)
			{
				luminanceMat.at<double>(x, y) = luminanceMat.at<double>(x, y) / (1 + luminanceMat.at<double>(x, y));
			}
		}
	}

	if (local)
	{
		mipmapMat = mipmap(luminanceMat, -1);
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

		for (int x = 0; x < luminanceMat.rows; x++)
		{
			for (int y = 0; y < luminanceMat.cols; y++)
			{
				luminanceMat.at<double>(x, y) = luminanceMat.at<double>(x, y) - averageValue;
			}
		}

		cv::Mat SAT;
		cv::integral(luminanceMat, SAT, CV_64F);

		for (int y = 0; y < luminanceMat.cols; y++)
		{
			for (int x = 0; x < luminanceMat.rows; x++)
			{
				int maxScale = getMaxScale(&SAT, x, y);
				double convolution = boxFilter(&SAT, x, y, maxScale);
				luminanceMat.at<double>(x, y) = luminanceMat.at<double>(x, y) / (1 + convolution) + averageValue;
			}
		}
	}

	double mesopicLightness = redResponseValue(10);

	if (mesopic && !varying) // Spatially-uniform mesopic vision reproduction operator
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
	}

	// apply the luminance compression, tone mapping
	pDst->Convert(TMO_RGB);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			pDst->GetPixel(x, y)[0] = luminanceMat.at<double>(x, y);
		}
	}
	pDst->Convert(TMO_RGB);

	return 0;
}
