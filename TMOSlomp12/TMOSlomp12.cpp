/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   Fast Local Tone Mapping, Summed-Area Tables                *
 *                          and Mesopic Vision Simulation                       *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

/**
 * @file TMOSlomp12.cpp
 * @brief Fast Local Tone Mapping, Summed-Area Tables and Mesopic Vision Simulation
 * @author Jan Findra
 * @class TMOSlomp12.cpp
 */

#include "TMOSlomp12.h"

/**
 * @brief Constructor
 */
TMOSlomp12::TMOSlomp12()
{
	SetName(L"Slomp12");
	SetDescription(L"Fast Local Tone Mapping, Summed-Area Tables and Mesopic Vision Simulation");

	local.SetName(L"local");
	local.SetDescription(L"Turn on or off local or global luminance compression (tone mapping).");
	local.SetDefault(true);
	local = true;
	this->Register(local);

	mesopic.SetName(L"mesopic");
	mesopic.SetDescription(L"Turn on or off the filter which simulates mesopic vision color change.");
	mesopic.SetDefault(true);
	mesopic = true;
	this->Register(mesopic);

	varying.SetName(L"varying");
	varying.SetDescription(L"Choose between spatially-varying and spatially-uniform mesopic vision reproduction operator.");
	varying.SetDefault(true);
	varying = true;
	this->Register(varying);
}

/**
 * @brief Destructor
 */
TMOSlomp12::~TMOSlomp12()
{
}

/**
 * @brief Function to convert the image to logarithmic luminance matrix
 * @return cv::Mat: logarithmic luminance matrix
 */
cv::Mat TMOSlomp12::TMOImageToLogLuminanceMat()
{
	cv::Mat logLuminanceMat(pSrc->GetWidth(), pSrc->GetHeight(), CV_64FC1);
	double delta = 0.00001; // small constant to avoid log(0)

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

/**
 * @brief Function to log the luminance image
 * @param luminanceMat pointer to the luminance matrix
 */
void TMOSlomp12::logLuminanceImage(cv::Mat *luminanceMat)
{
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pSrcPixel = pSrc->GetPixel(x, y);
			double *pDstPixel = pDst->GetPixel(x, y);
			double value = luminanceMat->at<double>(x, y);
			// luminance is changed, chromaticity is preserved
			pDstPixel[0] = value;
			pDstPixel[1] = pSrcPixel[1];
			pDstPixel[2] = pSrcPixel[2];
		}
	}
	// convert to black and white image
	pDst->Convert(TMO_RGB);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pPixel = pDst->GetPixel(x, y);
			// simple conversion to black and white
			double value = 0.299 * pPixel[0] + 0.587 * pPixel[1] + 0.114 * pPixel[2];
			pPixel[0] = value;
			pPixel[1] = value;
			pPixel[2] = value;
		}
	}
}

double TMOSlomp12::fullMipmap(cv::Mat *mat)
{
	// compute the number of levels
	int levels = (int)log2(std::min(mat->cols, mat->rows));
	cv::Mat mipmapMat = mat->clone();

	for (int i = 0; i < levels; i++)
	{
		// downsample the image
		cv::pyrDown(mipmapMat, mipmapMat);
	}

	// check if the mipmap matrix is a single cell, if multiple cells, compute the average value
	double avg;
	int numberOfCells = mipmapMat.rows * mipmapMat.cols;

	if (numberOfCells == 1) // single cell
	{
		avg = mipmapMat.at<double>(0, 0);
	}
	else // multiple cells
	{
		double sum = 0.;
		for (int y = 0; y < mipmapMat.cols; y++)
		{
			for (int x = 0; x < mipmapMat.rows; x++)
			{
				sum += mipmapMat.at<double>(x, y);
			}
		}
		avg = sum / numberOfCells;
	}

	return avg;
}

/**
 * @brief Function to scale the luminance matrix
 * @param luminanceMat pointer to the luminance matrix
 * @param keyValue key value
 */
void TMOSlomp12::scaleLuminance(cv::Mat *luminanceMat, double keyValue)
{
	for (int y = 0; y < luminanceMat->cols; y++)
	{
		for (int x = 0; x < luminanceMat->rows; x++)
		{
			luminanceMat->at<double>(x, y) = pSrc->GetPixel(x, y)[1] * (alpha / keyValue);
		}
	}
}

/**
 * @brief Function to scale the luminance image
 * @param luminanceMat pointer to the luminance matrix
 */
void TMOSlomp12::scaledLuminanceImage(cv::Mat *luminanceMat)
{
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pSrcPixel = pSrc->GetPixel(x, y);
			double *pDstPixel = pDst->GetPixel(x, y);
			// map relative luminance to normalized displayable range <0, 1)
			double value = luminanceMat->at<double>(x, y) / (1 + luminanceMat->at<double>(x, y));
			// luminance is changed, chromaticity is preserved
			pDstPixel[0] = value;
			pDstPixel[1] = pSrcPixel[1];
			pDstPixel[2] = pSrcPixel[2];
		}
	}
	// convert to black and white image
	pDst->Convert(TMO_RGB);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pPixel = pDst->GetPixel(x, y);
			// simple conversion to black and white
			double value = 0.299 * pPixel[0] + 0.587 * pPixel[1] + 0.114 * pPixel[2];
			pPixel[0] = value;
			pPixel[1] = value;
			pPixel[2] = value;
		}
	}
}

/**
 * @brief Function to compute the box filter
 * @param SAT pointer to the summed-area table
 * @param x x coordinate
 * @param y y coordinate
 * @param s scale
 * @return double: box filter value
 */
double TMOSlomp12::boxFilter(cv::Mat *SAT, int x, int y, int s)
{
	// compute the corners of the box

	int x0 = x - s - 1;
	int y0 = y - s - 1;
	// the elements out of the bottom or right edge of the image are set to the closest edge element, clamp-to-edge
	int x1 = std::min(SAT->rows - 1, x + s);
	int y1 = std::min(SAT->cols - 1, y + s);
	double A, B, C, D;

	A = SAT->at<double>(x1, y1);

	// the elements out of the top or left edge of the image are set to 0
	// first row and first column are skipped because they are 0 by the function cv::integral
	if (x0 < 1)
		B = 0;
	else
		B = SAT->at<double>(x0, y1);

	if (y0 < 1)
		C = 0;
	else
		C = SAT->at<double>(x1, y0);

	if (x0 < 1 || y0 < 1)
		D = 0;
	else
		D = SAT->at<double>(x0, y0);

	// compute the area of the box
	x0 = std::max(1, x0);
	y0 = std::max(1, y0);
	double area = (x1 - x0) * (y1 - y0);

	// compute the mean value of the box
	return (A - B - C + D) / area;
}

/**
 * @brief Function to compute the normalized difference
 * @param conv0 convolution 0
 * @param conv1 convolution 1
 * @param s scale
 * @return double: normalized difference
 */
double TMOSlomp12::getNormalizedDifference(double conv0, double conv1, int s)
{
	return std::abs((conv0 - conv1) / (std::pow(2, phi) * (alpha / std::pow(s, 2)) + conv0));
}

/**
 * @brief Function to get the maximum scale
 * @param SAT pointer to the summed-area table
 * @param x x coordinate
 * @param y y coordinate
 * @return int: maximum scale
 */
int TMOSlomp12::getMaxScale(cv::Mat *SAT, int x, int y)
{
	// default values
	int s = 0;
	double convolution0 = SAT->at<double>(x, y);
	double normalizedDifference = 0.;
	int maxScale = 0;

	int sLimit = std::max(std::max(x, y), std::max(SAT->rows - x, SAT->cols - y));

	// find the maximum scale while the normalized difference is smaller than epsilon
	while ((normalizedDifference < epsilon) && (s < sLimit))
	{
		maxScale = s;
		s++;
		double convolution1 = boxFilter(SAT, x, y, s);
		normalizedDifference = getNormalizedDifference(convolution0, convolution1, s);
		convolution0 = convolution1;
	}

	return maxScale;
}

/**
 * @brief Function to compute the red response value
 * @param illuminance illuminance
 * @return double: red response value
 */
double TMOSlomp12::redResponseValue(double illuminance)
{
	return 70 / (1 + pow(10 / illuminance, 0.383)) + 22;
}

/**
 * @brief Function to compute the arithmetic luminance average
 * @return double: arithmetic luminance average
 */
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

/**
 * @brief Function to apply the tone mapping operator
 * @return int: 0 = success, 1 = error
 */
int TMOSlomp12::Transform()
{
	pSrc->Convert(TMO_XYZ);
	// get the logarithmic luminance matrix
	cv::Mat luminanceMat = TMOImageToLogLuminanceMat();

	// uncomment the following lines to see the logarithmic luminance image
	// logLuminanceImage(&luminanceMat);
	// return 0;

	// full mipmap of the logarithmic luminance matrix
	double keyValue = fullMipmap(&luminanceMat);

	// exponentiate the key value
	keyValue = exp(keyValue);

	// scale the luminance matrix
	scaleLuminance(&luminanceMat, keyValue);

	// uncomment the following lines to see the scaled luminance image
	// scaledLuminanceImage(&luminanceMat);
	// return 0;

	if (!local) // global luminance compression
	{
		for (int y = 0; y < luminanceMat.cols; y++)
		{
			for (int x = 0; x < luminanceMat.rows; x++)
			{
				luminanceMat.at<double>(x, y) = luminanceMat.at<double>(x, y) / (1 + luminanceMat.at<double>(x, y));
			}
		}
	}
	else // local luminance compression
	{
		// find the average value of the luminance matrix to be subtracted to make the SAT non-monotonic

		double averageValue = fullMipmap(&luminanceMat);

		// subtract the average value from the luminance matrix
		for (int x = 0; x < luminanceMat.rows; x++)
		{
			for (int y = 0; y < luminanceMat.cols; y++)
			{
				luminanceMat.at<double>(x, y) = luminanceMat.at<double>(x, y) - averageValue;
			}
		}

		cv::Mat SAT;
		// generate the summed-area table (SAT), this function adds first column and row of zeros
		cv::integral(luminanceMat, SAT, CV_64F);

		for (int y = 0; y < luminanceMat.cols; y++)
		{
			for (int x = 0; x < luminanceMat.rows; x++)
			{
				// find the maximum scale and apply the box filter
				int maxScale = getMaxScale(&SAT, x + 1, y + 1);
				double convolution = boxFilter(&SAT, x + 1, y + 1, maxScale);
				luminanceMat.at<double>(x, y) = (luminanceMat.at<double>(x, y) + averageValue) / (1 + convolution + averageValue);
			}
		}
	}

	if (!mesopic) // no mesopic vision simulation
	{
		// copy the source image to the destination image
		pSrc->Convert(TMO_RGB);
		pDst->Convert(TMO_RGB);
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			for (int x = 0; x < pSrc->GetWidth(); x++)
			{
				double *srcPixel = pSrc->GetPixel(x, y);
				double *dstPixel = pDst->GetPixel(x, y);
				dstPixel[0] = srcPixel[0];
				dstPixel[1] = srcPixel[1];
				dstPixel[2] = srcPixel[2];
			}
		}
	}
	else // mesopic vision simulation
	{
		pSrc->Convert(TMO_RGB);
		pSrc->Convert(TMO_LAB);
		pDst->Convert(TMO_LAB);
		// value for normalization of the red response
		double mesopicLightness = redResponseValue(10);

		if (!varying) // spatially-uniform mesopic vision reproduction operator
		{
			// arithmetic luminance average is used as the absolute local area luminance
			double arithLuminanceAvg = arithLuminanceAverage();
			double coefficientRo = redResponseValue(arithLuminanceAvg) / mesopicLightness;

			for (int y = 0; y < pSrc->GetHeight(); y++)
			{
				for (int x = 0; x < pSrc->GetWidth(); x++)
				{
					// change the a value of the LAB color space, rest is preserved
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
		else // spatially-varying mesopic vision reproduction operator
		{
			for (int y = 0; y < pSrc->GetHeight(); y++)
			{
				for (int x = 0; x < pSrc->GetWidth(); x++)
				{
					// absolute local area luminance is used as the absolute local area luminance
					double absoluteLocalAreaLuminance = luminanceMat.at<double>(x, y) * (keyValue / alpha);
					double coefficientRo = redResponseValue(absoluteLocalAreaLuminance) / mesopicLightness;
					// change the a value of the LAB color space, rest is preserved
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
	}

	// apply the luminance compression, tone mapping
	pDst->Convert(TMO_RGB);
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			// global or local luminance compression is saved in the luminance matrix
			pDst->GetPixel(x, y)[0] = luminanceMat.at<double>(x, y);
		}
	}
	pDst->Convert(TMO_RGB);

	return 0;
}
