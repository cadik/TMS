/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                               *
*                                                                              *
*                       Bachelor thesis                                        *
*                       Author: Martin Molek [xmolek00 AT stud.fit.vutbr.cz]   *
*                       Brno 2017                                              *
*                                                                              *
*                       Fast Multispectral2Gray                                *
*                       Implementation of the Ali Alsam, Mark S. Drew          *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOAlsam09.cpp
 * @brief Fast Multispectral2Gray
 * @author Martin Molek
 * @class TMOALsam09 
 */

#include "TMOAlsam09.h"
#include <fftw3.h>

/**
  *  @brief Constructor
  */
TMOAlsam09::TMOAlsam09()
{
	SetName(L"Alsam09");
	SetDescription(L"Add your TMO description here");

	dParameter.SetName(L"ParameterName");
	dParameter.SetDescription(L"ParameterDescription");
	dParameter.SetDefault(1);
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0);
	this->Register(dParameter);
}

/**
  *  @brief Destructor
  */
TMOAlsam09::~TMOAlsam09()
{
}

/**
  *  @brief Max channel (R=>0, G=>1, B=>2)
  *  @param x Actual position (width)
  *  @param  y Actual position (height)
  *  @param  maxy Max height
  *  @param pStartSourceData Source data
  *
  *  @return Returns 0 if R>G and R>B, return 1 if G>R and G>B, else returns 2
  */
int maxChannel(int x, int y, int maxy, double *pStartSourceData)
{
	if ((pStartSourceData[3 * (y * maxy + x)] > pStartSourceData[3 * (y * maxy + x) + 1]) && (pStartSourceData[3 * (y * maxy + x)] > pStartSourceData[3 * (y * maxy + x) + 2])) //R>G && R>B
		return 0;
	if ((pStartSourceData[3 * (y * maxy + x) + 1] > pStartSourceData[3 * (y * maxy + x)]) && (pStartSourceData[3 * (y * maxy + x) + 1] > pStartSourceData[3 * (y * maxy + x) + 2])) //G>R && G>B
		return 1;
	return 2;
}

/**
  *  @brief Max pixel
  *  @param x Actual position (width)
  *  @param y Actual position (height)
  *  @param maxy Max height
  *  @param pStartSourceData Source data
  *
  *  @return Returns biggest value of pixel
  */
double pixelMax(int x, int y, int maxy, double *pStartSourceData)
{
	return std::max(std::max(pStartSourceData[3 * (y * maxy + x)], pStartSourceData[3 * (y * maxy + x) + 1]), pStartSourceData[3 * (y * maxy + x) + 2]);
}

/**
  *  @brief Sum of max pixels
  *  @param x Actual position (width)
  *  @param y Actual position (height)
  *  @param maxx Max width
  *  @param maxy Max height
  *  @param pStartSourceData Source data
  *
  *  @return Returns sum of max pixels
  */
double pixelMaxSum(int x, int y, int maxx, int maxy, double *pStartSourceData)
{
	double sum = 0.0;
	int outOfPicture = 0;

	if (x - 1 > 0)
		sum += std::max(std::max(pStartSourceData[3 * (y * maxy + x - 1)], pStartSourceData[3 * (y * maxy + x - 1) + 1]), pStartSourceData[3 * (y * maxy + x - 1) + 2]);
	else
		outOfPicture++;
	if (y - 1 > 0)
		sum += std::max(std::max(pStartSourceData[3 * ((y - 1) * maxy + x)], pStartSourceData[3 * ((y - 1) * maxy + x) + 1]), pStartSourceData[3 * ((y - 1) * maxy + x) + 2]);
	else
		outOfPicture++;
	if (x + 1 < maxx)
		sum += std::max(std::max(pStartSourceData[3 * (y * maxy + x + 1)], pStartSourceData[3 * (y * maxy + x + 1) + 1]), pStartSourceData[3 * (y * maxy + x + 1) + 2]);
	else
		outOfPicture++;
	if (y + 1 < maxy)
		sum += std::max(std::max(pStartSourceData[3 * ((y + 1) * maxy + x)], pStartSourceData[3 * ((y + 1) * maxy + x) + 1]), pStartSourceData[3 * ((y + 1) * maxy + x) + 2]);
	else
		outOfPicture++;
	sum += outOfPicture * std::max(std::max(pStartSourceData[3 * (y * maxy + x)], pStartSourceData[3 * (y * maxy + x) + 1]), pStartSourceData[3 * (y * maxy + x) + 2]);
	return sum;
}

/**
  *  @brief Sum of max pixels minus median
  *  @param x Actual position (width)
  *  @param y Actual position (height)
  *  @param maxx Max width
  *  @param maxy Max height
  *  @param channel Channel
  *  @param pStartSourceData Source data
  *
  *  @return Returns sum of gradients
  */
double gradSum(int x, int y, int maxx, int maxy, int channel, double *pStartSourceData)
{
	double sum = 0.0;
	double me = pStartSourceData[3 * (y * maxy + x) + channel];
	int outOfPicture = 0;
	if (x - 1 > 0)
		sum += pStartSourceData[3 * (y * maxy + x - 1) + channel] - me;
	if (y - 1 > 0)
		sum += pStartSourceData[3 * ((y - 1) * maxy + x) + channel] - me;
	if (x + 1 < maxx)
		sum += pStartSourceData[3 * (y * maxy + x + 1) + channel] - me;
	if (y + 1 < maxy)
		sum += pStartSourceData[3 * ((y + 1) * maxy + x) + channel] - me;
	return sum;
}

/**
  *  @brief Fast Multispectral2gray
  */
int TMOAlsam09::Transform()
{

	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();
	double r, g, b, shade;
	int mC, i;
	int size = pSrc->GetHeight() * pSrc->GetWidth() * 3;
	double *tmp1 = fftw_alloc_real(size);
	double *tmp2 = fftw_alloc_real(size);
	for (i = 0; i < size; i++)
		tmp2[i] = pSourceData[i];

	for (int omega = 0; omega < 5; omega++)
	{
		for (i = 0; i < size; i++)
			tmp1[i] = tmp2[i];
		i = 0;
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			for (int x = 0; x < pSrc->GetWidth(); x++)
			{

				mC = maxChannel(x, y, pSrc->GetWidth(), pSrc->GetData());
				tmp2[i + mC] = tmp1[i + mC] + gradSum(x, y, pSrc->GetHeight(), pSrc->GetWidth(), mC, tmp1) / 4.0;

				// if (x >	303 && y == 303)
				// std::cerr << "GS: " << gradSum(x, y, pSrc->GetHeight(), pSrc->GetWidth(), mC, tmp1) << std::endl;
				// shade = r * 0.299 + g * 0.587 + b * 0.114;

				if (tmp2[i + mC] > 1.0)
				{
					tmp2[i + mC] = 1.0;
				}
				else if (tmp2[i + mC] < 0.0)
				{
					tmp2[i + mC] = 0.0;
				}
				i += 3;
			}
			pSrc->ProgressBar(y, pSrc->GetHeight());
		}
	}

	for (i = 0; i < size; i += 3)
	{
		shade = tmp2[i] * 0.299 + tmp2[i + 1] * 0.587 + tmp2[i + 2] * 0.114;
		pDestinationData[i] = shade;
		pDestinationData[i + 1] = shade;
		pDestinationData[i + 2] = shade;
	}

	//pDst->Convert(TMO_RGB);
	fftw_free(tmp1);
	fftw_free(tmp2);
	return 0;
}
