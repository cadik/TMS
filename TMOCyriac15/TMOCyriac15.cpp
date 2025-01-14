/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   A Tone Mapping Operator Based on Neural and                *
 *                    Psychophysical Models of Visual Perception                *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMOCyriac15.h"

TMOCyriac15::TMOCyriac15()
{
	SetName(L"Cyriac15");							  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOCyriac15::~TMOCyriac15()
{
}

cv::Mat TMOCyriac15::normalizeAndLuminance()
{
	cv::Mat luminanceMat(pSrc->GetWidth(), pSrc->GetHeight(), CV_64FC1);

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double *pixel = pSrc->GetPixel(x, y);
			double R = pixel[0];
			double G = pixel[1];
			double B = pixel[2];
			pixel[0] = R / 255.0;
			pixel[1] = G / 255.0;
			pixel[2] = B / 255.0;
			luminanceMat.at<double>(x, y) = 0.2126 * R + 0.7152 * G + 0.0722 * B;
		}
	}

	return luminanceMat;
}

void TMOCyriac15::addToCumulativeHistogram(std::vector<double> *cumulativeHistogram, double value)
{
	for (int i = std::round(value * (bins - 1)); i < cumulativeHistogram->size(); i++)
	{
		cumulativeHistogram->at(i) += 1.0;
	}
}

int TMOCyriac15::Transform()
{
	cv::Mat luminanceMat = normalizeAndLuminance();

	double gammaSys = 0.4;
	double gammaPsy = 0.4;
	double gammaDec = 2.2;
	double previousF = 0.0;
	double difference = 1.0;

	while ((difference > 0.0) && (gammaSys < 2.2))
	{
		gammaSys += 0.1;
		std::vector<double> cumulativeHistogram(bins, 0.0);

		// compute cumulative histogram
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			for (int x = 0; x < pSrc->GetWidth(); x++)
			{
				double L = luminanceMat.at<double>(x, y);
				double Lstar = std::pow(L, gammaSys * gammaPsy);
				addToCumulativeHistogram(&cumulativeHistogram, Lstar);
			}
		}

		// compute sum
		double sum = 0.0;
		for (int i = 0; i < cumulativeHistogram.size(); i++)
		{
			sum += std::pow(cumulativeHistogram[i] - i, 2);
		}

		double F = 1 - std::sqrt(sum) / bins;
		difference = F - previousF;
		previousF = F;
	}

	double gammaEnc = gammaSys / gammaDec;

	return 0;
}
