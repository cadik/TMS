/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                      A Spatial Post-Processing Algorithm                     *
 *                         for Images of Night Scenes                           *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2025                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMOThompson02.h"

TMOThompson02::TMOThompson02()
{
	SetName(L"Thompson02");
	SetDescription(L"A Spatial Post-Processing Algorithm for Images of Night Scenes");

	sigmaBlur.SetName(L"SigmaBlur");
	sigmaBlur.SetDescription(L"Blurring parameter");
	sigmaBlur.SetDefault(1.6);
	sigmaBlur.SetRange(0.0, 10.0);
	this->Register(sigmaBlur);

	gammaEdge.SetName(L"GammaEdge");
	gammaEdge.SetDescription(L"Edge enhancement parameter");
	gammaEdge.SetDefault(1.25);
	gammaEdge.SetRange(1.0, 10.0);
	this->Register(gammaEdge);

	sigmaNoise.SetName(L"SigmaNoise");
	sigmaNoise.SetDescription(L"Additive noise parameter");
	sigmaNoise.SetDefault(0.0125);
	sigmaNoise.SetRange(0.0, 1.0);
	this->Register(sigmaNoise);
}

TMOThompson02::~TMOThompson02()
{
}

TMOThompson02::matAndDouble TMOThompson02::getLuminanceMat()
{
	matAndDouble result;
	// maximum luminance
	double max = 0.0;
	// create luminance matrix
	cv::Mat luminanceMat(pSrc->GetWidth(), pSrc->GetHeight(), CV_64F);

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double *pixel = pSrc->GetPixel(x, y);
			// calculate luminance
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

cv::Mat TMOThompson02::getScotopicLuminanceMat()
{
	cv::Mat scotopicLuminanceMat(pDst->GetWidth(), pDst->GetHeight(), CV_64F);

	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pixel = pDst->GetPixel(x, y);
			double X = pixel[0];
			double Y = pixel[1];
			double Z = pixel[2];
			// compute scotopic luminance
			double V = Y * (1.33 * (1 + (Y + Z) / X) - 1.68);
			scotopicLuminanceMat.at<double>(x, y) = V;
		}
	}

	return scotopicLuminanceMat;
}

double TMOThompson02::getMesopicFactor(double L)
{
	double sigma = 100.0;
	return (sigma - 0.25 * L) / (sigma + L);
}

cv::Mat TMOThompson02::TMO2mat()
{
	cv::Mat mat(pDst->GetWidth(), pDst->GetHeight(), CV_64FC3);

	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pixel = pDst->GetPixel(x, y);
			mat.at<cv::Vec3d>(x, y) = cv::Vec3d(pixel[0], pixel[1], pixel[2]);
		}
	}

	return mat;
}

void TMOThompson02::mat2TMO(cv::Mat &input)
{
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *pixel = pDst->GetPixel(x, y);
			cv::Vec3d vec = input.at<cv::Vec3d>(x, y);
			pixel[0] = vec[0];
			pixel[1] = vec[1];
			pixel[2] = vec[2];
		}
	}
}

cv::Mat TMOThompson02::applyGaussianBlur(cv::Mat &input, double sigma)
{
	cv::Mat blurred;
	int kernelSize = static_cast<int>(2 * std::ceil(2 * sigma) + 1);
	cv::GaussianBlur(input, blurred, cv::Size(kernelSize, kernelSize), sigma);

	return blurred;
}

cv::Mat TMOThompson02::applyNightFilter(cv::Mat &input)
{
	cv::Mat blurred = applyGaussianBlur(input, sigmaBlur);
	cv::Mat blurred2 = applyGaussianBlur(input, 1.6 * sigmaBlur);
	cv::Mat diff = blurred - blurred2;

	cv::Mat nightFiltered;
	cv::pow(diff, 1.0 / gammaEdge, diff);
	nightFiltered = blurred2 + diff;

	return nightFiltered;
}

cv::Mat TMOThompson02::addGaussianNoise(cv::Mat &input)
{
	cv::Mat noise = cv::Mat::zeros(input.size(), input.type());
	cv::randn(noise, 0, double(sigmaNoise));
	cv::Mat noisyImage = input + noise;

	return noisyImage;
}

int TMOThompson02::Transform()
{
	// get luminance matrix and maximum luminance
	TMOThompson02::matAndDouble result = getLuminanceMat();
	cv::Mat luminanceMat = result.mat;
	double maxLuminance = result.d;

	// map luminance
	mapLuminance(luminanceMat, maxLuminance);

	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);

	// change luminance in the destination image
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
	pDst->Convert(TMO_XYZ);

	cv::Mat scotopicLuminanceMat = getScotopicLuminanceMat();

	pDst->Convert(TMO_RGB);

	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double *dstPixel = pDst->GetPixel(x, y);
			double V = scotopicLuminanceMat.at<double>(x, y);
			double k = getMesopicFactor(V);

			// transform RGB values to bluish grey
			dstPixel[0] = k * V * bluishGreyRGB[0] + (1 - k) * dstPixel[0];
			dstPixel[1] = k * V * bluishGreyRGB[1] + (1 - k) * dstPixel[1];
			dstPixel[2] = k * V * bluishGreyRGB[2] + (1 - k) * dstPixel[2];
		}
	}

	// convert TMO data to matrix
	cv::Mat imageMat = TMO2mat();
	imageMat = applyNightFilter(imageMat);
	imageMat = addGaussianNoise(imageMat);
	// convert matrix back to TMO data
	mat2TMO(imageMat);

	return 0;
}
