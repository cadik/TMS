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

	sigmaBlur.SetName(L"sigmaBlur");
	sigmaBlur.SetDescription(L"Blurring parameter");
	sigmaBlur.SetDefault(1.6);
	sigmaBlur.SetRange(0.0, 10.0);
	this->Register(sigmaBlur);

	gammaEdge.SetName(L"gammaEdge");
	gammaEdge.SetDescription(L"Edge enhancement parameter");
	gammaEdge.SetDefault(1.25);
	gammaEdge.SetRange(1.0, 10.0);
	this->Register(gammaEdge);

	sigmaNoise.SetName(L"sigmaNoise");
	sigmaNoise.SetDescription(L"Additive noise parameter");
	sigmaNoise.SetDefault(0.0125);
	sigmaNoise.SetRange(0.0, 1.0);
	this->Register(sigmaNoise);

	rcf.SetName(L"rcf");
	rcf.SetDescription(L"Range compression factor. Just for manual day-for-night conversion.");
	rcf.SetDefault(1.0);
	rcf.SetRange(0.0, 1.0);
	this->Register(rcf);

	mf.SetName(L"mf");
	mf.SetDescription(L"Mesopic factor. If set, it is used instead of the computed one.");
	mf.SetDefault(0.0);
	mf.SetRange(0.0, 1.0);
	this->Register(mf);
}

TMOThompson02::~TMOThompson02()
{
}

double TMOThompson02::luminanceReduction(double Y, double YLogAvg, double Ymax)
{
	// get key value for luminance reduction
	double alpha = 1.03 - 2 / (2 + std::log10(YLogAvg + 1));
	// compute reduced luminance
	double Yr = (alpha * Y) / YLogAvg;
	// compute final, normalized luminance
	double Yn = (Yr * (1 + (Yr / std::pow(Ymax, 2)))) / (1 + Yr);
	return Yn;
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
	// if mesopic factor is set, return it
	if (mf != 0.0)
	{
		return mf;
	}
	// otherwise compute it based on scotopic luminance
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
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);

	// luminance reduction
	double epsilon = 1e-6;
	double sumLogY = 0.0;
	int pixelCount = pSrc->GetHeight() * pSrc->GetWidth();
	double Ymax = 0.0;

	// compute sum of logarithms of luminance and maximum luminance
	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double Y = pSrc->GetPixel(x, y)[0];
			sumLogY += std::log(Y + epsilon);
			if (Y > Ymax)
			{
				Ymax = Y;
			}
		}
	}

	// compute average luminance
	double YLogAvg = std::exp(sumLogY / pixelCount);

	// go through the image and apply luminance reduction
	for (int y = 0; y < pDst->GetHeight(); y++)
	{
		for (int x = 0; x < pDst->GetWidth(); x++)
		{
			double Y = pSrc->GetPixel(x, y)[0];
			double Yr = luminanceReduction(Y, YLogAvg, Ymax);
			pDst->GetPixel(x, y)[0] = rcf * Yr;
			// copy chromaticity values
			pDst->GetPixel(x, y)[1] = pSrc->GetPixel(x, y)[1];
			pDst->GetPixel(x, y)[2] = pSrc->GetPixel(x, y)[2];
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
			dstPixel[0] = k * V * bluishGreyRGB[0];
			dstPixel[1] = k * V * bluishGreyRGB[1];
			dstPixel[2] = k * V * bluishGreyRGB[2];
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
