/* --------------------------------------------------------------------------- *
 * TMOAubry14.cpp: implementation of the TMOAubry14 class.   *
 * Computational Photography - term project
 * Author: Tomas Hudziec, Brno 2018
 * --------------------------------------------------------------------------- */

#include "TMOAubry14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAubry14::TMOAubry14()
{
	SetName(L"Aubry14");
	SetDescription(L"Tone mapping and detail manipulation using fast local Laplacian filters");

	sigmaParameter.SetName(L"sigma");
	sigmaParameter.SetDescription(L"balance between local and global contrast");
	sigmaParameter.SetDefault(0.1);
	sigmaParameter.SetRange(0.01, 1.0);

	NParameter.SetName(L"N");
	NParameter.SetDescription(L"discretisation level");
	NParameter.SetDefault(10);
	NParameter.SetRange(2, 15);

	factParameter.SetName(L"factor");
	factParameter.SetDescription(L"multiply factor");
	factParameter.SetDefault(5);
	factParameter.SetRange(-10, 10);
	
	HDRParameter.SetName(L"HDR");
	HDRParameter.SetDescription(L"check when input image is HDR");
	HDRParameter.SetDefault(false);
	HDRParameter = false;	

	this->Register(HDRParameter);	
	this->Register(factParameter);
	this->Register(NParameter);
	this->Register(sigmaParameter);
}

TMOAubry14::~TMOAubry14()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of the tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAubry14::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	pSrc->Convert(TMO_RGB);

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	int height = pSrc->GetHeight();
	int width  = pSrc->GetWidth();

	cv::Mat I_RGB(height, width, CV_64FC3);
	cv::Mat I_Gray(height, width, CV_64FC1);

	double R, G, B;

	// Convert to grayscale
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			// need to store rgb in mat to calculate colour ratio later
			I_RGB.at<cv::Vec3d>(j,i)[0] = R = *pSourceData++;
			I_RGB.at<cv::Vec3d>(j,i)[1] = G = *pSourceData++;
			I_RGB.at<cv::Vec3d>(j,i)[2] = B = *pSourceData++;
			// convert to grayscale
			I_Gray.at<double>(j,i) = (0.2989*R + 0.5870*G + 0.1140*B);
		}
	}

	// make image of size of 2^n for pyramid functions
	int correctionWidth, correctionHeight;
	correctionWidth = std::ceil(log2(width));
	correctionHeight = std::ceil(log2(height));
	correctionHeight = std::pow(2, correctionHeight);
	correctionWidth = std::pow(2, correctionWidth);
	// resizing image to 2^n dimensions, borders are interpolated
	cv::copyMakeBorder(I_RGB, I_RGB,
		0, correctionHeight-height,
		0, correctionWidth-width,
		cv::BORDER_DEFAULT);
	cv::copyMakeBorder(I_Gray, I_Gray,
		0, correctionHeight-height,
		0, correctionWidth-width,
		cv::BORDER_DEFAULT);

	// calculate ratio for converting to rgb at the end
	cv::Mat I_ratio, I_gray_3c;
	cv::Mat grayChannels[] = {I_Gray, I_Gray, I_Gray};
	cv::merge(grayChannels, 3, I_gray_3c);
	cv::divide(I_RGB, I_gray_3c, I_ratio, 1, -1);

	// the method works with luminance part of image
	cv::Mat I = I_Gray;
	
	double eps = 1e-10;

	// convert HDR image to logarithmic domain
	if (HDRParameter) {
		cv::log(I + eps, I);
	}

	cv::normalize(I, I, 0.0, 1.0, cv::NORM_MINMAX, CV_64FC1);

	// The algorithm of Local Laplacian Filters follows

	// Build Gaussian pyramid
	int pyrLevels = std::ceil(log(std::min(height, width))-log(2))+2;
	std::vector<cv::Mat> inGaussianPyr;
	// 1.level is the image itself
	inGaussianPyr.push_back(I);
	cv::Mat GaussImg;
	for (size_t n = 1; n < pyrLevels; n++) {
		cv::pyrDown(inGaussianPyr[n-1], GaussImg);
		inGaussianPyr.push_back(GaussImg);
	}
	GaussImg.release();

	// Build Laplacian pyramid from Gaussian one
	// the last level is the same as last level of gaussian pyramid
	std::vector<cv::Mat> outLaplacePyr;
	outLaplacePyr.push_back(inGaussianPyr.back());
	cv::Mat smallerUpsampledGauss, LaplaceImg;
	for (size_t n = pyrLevels - 1; n > 0; n--) {
		cv::pyrUp(inGaussianPyr[n], smallerUpsampledGauss);
		cv::subtract(inGaussianPyr[n-1], smallerUpsampledGauss, LaplaceImg);
		outLaplacePyr.insert(outLaplacePyr.begin(), LaplaceImg);
	}
	LaplaceImg.release();  // necessary for later usage of LaplaceImg!

	double sigma = sigmaParameter.GetDouble();
	double fact = factParameter.GetDouble();
	int N = NParameter.GetInt();

	std::vector<double> discretisation = this->linspace(0, 1, N);
	double discretisationStep = discretisation[1];

	cv::Mat I_remap(I.size(), CV_64FC1);


	// main loop of the algorithm
	for (auto ref : discretisation) {
		// calculate I_remap
		for (j = 0; j < I_remap.rows; j++) {
			pSrc->ProgressBar(j, I_remap.rows);	// provide progress bar
			for (int i = 0; i < I_remap.cols; i++) {
				double pixI = I.at<double>(j,i);
				I_remap.at<double>(j,i) =
				fact*(pixI-ref)*exp(-(pixI-ref)*(pixI-ref)/(2.0*sigma*sigma));
			}
		}

		// Build temporary Laplacian pyramid
		std::vector<cv::Mat> tmpLaplacePyr;
		cv::Mat down, up;
		cv::Mat current = I_remap.clone();
		for (size_t n = 0; n < pyrLevels - 1; n++) {
			// apply low pass filter, and downsample
			cv::pyrDown(current, down);
			// in each level, store difference between image and upsampled low pass version
			cv::pyrUp(down, up);
			cv::subtract(current, up, LaplaceImg);
			tmpLaplacePyr.push_back(LaplaceImg);
			// continue with low pass image
			current = down;
		}
		// the coarest level contains the residual low pass image
		tmpLaplacePyr.push_back(current);

		down.release(); up.release(); current.release();
		LaplaceImg.release();

		// compute output Laplace pyramid
		for (size_t level = 0; level < pyrLevels - 1; level++) {
			for (j = 0; j < outLaplacePyr[level].rows; j++) {
				pSrc->ProgressBar(j, outLaplacePyr[level].rows);	// provide progress bar
				for (int i = 0; i < outLaplacePyr[level].cols; i++) {
					double pixInGaussPyr = inGaussianPyr[level].at<double>(j,i);
					double absDiff = abs(pixInGaussPyr - ref);
					if (absDiff < discretisationStep) {
						outLaplacePyr[level].at<double>(j,i) +=
						tmpLaplacePyr[level].at<double>(j,i)*
						(1-absDiff/discretisationStep);
					}
				}
			}
		}

	}// main loop of the algorithm

	// Reconstruct laplacian pyramid
	// start with low pass residual
	cv::Mat I_result_gray = outLaplacePyr.back();
	for (int lev = pyrLevels - 2; lev >= 0; --lev) {
		// upsample, and add to current level
		cv::pyrUp(I_result_gray, I_result_gray);
		I_result_gray += outLaplacePyr[lev];
	}

	// get HDR image from logarithmic domain
	if (HDRParameter) {
		cv::exp(I_result_gray, I_result_gray);
		I_result_gray -= eps;
	}

	// shift image values to positive
	cv::normalize(I_result_gray, I_result_gray, 0, 1, cv::NORM_MINMAX, CV_64FC1);

	// multiply result with ratio to get colours back
	cv::Mat resultChannels[] = {I_result_gray, I_result_gray, I_result_gray};
	cv::Mat I_result_gray_3c, I_result_RGB;
	cv::merge(resultChannels, 3, I_result_gray_3c);
	cv::multiply(I_result_gray_3c, I_ratio, I_result_RGB);

	// normalize to 0-255 range for display
	cv::normalize(I_result_RGB, I_result_RGB, 0, 255, cv::NORM_MINMAX, CV_64FC3);

	// for tone mapping, gamma correct linear intensities for display
	if (HDRParameter) {
		cv::pow(I_result_RGB, 1/2.2, I_result_RGB);
	}

	// output result
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			// put result to output, taking only the image itself, size correction is discarded
			*pDestinationData++ = I_result_RGB.at<cv::Vec3d>(j,i)[0];
			*pDestinationData++ = I_result_RGB.at<cv::Vec3d>(j,i)[1];
			*pDestinationData++ = I_result_RGB.at<cv::Vec3d>(j,i)[2];
		}
	}
	pDst->Convert(TMO_RGB);

	pSrc->ProgressBar(j, pSrc->GetHeight());
	return 0;
}

// helper function linspace
// author: Damith Suranga Jinasena
// https://dsj23.me/2013/02/13/matlab-linspace-function-written-in-c/
std::vector<double> TMOAubry14::linspace(double min, double max, int n)
{
    std::vector<double> result;
    for (int i = 0; i <= n-2; i++) {
        double temp = min + i*(max-min)/(std::floor((double)n) - 1);
        result.push_back(temp);
    }
    result.push_back(max);
    return result;
}
