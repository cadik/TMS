/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*             Fast Local Laplacian Filters: Theory and Applications (2014)     *
* by Mathieu Aubry, Sylvain Paris, Samuel W. Hasinoff, Jan Kautz, Fredo Durand *
*                         ACM Transactions on Graphics                         *
*                                                                              *
*             Author: Tomas Hudziec [xhudzi01 AT stud.fit.vutbr.cz]            *
*         Term project for Computational Photography course - 2018             *
*       Part of master thesis (HDR support, code reorganization) - 2019        *
*                                                                              *
*******************************************************************************/

#include "FastLocalLaplFilt.h"

// The algorithm of Fast Local Laplacian Filters
cv::Mat FastLocalLaplFilt(cv::Mat I, double sigma, double fact, int N, TMOImage *pSrc)
{
	int height = I.rows;
	int width = I.cols;
	
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

	std::vector<double> discretisation = linspace(0, 1, N);
	double discretisationStep = discretisation[1];

	cv::Mat I_remap(I.size(), CV_64FC1);

	// main loop of the algorithm
	for (auto ref : discretisation) {
		// calculate I_remap
		for (int j = 0; j < I_remap.rows; j++) {
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
			for (int j = 0; j < outLaplacePyr[level].rows; j++) {
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
	
	return I_result_gray;
}

// helper function linspace
// author: Damith Suranga Jinasena
// https://dsj23.me/2013/02/13/matlab-linspace-function-written-in-c/
static std::vector<double> linspace(double min, double max, int n)
{
	std::vector<double> result;
	for (int i = 0; i <= n-2; i++)
		result.push_back(min + i*(max-min)/(n - 1));
	result.push_back(max);
	return result;
}
