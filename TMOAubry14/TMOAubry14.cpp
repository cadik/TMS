/* --------------------------------------------------------------------------- *
 * TMOAubry14.cpp: implementation of the TMOAubry14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAubry14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAubry14::TMOAubry14()
{
	SetName(L"Aubry14");						// TODO - Insert operator name
	SetDescription(L"Tone mapping using fast local Laplacian filters");	// TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(0.0,1.0);				// TODO - Add acceptable range if needed
	// this->Register(dParameter);
}

TMOAubry14::~TMOAubry14()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAubry14::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	// pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	// pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	int height = pSrc->GetHeight();
	int width  = pSrc->GetWidth();

	// TODO comment properly variables
	cv::Mat I_RGB(height, width, CV_64FC3);
	cv::Mat I_Gray(height, width, CV_64FC1);

	double R, G, B;

	// Convert to grayscale
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// need to store rgb in mat to calculate ratio later
			I_RGB.at<cv::Vec3f>(j,i)[0] = R = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j,i)[1] = G = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j,i)[2] = B = *pSourceData++;
			// convert to grayscale
			I_Gray.at<double>(j,i) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
		}
	}

	// if I need 64b precision, I must convert to gray manually,
	// because cvtColor handles max 32b floats
	// cv::cvtColor(I_RGB, I_Gray, CV_RGB2GRAY);

	// this method works on grayscale image
	cv::Mat I = I_Gray;
	// normalize from 0-255 to 0-1
	// cv::normalize(I, I, 0.0, 1.0, cv::NORM_MINMAX, CV_64F);
	// std::cout << "gray input img normalized: " << std::endl << I << std::endl << std::endl;

	// calculate ratio for converting to rgb at the end
	// cv::Mat ratioMat, grayMat3;
	// cv::Mat grayChannels[] = {I_Gray, I_Gray, I_Gray};
	// cv::merge(grayChannels, 3, grayMat3);
	// cv::Mat dividendMat(height, width, CV_8UC3, cv::Scalar::all(255));
	// cv::divide(I_RGB, grayMat3, ratioMat, 1/255.0, -1);
	// std::cout << "ratioMat = " << std::endl << " " << ratioMat << std::endl << std::endl;
	// ratioMat is not the same with that from matlab
	// TODO dump out all mats and compare with matlab ones
	// ...

	// The algorithm of Local Laplacian Filters follows

	// Build Gaussian pyramid
	// FIXME should be pyrLevels int?
	double pyrLevels = std::ceil(log(std::min(height, width))-log(2))+2;
	// 1.level is the image itself
	std::vector<cv::Mat> inGaussianPyr;
	// inGaussianPyr.push_back(I_Gray);	// 1.level is the image itself
	inGaussianPyr.push_back(I);	// 1.level is the image itself
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
	// std::cout << '\n' << "init outLaplacePyr: " << '\n';
	// for (auto l : outLaplacePyr) {
	// 	std::cout << '\n' << l << "\n\n";
	// }

	// TODO make these parameters of the method
	double sigma = 0.1;
	double fact = 5;
	int N = 10;

	std::vector<double> discretisation = this->linspace(0, 1, N);
	double discretisationStep = discretisation[1];

	cv::Mat I_remap(I.size(), CV_64FC1);

	// std::cout << '\n' << "I " << '\n' << I << "\n\n";

	// main loop of the algorithm
	for (auto ref : discretisation) {
	// for (int k = 0; k <= 1; k++) {
		// double ref = discretisation[k];
		// calculate I_remap
		for (j = 0; j < I_remap.rows; j++) {
			pSrc->ProgressBar(j, I_remap.rows);	// provide progress bar
			for (int i = 0; i < I_remap.cols; i++) {
				double pixI = I.at<double>(j,i);
				I_remap.at<double>(j,i) =
				fact*(pixI-ref)*exp(-(pixI-ref)*(pixI-ref)/(2.0*sigma*sigma));
			}
		}
		// std::cout << "ref " << ref << '\n' << I_remap << "\n\n";

		// std::cout << "\noutLaplacePyr before tmpLaplacePyr:\n" << outLaplacePyr[0] << "\n\n";
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
		// std::cout << '\n' << "ref " << ref << '\n';
		// for (auto l : tmpLaplacePyr) {
		// 	std::cout << '\n' << l << "\n\n";
		// }

		// std::cout << "\ntmpLaplacePyr before colp():\n" << tmpLaplacePyr[0] << "\n\n";
		// std::cout << "\noutLaplacePyr before colp():\n" << outLaplacePyr[0] << "\n\n";
		// computation of output Laplacian pyramid
		colp(pyrLevels,
					discretisationStep,
					ref,
					inGaussianPyr,
					tmpLaplacePyr,
					outLaplacePyr
		);
		// for (size_t level = 0; level < pyrLevels - 1; level++) {
		// 	for (j = 0; j < outLaplacePyr[level].rows; j++) {
		// 		pSrc->ProgressBar(j, outLaplacePyr[level].rows);	// provide progress bar
		// 		for (int i = 0; i < outLaplacePyr[level].cols; i++) {
		// 			double pixInGaussPyr = inGaussianPyr[level].at<double>(j,i);
		// 			double absDiff = abs(pixInGaussPyr - ref);
		// 			if (absDiff < discretisationStep) {
		// 				outLaplacePyr[level].at<double>(j,i) +=
		// 				tmpLaplacePyr[level].at<double>(j,i)*
		// 				(1-absDiff/discretisationStep);
		// 			}
		// 		}
		// 	}
		// }
		// outLaplacePyr seems to be better now
		// std::cout << '\n' << "ref " << ref << '\n';
		// for (auto l : outLaplacePyr) {
		// 	std::cout << '\n' << l << "\n\n";
		// }

	}// main loop of the algorithm

	// std::cout << '\n' << "final outLaplacePyr:" << '\n';
	// for (auto l : outLaplacePyr) {
	// 	std::cout << '\n' << l << "\n\n";
	// }
	// Reconstruct laplacian pyramid
	// start with low pass residual
	cv::Mat result = outLaplacePyr.back();
	for (int lev = pyrLevels - 2; lev >= 0; --lev) {
		// upsample, and add to current level
		cv::pyrUp(result, result);
		result += outLaplacePyr[lev];
	}
	// FIXME result visually does not look like it should
	// std::cout << "\nresult: " << result << "\n\n";

	// show grayscale result
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// store results to the destination image
			*pDestinationData++ = result.at<double>(j,i);
			*pDestinationData++ = result.at<double>(j,i);
			*pDestinationData++ = result.at<double>(j,i);
		}
	}

	// ...

	// TODO multiply result with ratio to get colours back
	// ...

	pSrc->ProgressBar(j, pSrc->GetHeight());
	// pDst->Convert(TMO_RGB);
	return 0;
}

// author: Damith Suranga Jinasena
// https://dsj23.me/2013/02/13/matlab-linspace-function-written-in-c/
std::vector<double> TMOAubry14::linspace(double min, double max, int n)
{
    std::vector<double> result;
    // vector iterator
	// TODO remove this and use i instead
    int iterator = 0;

    for (int i = 0; i <= n-2; i++)
    {
        double temp = min + i*(max-min)/(std::floor((double)n) - 1);
        result.insert(result.begin() + iterator, temp);
        iterator += 1;
    }

    result.insert(result.begin() + iterator, max);
    return result;
}

// compute output Laplacian pyramid for debugging
void TMOAubry14::colp(double pyrLevels,
	double discretisationStep,
	double ref,
	const std::vector<cv::Mat> &inGaussianPyr,
	const std::vector<cv::Mat> &tmpLaplacePyr,
	std::vector<cv::Mat> &outLaplacePyr
)
{
	// std::cout << '\n' << "ref: " << ref << "\noutLaplacePyr: " << '\n';
	// std::cout << outLaplacePyr[0] << "\n\n";

	for (size_t level = 0; level < pyrLevels - 1; level++) {
		for (int j = 0; j < outLaplacePyr[level].rows; j++) {
			// pSrc->ProgressBar(j, outLaplacePyr[level].rows);	// provide progress bar
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
}
