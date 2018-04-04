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

	// double red, green, blue, gray;	// for manual conversion to gray

	// TODO comment properly variables
	cv::Mat I_RGB(height, width, CV_32FC3);
	cv::Mat I_Gray(height, width, CV_32FC1);


	// convert to grayscale
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// need to store rgb in mat to calculate ratio later
			I_RGB.at<cv::Vec3f>(j,i)[0] = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j,i)[1] = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j,i)[2] = *pSourceData++;
		}
	}

	// TODO if I need 64b precision, convert to gray manually
	// cvtColor handles max 32b floats
	cv::cvtColor(I_RGB, I_Gray, CV_RGB2GRAY);

	// TODO temporary, remove later
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// store results to the destination image
			*pDestinationData++ = I_Gray.at<float>(j,i);
			*pDestinationData++ = I_Gray.at<float>(j,i);
			*pDestinationData++ = I_Gray.at<float>(j,i);
		}
	}

	// this method works on grayscale image
	cv::Mat I = I_Gray;

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

	// calculate LLF
	// build Gaussian pyramid
	double pyrLevels = std::ceil(log(std::min(height, width))-log(2))+2;
	// 1.level is the image itself
	std::vector<cv::Mat> inGaussianPyr;
	inGaussianPyr.push_back(I_Gray);	// 1.level is the image itself
	cv::Mat GaussImg;
	for (size_t n = 1; n <= pyrLevels; n++) {
		cv::pyrDown(inGaussianPyr[n-1], GaussImg);
		inGaussianPyr.push_back(GaussImg);
	}

	// build Laplacian pyramid from Gaussian one
	// the last level is the same as last level of gaussian pyramid
	std::vector<cv::Mat> outLaplacePyr;
	outLaplacePyr.push_back(inGaussianPyr.back());
	cv::Mat smallerUpsampledGauss, LaplaceImg;
	for (size_t n = pyrLevels - 1; n > 0; n--) {
		cv::pyrUp(inGaussianPyr[n], smallerUpsampledGauss);
		cv::subtract(inGaussianPyr[n-1], smallerUpsampledGauss, LaplaceImg);
		outLaplacePyr.insert(outLaplacePyr.begin(), LaplaceImg);
	}

	// ...

	// multiply result with ratio
	// ...

	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
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
