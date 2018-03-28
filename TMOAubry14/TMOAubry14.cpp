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

	double red, green, blue, gray;

	cv::Mat inRGB, inGray;
	inRGB = cv::Mat::zeros (height, width, CV_32FC3);
	inGray = cv::Mat::zeros (height, width, CV_32FC1);


	// convert to grayscale
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// need to store rgb in mat to calculate ratio later
			inRGB.at<cv::Vec3f>(j,i)[0] = *pSourceData++;
			inRGB.at<cv::Vec3f>(j,i)[1] = *pSourceData++;
			inRGB.at<cv::Vec3f>(j,i)[2] = *pSourceData++;
		}
	}

	// cvtColor handles max 32b floats
	cv::cvtColor(inRGB, inGray, CV_RGB2GRAY);

	// temporary, remove later
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			// store results to the destination image
			*pDestinationData++ = inGray.at<float>(j,i);
			*pDestinationData++ = inGray.at<float>(j,i);
			*pDestinationData++ = inGray.at<float>(j,i);
		}
	}

	// calculate ratio for converting to rgb at the end
	// cv::Mat ratioMat, grayMat3;
	// cv::Mat grayChannels[] = {inGray, inGray, inGray};
	// cv::merge(grayChannels, 3, grayMat3);
	// cv::Mat dividendMat(height, width, CV_8UC3, cv::Scalar::all(255));
	// cv::divide(inRGB, grayMat3, ratioMat, 1/255.0, -1);
	// std::cout << "ratioMat = " << std::endl << " " << ratioMat << std::endl << std::endl;
	// ratioMat is not the same with that from matlab
	// TODO dump out all mats and compare with matlab ones
	// ...

	// calculate LLF
	// build Gaussian pyramid
	double numOfInGaussLevels = std::ceil(log(std::min(height, width))-log(2))+2;
	// 1.level is the image itself
	std::vector<cv::Mat> inGaussianPyr;
	inGaussianPyr.push_back(inGray);	// 1.level is the image itself
	cv::Mat GaussImg;
	for (size_t n = 1; n <= numOfInGaussLevels; n++) {
		cv::pyrDown(inGaussianPyr[n-1], GaussImg);
		inGaussianPyr.push_back(GaussImg);
	}

	// build Laplacian pyramid from Gaussian one
	// the last level is the same as last level of gaussian pyramid
	std::vector<cv::Mat> outLaplacePyr;
	outLaplacePyr.push_back(inGaussianPyr.back());
	cv::Mat smallerUpsampledGauss, LaplaceImg;
	for (size_t n = numOfInGaussLevels - 1; n > 0; n--) {
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
