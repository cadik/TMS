/* --------------------------------------------------------------------------- *
 * TMOFFLS08.cpp: implementation of the TMOFFLS08 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOFFLS08.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOFFLS08::TMOFFLS08()
{
	SetName(L"FFLS08");						// TODO - Insert operator name
	SetDescription(L"Edge-Preserving Decompositions for Multi-Scale Tone and Detail Manipulation");	// TODO - Insert description

	sigmaParameter.SetName(L"sigma");
	sigmaParameter.SetDescription(L"sigma color parameter");
	sigmaParameter.SetDefault(0.02);
	sigmaParameter.SetRange(0.01, 1.0);

	lambdaParameter.SetName(L"lambda");
	lambdaParameter.SetDescription(L"defining the amount of regularization");
	lambdaParameter.SetDefault(900);
	lambdaParameter.SetRange(100, 1000);

	multiplyParameter.SetName(L"scale");
	multiplyParameter.SetDescription(L"detail scale parameter");
	multiplyParameter.SetDefault(0.1);
	multiplyParameter.SetRange(0.1, 5.0);

	this->Register(sigmaParameter);
	this->Register(lambdaParameter);
	this->Register(multiplyParameter);
}

TMOFFLS08::~TMOFFLS08()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOFFLS08::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// convert input image to LAB color model
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array

	int height = pSrc->GetHeight();
	int width  = pSrc->GetWidth();

	cv::Mat I(height, width, CV_32FC3);		// INPUT IMAGE

	double r, g, b;

	// Convert to grayscale
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			I.at<cv::Vec3f>(j,i)[0] = r = *pSourceData++;
			I.at<cv::Vec3f>(j,i)[1] = g = *pSourceData++;
			I.at<cv::Vec3f>(j,i)[2] = b = *pSourceData++;
		}
	}

	// smoothing
  // fastGlobalSmootherFilter(InputArray guide,
	// 	InputArray src, OutputArray dst,
	// 	double lambda, double sigma_color,
	// 	double lambda_attenuation=0.25, int num_iter=3);

	cv::Mat B;		// BASE LAYER - edge-preserving smoothed input image
	cv::Mat guide;		// guide image for smoothing
	I.convertTo(guide, CV_8U);
	std::cout << "performing fastGlobalSmootherFilter with" << '\n';
	std::cout << "lambdaParameter: " << lambdaParameter << '\n';
	std::cout << "sigmaParameter: " << sigmaParameter << '\n';
	// lambda = 20^2 -- 30^2, sigma = 7e-2 -- 1e-1
	cv::ximgproc::fastGlobalSmootherFilter(guide, I, B, lambdaParameter, sigmaParameter*255.0);
	// cv::ximgproc::l0Smooth(I, B);
	std::cout << "fast global smoothing done" << '\n';
	guide.release();


	cv::Mat D;		// DETAIL LAYER
	D = I - B;

	double s = multiplyParameter;		// scale parameter

	cv::Mat E;		// ENHANCED IMAGE
	E = B + s*D;

	I.release(); B.release(); D.release();

	// output result
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			*pDestinationData++ = E.at<cv::Vec3f>(j,i)[0];
			*pDestinationData++ = E.at<cv::Vec3f>(j,i)[1];
			*pDestinationData++ = E.at<cv::Vec3f>(j,i)[2];
		}
	}

	E.release();

	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}
