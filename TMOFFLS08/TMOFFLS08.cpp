/* --------------------------------------------------------------------------- *
 * TMOFFLS08.cpp: implementation of the TMOFFLS08 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOFFLS08.h"

// util functions
// TODO remove for final version
void showImg(std::string name, const cv::Mat &img)
{
	cv::Mat show_img = img.clone();
	cv::normalize(show_img, show_img, 0, 1, cv::NORM_MINMAX, show_img.type());
  cv::cvtColor(show_img, show_img, cv::COLOR_RGB2BGR);
  // cv::cvtColor(show_img, show_img, cv::COLOR_Lab2BGR);
	cv::imshow(name, show_img);
	cv::waitKey(0);
	cv::destroyAllWindows();
}

void printMatRange(std::string name, const cv::Mat &mat)
{
	double min, max;
	cv::minMaxLoc(mat, &min, &max);
	std::cout << "name: " << name << '\n';
	std::cout << "min: " << min << '\n';
	std::cout << "max: " << max << "\n\n";
}

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
	lambdaParameter.SetRange(1, 1000);

	this->Register(sigmaParameter);
	this->Register(lambdaParameter);
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

	// convert input/output image data to RGB color model for smoothing
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array

	int height = pSrc->GetHeight();
	int width  = pSrc->GetWidth();

	cv::Mat I_RGB(height, width, CV_32FC3);		// RGB INPUT IMAGE
	cv::Mat I_Gray(height, width, CV_32FC1);		// GRAYSCALE INPUT IMAGE

	float r, g, b;

	// put source data into opencv matrices
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			// L.at<float>(j,i) = *pSourceData;
			I_RGB.at<cv::Vec3f>(j,i)[0] = r = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j,i)[1] = g = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j,i)[2] = b = *pSourceData++;
			// convert to grayscale
			// I_Gray.at<float>(j,i) = (0.2989*r + 0.5870*g + 0.1140*b);
		}
	}

	printMatRange("input rgb", I_RGB);
	// std::cout << "I_RGB" << I_RGB << '\n';
	// normalize to <0,255> range for smoothing
	cv::normalize(I_RGB, I_RGB, 0, 255, cv::NORM_MINMAX, I_RGB.type());
	printMatRange("input rgb 0-255", I_RGB);
	// convert to grayscale
  cv::cvtColor(I_RGB, I_Gray, cv::COLOR_RGB2GRAY);
	printMatRange("grayscale", I_Gray);
	// showImg("grayscale", I_Gray);

	// cv::Mat L;		// Lightness/Luminance from Lab
	// cv::extractChannel(I_Lab, L, 0);
	// shift minimum to zero
	// double L_min, L_max;
	// cv::minMaxLoc(L, &L_min, &L_max);
	// L = L - L_min;
	// TODO maybe scale L to 0-100?
	// cv::normalize(L, L, 0, 100, cv::NORM_MINMAX, L.type());

	// WARNING: L != I_Gray
	// printMatRange("Luminance from Lab", L);
	// showImg("Luminance from Lab", L);

	cv::Mat RGB_Smooth0, RGB_Smooth1;		// smoothed versions of grayscale input image
	cv::Mat guide;		// guide image for smoothing
	I_RGB.convertTo(guide, CV_8UC3);
	// L.convertTo(guide, CV_8U);
	// Lab.convertTo(guide, CV_8UC3);
	// showImg("I_RGB", I_RGB);

	// SMOOTHING
  // fastGlobalSmootherFilter(InputArray guide,
	// 	InputArray src, OutputArray dst,
	// 	double lambda, double sigma_color,
	// 	double lambda_attenuation=0.25, int num_iter=3);
	// in paper: lambda = 20^2 -- 30^2, sigma = 7e-2 -- 1e-1
	// create 2 versions of smoothed images: L0 and L1
	// TODO smooth image in grayscale
	// TODO try also with I_Gray instead of L (they are not the same!)
	// WARNING: This smoothing works on rgb or gray. I don't know what it does with L from Lab.
	cv::ximgproc::fastGlobalSmootherFilter(guide, I_RGB, RGB_Smooth0, 20.0, 0.02*255.0);
	cv::ximgproc::fastGlobalSmootherFilter(guide, I_RGB, RGB_Smooth1, 40.0, 0.03*255.0);
	// cv::ximgproc::fastGlobalSmootherFilter(guide, L, L0, lambdaParameter, sigmaParameter*255.0);
	// cv::ximgproc::fastGlobalSmootherFilter(guide, L, L1, lambdaParameter, sigmaParameter*255.0);
	// cv::ximgproc::l0Smooth(I, B);	// does not work for me
	guide.release();
	printMatRange("RGB_Smooth0", RGB_Smooth0);
	// showImg("RGB_Smooth1", RGB_Smooth1);
	printMatRange("RGB_Smooth1", RGB_Smooth1);

	// convert I_RGB to LAB
	cv::Mat I_Lab;
	cv::normalize(I_RGB, I_RGB, 0, 1, cv::NORM_MINMAX, I_RGB.type());
  cv::cvtColor(I_RGB, I_Lab, cv::COLOR_RGB2Lab);

	// convert smoothed versions to LAB and get Luminance channel
	// TODO try later with only grayscale conversion, if it will work also
	cv::Mat LAB_Smooth0, LAB_Smooth1;
	cv::Mat L0, L1 = cv::Mat::ones(height, width, CV_32FC1);
	cv::normalize(RGB_Smooth0, RGB_Smooth0, 0, 1, cv::NORM_MINMAX, RGB_Smooth0.type());
  cv::cvtColor(RGB_Smooth0, LAB_Smooth0, cv::COLOR_RGB2Lab);
	// std::cout << "RGB_Smooth0" << RGB_Smooth0 << '\n';
	// std::cout << "LAB_Smooth0" << LAB_Smooth0 << '\n';
	cv::extractChannel(LAB_Smooth0, L0, 0);
	// std::cout << "L0" << L0 << '\n';

	cv::normalize(RGB_Smooth1, RGB_Smooth1, 0, 1, cv::NORM_MINMAX, RGB_Smooth1.type());
	cv::cvtColor(RGB_Smooth1, LAB_Smooth1, cv::COLOR_RGB2Lab);
	// std::cout << "RGB_Smooth1" << RGB_Smooth1 << '\n';
	// std::cout << "LAB_Smooth1" << LAB_Smooth1 << '\n';
	// printMatRange("RGB_Smooth1", RGB_Smooth1);
	// printMatRange("LAB_Smooth1", LAB_Smooth1);
	cv::extractChannel(LAB_Smooth1, L1, 0);
	// std::cout << "L1" << L1 << '\n';
	// showImg("L0", L0);
	// showImg("L1", L1);
	// RGB_Smooth0.release();
	// RGB_Smooth1.release();
	// LAB_Smooth0.release();
	// LAB_Smooth1.release();
	// return 0;

	// TODO shift minimum of L's to 0
	double L0_min, L0_max;
	cv::minMaxLoc(L0, &L0_min, &L0_max);
	L0 = L0 - L0_min;
	double L1_min, L1_max;
	cv::minMaxLoc(L1, &L1_min, &L1_max);
	L1 = L1 - L1_min;


	// normalize to <0,1> range for tone mapping
	// cv::normalize(I_Lab, I_Lab, 0, 100, cv::NORM_MINMAX, I_Lab.type());
	// cv::normalize(L0, L0, 0, 100, cv::NORM_MINMAX, L0.type());
	// cv::normalize(L1, L1, 0, 100, cv::NORM_MINMAX, L1.type());

	// apply tonemapLAB algorithm from EPD method (FFLS08)
	double val0, val1, val2;
	double exposure, saturation, gamma;
	// fine details
	val0 = 25;
	val1 = 1;
	val2 = 1;
	exposure = 1.0;
	saturation = 1.1;
	gamma = 1.0;
	cv::Mat fine;
	cv::Mat Lab_channels[3];
	cv::split(I_Lab, Lab_channels);
	fine = Lab_channels[0];

	fine = tonemapLAB(I_Lab, L0, L1,
										val0, val1, val2,
										exposure, gamma, saturation);

	// normalize to 0-255 range for display
	// cv::normalize(fine, fine, 0, 255, cv::NORM_MINMAX, CV_32FC3);
	cv::normalize(fine, fine, 0, 255, cv::NORM_MINMAX, CV_32FC1);
	// cv::normalize(L0, L0, 0, 255, cv::NORM_MINMAX, CV_32FC1);
	// cv::normalize(L, L, 0, 255, cv::NORM_MINMAX, CV_32FC1);

	// output result
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			*pDestinationData++ = fine.at<cv::Vec3f>(j,i)[0];
			*pDestinationData++ = fine.at<cv::Vec3f>(j,i)[1];
			*pDestinationData++ = fine.at<cv::Vec3f>(j,i)[2];
			// *pDestinationData++ = fine.at<float>(j,i);
			// *pDestinationData++ = fine.at<float>(j,i);
			// *pDestinationData++ = fine.at<float>(j,i);
		}
	}

	// L.release();
	L0.release();
	L1.release();
	I_Lab.release();
	fine.release();

	pSrc->ProgressBar(j, pSrc->GetHeight());
	// pDst->Convert(TMO_RGB);

	return 0;
}

// TODO test it
// Applies a sigmoid function on the data X in [0-1] range.
// Then rescales the result so 0.5 will be mapped to itself.
cv::Mat TMOFFLS08::sigmoid(cv::Mat X, double a)
{
	double x, y, y05;
	cv::Mat Y(X.rows, X.cols, X.type());
	for (int j = 0; j < X.rows; j++) {
		pSrc->ProgressBar(j, X.rows);	// provide progress bar
		for (int i = 0; i < X.cols; i++) {
			// TODO maybe recognize data type dynamically, not just hard-coded 'float'
			x = X.at<float>(j,i);
			// apply sigmoid
			y = 1.0 / (1 + exp(-a*x)) - 0.5;
			// re-scale
			y05 = 1.0 / (1 + exp(-a*0.5)) - 0.5;
			Y.at<float>(j,i) = y * (0.5/y05);
		}
	}
	// std::cout << "X.type(): " << X.type() << '\n';
	// std::cout << "Y.type(): " << Y.type() << '\n';
	return Y;
}


// DESCR:
// This function gets an image in the CIELAB color
// space and tone maps it according to the parameters.
//
// PARAMS:
// lab is the image in CIELAB color space
// L0, L1 are smoothed versions of L of LAB
//
// val0-val3 compression/expansion params in [-1, 1] range
// exposure is in [0,inf) range
// gamma is in (0,1] range
// saturation is in [0,inf) range
//
// returns rgb image
cv::Mat TMOFFLS08::tonemapLAB(cv::Mat Lab, cv::Mat L0, cv::Mat L1,
															double val0, double val1, double val2,
															double exposure, double gamma, double saturation)
{
	cv::Mat L, a, b;

	cv::Mat Lab_channels[3];
	cv::split(Lab, Lab_channels);
	L = Lab_channels[0];
	// std::cout << "L" << L << '\n';
	a = Lab_channels[1];
	b = Lab_channels[2];
	// cv::extractChannel(Lab, L, 0);

	// shift minimum of L to zero
	double L_min, L_max;
	cv::minMaxLoc(L, &L_min, &L_max);
	L = L - L_min;

	printMatRange("L", L);
	printMatRange("L0", L0);
	printMatRange("L1", L1);

	// L's are in range 0-100
	cv::Mat diff0 = this->sigmoid((L-L0)/100, val0)*100;

	cv::Mat diff1 = this->sigmoid((L0-L1)/100, val0)*100;

	cv::Mat base = (this->sigmoid((exposure*L1-56)/100, val2)*100)+56;

	cv::Mat L_res = base + diff1 + diff0;

	// TODO maybe shift minimum of L_res back (from zero up)
	// ...
	Lab_channels[0] = L_res;
	// TODO multiply a and b channels with saturation
	// ...
	cv::Mat Lab_res;
	cv::merge(Lab_channels, 3, Lab_res);

	// convert to rgb
	cv::Mat RGB_res;
  cv::cvtColor(Lab_res, RGB_res, cv::COLOR_Lab2RGB);

	return RGB_res;
}

	// std::cout << "Lab.type(): " << Lab.type() << '\n';
	// std::cout << "L.type(): " << L.type() << '\n';
	// std::cout << "L0.type(): " << L0.type() << '\n';
	// std::cout << "L1.type(): " << L1.type() << '\n';
	// cv::imshow("L", L);
	// cv::waitKey(0);
