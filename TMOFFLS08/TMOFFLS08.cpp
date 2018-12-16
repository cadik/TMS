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
  // cv::cvtColor(show_img, show_img, cv::COLOR_RGB2BGR);
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

	bFineP.SetName(L"fine details");
	bFineP.SetDescription(L"include fine detail enhancement into result");
	bFineP.SetDefault(true);
	bFineP=true;

	bMediumP.SetName(L"medium details");
	bMediumP.SetDescription(L"include medium detail enhancement into result");
	bMediumP.SetDefault(true);
	bMediumP=true;

	bCoarseP.SetName(L"coarse details");
	bCoarseP.SetDescription(L"include coarse detail enhancement into result");
	bCoarseP.SetDefault(true);
	bCoarseP=true;

	// dExposureP.SetName(L"exposure");
	// dExposureP.SetDescription(L"exposure parameter");
	// dExposureP.SetDefault(1.0);
	// dExposureP.SetRange(0.1, 2);

	// this->Register(dExposureP);
	this->Register(bCoarseP);
	this->Register(bMediumP);
	this->Register(bFineP);
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

	// normalize to <0,255> range for smoothing
	cv::normalize(I_RGB, I_RGB, 0, 255, cv::NORM_MINMAX, I_RGB.type());
	// convert to grayscale
  // cv::cvtColor(I_RGB, I_Gray, cv::COLOR_RGB2GRAY);

	cv::Mat RGB_Smooth0, RGB_Smooth1;		// smoothed versions of grayscale input image
	cv::Mat guide;		// guide image for smoothing
	I_RGB.convertTo(guide, CV_8UC3);

	// SMOOTHING
  // fastGlobalSmootherFilter(InputArray guide,
	// 	InputArray src, OutputArray dst,
	// 	double lambda, double sigma_color,
	// 	double lambda_attenuation=0.25, int num_iter=3);
	// in paper: lambda = 20^2 -- 30^2, sigma = 7e-2 -- 1e-1
	// create 2 versions of smoothed images: L0 and L1
	cv::ximgproc::fastGlobalSmootherFilter(guide, I_RGB, RGB_Smooth0, 20.0, 0.02*255.0);
	cv::ximgproc::fastGlobalSmootherFilter(guide, I_RGB, RGB_Smooth1, 40.0, 0.03*255.0);
	// cv::ximgproc::l0Smooth(I, B);	// does not work for me
	guide.release();

	// convert I_RGB to LAB
	cv::Mat I_Lab;
	cv::normalize(I_RGB, I_RGB, 0, 1, cv::NORM_MINMAX, I_RGB.type());
  cv::cvtColor(I_RGB, I_Lab, cv::COLOR_RGB2Lab);

	// convert smoothed versions to LAB and get Luminance channel
	// TODO try later with only grayscale conversion, if it will work also
	cv::Mat LAB_Smooth0, LAB_Smooth1;
	cv::Mat L0, L1;

	cv::normalize(RGB_Smooth0, RGB_Smooth0, 0, 1, cv::NORM_MINMAX, RGB_Smooth0.type());
  cv::cvtColor(RGB_Smooth0, LAB_Smooth0, cv::COLOR_RGB2Lab);
	cv::extractChannel(LAB_Smooth0, L0, 0);

	cv::normalize(RGB_Smooth1, RGB_Smooth1, 0, 1, cv::NORM_MINMAX, RGB_Smooth1.type());
	cv::cvtColor(RGB_Smooth1, LAB_Smooth1, cv::COLOR_RGB2Lab);
	cv::extractChannel(LAB_Smooth1, L1, 0);

	RGB_Smooth0.release(); RGB_Smooth1.release();
	LAB_Smooth0.release(); LAB_Smooth1.release();

	// shift minimum of L's to 0
	// FIXME is it really necessary?
	double L0_min, L0_max;
	cv::minMaxLoc(L0, &L0_min, &L0_max);
	L0 = L0 - L0_min;
	double L1_min, L1_max;
	cv::minMaxLoc(L1, &L1_min, &L1_max);
	L1 = L1 - L1_min;


	// apply tonemapLAB algorithm from EPD method (FFLS08)
	// on different detail scales
	double val0, val1, val2;
	double exposure, gamma, saturation;
	cv::Mat fine = cv::Mat::zeros(height, width, CV_32FC3);
	cv::Mat medium = cv::Mat::zeros(height, width, CV_32FC3);
	cv::Mat coarse = cv::Mat::zeros(height, width, CV_32FC3);
	double count = 0.0;
	if(!(bFineP || bMediumP || bCoarseP)) {
		std::cerr << "none of detail check boxes were set, setting them all" << '\n';
		bFineP = true;
		bMediumP = true;
		bCoarseP = true;
		// TODO set to true also check boxes if possible
		count = 3.0;
	}

	// fine details
	if(bFineP) {
		val0 = 25;
		val1 = 1;
		val2 = 1;
		exposure = 1.0;
		saturation = 1.1;
		gamma = 1.0;
		fine = tonemapLAB(I_Lab, L0, L1,
											val0, val1, val2,
											exposure, gamma, saturation);
		count++;
	}

	// medium details
	if(bMediumP) {
		val0 = 1;
		val1 = 40;
		val2 = 1;
		exposure = 1.0;
		saturation = 1.1;
		gamma = 1.0;
		medium = tonemapLAB(I_Lab, L0, L1,
												val0, val1, val2,
												exposure, gamma, saturation);
		count++;
	}

	// coarse details
	if(bCoarseP) {
		val0 = 4;
		val1 = 1;
		val2 = 15;
		exposure = 1.10;
		saturation = 1.1;
		gamma = 1.0;
		coarse = tonemapLAB(I_Lab, L0, L1,
												val0, val1, val2,
												exposure, gamma, saturation);
		count++;
	}

	L0.release(); L1.release(); I_Lab.release();

	if(count < 1.0) count = 1.0;
	cv::Mat combined = (fine + medium + coarse) / count;
	fine.release(); medium.release(); coarse.release();

	// normalize to 0-255 range for display
	cv::normalize(combined, combined, 0, 255, cv::NORM_MINMAX, CV_32FC3);

	// output result
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// provide progress bar
		for (int i = 0; i < width; i++)
		{
			*pDestinationData++ = combined.at<cv::Vec3f>(j,i)[0];
			*pDestinationData++ = combined.at<cv::Vec3f>(j,i)[1];
			*pDestinationData++ = combined.at<cv::Vec3f>(j,i)[2];
		}
	}

	pSrc->ProgressBar(j, pSrc->GetHeight());
	// pDst->Convert(TMO_RGB);

	return 0;
}

// Applies a sigmoid function on the data X in [0-1] range.
// Then rescales the result so 0.5 will be mapped to itself.
cv::Mat TMOFFLS08::sigmoid(cv::Mat X, double a)
{
	double x, y, y05;
	cv::Mat Y(X.rows, X.cols, X.type());
	y05 = 1.0 / (1 + exp(-a*0.5)) - 0.5;
	for (int j = 0; j < X.rows; j++) {
		pSrc->ProgressBar(j, X.rows);	// provide progress bar
		for (int i = 0; i < X.cols; i++) {
			// TODO maybe recognize data type dynamically, not just hard-coded 'float'
			// or check type of input matrix that it must be float
			x = X.at<float>(j,i);
			// apply sigmoid
			y = 1.0 / (1 + exp(-a*x)) - 0.5;
			// re-scale
			Y.at<float>(j,i) = y * (0.5/y05);
		}
	}
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
// returns RGB image
cv::Mat TMOFFLS08::tonemapLAB(cv::Mat Lab, cv::Mat L0, cv::Mat L1,
															double val0, double val1, double val2,
															double exposure, double gamma, double saturation)
{
	cv::Mat L, a, b;

	cv::Mat Lab_channels[3];
	cv::split(Lab, Lab_channels);
	L = Lab_channels[0];
	a = Lab_channels[1];
	b = Lab_channels[2];

	// shift minimum of L to zero
	double L_min, L_max;
	cv::minMaxLoc(L, &L_min, &L_max);
	L = L - L_min;

	// L's are in range 0-100
	cv::Mat diff0 = this->sigmoid((L-L0)/100, val0)*100;

	cv::Mat diff1 = this->sigmoid((L0-L1)/100, val1)*100;

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
