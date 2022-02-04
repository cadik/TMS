/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                               *
*                                                                              *
*                       Diploma thesis                                         *
*                       Author: Tomas Hudziec [xhudzi01 AT stud.fit.vutbr.cz]  *
*                       Brno 2019                                              *
*                                                                              *
*                       Edge-Preserving Decomposition for Multi-Scale          *
*                       Tone and Detail Manipulation (2018)                    *
*                       by Z. Farbman, R. Fattal, D.Lischinski and R. Szeliski *
*                       http://www.cs.huji.ac.il/~danix/epd/                   *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOFarbman08.cpp
 * @brief Edge-Preserving Decomposition for Multi-Scale Tone and Detail Manipulation (2018) by Z. Farbman, R. Fattal, D.Lischinski and R. Szeliski http://www.cs.huji.ac.il/~danix/epd/ 
 * @author Tomas Hudziec
 * @class TMOFarbman08
 * 
 * WARNING: uses Fast Global Image Smoothing from OpenCV based on Weighted Least
 * Squares instead of original Weighted Least Squares smoothing itself
 */

#include "TMOFarbman08.h"

/**
  *  @brief Constructor
  */
TMOFarbman08::TMOFarbman08()
{
	SetName(L"Farbman08");
	SetDescription(L"Edge-Preserving Decompositions for Multi-Scale Tone and Detail Manipulation");

	bFineP.SetName(L"fine details");
	bFineP.SetDescription(L"include fine detail enhancement into result");
	bFineP.SetDefault(true);
	bFineP = true;
	dFineVal0P.SetName(L"val0");
	dFineVal0P.SetDescription(L"fine detail compression/expansion parameter for fine detail layer");
	dFineVal0P.SetRange(-100, 100);
	dFineVal0P.SetDefault(25);
	dFineVal0P = 25;
	dFineVal1P.SetName(L"val1");
	dFineVal1P.SetDescription(L"medium detail compression/expansion parameter for fine detail layer");
	dFineVal1P.SetRange(-100, 100);
	dFineVal1P.SetDefault(1);
	dFineVal1P = 1;
	dFineVal2P.SetName(L"val2");
	dFineVal2P.SetDescription(L"coarse detail compression/expansion parameter for fine detail layer");
	dFineVal2P.SetRange(-100, 100);
	dFineVal2P.SetDefault(1);
	dFineVal2P = 1;
	dFineExposureP.SetName(L"exposure");
	dFineExposureP.SetDescription(L"exposure for fine detail layer");
	dFineExposureP.SetDefault(1.0);
	dFineExposureP.SetRange(0.0, 100.0); // [0,inf)
	dFineExposureP = 1.0;
	dFineGammaP.SetName(L"gamma");
	dFineGammaP.SetDescription(L"gamma for fine detail layer");
	dFineGammaP.SetDefault(1.0);
	dFineGammaP.SetRange(0.1, 1.0); // (0,1]
	dFineGammaP = 1.0;
	dFineSaturationP.SetName(L"saturation");
	dFineSaturationP.SetDescription(L"saturation for fine detail layer");
	dFineSaturationP.SetDefault(1.1);
	dFineSaturationP.SetRange(0.0, 100.0); // [0,inf)
	dFineSaturationP = 1.1;

	bMediumP.SetName(L"medium details");
	bMediumP.SetDescription(L"include medium detail enhancement into result");
	bMediumP.SetDefault(true);
	bMediumP = true;
	dMediumVal0P.SetName(L"val0");
	dMediumVal0P.SetDescription(L"fine detail compression/expansion parameter for medium detail layer");
	dMediumVal0P.SetRange(-100, 100);
	dMediumVal0P.SetDefault(1);
	dMediumVal0P = 1;
	dMediumVal1P.SetName(L"val1");
	dMediumVal1P.SetDescription(L"medium detail compression/expansion parameter for medium detail layer");
	dMediumVal1P.SetRange(-100, 100);
	dMediumVal1P.SetDefault(40);
	dMediumVal1P = 40;
	dMediumVal2P.SetName(L"val2");
	dMediumVal2P.SetDescription(L"coarse detail compression/expansion parameter for medium detail layer");
	dMediumVal2P.SetRange(-100, 100);
	dMediumVal2P.SetDefault(1);
	dMediumVal2P = 1;
	dMediumExposureP.SetName(L"exposure");
	dMediumExposureP.SetDescription(L"exposure for medium detail layer");
	dMediumExposureP.SetDefault(1.0);
	dMediumExposureP.SetRange(0.0, 100.0); // [0,inf)
	dMediumExposureP = 1.0;
	dMediumGammaP.SetName(L"gamma");
	dMediumGammaP.SetDescription(L"gamma for medium detail layer");
	dMediumGammaP.SetDefault(1.0);
	dMediumGammaP.SetRange(0.1, 1.0); // (0,1]
	dMediumGammaP = 1.0;
	dMediumSaturationP.SetName(L"saturation");
	dMediumSaturationP.SetDescription(L"saturation for medium detail layer");
	dMediumSaturationP.SetDefault(1.1);
	dMediumSaturationP.SetRange(0.0, 100.0); // [0,inf)
	dMediumSaturationP = 1.1;

	bCoarseP.SetName(L"coarse details");
	bCoarseP.SetDescription(L"include coarse detail enhancement into result");
	bCoarseP.SetDefault(true);
	bCoarseP = true;
	dCoarseVal0P.SetName(L"val0");
	dCoarseVal0P.SetDescription(L"fine detail compression/expansion parameter for coarse detail layer");
	dCoarseVal0P.SetRange(-100, 100);
	dCoarseVal0P.SetDefault(4);
	dCoarseVal0P = 4;
	dCoarseVal1P.SetName(L"val1");
	dCoarseVal1P.SetDescription(L"medium detail compression/expansion parameter for coarse detail layer");
	dCoarseVal1P.SetRange(-100, 100);
	dCoarseVal1P.SetDefault(1);
	dCoarseVal1P = 1;
	dCoarseVal2P.SetName(L"val2");
	dCoarseVal2P.SetDescription(L"coarse detail compression/expansion parameter for coarse detail layer");
	dCoarseVal2P.SetRange(-100, 100);
	dCoarseVal2P.SetDefault(15);
	dCoarseVal2P = 15;
	dCoarseExposureP.SetName(L"exposure");
	dCoarseExposureP.SetDescription(L"exposure for coarse detail layer");
	dCoarseExposureP.SetDefault(1.1);
	dCoarseExposureP.SetRange(0.0, 100.0); // [0,inf)
	dCoarseExposureP = 1.1;
	dCoarseGammaP.SetName(L"gamma");
	dCoarseGammaP.SetDescription(L"gamma for coarse detail layer");
	dCoarseGammaP.SetDefault(1.0);
	dCoarseGammaP.SetRange(0.1, 1.0); // (0,1]
	dCoarseGammaP = 1.0;
	dCoarseSaturationP.SetName(L"saturation");
	dCoarseSaturationP.SetDescription(L"saturation for coarse detail layer");
	dCoarseSaturationP.SetDefault(1.1);
	dCoarseSaturationP.SetRange(0.0, 100.0); // [0,inf)
	dCoarseSaturationP = 1.1;

	this->Register(dCoarseSaturationP);
	this->Register(dCoarseGammaP);
	this->Register(dCoarseExposureP);
	this->Register(dCoarseVal2P);
	this->Register(dCoarseVal1P);
	this->Register(dCoarseVal0P);
	this->Register(bCoarseP);

	this->Register(dMediumSaturationP);
	this->Register(dMediumGammaP);
	this->Register(dMediumExposureP);
	this->Register(dMediumVal2P);
	this->Register(dMediumVal1P);
	this->Register(dMediumVal0P);
	this->Register(bMediumP);

	this->Register(dFineSaturationP);
	this->Register(dFineGammaP);
	this->Register(dFineExposureP);
	this->Register(dFineVal2P);
	this->Register(dFineVal1P);
	this->Register(dFineVal0P);
	this->Register(bFineP);
}

/**
  *  @brief Destructor
  */
TMOFarbman08::~TMOFarbman08()
{
}

/**
 * @brief Edge-Preserving Decomposition for Multi-Scale Tone and Detail Manipulation (2018)
 */
int TMOFarbman08::Transform()
{
	/** convert input/output image data to RGB color model for smoothing */
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double *pSourceData = pSrc->GetData();		/** You can work at low level data */
	double *pDestinationData = pDst->GetData(); /** Data are stored in form of array */

	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();

	cv::Mat I_RGB(height, width, CV_32FC3); /** RGB INPUT IMAGE */

	float r, g, b;

	/** put source data into opencv matrices */
	int j = 0;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height); /** provide progress bar */
		for (int i = 0; i < width; i++)
		{
			// L.at<float>(j,i) = *pSourceData;
			I_RGB.at<cv::Vec3f>(j, i)[0] = r = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j, i)[1] = g = *pSourceData++;
			I_RGB.at<cv::Vec3f>(j, i)[2] = b = *pSourceData++;
		}
	}

	/** normalize to <0,255> range for smoothing */
	cv::normalize(I_RGB, I_RGB, 0, 255, cv::NORM_MINMAX, I_RGB.type());

	cv::Mat RGB_Smooth0, RGB_Smooth1; /** smoothed versions of grayscale input image */
	cv::Mat guide;					  /** guide image for smoothing */
	I_RGB.convertTo(guide, CV_8UC3);

	/** SMOOTHING
	 *using fastGlobalSmootherFilter from paper
	 * Fast Global Image Smoothing Based on Weighted Least Squares
	 * (https://sites.google.com/site/globalsmoothing/)
	 * which is ~30x faster than original WLS smoothing used in EPD method

	 * fastGlobalSmootherFilter(InputArray guide,
	 * 	InputArray src, OutputArray dst,
	 * 	double lambda, double sigma_color,
	 * 	double lambda_attenuation=0.25, int num_iter=3);
	 * in paper: lambda = 20^2 -- 30^2, sigma = 7e-2 -- 1e-1
	 * create 2 versions of smoothed images: L0 and L1 */
	cv::ximgproc::fastGlobalSmootherFilter(guide, I_RGB, RGB_Smooth0, 20.0, 0.02 * 255.0);
	cv::ximgproc::fastGlobalSmootherFilter(guide, I_RGB, RGB_Smooth1, 40.0, 0.03 * 255.0);
	guide.release();

	/** convert I_RGB to LAB */
	cv::Mat I_Lab;
	cv::normalize(I_RGB, I_RGB, 0, 1, cv::NORM_MINMAX, I_RGB.type());
	cv::cvtColor(I_RGB, I_Lab, cv::COLOR_RGB2Lab);

	/** convert smoothed versions to LAB and get Luminance channel */
	cv::Mat LAB_Smooth0, LAB_Smooth1;
	cv::Mat L0, L1;

	cv::normalize(RGB_Smooth0, RGB_Smooth0, 0, 1, cv::NORM_MINMAX, RGB_Smooth0.type());
	cv::cvtColor(RGB_Smooth0, LAB_Smooth0, cv::COLOR_RGB2Lab);
	cv::extractChannel(LAB_Smooth0, L0, 0);

	cv::normalize(RGB_Smooth1, RGB_Smooth1, 0, 1, cv::NORM_MINMAX, RGB_Smooth1.type());
	cv::cvtColor(RGB_Smooth1, LAB_Smooth1, cv::COLOR_RGB2Lab);
	cv::extractChannel(LAB_Smooth1, L1, 0);

	RGB_Smooth0.release();
	RGB_Smooth1.release();
	LAB_Smooth0.release();
	LAB_Smooth1.release();

	/** apply tonemapLAB algorithm from EPD method (Farbman08)
	 * on different detail scales */
	double val0, val1, val2;
	double exposure, gamma, saturation;
	cv::Mat fine = cv::Mat::zeros(height, width, CV_32FC3);
	cv::Mat medium = cv::Mat::zeros(height, width, CV_32FC3);
	cv::Mat coarse = cv::Mat::zeros(height, width, CV_32FC3);
	double count = 0.0;
	if (!(bFineP || bMediumP || bCoarseP))
	{
		std::cout << "none of detail check boxes were set, setting them all" << '\n';
		bFineP = true;
		bMediumP = true;
		bCoarseP = true;
		count = 3.0;
	}

	/** fine details */
	if (bFineP)
	{
		val0 = dFineVal0P;
		val1 = dFineVal1P;
		val2 = dFineVal2P;
		exposure = dFineExposureP;
		saturation = dFineSaturationP;
		gamma = dFineGammaP;
		fine = tonemapLAB(I_Lab, L0, L1,
						  val0, val1, val2,
						  exposure, gamma, saturation);
		count++;
	}

	/** medium details */
	if (bMediumP)
	{
		val0 = dMediumVal0P;
		val1 = dMediumVal1P;
		val2 = dMediumVal2P;
		exposure = dMediumExposureP;
		saturation = dMediumSaturationP;
		gamma = dMediumGammaP;
		medium = tonemapLAB(I_Lab, L0, L1,
							val0, val1, val2,
							exposure, gamma, saturation);
		count++;
	}

	/** coarse details */
	if (bCoarseP)
	{
		val0 = dCoarseVal0P;
		val1 = dCoarseVal1P;
		val2 = dCoarseVal2P;
		exposure = dCoarseExposureP;
		saturation = dCoarseSaturationP;
		gamma = dCoarseGammaP;
		coarse = tonemapLAB(I_Lab, L0, L1,
							val0, val1, val2,
							exposure, gamma, saturation);
		count++;
	}

	L0.release();
	L1.release();
	I_Lab.release();

	if (count < 1.0)
		count = 1.0;
	cv::Mat combined = (fine + medium + coarse) / count;
	fine.release();
	medium.release();
	coarse.release();

	/** normalize to 0-1 range */
	cv::normalize(combined, combined, 0, 1, cv::NORM_MINMAX, CV_32FC3);

	/** output result */
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height); /** provide progress bar */
		for (int i = 0; i < width; i++)
		{
			*pDestinationData++ = combined.at<cv::Vec3f>(j, i)[0];
			*pDestinationData++ = combined.at<cv::Vec3f>(j, i)[1];
			*pDestinationData++ = combined.at<cv::Vec3f>(j, i)[2];
		}
	}

	combined.release();

	pSrc->ProgressBar(j, pSrc->GetHeight());

	return 0;
}

/**
 * @brief Edge-Preserving Decompositions for Multi-Scale Tone and Detail Manipulation
 * function sigmoid, rewritten from example matlab code at
 * http://www.cs.huji.ac.il/~danix/epd/msdm-example.zip
 *
 * Applies a sigmoid function on the data X in [0-1] range.
 * Then rescales the result so 0.5 will be mapped to itself. 
 * 
 * @param X 
 * @param a
 */
cv::Mat TMOFarbman08::sigmoid(cv::Mat X, double a)
{
	double x, y, y05;
	cv::Mat Y(X.size(), X.type());
	y05 = 1.0 / (1 + exp(-a * 0.5)) - 0.5;
	for (int j = 0; j < X.rows; j++)
	{
		pSrc->ProgressBar(j, X.rows); /** provide progress bar */
		for (int i = 0; i < X.cols; i++)
		{
			x = X.at<float>(j, i);
			/** apply sigmoid */
			y = 1.0 / (1 + exp(-a * x)) - 0.5;
			/** re-scale */
			Y.at<float>(j, i) = y * (0.5 / y05);
		}
	}
	return Y;
}

/**
 * @brief Edge-Preserving Decompositions for Multi-Scale Tone and Detail Manipulation
 * function tonemapLAB, rewritten from example matlab code at
 * http://www.cs.huji.ac.il/~danix/epd/msdm-example.zip
 * 
 * This function gets an image in the CIELAB color
 * space and tone maps it according to the parameters.
 * 
 * @param Lab is the image in CIELAB color space
 * @param L0 smoothed version of L of LAB
 * @param L1 smoothed version of L of LAB
 * @param val0 compression/expansion params in [-1, 1] range
 * @param val1 compression/expansion params in [-1, 1] range
 * @param val2 compression/expansion params in [-1, 1] range
 * @param exposure is in [0,inf) range
 * @param gamma is in (0,1] range
 * @param saturation is in [0,inf) range
 * @return returns RGB image
 */
cv::Mat TMOFarbman08::tonemapLAB(cv::Mat Lab, cv::Mat L0, cv::Mat L1,
								 double val0, double val1, double val2,
								 double exposure, double gamma, double saturation)
{
	cv::Mat Lab_channels[3];
	cv::split(Lab, Lab_channels);
	cv::Mat &L = Lab_channels[0];
	cv::Mat &a = Lab_channels[1];
	cv::Mat &b = Lab_channels[2];

	cv::Mat diff0, diff1, base;

	/** L's are in range 0-100 */
	if (val0 == 0)
		diff0 = L - L0;
	else if (val0 > 0)
		diff0 = this->sigmoid((L - L0) / 100, val0) * 100;
	else if (val0 < 0)
		diff0 = (1 + val0) * (L - L0);

	if (val1 == 0)
		diff1 = L0 - L1;
	else if (val1 > 0)
		diff1 = this->sigmoid((L0 - L1) / 100, val1) * 100;
	else if (val1 < 0)
		diff1 = (1 + val1) * (L0 - L1);

	if (val2 == 0)
		base = exposure * L1;
	else if (val2 > 0)
		base = (this->sigmoid((exposure * L1 - 56) / 100, val2) * 100) + 56;
	else if (val2 < 0)
		base = (1 + val2) * (exposure * L1 - 56) + 56;

	if (gamma == 1.0)
		L = base + diff1 + diff0;
	else
	{
		double minBase, maxBase;
		cv::minMaxLoc(base, &minBase, &maxBase);
		cv::normalize(base, base, 0, 1, cv::NORM_MINMAX, base.type());
		cv::pow(base, gamma, base);
		L = base * maxBase + diff1 + diff0;
	}
	diff0.release();
	diff1.release();
	base.release();

	if (saturation > 0.0)
	{
		a *= saturation;
		b *= saturation;
	}

	cv::Mat Lab_res;
	cv::merge(Lab_channels, 3, Lab_res);
	Lab_channels[0].release();
	Lab_channels[1].release();
	Lab_channels[2].release();

	cv::Mat RGB_res;
	cv::cvtColor(Lab_res, RGB_res, cv::COLOR_Lab2RGB);
	Lab_res.release();

	return RGB_res;
}
