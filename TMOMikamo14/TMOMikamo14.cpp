/* --------------------------------------------------------------------------- *
 * TMOMikamo14.cpp: implementation of the TMOMikamo14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMikamo14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMikamo14::TMOMikamo14()
{
	SetName(L"Mikamo14");
	SetDescription(L"Add your TMO description here");

	ri.SetName(L"ri");
	ri.SetDescription(L"Retinal illuminance in Trolands (Td)\n"
					  "Photopic vision: RI >= 10\n"
					  "Mesopic vision: 10 >= RI >= 0.01\n"
					  "Scotopic vision: 0.01 >= RI");
	ri.SetDefault(150.0);
	ri = 150.0;
	ri.SetRange(0.0, 1000.0);
	this->Register(ri);
}

TMOMikamo14::~TMOMikamo14()
{
}

// Function to compute adapted luminance (I)
double TMOMikamo14::computeAdaptedLuminance(int x, int y)
{
	// Compute the scene luminance from the RGB values of the pixel
	double L = pSrc->GetLuminance(x, y);

	// Compute the pupil area based on the luminance
	// Avoid very small or negative luminance values
	if (L < 0.001)
		L = 0.001;

	// Empirical formula to compute pupil diameter
	double pupilDiameter = (5.0 / 3.0) * std::exp(-0.1 * std::log(L));

	// Convert diameter to area (A = π * r^2)
	double A = M_PI * std::pow(pupilDiameter / 2.0, 2);

	// Adapted luminance I = A * L
	return A * L;
}

std::vector<double> TMOMikamo14::computeSigmoidParams(double I)
{
	// Compute parameter functions for different levels of adaptation
	std::vector<double> params(9);

	params[0] = -18.3 / (1 + 7.2 * std::pow(I, -0.7)) - 0.9;   // λl(I)
	params[1] = -44.6 / (1 + 35.4 * std::pow(I, -1.2)) + 22.0; // λm(I)
	params[2] = 43.0 / (1 + 9.0 * std::pow(I, -1.5)) + 28.0;   // λs(I)

	params[3] = 6.69 / (1 + 2500 * std::pow(I, -2.65)) + 0.80; // k1(I)
	params[4] = -6.24 / (1 + 2500 * std::pow(I, -2.5)) - 0.77; // k2(I)
	params[5] = 0.36 / (1 + 50.02 * std::pow(I, -1.5)) + 0.04; // k3(I)

	params[6] = 0.24 / (1 + 50.04 * std::pow(I, -1.7)) + 0.03; // k4(I)
	params[7] = 0.42 / (1 + 1.76 * std::pow(I, -0.02)) + 0.14; // k5(I)
	params[8] = 0.15 / (1 + 2.80 * std::pow(I, -0.46)) - 0.27; // k6(I)

	return params;
}

int TMOMikamo14::findBinFromI(std::vector<double>spectrum, double I)
{
	int bin = 0;
	double max = 0.0;

	for (int i = 0; i < bins; i++)
	{
		if (spectrum[i] > max)
		{
			max = spectrum[i];
			bin = i;
		}
	}

	return bin;
}

int TMOMikamo14::findNewBin(int oldBin, double diff)
{
	if (diff > 0.0)
	{
		int binTop = (oldBin + 1) * binWidth;
		int move = std::round((binTop - diff) / binWidth);
		return oldBin - move;
	}
	else if (diff < 0.0)
	{
		int binBottom = oldBin * binWidth;
		int move = std::round((binBottom - diff) / binWidth);
		return oldBin + move;
	}
	
	return oldBin;
}

cv::Mat TMOMikamo14::applyTwoStageModel(LMS lms, double I)
{
	std::vector<double> params = computeSigmoidParams(I);

	// find wavelength of original LMS values
	int binL = findBinFromI(lms.L, I);
	int binM = findBinFromI(lms.M, I);
	int binS = findBinFromI(lms.S, I);
	int newBinL = findNewBin(binL, params[0]);
	int newBinM = findNewBin(binM, params[1]);
	int newBinS = findNewBin(binS, params[2]);
	double newL = lms.L[newBinL];
	double newM = lms.M[newBinM];
	double newS = lms.S[newBinS];

	cv::Mat opponentColor(3, 1, CV_64F);
	cv::Mat LMS = (cv::Mat_<double>(3, 1) << newL, newM, newS);
	cv::Mat paramsMat = (cv::Mat_<double>(3, 3) << 0.6, 0.4, 0.0, params[3], params[4], params[5], params[6], params[7], params[8]);

	opponentColor = paramsMat * LMS;

	std::vector<double> newParams = computeSigmoidParams(ri);
	cv::Mat newParamsMat = (cv::Mat_<double>(3, 3) << 0.6, 0.4, 0.0, newParams[3], newParams[4], newParams[5], newParams[6], newParams[7], newParams[8]);
	newParamsMat = newParamsMat.inv();
	
	opponentColor = newParamsMat * opponentColor;

	return opponentColor;
}

std::vector<double> TMOMikamo14::RGBtoSpectrum(double red, double green, double blue)
{
	// Return vector
	std::vector<double> spectrum;

	for (int binIndex = 0; binIndex < bins; binIndex++)
	{
		// Access SPDs from color_data for the specific bin (bin_index)
		double white_spd = color_data[binIndex][White];
		double cyan_spd = color_data[binIndex][Cyan];
		double magenta_spd = color_data[binIndex][Magenta];
		double yellow_spd = color_data[binIndex][Yellow];
		double red_spd = color_data[binIndex][Red];
		double green_spd = color_data[binIndex][Green];
		double blue_spd = color_data[binIndex][Blue];
		double spd = 0.0;
		if (red <= green && red <= blue)
		{
			spd += white_spd * red;
			if (green <= blue)
			{
				spd += cyan_spd * (green - red);
				spd += blue_spd * (blue - green);
			}
			else
			{
				spd += cyan_spd * (blue - red);
				spd += green_spd * (green - blue);
			}
		}
		else if (green <= red && green <= blue)
		{
			spd += white_spd * green;
			if (red <= blue)
			{
				spd += magenta_spd * (red - green);
				spd += blue_spd * (blue - red);
			}
			else
			{
				spd += magenta_spd * (blue - green);
				spd += red_spd * (red - blue);
			}
		}
		else
		{
			spd += white_spd * blue;
			if (red <= green)
			{
				spd += yellow_spd * (red - blue);
				spd += green_spd * (green - red);
			}
			else
			{
				spd += yellow_spd * (green - blue);
				spd += red_spd * (red - green);
			}
		}
		spectrum.push_back(spd);
	}
	return spectrum;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMikamo14::Transform()
{
	pDst->Convert(TMO_Yxy);
	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double *pixel = pSrc->GetPixel(x, y);
			double I = computeAdaptedLuminance(x, y);
			LMS lms;
			lms.L = RGBtoSpectrum(pixel[0], 0, 0);
			lms.M = RGBtoSpectrum(0, pixel[1], 0);
			lms.S = RGBtoSpectrum(0, 0, pixel[2]);
			cv::Mat opponentColor = applyTwoStageModel(lms, I);

			double *outputPixel = pDst->GetPixel(x, y);
			outputPixel[0] = opponentColor.at<double>(0, 0);
			outputPixel[1] = opponentColor.at<double>(1, 0);
			outputPixel[2] = opponentColor.at<double>(2, 0);
		}
	}
	pDst->Convert(TMO_RGB);

	return 0;
}
