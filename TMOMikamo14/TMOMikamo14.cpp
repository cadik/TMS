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
	ri.SetDescription(L"Adapted retinal illuminance (ri) in Trolands; <0.0, 1000.0>");
	ri.SetDefault(0.0);
	ri = 0.0;
	ri.SetRange(0.0, 1000.0);
	this->Register(ri);

	al.SetName(L"al");
	al.SetDescription(L"Adapted luminance (al) in cd/m^2; <0.0, 1000.0>");
	al.SetDefault(0.0);
	al = 0.0;
	al.SetRange(0.0, 1000.0);
	this->Register(al);
}

TMOMikamo14::~TMOMikamo14()
{
}

// Function to compute adapted luminance (I)
double TMOMikamo14::computeAdaptedLuminance()
{
	if (ri != 0.0)
	{
		return ri;
	}

	double A = M_PI * std::pow(5.0 / 2.0, 2); // Assume a pupil diameter of 5mm

	if (al != 0.0)
	{
		return al * A;
	}

	double luminanceSum = 0.0;

	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double L = pSrc->GetLuminance(x, y);
			luminanceSum += L;
		}
	}

	double averageLuminance = luminanceSum / (pSrc->GetHeight() * pSrc->GetWidth());

	return averageLuminance * A;
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

std::vector<double> TMOMikamo14::lambdaAdjust(int cone, double lambdaDiff)
{
	int step = lambdaDiff / binWidth;
	std::vector<double> newSpectralSensitivity(bins, 0.0);
	for (int i = 0; i < bins; i++)
	{
		if (i - step >= 0 && i - step < bins)
		{
			newSpectralSensitivity[i] = LMSsensitivities[i - step][cone];
		}
	}

	return newSpectralSensitivity;
}

cv::Mat TMOMikamo14::applyTwoStageModel(std::vector<double> spd, double I)
{
	std::vector<double> params = computeSigmoidParams(I);

	double V = 0.0;
	double Org = 0.0;
	double Oyb = 0.0;
	std::vector<double> Cl = lambdaAdjust(0, params[0]);
	std::vector<double> Cm = lambdaAdjust(1, params[1]);
	std::vector<double> Cs = lambdaAdjust(2, params[2]);

	for (int i = 0; i < bins; i++)
	{
		cv::Mat CmClCs = (cv::Mat_<double>(3, 1) << Cl[i], Cm[i], Cs[i]);
		cv::Mat M = (cv::Mat_<double>(3, 3) << 0.6, 0.4, 0.0, params[3], params[4], params[5], params[6], params[7], params[8]);
		cv::Mat z = M * CmClCs;
		V += spd[i] * z.at<double>(0, 0) * binWidth;
		Org += spd[i] * z.at<double>(1, 0) * binWidth;
		Oyb += spd[i] * z.at<double>(2, 0) * binWidth;
	}

	cv::Mat opponentColor = (cv::Mat_<double>(3, 1) << V, Org, Oyb);
	std::vector<double> newParams = computeSigmoidParams(150.0);
	cv::Mat invM = (cv::Mat_<double>(3, 3) << 0.6, 0.4, 0.0, newParams[3], newParams[4], newParams[5], newParams[6], newParams[7], newParams[8]);
	invM = invM.inv();
	opponentColor = invM * opponentColor;

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
	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();
	double I = computeAdaptedLuminance();
	for (int y = 0; y < pSrc->GetHeight(); y++)
	{
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
			double *pixel = pSrc->GetPixel(x, y);
			std::vector<double> spd = RGBtoSpectrum(*pSourceData++, *pSourceData++, *pSourceData++);
			cv::Mat opponentColor = applyTwoStageModel(spd, I);

			*pDestinationData++ = opponentColor.at<double>(0, 0) / 100;
			*pDestinationData++ = opponentColor.at<double>(1, 0);
			*pDestinationData++ = opponentColor.at<double>(2, 0);
		}
	}
	pDst->Convert(TMO_RGB);

	return 0;
}
