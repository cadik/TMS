/* --------------------------------------------------------------------------- *
 * TMOMikamo14.cpp: implementation of the TMOMikamo14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMikamo14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMikamo14::TMOMikamo14()
{
	SetName(L"Mikamo14");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOMikamo14::~TMOMikamo14()
{
}

std::vector<double> TMOMikamo14::RGBtoSpectrum(int x, int y)
{
	// Return vector
	std::vector<double> spectrum;
	// Get the RGB values of the pixel
	double red = pSrc->GetPixel(x, y)[0];
	double green = pSrc->GetPixel(x, y)[1];
	double blue = pSrc->GetPixel(x, y)[2];

	for (int binIndex = 0; binIndex < bins; binIndex++)
	{
		// Access SPDs from color_data for the specific bin (bin_index)
		double white_spd   = color_data[binIndex][White];
		double cyan_spd    = color_data[binIndex][Cyan];
		double magenta_spd = color_data[binIndex][Magenta];
		double yellow_spd  = color_data[binIndex][Yellow];
		double red_spd     = color_data[binIndex][Red];
		double green_spd   = color_data[binIndex][Green];
		double blue_spd    = color_data[binIndex][Blue];
		int spd = 0;
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
	// Create a 3D array to store the spectral representation of the image
	double imageSpectrum[pSrc->GetWidth()][pSrc->GetHeight()][bins];
	// transform the image to spectral representation
	for (int x = 0; x < pSrc->GetWidth(); x++)
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			std::vector<double> pixelSpectrum = RGBtoSpectrum(x, y);
			for (int binIndex = 0; binIndex < bins; binIndex++)
			{
				imageSpectrum[x][y][binIndex] = pixelSpectrum[binIndex];
			}
		}

	// Calculate the average luminance of the image which is also the intensity
	double luminanceAverage = 0;

	for (int x = 0; x < pSrc->GetWidth(); x++)
		for (int y = 0; y < pSrc->GetHeight(); y++)
		{
			luminanceAverage += pSrc->GetLuminance(x, y);
		}

	luminanceAverage /= pSrc->GetHeight() * pSrc->GetWidth();
	double intensity = luminanceAverage;

	// compute the sets of parameters that satisfy the wavelength discrimination
	double lambda_l, lambda_m, lambda_s, k1, k2, k3, k4, k5, k6;

	lambda_l = -18.3 / (1 + 7.2 * pow(intensity, -0.7)) - 0.9;
	lambda_m = -44.6 / (1 + 35.4 * pow(intensity, -1.2)) + 22.0;
	lambda_s = 43.0 / (1 + 9.0 * pow(intensity, -1.5)) + 28.0;
	k1 = 6.69 / (1 + 2500 * pow(intensity, -2.65)) + 0.80;
	k2 = -6.24 / (1 + 2500 * pow(intensity, -2.5)) - 0.77;
	k3 = 0.36 / (1 + 50.02 * pow(intensity, -1.5)) + 0.04;
	k4 = 0.24 / (1 + 50.04 * pow(intensity, -1.7)) + 0.03;
	k5 = 0.42 / (1 + 1.76 * pow(intensity, -0.02)) + 0.14;
	k6 = 0.15 / (1 + 2.8 * pow(intensity, -0.46)) - 0.27;

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py;

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			pY *= dParameter; // Parameters can be used like
							  // simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	return 0;
}
