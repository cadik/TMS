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

	ri.SetName(L"ri");				// TODO - Insert parameters names
	ri.SetDescription(L"Retinal illuminance in Trolands (Td)\n"
										"Photopic vision: RI >= 10\n"
										"Mesopic vision: 10 >= RI >= 0.01\n"
										"Scotopic vision: 0.01 >= RI"); // TODO - Insert parameter descriptions
	ri.SetDefault(150.0);							// TODO - Add default values
	ri = 150.0;
	ri.SetRange(0.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(ri);
}

TMOMikamo14::~TMOMikamo14()
{
}

// Function to compute adapted luminance (I)
double TMOMikamo14::computeAdaptedLuminance(int x, int y) {
    // Compute the scene luminance from the RGB values of the pixel
    double L = pSrc->GetLuminance(x, y);
    
    // Compute the pupil area based on the luminance
    // Avoid very small or negative luminance values
    if (L < 0.001) L = 0.001;

    // Empirical formula to compute pupil diameter
    double pupilDiameter = (5.0 / 3.0) * std::exp(-0.1 * std::log(L));
    
    // Convert diameter to area (A = π * r^2)
    double A = M_PI * std::pow(pupilDiameter / 2.0, 2);
    
    // Adapted luminance I = A * L
    return A * L;
}

std::vector<double> TMOMikamo14::computeSigmoidParams(double I) {
    // Compute parameter functions for different levels of adaptation
    std::vector<double> params(9);

    params[0] = -18.3 / (1 + 7.2 * std::pow(I, -0.7)) - 0.9; // λl(I)
    params[1] = -44.6 / (1 + 35.4 * std::pow(I, -1.2)) + 22.0; // λm(I)
    params[2] = 43.0 / (1 + 9.0 * std::pow(I, -1.5)) + 28.0;  // λs(I)

    params[3] = 6.69 / (1 + 2500 * std::pow(I, -2.65)) + 0.80; // k1(I)
    params[4] = -6.24 / (1 + 2500 * std::pow(I, -2.5)) - 0.77; // k2(I)
    params[5] = 0.36 / (1 + 50.02 * std::pow(I, -1.5)) + 0.04; // k3(I)

    params[6] = 0.24 / (1 + 50.04 * std::pow(I, -1.7)) + 0.03; // k4(I)
    params[7] = 0.42 / (1 + 1.76 * std::pow(I, -0.02)) + 0.14; // k5(I)
    params[8] = 0.15 / (1 + 2.80 * std::pow(I, -0.46)) - 0.27; // k6(I)

    return params;
}

double TMOMikamo14::getIntensityFromWavelength(std::vector<double> spectrum, double lambda) {
	// spectrum is a vector of 10 values, each representing the SPD at a specific wavelength (380nm to 720nm)
	double position = lambda - 380.0 / 340.0; // - first wavelength / wavelength range

	if (position < 0) {
		return spectrum[0];
	}

	if (position >= 9) {
		return spectrum[9];
	}

	int i = static_cast<int>(std::floor(position));
	double fraction = position - i;

	return spectrum[i] + (spectrum[i + 1] - spectrum[i]) * fraction;
}

double TMOMikamo14::getWavelengthFromIntensity(std::vector<double> spectrum, double intensity) {
	// spectrum is a vector of 10 values, each representing the SPD at a specific wavelength (380nm to 720nm)
	for (int i = 0; i < 9; i++) {
		if (spectrum[i] <= intensity && spectrum[i + 1] >= intensity) {
			double fraction = (intensity - spectrum[i]) / (spectrum[i + 1] - spectrum[i]);
			return 380.0 + 340.0 * (i + fraction);
		}
	}

	return 720.0;
}

TMOMikamo14::OpponentColor TMOMikamo14::applyTwoStageModel(std::vector<double> spectrum, double I) {
	std::vector<double> paramsOld = computeSigmoidParams(I);
	std::vector<double> paramsNew = computeSigmoidParams(ri);

	// TODO: How to calculate Cl, Cm, Cs?
	double L_wavelength = getWavelengthFromIntensity(spectrum, I);
	double M_wavelength = getWavelengthFromIntensity(spectrum, I);
	double S_wavelength = getWavelengthFromIntensity(spectrum, I);


    double L_adj = getIntensityFromWavelength(spectrum, L_wavelength - paramsOld[0]);
    double M_adj = getIntensityFromWavelength(spectrum, M_wavelength - paramsOld[1]);
    double S_adj = getIntensityFromWavelength(spectrum, S_wavelength - paramsOld[2]);
    
    OpponentColor color;
    color.achromatic = 0.6 * L_adj + 0.4 * M_adj;
    color.redGreen = paramsNew[3] * L_adj + paramsNew[4] * M_adj + paramsNew[5] * S_adj;
    color.yellowBlue = paramsNew[6] * L_adj + paramsNew[7] * M_adj + paramsNew[8] * S_adj;

    return color;
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
			double* pixel = pSrc->GetPixel(x, y);
			double I = computeAdaptedLuminance(x, y);
			std::vector<double> spectrum = RGBtoSpectrum(x, y);
			OpponentColor oColor = applyTwoStageModel(spectrum, I);

			double* outputPixel = pDst->GetPixel(x, y);
			outputPixel[0] = oColor.achromatic;
			outputPixel[1] = oColor.redGreen;
			outputPixel[2] = oColor.yellowBlue;
		}
	}
	pDst->Convert(TMO_RGB);

	return 0;
}
