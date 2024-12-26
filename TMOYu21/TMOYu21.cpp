/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOYu21.h"
#include <fstream>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOYu21::TMOYu21()
{
	SetName(L"Yu21");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOYu21::~TMOYu21()
{
}
/*
	Calculate all three correlation coefficients. In the result
	array are: corrRG, corrGB, corrBR
*/
TMOYu21::SImageStats TMOYu21::computeCorrelationCoefficient()
{
	double *pSourceData(pSrc->GetData());

	SImageStats result;

	// Compute mean values
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			result.meanR += *pSourceData++;
			result.meanG += *pSourceData++;
			result.meanB += *pSourceData++;
		}
	}
	double invNumValues(1.0 / double(pSrc->GetWidth() * pSrc->GetHeight()));
	result.meanR *= invNumValues;
	result.meanG *= invNumValues;
	result.meanB *= invNumValues;

	pSourceData = pSrc->GetData();

	// Finalize correlation coefficients
	double numeratorRG(0), numeratorGB(0), numeratorBR(0);
	double denominatorR(0), denominatorG(0), denominatorB(0);
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

			double diffR = pR - result.meanR;
			double diffG = pG - result.meanG;
			double diffB = pB - result.meanB;

			numeratorRG += diffR * diffG;
			numeratorGB += diffG * diffB;
			numeratorBR += diffB * diffR;

			denominatorR += diffR * diffR;
			denominatorG += diffG * diffG;
			denominatorB += diffB * diffB;
		}
	}

	//Just in that highly unlike case, that denominator would be exactly 0
	if(denominatorR == 0.0 || denominatorG == 0.0 || denominatorB == 0.0)
	{
		throw std::runtime_error("Standard deviation is zero.");
	}

	result.stddevR = std::sqrt(denominatorR * invNumValues);
	result.stddevG = std::sqrt(denominatorG * invNumValues);
	result.stddevB = std::sqrt(denominatorB * invNumValues);

	result.covRG = numeratorRG * invNumValues;
	result.covGB = numeratorGB * invNumValues;
	result.covBR = numeratorBR * invNumValues;

	result.Krg = numeratorRG / std::sqrt(denominatorR * denominatorG);
	result.Kgb = numeratorGB / std::sqrt(denominatorG * denominatorB);
	result.Kbr = numeratorBR / std::sqrt(denominatorB * denominatorR);

	/* Krg/gb/br can be negative number. But author does not say how to work with that
	This is attempt to get their contrast picture...*/
	result.Krg = (result.Krg+1)/2;
	result.Kgb = (result.Kgb+1)/2;
	result.Kbr = (result.Kbr+1)/2;
	
	/*result.Krg = std::abs(result.Krg);
	result.Kgb = std::abs(result.Kgb);
	result.Kbr = std::abs(result.Kbr);*/

	return result;
}

//Function for computing contrast image
TMOYu21::CImagePlusStats TMOYu21::createContrastImage(const SImageStats &imageStatistics)
{
	// Create data output
	CImagePlusStats result;

	result.contrastPicture = std::make_unique<std::vector<double>>(pSrc->GetWidth() * pSrc->GetHeight());
	result.meanC = 0;
	result.stddevC = 0;
	double *pSourceData(pSrc->GetData());

	double Krg(imageStatistics.Krg);
	double Kgb(imageStatistics.Kgb);
	double Kbr(imageStatistics.Kbr);

	// Fill contrast picture and compute mean
	auto itOut = result.contrastPicture->begin();

	double min(1000), max(-1000);
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

			*itOut = 0.5*(1+(0.5 * (Krg * (pR + pG) + Kgb * (pG + pB) + Kbr * (pB + pR)))); 
			min = std::min(min, *itOut);
			max = std::max(max, *itOut);

			result.meanC += *itOut;
			++itOut;
		}
	}

	// Normalize image, so we get rid of negative
	/*if(max > min)
	{
		result.meanC = 0.0;
		itOut = result.contrastPicture->begin();
		double divider(1.0 / (max - min));
		for (int j = 0; j < pSrc->GetHeight(); j++)
		{
			for (int i = 0; i < pSrc->GetWidth(); i++)
			{
				*itOut = (*itOut - min) * divider;
				result.meanC += *itOut;
				++itOut; 
			}
		}
	}*/

	//Computing mean
	double invNumValues(1.0 / double(pSrc->GetWidth() * pSrc->GetHeight()));
	result.meanC *= invNumValues;

	//Computing standard deviation
	double denominator(0);

	auto iteContrast = result.contrastPicture->begin();

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pixel = *iteContrast;
			double diff = pixel - result.meanC;
			denominator += diff * diff;
			++iteContrast;
		}
	}

	//Just in that highly unlike case, that denominator would be exactly 0
	if(denominator == 0.0)
	{
		throw std::runtime_error("Standard deviation is zero.");
	}

	result.stddevC = std::sqrt(denominator * invNumValues);

	return result;
}

//Computing cov(RC, GC, BC)
std::array<double, 3> TMOYu21::computeCovContrastRGB(const SImageStats &imageStatistics, const CImagePlusStats &contrastImageStat)
{
	// Create data output
	std::array<double, 3> result;

	double *pSourceData(pSrc->GetData());
	double numeratorRC(0), numeratorGC(0), numeratorBC(0);

	auto iteContrast = contrastImageStat.contrastPicture->begin();
	
	//Computing cov(c,r)(c,g)(c,b)
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pC = *iteContrast;

			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

			double diffR = pR - imageStatistics.meanR;
			double diffG = pG - imageStatistics.meanG;
			double diffB = pB - imageStatistics.meanB;
			double diffC = pC - contrastImageStat.meanC;
			

			numeratorRC += diffR * diffC;
			numeratorGC += diffG * diffC;
			numeratorBC += diffB * diffC;
			++iteContrast;
		}
	}
	double invNumValues(1.0 / double(pSrc->GetWidth() * pSrc->GetHeight()));

	result[0] = numeratorRC * invNumValues;
	result[1]= numeratorGC * invNumValues;
	result[2] = numeratorBC * invNumValues;

	return result;
}

//Computing SSIM(R,C)(G,C)(B,C) for this picture
std::array<double, 3> TMOYu21::computeSSIM(const SImageStats &imageStatistics, const CImagePlusStats &contrastImageStat)
{
	// Create data output
	std::array<double, 3> resultSSIM;

	std::array<double, 3> covRC_GC_BC = computeCovContrastRGB(imageStatistics, contrastImageStat);
	//Constant to prevent dividing by zero, small enough to be ignored
	double C1 = 0.0001; //Typical constant for interval [0,1] is 0.0001 or 0.0003

	//Compute l, d and s for SSIM for each combination of RGB pictures and Contrast picture

	//Combination R from RGB and contrast
	double l = (2 * imageStatistics.meanR * contrastImageStat.meanC + C1) /
		(imageStatistics.meanR*imageStatistics.meanR + contrastImageStat.meanC* contrastImageStat.meanC + C1);
	double d = (2 * imageStatistics.stddevR * contrastImageStat.stddevC + C1) / 
		(imageStatistics.stddevR*imageStatistics.stddevR + contrastImageStat.stddevC* contrastImageStat.stddevC + C1);
	double s = (covRC_GC_BC[0] + C1) / (2 * imageStatistics.stddevR * contrastImageStat.stddevC + C1);

	resultSSIM[0] = l*d*s;

	//Combination G from RGB and contrast
	l = (2 * imageStatistics.meanG * contrastImageStat.meanC + C1) /
		(imageStatistics.meanG*imageStatistics.meanG + contrastImageStat.meanC* contrastImageStat.meanC + C1);
	d = (2 * imageStatistics.stddevG * contrastImageStat.stddevC + C1) / 
		(imageStatistics.stddevG*imageStatistics.stddevG + contrastImageStat.stddevC* contrastImageStat.stddevC + C1);
	s = (covRC_GC_BC[1] + C1) / (2 * imageStatistics.stddevG * contrastImageStat.stddevC + C1);

	resultSSIM[1] = l*d*s;

	//Combination B from RGB and contrast
	l = (2 * imageStatistics.meanB * contrastImageStat.meanC + C1) /
		(imageStatistics.meanB*imageStatistics.meanB + contrastImageStat.meanC* contrastImageStat.meanC + C1);
	d = (2 * imageStatistics.stddevB * contrastImageStat.stddevC + C1) / 
		(imageStatistics.stddevB*imageStatistics.stddevB + contrastImageStat.stddevC* contrastImageStat.stddevC + C1);
	s = (covRC_GC_BC[2] + C1) / (2 * imageStatistics.stddevB * contrastImageStat.stddevC + C1);

	resultSSIM[2] = l*d*s;


	
	// Attempt to get rid of negative values by giving negative SSIM to zero 
	// This makes SSIM sometimes zero, that makes some of kc zero, co sum of kc's is not 1
	/*for(auto &v : resultSSIM)
		v = std::clamp(v, 0.0, 1.0);	
	*/
	
	//Second attempt to solve problem with negative SSIM - shift all values form [-1,1] to [0,1] 
	for(auto &v : resultSSIM)
		v = (v + 1) * 0.5;

	
	return resultSSIM;
}

//Computing constants kr, kg and kb for this picture
std::array<double, 3> TMOYu21::computeK(const SImageStats &imageStatistics)
{
	CImagePlusStats contrastImageStat = createContrastImage(imageStatistics);
	std::array<double, 3> SSIM_RC_GC_BC =computeSSIM(imageStatistics, contrastImageStat);

	std::array<double, 3> result;

	double invSSIMadition = 1/(SSIM_RC_GC_BC[0] + SSIM_RC_GC_BC[1] + SSIM_RC_GC_BC[2]);

	// This is correct, but I am testing pixture with kc = 1
	result[0] = SSIM_RC_GC_BC[0] * invSSIMadition;
	result[1] = SSIM_RC_GC_BC[1] * invSSIMadition;
	result[2] = SSIM_RC_GC_BC[2] * invSSIMadition;

	/*result[0] = 1,
	result[1] = 1;
	result[2] = 1;*/
	return result;
}

//Function for storing and getting pixels in creating resized picture
inline double TMOYu21::getPixel(const double* data, int width, int x, int y, int channel) {
    return data[(y * width + x) * 3 + channel];
}

inline void TMOYu21::setPixel(double* data, int width, int x, int y, int channel, double value) {
    data[(y * width + x) * 3 + channel] = value;
}

//Create resized picture
std::unique_ptr<double[]> TMOYu21::resizeImage(const double* input, int srcWidth, int srcHeight, int destWidth, int destHeight)
{
    //Allocating new picture
    std::unique_ptr<double[]> output(new double[destWidth * destHeight * 3]);

    for (int y = 0; y < destHeight; ++y) {
        for (int x = 0; x < destWidth; ++x) {
            // Recalculating pixels from original to output picture
            double srcX = x * (static_cast<double>(srcWidth) / destWidth);
            double srcY = y * (static_cast<double>(srcHeight) / destHeight);

            // Surrounding pixels
            int x0 = std::clamp<double>(static_cast<int>(std::floor(srcX)), 0, srcWidth - 1);
            int x1 = std::clamp<double>(std::min(x0 + 1, srcWidth - 1), 0, srcWidth - 1);
            int y0 = std::clamp<double>(static_cast<int>(std::floor(srcY)), 0, srcHeight - 1);
            int y1 = std::clamp<double>(std::min(y0 + 1, srcHeight - 1), 0, srcHeight - 1);

            // Interpolation
            double dx = srcX - x0;
            double dy = srcY - y0;

            for (int c = 0; c < 3; ++c) { // R, G, B canals
                double value =
                    (1 - dx) * (1 - dy) * getPixel(input, srcWidth, x0, y0, c) +
                    dx * (1 - dy) * getPixel(input, srcWidth, x1, y0, c) +
                    (1 - dx) * dy * getPixel(input, srcWidth, x0, y1, c) +
                    dx * dy * getPixel(input, srcWidth, x1, y1, c);

                setPixel(output.get(), destWidth, x, y, c, value);
            }
        }
    }

    return output;
}

// Compuing intenstity difference from random and neighboring pairs
std::shared_ptr<std::vector<double>> TMOYu21::computeContrastDifferences(const double* image64, const double* image32, int channel) 
{
   
   auto contrastDifferences = std::make_shared<std::vector<double>>();

    // Random pairs from 64*64 picture
    for (int y = 0; y < 64; ++y) {
        for (int x = 0; x < 64; ++x) {
            // Random pixel
            int randomX = rand() % 64;
            int randomY = rand() % 64;

            // Contrast difference
			double pix1  = getPixel(image64, 64, x, y, channel);
			double pix2 = getPixel(image64, 64, randomX, randomY, channel);
			// This can be negative MAYBE ABS
            double diff = getPixel(image64, 64, x, y, channel) - getPixel(image64, 64, randomX, randomY, channel);
            contrastDifferences->push_back(diff);
        }
    }
    // Neighboring pairs from 32*32
    for (int y = 0; y < 32; ++y) {
        for (int x = 0; x < 32; ++x) {
			// For each neighbor
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
					// If not the same pixel, copute contr. difference
					if ((dx == 0 || dy == 0) && (dx != 0 || dy != 0) &&
            				x + dx >= 0 && x + dx < 32 && y + dy >= 0 && y + dy < 32) {
						// This can be negative MAYBE ABS
                        double diff = getPixel(image32, 32, x, y, channel) - getPixel(image32, 32, x + dx, y + dy, channel);
                        contrastDifferences->push_back(diff);
                    }
                }
            }
        }
    }
    return contrastDifferences;
}


//Pomocná funkce pro ukládání větších výpisů do souboru
void TMOYu21::logDebug(const std::string& message) {
	static bool firstRun = true;
    std::ofstream logFile("chyba.txt", firstRun ? (std::ios::out | std::ios::trunc) : std::ios::app);
    firstRun = false;
    if (logFile.is_open()) {
        logFile << message << std::endl;
    }
}

//Computing weights wr, wg, wb
std::array<double, 3> TMOYu21::computeWeights(const std::vector<double> &allIr, const std::vector<double> &allIg,
				const std::vector<double> &allIb, const std::array<double, 3> &kr_kg_kb)
{
	std::array<double, 3> result_wr_wg_wb;


	//For each from 66 w-rgb combination ()
	const double step = 0.1;
    const int precision = 1; // One decimal place

	const double kr = kr_kg_kb[0];
	const double kg = kr_kg_kb[1];
	const double kb = kr_kg_kb[2];

	double countR = 0;
	double minEnergy = std::numeric_limits<double>::max();
    for (double wr = 0.0; wr <= 1.0; wr += step) {
        for (double wg = 0.0; wg <= 1.0 - wr; wg += step) {
            double wb = 1.0 - wr - wg;

            // Check if wb is a valid value (multiple of step and within range)
            if (wb >= 0.0 && wb <= 1.0) 
			{
				double energyR = computeColorEnergy({wr, wg, wb}, kr, {allIr, allIg, allIb}, 0);
				double energyG = computeColorEnergy({wr, wg, wb}, kg, {allIr, allIg, allIb}, 1);
				double energyB = computeColorEnergy({wr, wg, wb}, kb, {allIr, allIg, allIb}, 2);
				
				// Sum all 3 colors weight and contrast difference
				double totalEnergy = energyR + energyG + energyB;

				if (totalEnergy < minEnergy)
				{
					minEnergy = totalEnergy;
					result_wr_wg_wb[0] = wr;
					result_wr_wg_wb[1] = wg;
					result_wr_wg_wb[2] = wb;
				}
  			}
		}
    }

	return result_wr_wg_wb;
}

double TMOYu21::computeColorEnergy(const std::array<double, 3> &w, double k, const std::array<std::vector<double>, 3> &I, size_t colorIndex)
{
	const double epsilon = 0.15f;  //According the authors, the best constant for Cadik dataset
	double energy = 0.0;
	size_t numPairs = I[0].size();

	for (size_t i = 0; i < numPairs; ++i) 
	{
		const double Ir = I[0][i];
		const double Ig = I[1][i];
		const double Ib = I[2][i];

		const double wr = w[0];
		const double wg = w[1];
		const double wb = w[2];

		double weightedSum = wr * Ir + wg * Ig + wb * Ib;
		double absWeightedSum = std::abs(weightedSum);

		double numerator = absWeightedSum - k * std::abs(I[colorIndex][i])  - epsilon;
		double denominator = absWeightedSum + k * std::abs(I[colorIndex][i]) + epsilon;

		energy += std::abs(numerator) / denominator;
	}

	return energy;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOYu21::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py;

	SImageStats imageStatistics = computeCorrelationCoefficient();
	
	std::array<double, 3> kr_kg_kb = computeK(imageStatistics);

	std::unique_ptr<double[]> resized32 = resizeImage(pSourceData, pSrc->GetWidth(), pSrc->GetHeight(), 32, 32);
	std::unique_ptr<double[]> resized64 = resizeImage(pSourceData, pSrc->GetWidth(), pSrc->GetHeight(), 64, 64);

	{	
		auto img32 = createImage(resized32.get(), 32, 32);
		img32->SaveAs("../../resized32.png", TMO_PNG_8);
		auto img64 = createImage(resized64.get(), 64, 64);
		img64->SaveAs("../../resized64.png", TMO_PNG_8);
	}
	auto allIr = computeContrastDifferences(resized64.get(), resized32.get(), 0);
	auto allIg = computeContrastDifferences(resized64.get(), resized32.get(), 1);
	auto allIb = computeContrastDifferences(resized64.get(), resized32.get(), 2);

	auto wr_wg_wb = computeWeights(*allIr, *allIg, *allIb, kr_kg_kb);
	
	int j = 0;
	int k = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		//pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			
			auto R = *pSourceData++;
			auto G = *pSourceData++;
			auto B = *pSourceData++;

			auto intensity = wr_wg_wb[0] * R + wr_wg_wb[1] * G + wr_wg_wb[2] * B;
			//auto intensity =  0.5 * R + 0.5 * G + 0 * B; // Correct for picture with girl
			//auto intensity =  0.8 * R + 0.2 * G + 0 * B; // Correct for picture with girl - kc = 1
			*pDestinationData++ = intensity;
			*pDestinationData++ = intensity;
			*pDestinationData++ = intensity;


		}
	}


	//pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}

std::unique_ptr<TMOImage> TMOYu21::createImage(const double *data, int width, int height)
{
	auto pImage = std::make_unique<TMOImage>();
	pImage->New(width, height);
	auto dataCopy = new double[width * height * 3];
	memcpy(dataCopy, data, width * height * 3 * sizeof(double));

	pImage->SetData(dataCopy);

	return pImage;
}