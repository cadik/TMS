/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*             Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]         *
*                                    Brno 2025                                 *
*                                                                              *
*                      Implementation of the TMOYu21 class                     *
*                                                                              *
*******************************************************************************/


#include "TMOYu21.h"
#include <fstream>

TMOYu21::TMOYu21()
{
	SetName(L"Yu21");					  
	SetDescription(L"Add your TMO description here");
}

TMOYu21::~TMOYu21()
{
}

/*
*	This function calculates the mean, standard deviation, covariance, and correlation 
*  coefficients (Krg, Kgb, Kbr) for the red, green, and blue channels of the input image.
*/
TMOYu21::SImageStats TMOYu21::computeStats()
{
	double *pSourceData(pSrc->GetData());

	SImageStats result;

	// Compute mean values for each channel
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

   // Reset pointer for second pass
	pSourceData = pSrc->GetData();

	// Variables for covariance and variance computation
	double numeratorRG(0), numeratorGB(0), numeratorBR(0);
	double denominatorR(0), denominatorG(0), denominatorB(0);
	
   // Compute covariance and variance
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

	// Check for zero variance to prevent division by zero
	if(denominatorR == 0.0 || denominatorG == 0.0 || denominatorB == 0.0)
	{
		throw std::runtime_error("Standard deviation is zero.");
	}

   // Compute standard deviations
	result.stddevR = std::sqrt(denominatorR * invNumValues);
	result.stddevG = std::sqrt(denominatorG * invNumValues);
	result.stddevB = std::sqrt(denominatorB * invNumValues);

   // Compute covariance values
	result.covRG = numeratorRG * invNumValues;
	result.covGB = numeratorGB * invNumValues;
	result.covBR = numeratorBR * invNumValues;

   // Compute correlation coefficients
	result.Krg = result.covRG / (result.stddevR * result.stddevG);
	result.Kgb = result.covGB / (result.stddevG * result.stddevB);
	result.Kbr = result.covBR / (result.stddevB * result.stddevR);

   // --- EXPERIMENTAL APPROACHES BELOW ---
    
    /*
     * In the article was not specified, how to compute krg, so that I treid different approach.
     * Alternative computation of Krg, Kgb, and Kbr using square root normalization.
     * These versions use the geometric mean of variances instead of standard deviations.
     * The results are similar, but 
     */
    // result.Krg = numeratorRG / sqrt(denominatorR * denominatorG);
    // result.Kgb = numeratorGB / sqrt(denominatorG * denominatorB);
    // result.Kbr = numeratorBR / sqrt(denominatorB * denominatorR);

    /*
     * Handling negative correlation values.
     * The paper does not specify how to solve problem with negative values, so here are two different approaches:
     *
     * 1) Shift correlation values into the range [0, 1] to make them strictly positive.
     * 2) Take the absolute value to remove sign dependence.
     */
    
    // Attempt 1: Normalize to range [0, 1]
    // result.Krg = (result.Krg + 1) / 2;
    // result.Kgb = (result.Kgb + 1) / 2;
    // result.Kbr = (result.Kbr + 1) / 2;

    // Attempt 2: Use absolute values
    // result.Krg = std::abs(result.Krg);
    // result.Kgb = std::abs(result.Kgb);
    // result.Kbr = std::abs(result.Kbr);

	return result;
}

/**
 *
 * This function computes a contrast-enhanced image using correlation coefficients Krg, Kgb, and Kbr.
 * The contrast image is stored in `contrastPicture`, and its mean and standard deviation are also computed.
 *
 * The contrast image is computed as:
 * \f[
 * C(x, y) = 0.5 \cdot (K_{rg} \cdot (R + G) + K_{gb} \cdot (G + B) + K_{br} \cdot (B + R))
 * \f]
 */
TMOYu21::CImagePlusStats TMOYu21::createContrastImage(const SImageStats &imageStatistics)
{
	// Create data output
	CImagePlusStats result;

	result.contrastPicture = std::make_unique<std::vector<double>>(pSrc->GetWidth() * pSrc->GetHeight());
	result.meanC = 0;
	result.stddevC = 0;
	double *pSourceData(pSrc->GetData());

   // Load correlation coefficients
	double Krg(imageStatistics.Krg);
	double Kgb(imageStatistics.Kgb);
	double Kbr(imageStatistics.Kbr);

   /*
     * --- EXPERIMENTAL NORMALIZATION OF CORRELATION COEFFICIENTS ---
     *
     * This section explores adjusting the coefficients before applying them.
     * Two alternative approaches are tested:
     *
     * 1) Shift all coefficients so the smallest becomes zero.
     * 2) Normalize the coefficients by their sum to scale them proportionally.
     */

    {
        // Alternative 1: Shift coefficients so the smallest becomes zero
        // auto min = std::min(std::min(Krg, Kgb), Kbr);
        // Krg -= min;
        // Kgb -= min;
        // Kbr -= min;

        // Alternative 2: Normalize coefficients by their sum
        // auto invSum = Krg + Kgb + Kbr;
        // Krg *= invSum;
        // Kgb *= invSum;
        // Kbr *= invSum;
    }


   // --- COMPUTE CONTRAST IMAGE ---

	auto itOut = result.contrastPicture->begin();

	double min(std::numeric_limits<double>::max());
	double max(-min);

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

         // Compute contrast image based on the original formula from the paper
         *itOut = 0.5 * (Krg * (pR + pG) + Kgb * (pG + pB) + Kbr * (pB + pR)); 


         // --- NORMALIZATION EXPERIMENTS ---
         // The original equation can produce negative contrast values, witch leads to negative values in further processing. 
         // To shift the range into positive values, I experimented with an offset:
         //
         // *itOut = 0.5 * (1 + (Krg * (pR + pG) + Kgb * (pG + pB) + Kbr * (pB + pR)));
         //
         

         // Accumulate for mean computation
			result.meanC += *itOut;
			++itOut;
		}
	}

   // --- NORMALIZATION EXPERIMENTS ---
   // These functions attempt to remap the contrast image to improve visual representation.
   // remapContrastToInputRange(*result.contrastPicture);
   // auto iminmax = getImageMinMax(*pSrc);
	auto cminmax = getContrastImageMinMax(*result.contrastPicture);

   // --- COMPUTE FINAL MEAN ---
	result.meanC = 0.0;
	for(double v : *result.contrastPicture)
	{
		result.meanC += v;
	}

	//Finalize mean computation 
	double invNumValues(1.0 / double(pSrc->GetWidth() * pSrc->GetHeight()));
	result.meanC *= invNumValues;


	// --- COMPUTE STANDARD DEVIATION ---
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

	// Handle unlikely case where variance is zero
	if(denominator == 0.0)
	{
		throw std::runtime_error("Standard deviation is zero.");
	}

	result.stddevC = sqrt(denominator * invNumValues);

	return result;
}
/*
* Computes the covariance between contrast (C) and color channels (R, G, B)
*/
std::array<double, 3> TMOYu21::computeCovContrastRGB(const SImageStats &imageStatistics, const CImagePlusStats &contrastImageStat)
{
	// Output array for cov(C, R), cov(C, G), cov(C, B)
	std::array<double, 3> result;

	double *pSourceData(pSrc->GetData());
	double numeratorRC(0), numeratorGC(0), numeratorBC(0);

	auto iteContrast = contrastImageStat.contrastPicture->begin();
	
	// Compute covariance values for (C, R), (C, G), (C, B)
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

   // Normalize the covariance values
	result[0] = numeratorRC * invNumValues;
	result[1] = numeratorGC * invNumValues;
	result[2] = numeratorBC * invNumValues;
	{
   /*
      // Experimental adjustment: Normalize covariance values 
       
      double min = std::min({result[0], result[1], result[2]});

      // Shift values to prevent negative covariance influence
      result[0] -= min;
      result[1] -= min;
      result[2] -= min;

      // Normalize to ensure the values sum to 1
      double invSum = 1.0 / (result[0] + result[1] + result[2]);
      result[0] *= invSum;
      result[1] *= invSum;
      result[2] *= invSum;
   */
	}
	return result;
}

/*
* Computing SSIM(R,C)(G,C)(B,C) for this picture
*/
std::array<double, 3> TMOYu21::computeSSIM(const SImageStats &imageStatistics, const CImagePlusStats &contrastImageStat)
{
	// Output array for SSIM values
	std::array<double, 3> resultSSIM;

   // Compute covariance values between contrast (C) and color channels (R, G, B)
	std::array<double, 3> covRC_GC_BC = computeCovContrastRGB(imageStatistics, contrastImageStat);

	// Constants to stabilize division (small enough to be negligible in real images)
	double L = 0.5; // Assumed dynamic range of image values
	double K1 = 0.01, K2 = 0.03;
	double C1 = K1 * L; C1 = C1 * C1;
	double C2 = K2 * L; C2 = C2 * C2;
	double C3 = C2 * 0.5;

	// Compute SSIM components for each color channel (R, G, B) compared to contrast (C)

	// SSIM for Red (R) and Contrast (C)
	double l = (2.0 * imageStatistics.meanR * contrastImageStat.meanC + C1) /
		(imageStatistics.meanR*imageStatistics.meanR + contrastImageStat.meanC * contrastImageStat.meanC + C1);
	double d = (2.0 * imageStatistics.stddevR * contrastImageStat.stddevC + C2) / 
		(imageStatistics.stddevR*imageStatistics.stddevR + contrastImageStat.stddevC* contrastImageStat.stddevC + C2);
	double s = (covRC_GC_BC[0] + C3) / (2.0 * imageStatistics.stddevR * contrastImageStat.stddevC + C3);

	resultSSIM[0] = (l*d*s);

	// SSIM for Green (G) and Contrast (C)
	l = (2.0 * imageStatistics.meanG * contrastImageStat.meanC + C1) /
		(imageStatistics.meanG*imageStatistics.meanG + contrastImageStat.meanC* contrastImageStat.meanC + C1);
	d = (2.0 * imageStatistics.stddevG * contrastImageStat.stddevC + C2) / 
		(imageStatistics.stddevG*imageStatistics.stddevG + contrastImageStat.stddevC* contrastImageStat.stddevC + C2);
	s = (covRC_GC_BC[1] + C3) / (2.0 * imageStatistics.stddevG * contrastImageStat.stddevC + C3);

	resultSSIM[1] = (l*d*s);

	// SSIM for Blue (B) and Contrast (C)
	l = (2 * imageStatistics.meanB * contrastImageStat.meanC + C1) /
		(imageStatistics.meanB*imageStatistics.meanB + contrastImageStat.meanC* contrastImageStat.meanC + C1);
	d = (2 * imageStatistics.stddevB * contrastImageStat.stddevC + C2) / 
		(imageStatistics.stddevB*imageStatistics.stddevB + contrastImageStat.stddevC* contrastImageStat.stddevC + C2);
	s = (covRC_GC_BC[2] + C3) / (2.0 * imageStatistics.stddevB * contrastImageStat.stddevC + C3);

	resultSSIM[2] = (l*d*s);


   /*
      // Alternative adjustments to handle negative SSIM values

      // Attempt 1: Take the absolute value of SSIM
      // This may artificially boost SSIM values, leading to incorrect results.
      // for(auto &v : resultSSIM)
      //     v = abs(v);

      // Attempt 2: Clamp negative SSIM values to zero
      // This ensures SSIM is always non-negative, but it may cause some weights to be zero,
      // which could affect the normalization of coefficients.
      // for(auto &v : resultSSIM)
      //     v = std::clamp(v, 0.0, 1.0);

      // Attempt 3: Shift SSIM range from [-1,1] to [0,1]
      // This ensures all values are in a valid range but may distort the original SSIM interpretation.
      // for(auto &v : resultSSIM)
      //     v = (v + 1) * 0.5;
   */

	
	return resultSSIM;
}

/*
* Computing constants kr, kg, and kb for the given image based on SSIM (Structural Similarity Index)
*/
std::array<double, 3> TMOYu21::computeK(const SImageStats &imageStatistics)
{
   // Create a contrast image based on the input image statistics
	CImagePlusStats contrastImageStat = createContrastImage(imageStatistics);

	/*
      // Save the contrast image for analysis (debugging purposes)
		auto cimage = createImageFromIntenzities(&((*contrastImageStat.contrastPicture)[0]), pSrc->GetWidth(), pSrc->GetHeight());
		cimage->SaveAs("../../contrast.png", TMO_PNG_8);
	*/

   // Compute SSIM (Structural Similarity Index) for each color channel (Red, Green, Blue)
	std::array<double, 3> SSIM_RC_GC_BC = computeSSIM(imageStatistics, contrastImageStat);

	std::array<double, 3> result;

   // Calculate the sum of SSIM values for normalization
	double invSSIMSum = 1.0/(SSIM_RC_GC_BC[0] + SSIM_RC_GC_BC[1] + SSIM_RC_GC_BC[2]);

   // Use the SSIM values to compute the constants for each channel
	result[0] = SSIM_RC_GC_BC[0] * invSSIMSum;
	result[1] = SSIM_RC_GC_BC[1] * invSSIMSum;
	result[2] = SSIM_RC_GC_BC[2] * invSSIMSum;

   /*
   // Alternative approach for debugging purposes (not used in the final result)
	result[0] = 1,
	result[1] = 1;
	result[2] = 1;*/

	return result;
}


/*
* Function to retrieve pixel value for the specified channel (R, G, or B) from the image data
*/
inline double TMOYu21::getPixel(const double* data, int width, int x, int y, int channel) {
   // Calculate the correct index for the pixel in the 1D array
    return data[(y * width + x) * 3 + channel];
}

/*
* Function to set the pixel value for a specific channel (R, G, or B) in the image data
*/
inline void TMOYu21::setPixel(double* data, int width, int x, int y, int channel, double value) {
   // Assign the pixel value to the corresponding position in the 1D array
    data[(y * width + x) * 3 + channel] = value;
}

/*
* Function to resize an image using bilinear interpolation
*/
std::unique_ptr<double[]> TMOYu21::resizeImage(const double* input, int srcWidth, int srcHeight, int destWidth, int destHeight)
{
    // Allocate memory for the resized image (output)
    std::unique_ptr<double[]> output(new double[destWidth * destHeight * 3]);

    // Iterate over each pixel of the destination image
    for (int y = 0; y < destHeight; ++y) {
        for (int x = 0; x < destWidth; ++x) {
            // Calculate the corresponding pixel coordinates in the source image
            double srcX = x * (static_cast<double>(srcWidth) / destWidth);
            double srcY = y * (static_cast<double>(srcHeight) / destHeight);

            // Identify the surrounding pixels in the source image (clamping to ensure valid indices)
            int x0 = std::clamp<double>(static_cast<int>(std::floor(srcX)), 0, srcWidth - 1);
            int x1 = std::clamp<double>(std::min(x0 + 1, srcWidth - 1), 0, srcWidth - 1);
            int y0 = std::clamp<double>(static_cast<int>(std::floor(srcY)), 0, srcHeight - 1);
            int y1 = std::clamp<double>(std::min(y0 + 1, srcHeight - 1), 0, srcHeight - 1);

            // Compute interpolation factors for the current pixel
            double dx = srcX - x0;
            double dy = srcY - y0;

            // Interpolate for each color channel (R, G, B)
            for (int c = 0; c < 3; ++c) {
               double value =
                    (1 - dx) * (1 - dy) * getPixel(input, srcWidth, x0, y0, c) +
                    dx * (1 - dy) * getPixel(input, srcWidth, x1, y0, c) +
                    (1 - dx) * dy * getPixel(input, srcWidth, x0, y1, c) +
                    dx * dy * getPixel(input, srcWidth, x1, y1, c);

               // Set the interpolated value to the corresponding pixel in the output image
               setPixel(output.get(), destWidth, x, y, c, value);
            }
        }
    }

    return output;
}

// Function to generate random pairs of coordinates within the given size
std::vector<std::pair<int, int>> TMOYu21::findRandomPairs(int size) const
{
	std::vector<std::pair<int, int>> result;
	result.reserve(size * size);

   // Generate random pairs of coordinates within the specified size
	for(int i = 0; i < size * size; ++i)
	{
		result.push_back({rand() % size, rand() % size});
	}

	return result;
}

// Compuing intenstity difference from random and neighboring pairs of pixels
std::shared_ptr<std::vector<double>> TMOYu21::computeContrastDifferences(const std::vector<std::pair<int, int>> &pairs, const double* image64, const double* image32, int channel) 
{
   // Create a shared pointer to store the contrast differences
   auto contrastDifferences = std::make_shared<std::vector<double>>();

   // Process random pairs from a 64x64 imag
	auto itPair = pairs.begin();

    for (int y = 0; y < 64; ++y) {
        for (int x = 0; x < 64; ++x) {
            // Get random pixel coordinates from the pair list
            int randomX = itPair->first;
            int randomY = itPair->second;
			   ++itPair;

             // Calculate the contrast difference between the current pixel and the randomly selected pixel
			   double pix1  = getPixel(image64, 64, x, y, channel);  // Current pixel intensity
			   double pix2 = getPixel(image64, 64, randomX, randomY, channel); // Random pixel intensity
            double diff = getPixel(image64, 64, x, y, channel) - getPixel(image64, 64, randomX, randomY, channel);
            contrastDifferences->push_back(diff);
        }
    }
    // Process neighboring pairs from a 32x32 image
    for (int y = 0; y < 32; ++y) {
        for (int x = 0; x < 32; ++x) {
			// For each pixel, check its neighbors (within a 3x3 grid)
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
					 // Skip the pixel itself (dx == 0 && dy == 0)
					if ((dx == 0 || dy == 0) && (dx != 0 || dy != 0) &&
            				x + dx >= 0 && x + dx < 32 && y + dy >= 0 && y + dy < 32) {
                        // Calculate contrast difference between the current pixel and its neighbor
                        double diff = getPixel(image32, 32, x, y, channel) - getPixel(image32, 32, x + dx, y + dy, channel);
                        contrastDifferences->push_back(diff);
                    }
                }
            }
        }
    }

    return contrastDifferences;
}


/*
* Helper function to log debug messages into a file, appending new messages after the first run
*/
void TMOYu21::logDebug(const std::string& message) {
	static bool firstRun = true;
    std::ofstream logFile("chyba.txt", firstRun ? (std::ios::out | std::ios::trunc) : std::ios::app);
    firstRun = false;
    if (logFile.is_open()) {
        logFile << message << std::endl;
    }
}

/*
* Computes the optimal weights for red, green, and blue channels (wr, wg, wb) that minimize the color energy
*/
std::array<double, 3> TMOYu21::computeWeights(const std::vector<double> &allIr, const std::vector<double> &allIg,
				const std::vector<double> &allIb, const std::array<double, 3> &kr_kg_kb)
{
	std::array<double, 3> result_wr_wg_wb;


	// Step size for iterating over weight combinations
	const double step = 0.1;
 
	const double kr = kr_kg_kb[0];
	const double kg = kr_kg_kb[1];
	const double kb = kr_kg_kb[2];

	double minEnergy = std::numeric_limits<double>::max();

   // Iterate over all possible combinations of wr, wg, wb
   for (double wr = 0.0; wr <= 1.0; wr += step) {
      for (double wg = 0.0; wg <= 1.0 - wr; wg += step) {
         double wb = 1.0 - wr - wg;

         // Check if wb is a valid value (multiple of step and within range)
         if (wb >= 0.0 && wb <= 1.0) 
			{
            // Calculate the energy for each color channel (R, G, B)
				double energyR = computeColorEnergy({wr, wg, wb}, kr, {allIr, allIg, allIb}, 1);
				double energyG = computeColorEnergy({wr, wg, wb}, kg, {allIr, allIg, allIb}, 1);
				double energyB = computeColorEnergy({wr, wg, wb}, kb, {allIr, allIg, allIb}, 2);
				
				// Calculate the total energy for this combination of weights
				double totalEnergy = energyR + energyG + energyB;

            // Alternative energy computation using a modified function
            // double totalEnergy = computeColorEnergy2({wr, wg, wb}, {kr, kg, kb}, {allIr, allIg, allIb});

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

/* 
* Computes the color energy of an image based on given weights and a scaling factor.
* The function follows the formula proposed in the referenced research paper, but it does not always produce expected results, as observed in experiments.
* 
* - `w` : Array of three weights (for R, G, B channels).
* - `k` : Computed constant
* - `I` : Array of three vectors representing pixel contrasts between color (R, G or B) and contrast image.
* - `colorIndex` : Index (0 = R, 1 = G, 2 = B) of the color channel used for contrast computation.
*/
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

      // Compute the weighted sum of the RGB intensities
		double weightedSum = wr * Ir + wg * Ig + wb * Ib;
		double absWeightedSum = std::abs(weightedSum);

      // Compute the numerator and denominator of the energy formula
		double numerator = absWeightedSum - k * std::abs(I[colorIndex][i]) - epsilon;
		double denominator = absWeightedSum + k * std::abs(I[colorIndex][i]) + epsilon;

      // Accumulate the energy value using the absolute ratio
		energy += std::abs(numerator) / denominator;
	}

	return energy;
}


/*
* 
* Modified version of the original formula.
* 
* This function incorporates a different approach to computing contrast differences.
* Instead of using a single scaling factor `k` for one color, it applies a separate `k` for each color channel computes the weighted sum of the contrast values.
* However, this approach also does not always work as expected, as observed in experiments.
* 
* - `w` : Array of three weights (for R, G, and B channels).
* - `k` : Array of three scaling factors (one for each color channel).
* - `I` : Array of three vectors representing pixel contrasts between the color image (RGB) and the contrast image.
*
*/
double TMOYu21::computeColorEnergy2(const std::array<double, 3> &w, const std::array<double, 3> &k, const std::array<std::vector<double>, 3> &I)
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

      // Compute the weighted sum of the RGB intensities (grayscale approximation)
		double weightedSumGray = wr * Ir + wg * Ig + wb * Ib;
		double absWeightedSumGray = std::abs(weightedSumGray);

      // Compute the weighted sum for contrast adjustments
		double weightedSumContrast = std::abs(k[0] * Ir + k[1] * Ig + k[2] * Ib); 

      // Compute the numerator and denominator of the energy formula
		double numerator = absWeightedSumGray - weightedSumContrast  - epsilon;
		double denominator = absWeightedSumGray + weightedSumContrast + epsilon;

		energy += std::abs(numerator) / denominator;
	}

	return energy;
}

/* --------------------------------------------------------------------------- *
 * Applies the tone mapping operator to transform the image. 	                *
 * --------------------------------------------------------------------------- */
int TMOYu21::Transform()
{
   // Get the source image data (R, G, B components for each pixel)
	double *pSourceData = pSrc->GetData();	

   // Get the destination image data (to store the result after transformation)	
	double *pDestinationData = pDst->GetData(); 

   // Compute image statistics like mean and standard deviation
	SImageStats imageStatistics = computeStats();
	
   // Compute coefficients (kr, kg, kb) based on image statistics
	std::array<double, 3> kr_kg_kb = computeK(imageStatistics);

   // Resize the source image to 32x32 and 64x64 for further processing
	std::unique_ptr<double[]> resized32 = resizeImage(pSourceData, pSrc->GetWidth(), pSrc->GetHeight(), 32, 32);
	std::unique_ptr<double[]> resized64 = resizeImage(pSourceData, pSrc->GetWidth(), pSrc->GetHeight(), 64, 64);

   // Find random pairs of pixels for contrast calculation
	auto randomPairs = findRandomPairs(64);

   // Compute the contrast differences for each color channel (R, G, B)
	auto allIr = computeContrastDifferences(randomPairs, resized64.get(), resized32.get(), 0);
	auto allIg = computeContrastDifferences(randomPairs, resized64.get(), resized32.get(), 1);
	auto allIb = computeContrastDifferences(randomPairs, resized64.get(), resized32.get(), 2);

   // Compute the weights for each color channel based on the contrast differences and the coefficients
	auto wr_wg_wb = computeWeights(*allIr, *allIg, *allIb, kr_kg_kb);
	
	int j = 0;
	int k = 0;

	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			
			auto R = *pSourceData++;
			auto G = *pSourceData++;
			auto B = *pSourceData++;

			auto intensity = wr_wg_wb[0] * R + wr_wg_wb[1] * G + wr_wg_wb[2] * B;

			*pDestinationData++ = intensity;
			*pDestinationData++ = intensity;
			*pDestinationData++ = intensity;
		}
	}

	// Normalize the resulting grayscale image (optional, not in the article)
   // normalizeGrayscaleImage(*pDst);

	return 0;
}

/*
* Normalize the grayscale image by scaling pixel values to the range [0, 1]
*/
void TMOYu21::normalizeGrayscaleImage(TMOImage &image)
{
    // Get the minimum and maximum pixel values in the image
	auto minmax = getImageMinMax(image);
	double min(minmax.first);
	double max(minmax.second);

   
   // Only normalize if the image has more than one unique pixel value
	if(max > min)
	{
      // Get the image data (pixel values)
		auto data = image.GetData();
		double invRange(1.0 / (max - min));

      // Normalize the pixel values to [0, 1] by scaling based on the min/max values
		for(size_t i = 0; i < 3 * image.GetWidth() * image.GetHeight(); ++i)
		{
			*data++ = (*data - min) * invRange;
		}
	}
}


/*
* Create a new image object from raw pixel data, specifying width and height
*/
std::unique_ptr<TMOImage> TMOYu21::createImage(const double *data, int width, int height)
{
   // Create a new image object
	auto pImage = std::make_unique<TMOImage>();
	pImage->New(width, height);

   // Allocate memory for the image data and copy the raw data into it
	auto dataCopy = new double[width * height * 3];
	memcpy(dataCopy, data, width * height * 3 * sizeof(double));

   // Set the image data
	pImage->SetData(dataCopy);

	return pImage;
}

/*
* Create a new image from intensity values, normalizing them to [0, 1]
*/
std::unique_ptr<TMOImage> TMOYu21::createImageFromIntenzities(const double *data, int width, int height)
{
   // Create a new image object
	auto pImage = std::make_unique<TMOImage>();
	pImage->New(width, height);
	auto dataCopy = new double[width * height * 3];

   // Find the minimum and maximum intensity values
	double min(std::numeric_limits<double>::max()), max(-std::numeric_limits<double>::max());
	for(size_t i = 0; i < width * height; ++i)
	{
		min = std::min(min, data[i]);
		max = std::max(max, data[i]);
	}

	double range(max - min);
	if(range > 0.0)
	{
		min = 0;
		range = 1.0;
		double denom(1.0 / range);
		
      // Normalize the intensity values to the range [0, 1] and assign them to all color channels
		auto dest(dataCopy);
		for(size_t i = 0; i < width * height; ++i)
		{
			*dest = (data[i] - min) * denom; ++dest;
			*dest = (data[i] - min) * denom; ++dest;
			*dest = (data[i] - min) * denom; ++dest;
		}
	}
   // Set the image data
	pImage->SetData(dataCopy);

	return pImage;
}

/*
* Get the minimum and maximum pixel values from the image data
*/
std::pair<double, double> TMOYu21::getImageMinMax(TMOImage &image)
{
	// Compute stats
	double min(std::numeric_limits<double>::max());
	double max(-max); // Negative infinity

   // Get the image data
	auto data = image.GetData();

   // Traverse through all pixels (R, G, B) and find the min and max values
	for(size_t i = 0; i < 3 * image.GetWidth() * image.GetHeight(); ++i)
	{
		double value = *data++;
		min = std::min(min, value);
		max = std::max(max, value);
	}

	return std::make_pair(min, max);
}

/*
* Remap the contrast image to the input image's pixel range
*/
void TMOYu21::remapContrastToInputRange(std::vector<double> &contrastImage)
{
	// Get the min and max range of the input image
	auto inputRange = getImageMinMax(*pSrc);

	// Get the min and max range of the contrast image
	auto contrastRange = getContrastImageMinMax(contrastImage);

   // Only perform remapping if both ranges are valid
	if(inputRange.first < inputRange.second && contrastRange.second > contrastRange.first)
	{
		double scale((inputRange.second - inputRange.first)/(contrastRange.second - contrastRange.first));

		for(double &v : contrastImage)
		{
			v = (v - contrastRange.first) * scale + inputRange.first;
		}
	}
}

/*
* Get the minimum and maximum values from the contrast image data
*/
std::pair<double, double> TMOYu21::getContrastImageMinMax(const std::vector<double> &image)
{
   // Initialize min and max values to extreme bounds
	double cmin(std::numeric_limits<double>::max());
	double cmax(-cmin);

   // Traverse through all contrast values and find the min and max
	for(double v : image)
	{
		cmin = std::min(v, cmin);
		cmax = std::max(v, cmax);
	}

	return std::make_pair(cmin, cmax);
}