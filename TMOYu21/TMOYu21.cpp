/***********************************************************************************
*                                                                                  *
*                           Brno University of Technology                          *
*                         Faculty of Information Technology                        *
*                                                                                  *
*                           Color-to-Grayscale Conversions                         *
*                                                                                  *
*               Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]           *
*                                      Brno 2025                                   *
*                                                                                  *
*                        Implementation of the TMOYu21 class                       *
*   Contrast preserving decolorization based on the weighted normalized L1 norm    *
*                    https://doi.org/10.1007/s11042-021-11172-9                    *
*                                                                                  *
***********************************************************************************/


#include "TMOYu21.h"
#include <fstream>
#include <opencv2/opencv.hpp>

TMOYu21::TMOYu21()
{
	SetName(L"Yu21");					  
	SetDescription(L"Contrast preserving decolorization based on the weighted normalized L1 norm");

   HDRParameter.SetName(L"HDR");
	HDRParameter.SetDescription(L"is input image HDR");
	HDRParameter.SetDefault(false);
	HDRParameter = false;

   epsilonParameter.SetName(L"epsilon");
	epsilonParameter.SetDescription(L"epsilon - default 0.26, Cadik 0.15, Color250 0.2, CSDD 0.26");
	epsilonParameter.SetDefault(0.26);
	epsilonParameter = 0.26;

   this->Register(HDRParameter);
   this->Register(HDRParameter);
}

TMOYu21::~TMOYu21()
{
}

/*
*	This function calculates correlation coefficients (Krg, Kgb, Kbr) 
*  for the red, green, and blue channels of the input image.
*/
std::array<double, 3> TMOYu21::computeKrg_Kgb_Kbr()
{
	double *pSourceData(pSrc->GetData());

	std::array<double, 3> result;

   double meanR = 0.0, meanG = 0.0, meanB = 0.0;

	// Compute mean values for each channel
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
         double pR = *pSourceData++;
         double pG = *pSourceData++;
         double pB = *pSourceData++;
      
      
	      meanR += pR;
	      meanG += pG;
	      meanB += pB;
		}
	}

	double invNumValues(1.0 / double(pSrc->GetWidth() * pSrc->GetHeight()));
	meanR *= invNumValues;
	meanG *= invNumValues;
	meanB *= invNumValues;

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

			double diffR = pR - meanR;
			double diffG = pG - meanG;
			double diffB = pB - meanB;

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

   // Compute correlation coefficients - according to authors source code
   result[0] = numeratorRG / (denominatorR * denominatorG);
   result[1] = numeratorGB / (denominatorG * denominatorB);
   result[2] = numeratorBR / (denominatorB * denominatorR);

	return result;
}

/**
 *
 * This function computes a contrast-enhanced image using correlation coefficients Krg, Kgb, and Kbr.
 *
 * The contrast image is computed (according to the provided sourcecode, not article) as:
 * \f[
 * C(x, y) = \cdot (K_{rg} \cdot (R + G) + K_{gb} \cdot (G + B) + K_{br} \cdot (B + R))
 * \f]
 */
std::unique_ptr<std::vector<double>>  TMOYu21::createContrastImage(std::array<double, 3> &Krg_Kgb_Kbr)
{
	// Create data output
	std::unique_ptr<std::vector<double>> contrastPicture = std::make_unique<std::vector<double>>(pSrc->GetWidth() * pSrc->GetHeight());

	double *pSourceData(pSrc->GetData());

   // Load correlation coefficients
	double Krg(Krg_Kgb_Kbr[0]);
	double Kgb(Krg_Kgb_Kbr[1]);
	double Kbr(Krg_Kgb_Kbr[2]);


	auto itOut = contrastPicture->begin();

	double min = std::numeric_limits<double>::max();
   double max = std::numeric_limits<double>::lowest();


	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

      

         // Compute contrast image based on the provided sourcecode
         *itOut = Krg * (pR + pG) + Kgb * (pG + pB) + Kbr * (pB + pR); 
         if (min > *itOut)
            min = *itOut;

         if (max < *itOut)
            max = *itOut;

			++itOut;
		}
	}

   // Normalize values 
   for (auto it = contrastPicture->begin(); it != contrastPicture->end(); ++it)
   {
      *it = (*it - min) / (max - min);
   }

	return contrastPicture;
}

/*
* Computing SSIM(R,C)(G,C)(B,C) for this picture
*/
std::array<double, 3> TMOYu21::computeSSIM(std::unique_ptr<std::vector<double>> &contrastImage, cv::Mat imgR, cv::Mat imgG, cv::Mat imgB)
{
	// Output array for SSIM values
	std::array<double, 3> resultSSIM;

   // Create contrast cv::Mat img
   cv::Mat imgC(pSrc->GetHeight(), pSrc->GetWidth(), CV_64F, contrastImage->data());

	// Constants to stabilize division 
	const double K1 = 0.01, K2 = 0.03;
   const double L = 1.0; 
   const double C1 = (K1 * L) * (K1 * L);
   const double C2 = (K2 * L) * (K2 * L);
 
   // Convert images to to float 
   cv::Mat imgR_f, imgG_f, imgB_f, imgC_f;
   imgR.convertTo(imgR_f, CV_32F);
   imgG.convertTo(imgG_f, CV_32F);
   imgB.convertTo(imgB_f, CV_32F);
   imgC.convertTo(imgC_f, CV_32F);

   // Gaussian filter for means
   cv::Mat muR, muG, muB, muC;
   cv::GaussianBlur(imgR_f, muR, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgG_f, muG, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgB_f, muB, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgC_f, muC, cv::Size(11, 11), 1.5);

   // Squares of averages
   cv::Mat muR_sq = muR.mul(muR);
   cv::Mat muG_sq = muG.mul(muG);
   cv::Mat muB_sq = muB.mul(muB);
   cv::Mat muC_sq = muC.mul(muC);

   cv::Mat muR_muC = muR.mul(muC);
   cv::Mat muG_muC = muG.mul(muC);
   cv::Mat muB_muC = muB.mul(muC);

   // Gaussian filter for square pixels
   cv::Mat sigmaR_sq, sigmaG_sq, sigmaB_sq, sigmaC_sq, sigmaRC, sigmaGC, sigmaBC;
   cv::GaussianBlur(imgR_f.mul(imgR_f), sigmaR_sq, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgG_f.mul(imgG_f), sigmaG_sq, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgB_f.mul(imgB_f), sigmaB_sq, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgC_f.mul(imgC_f), sigmaC_sq, cv::Size(11, 11), 1.5);
   
   cv::GaussianBlur(imgR_f.mul(imgC_f), sigmaRC, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgG_f.mul(imgC_f), sigmaGC, cv::Size(11, 11), 1.5);
   cv::GaussianBlur(imgB_f.mul(imgC_f), sigmaBC, cv::Size(11, 11), 1.5);

   // Calculation of variances
   sigmaR_sq -= muR_sq;
   sigmaG_sq -= muG_sq;
   sigmaB_sq -= muB_sq;
   sigmaC_sq -= muC_sq;

   sigmaRC -= muR_muC;
   sigmaGC -= muG_muC;
   sigmaBC -= muB_muC;

   // Compute SSIM - RC
   cv::Mat numerator = (2 * muR_muC + C1).mul(2 * sigmaRC + C2);
   cv::Mat denominator = (muR_sq + muC_sq + C1).mul(sigmaR_sq + sigmaC_sq + C2);
   cv::Mat ssim_map;
   cv::divide(numerator, denominator, ssim_map);

   // Average SSIM for the channel
   resultSSIM[0] = cv::mean(ssim_map)[0];   

   // Compute SSIM - GC
   numerator = (2 * muG_muC + C1).mul(2 * sigmaGC + C2);
   denominator = (muG_sq + muC_sq + C1).mul(sigmaG_sq + sigmaC_sq + C2);
   cv::divide(numerator, denominator, ssim_map);

   // Average SSIM for the channel
   resultSSIM[1] = cv::mean(ssim_map)[0];   

   // Compute SSIM
   numerator = (2 * muB_muC + C1).mul(2 * sigmaBC + C2);
   denominator = (muB_sq + muC_sq + C1).mul(sigmaB_sq + sigmaC_sq + C2);
   cv::divide(numerator, denominator, ssim_map);

   // Average SSI for the channel
   resultSSIM[2] = cv::mean(ssim_map)[0];

   return resultSSIM;
}

/*
* Computing constants kr, kg, and kb for the given image based on SSIM (Structural Similarity Index)
*/
std::array<double, 3> TMOYu21::computeK(std::array<double, 3>  &Krg_Kgb_Kbr, cv::Mat imgR, cv::Mat imgG, cv::Mat imgB)
{
   // Create a contrast image based on the input image statistics
	std::unique_ptr<std::vector<double>> contrastImageStat = createContrastImage(Krg_Kgb_Kbr);	

   // Compute SSIM (Structural Similarity Index) for each color channel (Red, Green, Blue)
	std::array<double, 3> SSIM_RC_GC_BC = computeSSIM(contrastImageStat, imgR, imgG, imgB);

   auto SSIM_RC_ABS = std::abs(SSIM_RC_GC_BC[0]);
   auto SSIM_GC_ABS = std::abs(SSIM_RC_GC_BC[1]);
   auto SSIM_BC_ABS = std::abs(SSIM_RC_GC_BC[2]);

	std::array<double, 3> result;

   // Calculate the sum of SSIM values for normalization
	double invSSIMSum = 1.0/(SSIM_RC_ABS + SSIM_GC_ABS + SSIM_BC_ABS);

   // Use the SSIM values to compute the constants for each channel
	result[0] = SSIM_RC_ABS * invSSIMSum;
	result[1] = SSIM_GC_ABS * invSSIMSum;
	result[2] = SSIM_BC_ABS * invSSIMSum;


	return result;
}


/*
* Function to retrieve pixel value for the specified channel (R, G, or B) from the image data
*/
inline double TMOYu21::getPixel(const double* data, int width, int x, int y, int channel) {
   return data[(y * width + x) * 3 + channel];
}

/*
* Function to set the pixel value for a specific channel (R, G, or B) in the image data
*/
inline void TMOYu21::setPixel(double* data, int width, int x, int y, int channel, double value) {
   data[(y * width + x) * 3 + channel] = value;
}

/*
* Function to resize an image using bilinear interpolation
*/
std::unique_ptr<double[]> TMOYu21::resizeImage(double* input, int srcWidth, int srcHeight, int destWidth, int destHeight)
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
* Computes the optimal weights for red, green, and blue channels (wr, wg, wb) that minimize the color energy
*/
std::array<double, 3> TMOYu21::computeWeights(const std::vector<double> &allIr, const std::vector<double> &allIg,
				const std::vector<double> &allIb, const std::array<double, 3> &kr_kg_kb)
{
	std::array<double, 3> result_wr_wg_wb;


   // Filter out low-contrast values based on the provided source code
   std::vector<double> filteredIr;
   std::vector<double> filteredIg;
   std::vector<double> filteredIb;
   for (size_t i = 0; i < allIr.size(); ++i) {
      double det = std::sqrt(allIr[i] * allIr[i] + allIg[i] * allIg[i] + allIb[i] * allIb[i]) / 1.41;
      if (det >= 0.05) {
         filteredIr.push_back(allIr[i]);
         filteredIg.push_back(allIg[i]);
         filteredIb.push_back(allIb[i]);
      }
   }

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
				double energyR = computeColorEnergy({wr, wg, wb}, kr, {filteredIr, filteredIg, filteredIb}, 0);
				double energyG = computeColorEnergy({wr, wg, wb}, kg, {filteredIr, filteredIg, filteredIb}, 1);
				double energyB = computeColorEnergy({wr, wg, wb}, kb, {filteredIr, filteredIg, filteredIb}, 2);
				
				// Calculate the total energy for this combination of weights
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

/* 
* Computes the color energy of an image based on given weights and a scaling factor.
* The function's result is highly influenced by the epsilon variable,  
*  despite the authors claiming otherwise in their paper.
* 
* - `w` : Array of three weights (for R, G, B channels).
* - `k` : Computed constant
* - `I` : Array of three vectors representing pixel contrasts between color (R, G or B) and contrast image.
* - `colorIndex` : Index (0 = R, 1 = G, 2 = B) of the color channel used for contrast computation.
* - `epsilon` : key constant for conversion. Differes for datasets - Cadik 0.15, Color250 0.2, CSDD 0.26  
*/
double TMOYu21::computeColorEnergy(const std::array<double, 3> &w, double k, const std::array<std::vector<double>, 3> &I, size_t colorIndex)
{
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
		double numerator = absWeightedSum - k * std::abs(I[colorIndex][i]) - epsilonParameter;
		double denominator = absWeightedSum + k * std::abs(I[colorIndex][i]) + epsilonParameter;

      // Accumulate the energy value using the absolute ratio
		energy += std::abs(numerator) / denominator;
	}

	return energy;
}

// Finds if range is 0-1 or in 0-255
bool TMOYu21::isInRange0to1(double *pSourceData, int numPix)
{
   for (int i = 0; i < numPix * 3; i++)
   {
      if(pSourceData[i] > 1)
         return false;
   }
   return true;
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

   // If picture range is in range 255 instead of 0-1
   bool range0to1 = isInRange0to1(pSourceData, pSrc->GetHeight() * pSrc->GetWidth());

   if(!range0to1 && !HDRParameter)
   {
      for(int i = 0; i < pSrc->GetHeight() * pSrc->GetWidth() * 3; i++)
      {
         pSourceData[i] /= 255;
      }   
   }

   // Compute image statistics like mean and standard deviation
	std::array<double, 3> Krg_Kgb_Kbr = computeKrg_Kgb_Kbr();

   // Create cv::Mat pictures to compute ssim
   double* red = new double[pSrc->GetHeight() * pSrc->GetWidth()];   
   double* green = new double[pSrc->GetHeight() * pSrc->GetWidth()];
   double* blue = new double[pSrc->GetHeight() * pSrc->GetWidth()];
   double* pCurrentSource = pSourceData;

   for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{	
         red[j * pSrc->GetWidth() + i] = *pCurrentSource++;
         green[j * pSrc->GetWidth() + i] = *pCurrentSource++;
         blue[j * pSrc->GetWidth() + i] = *pCurrentSource++;
		}
	}

   cv::Mat imgR(pSrc->GetHeight(), pSrc->GetWidth(), CV_64F, red);
   cv::Mat imgG(pSrc->GetHeight(), pSrc->GetWidth(), CV_64F, green);
   cv::Mat imgB(pSrc->GetHeight(), pSrc->GetWidth(), CV_64F, blue);

	
   // Compute coefficients (kr, kg, kb) based on image statistics
	std::array<double, 3> kr_kg_kb = computeK(Krg_Kgb_Kbr, imgR, imgG, imgB);

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

   delete[] red;
   delete[] green;
   delete[] blue;

	return 0;
}