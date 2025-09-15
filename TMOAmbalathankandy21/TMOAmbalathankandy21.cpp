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
*                 Implementation of the TMOAmbalathankandy21 class             *
*                                                                              *
*               A color temperature-based high-speed decolorization:           *
*               an empirical approach for tone mapping applications            * 	
*                    https://doi.org/10.48550/arXiv.2108.13656                 *
*                                                                              *
*******************************************************************************/

#include "TMOAmbalathankandy21.h"
#include <fstream>


TMOAmbalathankandy21::TMOAmbalathankandy21()
{
	SetName(L"Ambalathankandy21");					 
	SetDescription(L"A color temperature-based high-speed decolorization"); 

   HDRParameter.SetName(L"HDR");
	HDRParameter.SetDescription(L"is input image HDR");
	HDRParameter.SetDefault(false);
	HDRParameter = false;

   this->Register(HDRParameter);
}

TMOAmbalathankandy21::~TMOAmbalathankandy21()
{
}

// Computes the minimum and maximum values in the image for normalization
std::pair<double, double> TMOAmbalathankandy21::getImageMinMax(TMOImage &image)
{
	double min(std::numeric_limits<double>::max());
	double max(-max);

	auto data = image.GetData();
	for(size_t i = 0; i < 3 * image.GetWidth() * image.GetHeight(); ++i)
	{
		double value = *data++;
		min = std::min(min, value);
		max = std::max(max, value);
	}

	return std::make_pair(min, max);
}

// Normalizes the grayscale image (not part of the original article, added for comparison)
void TMOAmbalathankandy21::normalizeGrayscaleImage(TMOImage &image)
{
	auto minmax = getImageMinMax(image);

	double min(minmax.first);
	double max(minmax.second);

	if(max > min)
	{
		auto data = image.GetData();
		double invRange(1.0 / (max - min));

		for(size_t i = 0; i < 3 * image.GetWidth() * image.GetHeight(); ++i)
		{
			*data++ = (*data - min) * invRange;
		}
	}
}

// Finds if range is 0-1 or in 0-255
bool TMOAmbalathankandy21::isInRange0to1(double *pSourceData, int numPix)
{
   for (int i = 0; i < numPix * 3; i++)
   {
      if(pSourceData[i] > 1)
         return false;
   }
   return true;
}


/* --------------------------------------------------------------------------- *
 * Applies the tone mapping operator to transform the image. 				       *
 * --------------------------------------------------------------------------- */
int TMOAmbalathankandy21::Transform()
{
	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData(); 

	double R, G, B, L_White, L_RG, L_Warm, L_Cool, L_B, L; 
	
	// Constants based on the reference article
	double betaR = 0.55;
	double betaK = 0.7;

   bool range0to1 = isInRange0to1(pSourceData, pSrc->GetHeight() * pSrc->GetWidth());

	// Iterate over each pixel in the image
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // Progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			// Read RGB values from the input image
			R = *pSourceData++;
			G = *pSourceData++;
			B = *pSourceData++;

         // If format is in range 0-255
         if (!range0to1 && !HDRParameter)
         {
            R /= 255;
            G /= 255;
            B /= 255;
         }

			// Compute luminance components
        	L_White = sqrt((R*R + G*G + B*B)/3);
        	L_B = B;
        	L_RG = sqrt(betaR * R * R + (1 - betaR) * G * G);
        	L_Warm = R / (R + G + B)*L_RG + (1 - R / (R + G + B)) * L_White;
        	L_Cool = (1 - B / (R + G + B))*L_B + B/(R + G + B) * L_White;
        	L = sqrt(betaK *L_Warm * L_Warm + (1 - betaK) * L_Cool * L_Cool);

         //L = 0.299 * R + 0.587 * G + 0.114 * B;

			// Store the computed luminance as grayscale output
			*pDestinationData++ = L;
			*pDestinationData++ = L;			
         *pDestinationData++ = L;		
      }
	}

	// Optional: Normalize the output image (not part of the original article)
	//normalizeGrayscaleImage(*pDst);
	return 0; 
}
