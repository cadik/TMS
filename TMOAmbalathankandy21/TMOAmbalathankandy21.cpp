/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAmbalathankandy21.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAmbalathankandy21::TMOAmbalathankandy21()
{
	SetName(L"Ambalathankandy21");					 
	SetDescription(L"Method ambalathankandy21"); 
}

TMOAmbalathankandy21::~TMOAmbalathankandy21()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAmbalathankandy21::Transform()
{
	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData(); 

	double R, G, B, L_White, L_RG, L_Warm, L_Cool, L_B, L; 
	double betaR = 0.55; // Contstants according to the article
	double betaK = 0.7;

	// For each pixel
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // Progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			// RGB pixel form input picture
			R = *pSourceData++;
			G = *pSourceData++;
			B = *pSourceData++;

			// Computing luminance
        	L_White = sqrt((R*R + G*G + B*B)/3);
        	L_B = B;
        	L_RG = sqrt(betaR * R * R + (1 - betaR) * G * G);
        	L_Warm = R / (R + G + B)*L_RG + (1 - R / (R + G + B)) * L_White;
        	L_Cool = (1 - B / (R + G + B))*L_B + B/(R + G + B) * L_White;
        	L = sqrt(betaK *L_Warm * L_Warm + (1 - betaK) * L_Cool * L_Cool);

			// Computed pixel to output
			*pDestinationData++ = L;
			*pDestinationData++ = L;
			*pDestinationData++ = L;
		}
	}
	return 0; 
}
