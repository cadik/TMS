/* --------------------------------------------------------------------------- *
 * TMOHafner16.cpp: implementation of the TMOHafner16 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOHafner16.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOHafner16::TMOHafner16()
{
	SetName(L"Hafner16");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOHafner16::~TMOHafner16()
{
}

void RGB2YCbCr(double* data, int numPix)
{
   double r, g, b;
   for(int i = 0; i < numPix*3; i += 3)
   {
      r = data[i];
      g = data[i+1];
      b = data[i+2];

      data[i] = 0.299 * r + 0.587 * g + 0.114 * b; // Y
      data[i+1] = 128 + (-0.168736 * r - 0.33126 * g + 0.5 * b); // Cb
      data[i+2] = 128 + (0.5 * r - 0.418688 * g - 0.081312 * b); // Cr
   }
}
/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOHafner16::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format


	double *pSourceData = pSrc->GetData();	
	double *pDestinationData = pDst->GetData(); 
   pDst->Convert(TMO_Yxy);

   int numPix = pSrc->GetWidth() * pSrc->GetHeight();
   RGB2YCbCr(pSourceData, numPix);
	double pY, pCb, pCr;

   int numIter = 100; // Number of optimization iterations
   double alpha = 1.0, beta = 1.0, gamma = 0.25; // Energy weights

   // Initialize weights evenly
   std::vector<double> wR(numPix, 1.0 / 3.0);
   std::vector<double> wG(numPix, 1.0 / 3.0);
   std::vector<double> wB(numPix, 1.0 / 3.0);

   for (int iter = 0; iter < numIter; ++iter) {
      for (int i = 0; i < numPix; ++i) {
          double Y = pSourceData[i * 3];
          double Cb = pSourceData[i * 3 + 1];
          double Cr = pSourceData[i * 3 + 2];
          
          // Compute contrast-enhanced term (simplified)
          double contrast = gamma * (wR[i] * Y + wG[i] * Cb + wB[i] * Cr);
          
          // Compute new weights (gradient step)
          wR[i] += alpha * (Y - contrast);
          wG[i] += alpha * (Cb - contrast);
          wB[i] += alpha * (Cr - contrast);
      }
      
      // Project weights onto simplex to ensure sum = 1
      for (int i = 0; i < numPix; ++i) {
          double sumW = wR[i] + wG[i] + wB[i];
          wR[i] /= sumW;
          wG[i] /= sumW;
          wB[i] /= sumW;
      }
  }


	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
         double Y = wR[i] * pSourceData[i * 3] + 
            wG[i] * pSourceData[i * 3 + 1] + 
            wB[i] * pSourceData[i * 3 + 2];

			// and store results to the destination image
			*pDestinationData++ = Y;
			*pDestinationData++ = 0;
			*pDestinationData++ = 0;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}
