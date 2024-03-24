/* --------------------------------------------------------------------------- *
 * TMOTMONafchi17.cpp: implementation of the TMOTMONafchi17 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOTMONafchi17.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOTMONafchi17::TMOTMONafchi17()
{
	SetName(L"TMONafchi17");
	SetDescription(L"Color to Gray Conversion by Correlation.");

	dParameter.SetName(L"Inverse correlation contribution");
	dParameter.SetDescription(L"Controls the level of contribution of the inverse correlations");
	dParameter.SetDefault(0.5);
	dParameter = 0.5;
	dParameter.SetRange(0.0, 1.0);
	this->Register(dParameter);

   bParameter.SetName(L"Complement standard deviation");
	bParameter.SetDescription(L"Standard deviation image will be subsituted for it`s complement");
	bParameter.SetDefault(false);
	bParameter = false;
	this->Register(bParameter);
}

TMOTMONafchi17::~TMOTMONafchi17()
{
}

double getCorr(std::vector<double> *X, double meanX, std::vector<double> *Y, double meanY)
{
   double sumOfDiffs = 0.0, sumOfDiffX = 0.0, sumOfDiffY = 0.0;
   double downPart;

   for (int i = 0; i < X->size(); i++) {
      sumOfDiffs += (X->at(i) - meanX) * (Y->at(i) - meanY);
      sumOfDiffX += std::pow(X->at(i) - meanX, 2);
      sumOfDiffY += std::pow(Y->at(i) - meanY, 2);
   }

   downPart = std::sqrt(sumOfDiffX * sumOfDiffY);

   return sumOfDiffs / downPart;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOTMONafchi17::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
   pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py, pR, pG, pB, pQ, u, d, diff, g;
   double meanR, meanG, meanB, meanQ;
   double sumR = 0.0, sumG = 0.0, sumB = 0.0, sumQ = 0.0;
   double pearsRQ, pearsGQ, pearsBQ;

   double betaR, betaG, betaB;
   double weightR, weightG, weightB;
   double lambdaR, lambdaG, lambdaB, lambdaSum;
   double pearsMin, pearsMax, pearsSum;

   int totalPixels = pSrc->GetHeight() * pSrc->GetWidth();
   std::vector<double> qVect(totalPixels);
   std::vector<double> rVect(totalPixels);
   std::vector<double> gVect(totalPixels);
   std::vector<double> bVect(totalPixels);

   const double dThrs = (147.2243f / 255.0f);

	int j = 0, index = 0;

   for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight() * 2);
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
         index = (j * pSrc->GetWidth()) + i;

			pR = *pSourceData++;
			pG = *pSourceData++;
			pB = *pSourceData++;
         sumR += pR;
         sumG += pG;
         sumB += pB;

         u = (pR + pG + pB) / 3.0f;
         d = std::sqrt((std::pow(pR -  u, 2) + std::pow(pG -  u, 2) + std::pow(pB -  u, 2)) / 2.0f) / dThrs;

         if (bParameter) {
            d = 1 - d;
         }

         pQ = u * d;

         sumQ += pQ;

         qVect[index] = pQ;
         rVect[index] = pR;
         gVect[index] = pG;
         bVect[index] = pB;
		}
	}

   meanR = sumR / totalPixels;
   meanG = sumG / totalPixels;
   meanB = sumB / totalPixels;
   meanQ = sumQ / totalPixels;

   pearsRQ = getCorr(&qVect, meanQ, &rVect, meanR);
   pearsGQ = getCorr(&qVect, meanQ, &gVect, meanG);
   pearsBQ = getCorr(&qVect, meanQ, &bVect, meanB);

   pearsMax = std::max(pearsRQ, std::max(pearsGQ, pearsBQ));
   pearsMin = std::min(pearsRQ, std::min(pearsGQ, pearsBQ));

   pearsSum = std::abs(pearsRQ) + std::abs(pearsGQ)  + std::abs(pearsBQ);
   betaR = std::abs(pearsRQ) / pearsSum;
   betaG = std::abs(pearsGQ) / pearsSum;
   betaB = std::abs(pearsBQ) / pearsSum;

   weightR = ((pearsRQ - pearsMin) / (pearsMax - pearsMin)) - dParameter;
   weightG = ((pearsGQ - pearsMin) / (pearsMax - pearsMin)) - dParameter;
   weightB = ((pearsBQ - pearsMin) / (pearsMax - pearsMin)) - dParameter;

   lambdaR = std::abs(betaR + std::min(betaR, weightR));
   lambdaG = std::abs(betaG + std::min(betaG, weightG));
   lambdaB = std::abs(betaB + std::min(betaB, weightB));
   lambdaSum = lambdaR + lambdaG + lambdaB;

   lambdaR /= lambdaSum;
   lambdaG /= lambdaSum;
   lambdaB /= lambdaSum;

   for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j + pSrc->GetHeight(), pSrc->GetHeight() * 2);
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
         index = (j * pSrc->GetWidth()) + i;
         g = ((lambdaR * rVect[index]) + (lambdaG * gVect[index]) + (lambdaB * bVect.at(index)));;


         std::cout << "G is: " << g << " on index: " << index << "\n";

			*pDestinationData++ = g;
         *pDestinationData++ = g;
         *pDestinationData++ = g;
		}
	}


   pSrc->ProgressBar(pSrc->GetHeight() * 2, pSrc->GetHeight() * 2);
	pDst->Convert(TMO_RGB);
	return 0;
}
