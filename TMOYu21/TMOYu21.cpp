/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOYu21.h"

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
	
	return result;
}

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

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

			*itOut = 0.5 * (Krg * (pR + pG) + Kgb * (pG + pB) + Kbr * (pB + pR)); 
			++result.meanC;
			++itOut;
		}
	}

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
	double C1 = 0.01;

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
	
	return resultSSIM;
}

//Computing constants k for this picture
std::array<double, 3> TMOYu21::computeK(const SImageStats &imageStatistics)
{
	CImagePlusStats contrastImageStat = createContrastImage(imageStatistics);
	std::array<double, 3> SSIM_RC_GC_BC =computeSSIM(imageStatistics, contrastImageStat);
	return {};
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
	pDst->Convert(TMO_RGB);
	return 0;
}
