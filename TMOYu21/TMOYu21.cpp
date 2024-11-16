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

std::unique_ptr<std::vector<double>> TMOYu21::createContrastImage(const SImageStats &correlationCoefficients)
{
	// Create data output
	auto result = std::make_unique<std::vector<double>>(pSrc->GetWidth() * pSrc->GetHeight());

	double *pSourceData(pSrc->GetData());

	double Krg(correlationCoefficients.Krg);
	double Kgb(correlationCoefficients.Kgb);
	double Kbr(correlationCoefficients.Kbr);

	// Fill
	auto itOut = result->begin();

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pR = *pSourceData++;
			double pG = *pSourceData++;
			double pB = *pSourceData++;

			*itOut = 0.5 * (Krg * (pR + pG) + Kgb * (pG + pB) + Kbr * (pB + pR)); 

			++itOut;
		}
	}

	return result;
}

//Computing constants k for this picture
std::array<double, 3> TMOYu21::computeK()
{
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
