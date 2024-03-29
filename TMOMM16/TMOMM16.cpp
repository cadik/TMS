/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio	                                    *
*                                                                                   *
*                       Bachelor thesis                                             *
*                       Author: Martin Molek [xmolek00 AT stud.fit.vutbr.cz]        *
*                       Brno 2017                                                   *
*                                                                                   *
*                       Implementation of TMOMM16 class by Martin Molek             *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOMM16.cpp
 * @brief Implementation of TMOMM16 class by Martin Molek
 * @author Martin Molek
 * @class TMOMM16.cpp
 */

#include "TMOMM16.h"

/**
  *  @brief Constructor
  */
TMOMM16::TMOMM16()
{
	SetName(L"MM16");
	SetDescription(L"Co jsem si napsal");

	red.SetName(L"red");
	red.SetDescription(L"RED");
	red.SetDefault(0.299);
	red = 0.299;
	red.SetRange(0.0, 1.0);
	this->Register(red);

	green.SetName(L"green");
	green.SetDescription(L"GREEN");
	green.SetDefault(0.587);
	green = 0.587;
	green.SetRange(0.0, 1.0);
	this->Register(green);

	blue.SetName(L"blue");
	blue.SetDescription(L"BLUE");
	blue.SetDefault(0.114);
	blue = 0.114;
	blue.SetRange(0.0, 1.0);
	this->Register(blue);

	negative.SetName(L"neg");
	negative.SetDescription(L"Negative");
	negative.SetDefault(false);
	negative = false;
	this->Register(negative);
}

/**
  *  @brief Destructor
  */
TMOMM16::~TMOMM16()
{
}

/**
 * @brief transformation function
 */
int TMOMM16::Transform()
{
	//pDst->Convert(TMO_RGB);
	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();

	double r, g, b, shade;
	//int inputShade;
	//int value;
	//int outputShade;

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			r = *pSourceData++;
			g = *pSourceData++;
			b = *pSourceData++;

			shade = r * red + g * green + b * blue;

			if (shade > 1.0)
			{
				shade = 1.0;
			}
			else if (shade < 0.0)
			{
				shade = 0.0;
			}

			if (negative)
			{
				*pDestinationData++ = 1.0 - shade;
				*pDestinationData++ = 1.0 - shade;
				*pDestinationData++ = 1.0 - shade;
			}
			else
			{
				*pDestinationData++ = shade;
				*pDestinationData++ = shade;
				*pDestinationData++ = shade;
			}
		}
	}

	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);

	return 0;
}
