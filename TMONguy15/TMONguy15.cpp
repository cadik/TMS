/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                       Author: Matej Valek                                         *
*                       Brno 2020                                                   *
*                                                                                   *
*                       Implementation of the TMONguy15 class                       *
*                                                                                   *
************************************************************************************/
/**
 * @file TMONguy15.cpp
 * @brief Implementation of the TMONguy15 class
 * @author Matej Valek
 * @class TMONguy15.cpp
 * 
 * @todo Insert operator name
 */ 

#include "TMONguy15.h"
#include "mainprepare.h"

/**
  *  @brief Constructor
  */
TMONguy15::TMONguy15()
{
	SetName(L"Nguy15");						// TODO - Insert operator name
	SetDescription(L"Color to grayscale with nonlinear quadratic programing");	// TODO - Insert description

	//this->Register(dParameter);
}

/**
  *  @brief Destructor
  */
TMONguy15::~TMONguy15()
{
}

/**
  *  @brief Converts image
  * 
  *  Source image is stored in local parameter pSrc
  *  Destination image is in pDst
  *  Initialy images are in RGB format, but you can convert it into other format
  */
 int TMONguy15::Transform()
{
    double* pSourceData = pSrc->GetData();				/** You can work at low level data */

    double* pDestinationData = pDst->GetData();			/** Data are stored in form of array */
    double* pom = mainprepare(pSourceData,pSrc->GetWidth(),pSrc->GetHeight());		/** entering main calculations */
    
	double pY, px, py; /** three colour components */

	int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pom++;
			px = *pom++;
			py = *pom++;

			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}

	return 0;
}
