/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                       Bachelor thesis                                             *
*                       Author: Martin Molek [xmolek00 AT stud.fit.vutbr.cz]        *
*                       Brno 2017                                                   *
*                                                                                   *
*                       Converting color images to grayscale by reducing dimensions *
*                       Tae-Hee Lee, Byoung-Kwang Kim, Woo-Jin Song                 *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOSong10.cpp
 * @brief Implementation of "Art-Photographic Detail Enhancement" Minjung Son, Yunjin Lee, Henry Kang, Seungyong Lee, Computer Graphics Forum 2014
 * @author Martin Molek
 * @class TMOSong10.cpp
 */

#include "TMOSong10.h"

/**
  *  @brief Constructor
  */
TMOSong10::TMOSong10()
{
	SetName(L"Song10");
	SetDescription(L"Reducing dimensions");

	k.SetName(L"k");
	k.SetDescription(L"Contrast");
	k.SetDefault(3.0);
	k = 3.0;
	k.SetRange(-10.0, 10.0);
	this->Register(k);

	alpha.SetName(L"alpha");
	alpha.SetDescription(L"Range");
	alpha.SetDefault(0.5);
	alpha = 0.5;
	alpha.SetRange(0.0, 1.0);
	this->Register(alpha);
}

/**
  *  @brief Destructor
  */
TMOSong10::~TMOSong10()
{
}

/**
 * @brief brief Converting color images to grayscale by reducing dimensions
 */
int TMOSong10::Transform()
{
	//pDst->Convert(TMO_RGB);
	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();

	double Y, Cb, Cr, Cd, r, g, b, k2;
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			r = *pSourceData++;
			g = *pSourceData++;
			b = *pSourceData++;

			Y = 0.299 * r + 0.587 * g + 0.114 * b;
			Cr = 0.500 * r - 0.419 * g - 0.081 * b;	 // + 128;
			Cb = -0.169 * r - 0.331 * g + 0.500 * b; // + 128;
			Cd = Cr - Cb;							 /** hence Cd =  0.669*r-0.088*g-0.581*b */

			if ((Cd) < 0.0)
			{
				k2 = -k;
				Cd = -Cd;
			}
			else
				k2 = k;

			Y = Y + k2 * Cd * pow(Cd, alpha);

			if (Y > 1.0)
				Y = 1.0;
			else if (Y < 0.0)
				Y = 0.0;

			*pDestinationData++ = Y;
			*pDestinationData++ = Y;
			*pDestinationData++ = Y;
		}
	}

	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);

	return 0;
}
