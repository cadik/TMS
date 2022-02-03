/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2021                                              *
*                                                                              *
*                       Implementation of the TMODrago03 class                 *
*                                                                              *
*******************************************************************************/
#include <math.h>

#include "TMODrago03.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMODrago03::TMODrago03()
{
	SetName(L"Drago03");
	SetDescription(L"Adaptive Logarithmic Mapping For Displaying High Contrast Scenes");	
	
	/* Kernel size multiplier */
	kernel.SetName(L"kernel");
	kernel.SetDescription(L"Kernel size multiplier: <0.1,1.0>");
	kernel.SetDefault(0.125);
	kernel=0.125;
	this->Register(kernel);
	kernel.SetRange(0.1, 1.0);

	/* y coordinate for center-weighting */
	centerY.SetName(L"centerY");
	centerY.SetDescription(L"y coordinate for center-weighting");
	centerY.SetDefault(false);
	centerY=false;
	this->Register(centerY);	

	/* x coordinate for center-weighting */
	centerX.SetName(L"centerX");
	centerX.SetDescription(L"x coordinate for center-weighting");
	centerX.SetDefault(0);
	centerX=0;
	this->Register(centerX);	

	/* Center-weighted scalefactor */
	center.SetName(L"center");
	center.SetDescription(L"Use a center-weighted scalefactor");
	center.SetDefault(false);
	center=false;
	this->Register(center);		

	/* Exposure */
	exposure.SetName(L"exposure");
	exposure.SetDescription(L"Exposure scale factor: <0,100>");
	exposure.SetDefault(0.0);
	exposure=0.0;
	exposure.SetRange(0,100);
	this->Register(exposure);

	/* Gamma */
	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma correction value: <1.0e-3,1.0e+2>");
	gamma.SetDefault(1.125);
	gamma=1.125;
	gamma.SetRange(1.0e-3,1.0e+2);
	this->Register(gamma);

	/* Bias parameter b */
	bias.SetName(L"bias");
	bias.SetDescription(L"Bias parameter b: <0.7,0.9>");
	bias.SetDefault(0.85);
	bias=0.85;
	bias.SetRange(0.7,0.9);
	this->Register(bias);	
}

TMODrago03::~TMODrago03()
{
}

double biasFunc (double t, double bias)
{
	const double LOG05 = -0.693147;

	return pow(t, log(bias)/LOG05);
}

void setExp (double* exp_d)
{
	*exp_d = pow(2, *exp_d);
}

int TMODrago03::Transform()
{
	double X, Y, Z;	
	double L_w, L_d, L_s;
	double interpol, divider;
	double exp_d;
	double L_max, L_av;
	double biasValue;

	double* pSourceData;
	double* pDestinationData;

	/* Set exposure */
	exp_d = exposure.GetDouble();
	setExp(&exp_d);
	
	pSrc->Convert(TMO_XYZ);
	pDst->Convert(TMO_XYZ);

	pSourceData = pSrc->GetData();
	pDestinationData = pDst->GetData();

	/* Set L_max and L_av */
	pSrc->CalculateLuminance(L_max, L_av);
	
	if (center.GetBool())
	{
		pSrc->CenterWeight(centerX.GetInt(), centerY.GetInt(), kernel.GetDouble(), &L_av);
	}

	/* Tone mapping */	
	L_max /= L_av;
	divider = log10(L_max+1.0);
	
	biasValue = bias.GetDouble();

	int j = 0;

	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			X = *pSourceData++;
			Y = *pSourceData++;
			Z = *pSourceData++;
			
			L_w = Y / L_av;
			
			if (exp_d != 1.0)
			{
				L_w *= exp_d;
			}
			
			interpol = log(2.0 + biasFunc(L_w / L_max, biasValue) * 8.0);			
			L_d = (log(L_w+1.0)/interpol) / divider;

			L_s = L_d / Y;

			*pDestinationData++ = X * L_s;			
			*pDestinationData++ = L_d;			
			*pDestinationData++ = Z * L_s;	
		}
	}

	pDst->Convert(TMO_RGB);
	/*
	 * Due to different data processing in pfstools
	 * and TMS it is necessary to use different gamma
	 * correction value compared to the original article
	 */
	pDst->CorrectGamma(gamma.GetDouble());
	pSrc->ProgressBar(j, pSrc->GetHeight());

	return 0;
} /* Transform */
