/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
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

	/* Gamma transfer function */
	gammaForm.SetName(L"gammaForm");
	gammaForm.SetDescription(L"Use a gamma transfer function");
	gammaForm.SetDefault(false);
	gammaForm=false;
	this->Register(gammaForm);

	/* Gamma */
	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma correction value: <1.0e-3,1.0e+2>");
	gamma.SetDefault(1.0);
	gamma=1.0;
	gamma.SetRange(1.0e-3,1.0e+2);
	this->Register(gamma);

	/* Exposure */
	exposure.SetName(L"exposure");
	exposure.SetDescription(L"Exposure scale factor: <0,100>");
	exposure.SetDefault(0.0);
	exposure=0.0;	
	exposure.SetRange(0,100);
	this->Register(exposure);

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

float BiasFunc (float t, float bias)
{
	const float LOG05 = -0.693147f;

	return pow(t, log(bias)/LOG05);
}

void SetExp (double* exp_d)
{
	*exp_d = pow(2, *exp_d);
}

int TMODrago03::Transform()
{
	double Y, x, y;	
	double L_av;
	double exp_d;
	double L_w, L_d, interpol;
	double L_min=0.;
	double L_max=0.;
	double L_world=0.;

	double* pSourceData;
	double* pDestinationData;

	pSrc->GetMinMaxAvgWorldAdapt(&L_min, &L_max, &L_world);	

	/* Set exposure */
	exp_d = exposure.GetDouble();
	SetExp(&exp_d);
	
	pSrc->Convert(TMO_Yxy);
	pDst->Convert(TMO_Yxy);

	pSourceData = pSrc->GetData();
	pDestinationData = pDst->GetData();

	if (center.GetBool())
	{
		pSrc->CenterWeight(centerX.GetInt(), centerY.GetInt(), (float)kernel.GetDouble(), &L_world);
	}
	
	L_av = exp(L_world)/1.0;
	L_max /= L_av;

	/* Tone mapping */
	int j = 0;

	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			Y = *pSourceData++;
			x = *pSourceData++;
			y = *pSourceData++;
						
			L_w = Y / L_av;			
			
			if (exp_d != 1.0)
			{
				L_w *= exp_d;
			}			
			
			interpol = log (2.0f + BiasFunc(L_w / L_max, bias.GetDouble()) * 8.0f);			
			L_d = (log(L_w+1.0f)/interpol) / log10(L_max+1.0f);

			*pDestinationData++ = L_d;
			*pDestinationData++ = x;
			*pDestinationData++ = y;
		}
	}

	pDst->Convert(TMO_RGB);

	/* Gamma */
	if (gamma.GetDouble() != 1.0)
	{
		if (gammaForm.GetBool())
		{		
			pDst->RecCorrectGamma(gamma.GetDouble());
		}
		else
		{
			pDst->CorrectGamma(gamma.GetDouble());
		}
	}
		
	pSrc->ProgressBar(j, pSrc->GetHeight());

	return 0;
} /* Transform */
