//(c)Martin Cadik
//see [tumblin99two], pp. 73
//11/2004

/* --------------------------------------------------------------------------- *
 * TMOWard94.cpp: implementation of the TMOWard94 class.                       *
 * --------------------------------------------------------------------------- */

#include "./TMOWard94.h"

#ifdef LINUX
 #define max(X,Y) (X<Y?Y:X)
#endif

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOWard94::TMOWard94()
{
	SetName(L"Ward94");
	SetDescription(L"A contrast-based scalefactor for luminance display");

	Ld_max.SetName(L"Ld_max");
	Ld_max.SetDescription(L"Maximum display luminance [cd/m^2]");
	Ld_max.SetDefault(100.);
	Ld_max=100.;
	Ld_max.SetRange(1.0e-3,1.0e+3);
	this->Register(Ld_max);

	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma Correction Factor [-]");
	gamma.SetDefault(2.2);
	gamma=2.2;
	gamma.SetRange(1.0e-3,1.0e+2);
	this->Register(gamma);
}

TMOWard94::~TMOWard94()
{
}

/* --------------------------------------------------------------------------- *
 * An implementation of tone mapping operator                                  *
 * --------------------------------------------------------------------------- */
int TMOWard94::Transform()
{
	int tmp_y=0, i=0, j=0;

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double stonits=pSrc->GetStonits();
	fprintf(stderr, "\nINPUT SCENE: \nStonits: %g\n", stonits);

	double m = 0.;                   // Ward's scale factor
	
	double L_max_log=0., L_max=0.,   // max luminance 
		   L_min_log=0., L_min=0.,				 // min luminance
		   L_wa_log = 0., L_wa=0.;   // adaptation luminance

	pSrc->GetMinMaxAvg(&L_min, &L_max, &L_wa);
	fprintf(stderr, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * stonits, L_max * stonits, L_wa * stonits);
    pSrc->GetMinMaxAvgLog10(&L_min_log, &L_max_log, &L_wa_log);
	fprintf(stderr, "Min luminance (log): %g[cd/m^2]\nMax luminance (log): %g[cd/m^2]\nAvg luminance (log): %g[cd/m^2]\n", L_min_log + TAKE_LOG10(stonits), L_max_log + TAKE_LOG10(stonits), (L_wa_log + TAKE_LOG10(stonits)));

	pSrc->Convert(TMO_Yxy);						
	pDst->Convert(TMO_Yxy);					

	//fprintf(stderr, "\nINPUT PARAMETERS: \n%s: %g\n", Ld_max.GetDescription().c_str(), Ld_max.GetDouble());
	fprintf(stderr, "\nINPUT PARAMETERS: \nLd_max: %g[cd/m^2]\n", (double)Ld_max);

	L_wa=pow((double)10, (L_wa_log + TAKE_LOG10(stonits)));  // log(L_wa) = mean (log(L_w)) - adaptation luminance
	fprintf(stderr, "\nOPERATOR RESULTS:\nAdaptation luminance: %g\n", L_wa);

	// Ward's formula
	m = 1/Ld_max * pow( (1.219+pow(Ld_max/2., 0.4))/(1.219+pow(L_wa, 0.4)), 2.5 );
	fprintf(stderr, "Scale factor (m): %g\n", m);

	double Ld_out_max=-HUGE_VAL;
//	double min=100;
//	struct timespec reg;
//	reg.tv_sec=0;
//	reg.tv_nsec=1;
	for (i = 0; i < pSrc->GetHeight(); i++)
	{
//	nanosleep(&reg,NULL);
		pSrc->ProgressBar(i, pSrc->GetHeight());
		for (int j = 0; j < pSrc->GetWidth(); j++)
		{
			double L_w = *pSourceData++;
			double x_w = *pSourceData++;
			double y_w = *pSourceData++;

//			if(min>(L_w*stonits) && (L_w*stonits>0)) min=L_w*stonits;

			Ld_out_max=max(Ld_out_max, m * L_w);
			*pDestinationData++ = m * L_w ;  //L_d = m * L_w
			*pDestinationData++ = x_w; //colors remain unchanged
			*pDestinationData++ = y_w; //colors remain unchanged
//		for(unsigned i=0;i<100000;i++);
		}
	}	

//	printf("mininum je: %g\n", min);

	////normalization using Ld_out_max
	//printf("\nMax output luminance: %e\n", Ld_out_max);
	//pDestinationData = pDst->GetData();
	//for (i = 0; i < pSrc->GetHeight(); i++)
	//{
	//	pSrc->ProgressBar(i, pSrc->GetHeight());
	//	for (int j = 0; j < pSrc->GetWidth(); j++)
	//	{
	//		*pDestinationData++ /= Ld_out_max;
	//		*pDestinationData++;
	//		*pDestinationData++;
	//	}
	//}	

	pDst->Convert(TMO_RGB);
	pDst->CorrectGamma(gamma);
	pSrc->ProgressBar(i, pSrc->GetHeight());

	return 0;
}//Transform
