/* --------------------------------------------------------------------------- *
 * TMOMantiuk08.cpp: implementation of the TMOMantiuk08 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMantiuk08.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMantiuk08::TMOMantiuk08()
{
	SetName(L"Mantiuk08");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

#define vm 7
#define vm_contrast 2

struct DisplaySize{
   float view_dist;
   float pp_dist;
};

struct DisplayFunc{
   float gamma;
   float L_black;
   float L_max;
   float L_refl;
   float k;
   float E_amb;

};

// Display size functions

void displaySize(int res, float vd_screen, float vd_size, DisplaySize& disp)
{
   disp.view_dist = vd_size;
   disp.pp_dist = res * M_PI / (360 * atan(0.5/vd_screen));
}

// Display function functions
void DisplayFuncInit(float gamma, float L_black, float L_max, float k, float E_amb, DisplayFunc& df)
{
   df.gamma = gamma;
   df.L_black = L_black;
   df.L_max = L_max;
   df.k = k;
   df.E_amb = E_amb;
   df.L_refl = k/M_PI * E_amb;
}
float calcDisplayFunc(float p, DisplayFunc& df)
{
   if(p >= 0 && p <= 1)
   {
      return pow(p, df.gamma) * (df.L_max - df.L_black) + df.L_black + df.L_refl;
   }
}
float calcInverseDisplayFunc(float p, DisplayFunc& df)
{
   if(p < df.L_refl){ p = df.L_refl;}
   if(p > df.L_refl + df.L_max){ p = df.L_refl + df.L_max;}
   return pow((p - df.L_refl)/(df.L_max - df.L_black), 1/df.gamma);

}

// Human visual system functions
double transducer(double G, double sens)
{
   if(vm & vm_contrast)
   {
      double W = pow(10, fabs(G)) - 1;
      double k = 0.2599, q = 3, a = 3.291, b = 3.433, e = 0.8;
      double SW = W * sens;
      int sign = 0;
      if(G < 0)
      {
         sign = -1;
      }
      else{
         sign = 1;
      }
      return sign * (a*(pow(1+pow(SW,q),1.0/3.0) - 1))/(k * pow(b + SW, e));
   }
}

// Daly's contrast sensitivity function
double cs_daly(double rho, double view_dist = 0.5, double img_size, double theta, double adapt_lum)
{
   double P = 250.f;
   double eps = 0.9;
   double A = 0.801*(pow(1.0 + 0.7 * 1/adapt_lum, -0.2));
   double B = 0.3 * (pow(1.0 + 100 * 1/adapt_lum, 0.15));
   double r_a = 0.856 * powf(view_dist, 0.14);
   double c = 0.0;
   double r_c = 1.0/(1.0 + 0.24 * c);
   double r_theta = 0.11 * cosf(4.0 * theta) + 0.89;
   double b_eps_rho = B * eps * rho;
   double S1 = pow(pow(3.23 * pow(rho*rho*img_size,-0.3),5)+1.0, -0.2) * A * eps * rho * exp(-b_eps_rho) * sqrt(1 + 0.06*exp(b_eps_rho));
   double new_rho = rho / (r_a * r_c * r_theta);
   double b_eps_newrho = B * eps * new_rho;
   double S2 = pow(pow(3.23 * pow(new_rho*new_rho*img_size,-0.3),5)+1.0, -0.2) * A * eps * new_rho * exp(-b_eps_newrho) * sqrt(1 + 0.06*exp(b_eps_newrho));
   if(S1 > S2)
   {
      return S2 * P;
   }
   else{
      return S1 * P;
   }
}
TMOMantiuk08::~TMOMantiuk08()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMantiuk08::Transform()
{
   DisplaySize ds;
   DisplayFunc df;
   float res = 1024;
   float vd_screen = 2;
   float vd_size = 0.5;
   displaySize(res, vd_screen, vd_size, ds);
   DisplayFuncInit(2.2f, 0.8, 200, 0.01, 60, df);
   
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	pDst->Convert(TMO_Yxy); // x, y as color information

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
