/* --------------------------------------------------------------------------- *
 * TMOYee03.cpp: implementation of the TMOYee03 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOYee03.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOYee03::TMOYee03()
{
	SetName(L"Yee03");					  // TODO - Insert operator name
	SetDescription(L"Segmentation and adaptive assimilation for detail-preserving"); // TODO - Insert description

	bin_size1.SetName(L"bin_size1");				// TODO - Insert parameters names
	bin_size1.SetDescription(L"bin size determines contrast captured"); // TODO - Insert parameter descriptions
	bin_size1.SetDefault(0.5);							// TODO - Add default values
	bin_size1 = 0.5;
	bin_size1.SetRange(0.5, 1.0); // TODO - Add acceptable range if needed
	this->Register(bin_size1);

   bin_size2.SetName(L"bin_size2");				
	bin_size2.SetDescription(L"bin size determines contrast captured"); 
	bin_size2.SetDefault(1.);							
	bin_size2 = 1.;
	bin_size2.SetRange(1., 2.0); 
	this->Register(bin_size2);

   small_threshold.SetName(L"small_threshold");				
	small_threshold.SetDescription(L"threshold determines when larger groups assimilate small groups"); 
	small_threshold.SetDefault(0.001);						
	small_threshold = 0.001;
	small_threshold.SetRange(0.001, 3.0); 
	this->Register(small_threshold);

   big_threshold.SetName(L"big_threshold");				
	big_threshold.SetDescription(L"threshold determines when larger groups assimilate small groups"); 
	big_threshold.SetDefault(1.);							
	big_threshold = 1.;
	big_threshold.SetRange(1., 10.0); 
	this->Register(big_threshold);   

   max_layers.SetName(L"max_layers");				
	max_layers.SetDescription(L"max_layers controls how many layers to generate, its used to smooth out the boundaries"); 
	max_layers.SetDefault(16.);							
	max_layers = 16.;
	max_layers.SetRange(16., 96.0); 
	this->Register(max_layers);

}

#define cdm2ToLambert(C) (C*0.001/3.18)
#define lambertToCmd2(L) (L*3.18*1000)
#define rgb2luminance(R,G,B) (R*0.21 + G*0.71 + B*0.07)
#define MAX_DISPLAY_LUMINANCE 125.0
#define DISPLAY_ADAPTATION_LUMINANCE 25.0

unsigned int IMAGE_HEIGHT;
unsigned int IMAGE_WIDTH;

struct Point{
   unsigned short x,y;
};

typedef std::vector<Point> PointVector;

typedef std::list<unsigned short> NeighbourContainer;
struct Group_Record{
   PointVector Memebers;
   NeighbourContainer Neighbours;
   double Sum;
   unsigned int Count;
};

typedef std::vector<Group_Record *> Groups;

typedef std::vector<Groups> Layers;

TMOYee03::~TMOYee03()
{
}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOYee03::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pR, pG, pB;
   IMAGE_HEIGHT = pSrc->GetHeight();
   IMAGE_WIDTH = pSrc->GetWidth();

   double stonits = pSrc->GetStonits();
   double MinimalImageLuminance = 0., MaximalImageLuminance = 0., AverageImageLuminance = 0.;
   pSrc->GetMinMaxAvg(&MinimalImageLuminance, &MaximalImageLuminance, &AverageImageLuminance);
   MinimalImageLuminance = MinimalImageLuminance*stonits;
   MaximalImageLuminance = MaximalImageLuminance*stonits;
   // Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information
   double AdaptationLuminance = 800.0;
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			//pR = *(pSourceData+((j*IMAGE_WIDTH*3)+(i*3)));
			//pG = *(pSourceData+((j*IMAGE_WIDTH*3)+(i*3)+1));
			//pB = *(pSourceData+((j*IMAGE_WIDTH*3)+(i*3)+2));
         pR = *pSourceData++;
			pG = *pSourceData++;
			pB = *pSourceData++;

         fprintf(stderr, "pixel %i,%i value R: %g G: %g B: %g",j,i,pR,pG,pB);

			// Here you can use your transform
			// expressions and techniques...
			//pR *= dParameter; // Parameters can be used like
							  // simple variables
         double L_wa = cdm2ToLambert(AdaptationLuminance);
         double L_w = cdm2ToLambert(rgb2luminance(pR, pG, pB));
         double f_r = pR/L_w, f_g = pG/L_w, f_b = pB/L_w;

         double S_w = 100.0 + 10.0*log10(L_wa);
         double R_w = 10.0*log10(L_wa/L_w);

         double L_dMax = cdm2ToLambert(MAX_DISPLAY_LUMINANCE);
         double L_da = cdm2ToLambert(DISPLAY_ADAPTATION_LUMINANCE);
         double S_d = 100.0 + 10.0*log10(L_da);
         double R_d = 8.4 - (S_w-27)*(8.4-R_w)/(S_d-27);

         double L_d = lambertToCmd2(L_da*pow(10,-0.1*R_d))/MAX_DISPLAY_LUMINANCE;
         
         fprintf(stderr, " cR: %g cG: %g cB: %g\n",L_d*f_r, L_d*f_g, L_d*f_b);
         // and store results to the destination image
			//*pDestinationData++ = pR;
			//*pDestinationData++ = pG;
			//*pDestinationData++ = pB;
         *pDestinationData++ = MIN(1.0,L_d*f_r);
			*pDestinationData++ = MIN(1.0,L_d*f_g);
			*pDestinationData++ = MIN(1.0,L_d*f_b);
		}
	}
   fprintf(stderr, "\nMinimal image luminance: %g\n", MinimalImageLuminance);
   fprintf(stderr, "Maximal image luminance: %g\n",MaximalImageLuminance);
	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0;
}
