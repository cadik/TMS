/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAncuti16.h"
#include <boost/concept_check.hpp>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAncuti16::TMOAncuti16()
{
	SetName(L"Ancuti16");						// TODO - Insert operator name
	SetDescription(L"Add your TMO description here");	// TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOAncuti16::~TMOAncuti16()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAncuti16::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	//pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	//pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
	int h=	pSrc->GetHeight();
	int w =pSrc->GetWidth();
										// three colour components
		
	double *red = (double*)malloc( h * w * sizeof(double));
	double *green = (double*)malloc( h * w * sizeof(double));
	double *blue = (double*)malloc( h * w * sizeof(double));
	
	double *redLap = (double*)malloc( h * w * sizeof(double));
	double *greenLap = (double*)malloc( h * w * sizeof(double));
	double *blueLap = (double*)malloc( h * w * sizeof(double));
	
	double *lapWeightMapR =(double*)malloc( h * w * sizeof(double));
	double *lapWeightMapG =(double*)malloc( h * w * sizeof(double));
	double *lapWeightMapB =(double*)malloc( h * w * sizeof(double));
	double *globWeightMapR =(double*)malloc( h * w * sizeof(double));
	double *globWeightMapG =(double*)malloc( h * w * sizeof(double));
	double *globWeightMapB =(double*)malloc( h * w * sizeof(double));
	
	float kernel[3][3] = {{0,-1,0},
			      {-1,4,-1},
			      {0,-1,0}};
	double sumRed=0.0;
	double sumGreen=0.0;
	double sumBlue=0.0;

        int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
 		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			red[i+j*w] = *pSourceData++;
			green[i+j*w] = *pSourceData++;
			blue [i+j*w]= *pSourceData++;

			// Here you can use your transform 
			// expressions and techniques...
			//pY *= dParameter;							// Parameters can be used like
														// simple variables

			// and store results to the destination image
			/**pDestinationData++ = 0;
			*pDestinationData++ = 0;
			*pDestinationData++ = b;*/

			


		}
	}
	
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		
	    pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
	    for (int i = 0; i < pSrc->GetWidth(); i++)  ////vynechvam okraje => treba sa spytat ako to riesit
	    {
	      sumRed = 0.0;
	      sumGreen = 0.0;
	      sumBlue = 0.0;
	    
	    
	      
		for(int k = -1; k <=1; k++)
		{
		  for(int l = -1; l <=1; l++)
		  {

		      sumRed = sumRed + getSum(i,j,kernel,red,l,k);
		      sumGreen = sumGreen + getSum(i,j,kernel,green,l,k);
		      sumBlue = sumBlue + getSum(i,j,kernel,blue,l,k);
		  } 
		}
		redLap[i+j*w] = sumRed;
		greenLap[i+j*w] =sumGreen;
		blueLap[i+j*w] =sumBlue;
	    
		blue++;
		red++;
		green++;
	    
	    
	    
	      }
	      
		
		
	}
	
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		
	  for (int i = 0; i < pSrc->GetWidth(); i++)
	  {
	    double mean;
	   
	   
	  
	    mean = getLaplacianMean(i,j,redLap,w);
	    lapWeightMapR[i+j*w]=mean+std::abs(*redLap);
	    globWeightMapR[i+j*w]=std::pow((*redLap)-mean,2);
	  
	    mean = getLaplacianMean(i,j,greenLap,w);
	    lapWeightMapG[i+j*w]=mean+std::abs(*greenLap);
	    globWeightMapG[i+j*w]=std::pow((*greenLap)-mean,2);
	  
	    mean = getLaplacianMean(i,j,blueLap,w);
	    lapWeightMapB[i+j*w]=mean+std::abs(*blueLap);
	    globWeightMapB[i+j*w]=std::pow((*blueLap)-mean,2);
	    *pDestinationData++ = globWeightMapB[i+j*w];
	    *pDestinationData++ = globWeightMapB[i+j*w];
	    *pDestinationData++ =globWeightMapB[i+j*w];
	      redLap++;
	      greenLap++;
	      blueLap++;
	  
	  } ///wrok i n progress waintg fot mr ancutis reasponse
	}
	//printf("%f\n",r[0]);
	
	
	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	//free(red);
	return 0;
}

////wrok in progresss gottamake it neater, was flustrated when doing it
double TMOAncuti16::getLaplacianMean(int i, int j, double* laplacianOfColor, int w)
{
  double mean;
  if(i == 0 && j==0)
  {
    mean = *laplacianOfColor + *(laplacianOfColor+1) + *(laplacianOfColor) + *(laplacianOfColor) + *(laplacianOfColor + w) +
	*(laplacianOfColor+1+w) + *(laplacianOfColor+w) + *(laplacianOfColor)+*(laplacianOfColor+1);

	return mean/9;
  }
  if(i==0)
  {
    mean = *laplacianOfColor + *(laplacianOfColor+1) + *(laplacianOfColor) + *(laplacianOfColor - w) + *(laplacianOfColor + w) +
	*(laplacianOfColor+1+w) + *(laplacianOfColor+w) + *(laplacianOfColor-w)+*(laplacianOfColor+1-w);

	return mean/9;
  }
  if(j==0)
  {
    mean = *laplacianOfColor + *(laplacianOfColor+1) + *(laplacianOfColor-1) + *(laplacianOfColor) + *(laplacianOfColor + w) +
	*(laplacianOfColor+1+w) + *(laplacianOfColor-1+w) + *(laplacianOfColor-1)+*(laplacianOfColor+1);
	
	return mean/9;
  }
  if(i==pSrc->GetWidth()-1 && j==pSrc->GetHeight()-1)
  {
    mean = *laplacianOfColor + *(laplacianOfColor) + *(laplacianOfColor-1) + *(laplacianOfColor - w) + *(laplacianOfColor) +
	*(laplacianOfColor) + *(laplacianOfColor-1) + *(laplacianOfColor-1-w)+*(laplacianOfColor-w);
	
	return mean/9;
  }
  if(i==pSrc->GetWidth()-1 )
  {
    mean = *laplacianOfColor + *(laplacianOfColor) + *(laplacianOfColor-1) + *(laplacianOfColor - w) + *(laplacianOfColor + w) +
	*(laplacianOfColor+w) + *(laplacianOfColor-1+w) + *(laplacianOfColor-1-w)+*(laplacianOfColor-w);
	
	return mean/9;
  }
  if(j==pSrc->GetHeight()-1)
  {
      mean = *laplacianOfColor + *(laplacianOfColor+1) + *(laplacianOfColor-1) + *(laplacianOfColor - w) + *(laplacianOfColor ) +
	*(laplacianOfColor+1) + *(laplacianOfColor-1) + *(laplacianOfColor-1-w)+*(laplacianOfColor+1-w);
	
	return mean/9;
  }
  mean = *laplacianOfColor + *(laplacianOfColor+1) + *(laplacianOfColor-1) + *(laplacianOfColor - w) + *(laplacianOfColor + w) +
		    *(laplacianOfColor+1+w) + *(laplacianOfColor-1+w) + *(laplacianOfColor-1-w)+*(laplacianOfColor+1-w);
		    
  return mean/9;
  
}
double TMOAncuti16::getSum(int i, int j, float kernel[3][3], double* colorChannel, int l, int k)
{
  double sum;
  if(i==0 || j==0 || i==pSrc->GetWidth()-1 || j==pSrc->GetHeight()-1)
  {
    if(i==0 && l==1)
    {
      if(j==0 && k==1)
      {
	sum =  kernel[l+1][k+1] * *(colorChannel );
	return sum;
      }
      else 
      {
	sum = kernel[l+1][k+1] * *((colorChannel ) + (pSrc->GetWidth() * (-k)));
	return sum;
      }
      
    }
    if(i==pSrc->GetWidth()-1 && l==-1)
    {
      if(j==pSrc->GetHeight()-1 && k==-1)
      {
	sum=kernel[l+1][k+1] * *(colorChannel );
	return sum;

      }
      else 
      {
	sum=kernel[l+1][k+1] * *((colorChannel  ) + (pSrc->GetWidth() * (-k)));
	return sum;

      }
    }
  }
  
  
    sum= kernel[l+1][k+1] * *((colorChannel - l) + (pSrc->GetWidth() * (-k)));
    int a = kernel[l+1][k+1];
    return sum;

  
}

