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
	;
	
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
		/*if (j==1) 
		{
		  blue=(blue+w);
		  red=(red+w);
		  green=(green+w);
		}*/
		
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
 		for (int i = 0; i < pSrc->GetWidth(); i++)  ////vynechvam okraje => treba sa spytat ako to riesit
		{
		  sumRed = 0.0;
		  sumGreen = 0.0;
		  sumBlue = 0.0;
		  blue++;
		  red++;
		  green++;
		  
		  //g=green[i+j*w];
		 // *green=*(green+i+j*w);
		  
		  /**pDestinationData++ = *red++;
			*pDestinationData++ = *green++;
			*pDestinationData++ = *blue++;*/
		
		
		  
		    for(int k = -1; k <=1; k++)
		    {
		      for(int l = -1; l <=1; l++)
		      {
			
			if(i==0 || j==0 || i==pSrc->GetWidth()-1 || j==pSrc->GetHeight()-1)
			{
			  if(i==0 && l==1)
			  {
			    if(j==0 && k==1)
			    {
			      sumRed = sumRed + kernel[l+1][k+1] * *(red );
			      sumGreen = sumGreen + kernel[l+1][k+1] * *(green );
			      sumBlue= sumBlue + kernel[l+1][k+1] * *(blue);
			    }
			    else 
			    {
			      sumRed = sumRed + kernel[l+1][k+1] * *((red ) + (w * (-k)));
			      sumGreen = sumGreen + kernel[l+1][k+1] * *((green ) + (w * (-k)));
			      sumBlue = sumBlue + kernel[l+1][k+1] * *((blue ) + (w * (-k)));
			    }
			    
			  }
			  if(i==pSrc->GetWidth()-1 && l==-1)
			  {
			    if(j==pSrc->GetHeight()-1 && k==-1)
			    {
			      sumRed = sumRed + kernel[l+1][k+1] * *(red );
			      sumGreen = sumGreen + kernel[l+1][k+1] * *(green );
			      sumBlue= sumBlue + kernel[l+1][k+1] * *(blue);
			    }
			    else 
			    {
			      sumRed = sumRed + kernel[l+1][k+1] * *((red ) + (w * (-k)));
			      sumGreen = sumGreen + kernel[l+1][k+1] * *((green ) + (w * (-k)));
			      sumBlue = sumBlue + kernel[l+1][k+1] * *((blue ) + (w * (-k)));
			    }
			  }
			}
			else
			{
			  sumRed = sumRed + kernel[l+1][k+1] * *((red - l) + (w * (-k)));
			  sumGreen = sumGreen + kernel[l+1][k+1] * *((green - l) + (w * (-k)));
			  sumBlue = sumBlue + kernel[l+1][k+1] * *((blue - l) + (w * (-k)));
			}
			
		  /* sumRed = sumRed + getSum(i,j,kernel,red,l,k);
		    sumGreen = sumGreen + getSum(i,j,kernel,green,l,k);
		   sumBlue = sumBlue + getSum(i,j,kernel,blue,l,k);*/
			 
			  
			
		      } ///rozdiel medzi tym ked to dam cez funkciu a ked to necham tuna treba dokoncit
		    }
		    *redLap++ = sumRed;
		    *greenLap++ =sumGreen;
		  *blueLap++ =sumBlue;
		*pDestinationData++ = 0;
		*pDestinationData++ = sumGreen;
		*pDestinationData++ = 0;
		
		}
		
		
		
	}
	//printf("%f\n",r[0]);
	
	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	//free(red);
	return 0;
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
    return sum;

  
}

