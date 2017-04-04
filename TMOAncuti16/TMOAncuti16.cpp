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
	SetDescription(L"Image decolorization using Laplacian operator and multi-scale fusion");	// TODO - Insert description

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
	
	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
	int h=	pSrc->GetHeight();
	int w =pSrc->GetWidth();
	//pDst->SetDimensions(w/2,h/2);
	
	// Separable kernel
    float kernelX[5] = { 1/16.0f,  4/16.0f,  6/16.0f,  4/16.0f, 1/16.0f };
    float kernelY[5] = { 1/16.0f,  4/16.0f,  6/16.0f,  4/16.0f, 1/16.0f };
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
			      {-1,4,-1}, ///laplacian kernel
			      {0,-1,0}};
	double sumRed=0.0;
	double sumGreen=0.0;
	double sumBlue=0.0;
	
	double lapMapMax=0;
	double lapMapMin=0;

        
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	//getting separate RGB channels
 		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			red[i+j*w] = *pSourceData++;
			green[i+j*w] = *pSourceData++;
			blue [i+j*w]= *pSourceData++;

		}
	}
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		
	    pSrc->ProgressBar(j, pSrc->GetHeight());	// get the laplacian for each channel
	    for (int i = 0; i < pSrc->GetWidth(); i++)
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
	red = red - w * h;
	green = green - w * h;  ///reseting the pointers
	blue = blue - w * h;
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		
	  for (int i = 0; i < pSrc->GetWidth(); i++) ///creating weight maps
	  {
	    double mean;
	   
	   
	  
	    mean = getLaplacianMean(i,j,redLap,w);
	    lapWeightMapR[i+j*w]=mean+std::abs(*redLap);
	    globWeightMapR[i+j*w]=std::pow((red[i+j*w]-mean),2);
	    double g = lapWeightMapR[i+j*w] / (lapWeightMapR[i+j*w]+globWeightMapR[i+j*w]);
	    
	 // if(i>=2 && i<pSrc->GetWidth()-2  &&j>=2 && j<pSrc->GetHeight()-2 )
	 //  {
	 //    
	     
	 //  }
	  
	    mean = getLaplacianMean(i,j,greenLap,w);
	    lapWeightMapG[i+j*w]=mean+std::abs(*greenLap);
	    globWeightMapG[i+j*w]=std::pow((green[i+j*w]-mean),2);
	  
	    mean = getLaplacianMean(i,j,blueLap,w);
	    lapWeightMapB[i+j*w]=mean+std::abs(*blueLap);
	    globWeightMapB[i+j*w]=std::pow((blue[i+j*w]-mean),2);
	    
	   // if(lapMapMin > lapWeightMapG[i+j*w] ) lapMapMin = lapWeightMapG[i+j*w];
	   // else if( lapMapMax < lapWeightMapG[i+j*w]) lapMapMax = lapWeightMapG[i+j*w]; 
	      
	  
	      redLap++;
	      greenLap++;
	      blueLap++;
	      
	      
	  
	  } 
	}
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		
	  for (int i = 0; i < pSrc->GetWidth(); i++) ///creating weight maps
	  {
	    globWeightMapG[i+j*w]=getGaussianBlurPix(i, j, kernelX, kernelY, red, w);
	    *pDestinationData++ = globWeightMapG[i+j*w];
	    *pDestinationData++ = 0;//globWeightMapG[i+j*w];
	    *pDestinationData++ =0;//globWeightMapG[i+j*w];
	  }
	}

	
	
	
	
	//pSrc->ProgressBar(j, pSrc->GetHeight());
	
	//free(red);
	return 0;
}
double TMOAncuti16::getGaussianBlurPix(int i, int j, float kernelX[5], float kernelY[5], double* map, int w)
{
  double tmp[5]={0,0,0,0,0};
  int p= 0;
  double res=0;
   double e =0;
  for(int k=-2; k<=2; k++)
  {
    int  kernelPos=0;
    for(int l = -2; l<=2; l++)
    {
      if(i<2 && j<2 && k<0 && l<0) tmp[p]+=*(map +((i)+(j)*w))* kernelY[kernelPos];
      else if(i<2 && k<0) tmp[p]+=*(map +((i)+(j+l)*w))* kernelY[kernelPos];
      else if(j<2 && l<0) tmp[p]+=*(map +((i+k)+(j)*w))* kernelY[kernelPos];
      else if(i> pSrc->GetWidth()-2 && k>0) tmp[p]+=*(map +((i)+(j+l)*w))* kernelY[kernelPos];
       else if(j> pSrc->GetHeight()-2 && l>0) tmp[p]+=*(map +((i+k)+(j)*w))* kernelY[kernelPos];
       else if(i> pSrc->GetWidth()-2 && j> pSrc->GetHeight()-2 && k>0 && l>0) tmp[p]+=*(map +((i)+(j)*w))* kernelY[kernelPos];
      else tmp[p]+=*(map +((i+k)+(j+l)*w))* kernelY[kernelPos];
      
     
      e=tmp[p];
      kernelPos++;
    }
    res += tmp[p] * kernelX[p];
    p++;
  }
  return res;
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

/// getting the sum for the computation of the laplacian
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

