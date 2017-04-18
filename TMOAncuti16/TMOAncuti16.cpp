/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAncuti16.h"
#include <boost/concept_check.hpp>
#include <complex>



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
	
	double *normWeightMapR =(double*)malloc( h * w * sizeof(double));
	double *normWeightMapG =(double*)malloc( h * w * sizeof(double));
	double *normWeightMapB =(double*)malloc( h * w * sizeof(double));
	
	float kernel[3][3] = {{0,-1,0},
			      {-1,4,-1}, ///laplacian kernel
			      {0,-1,0}};
	double sumRed=0.0;
	double sumGreen=0.0;
	double sumBlue=0.0;
	
	double max=0;
	double res=0;
        
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
	   
	   
	  ////for red channel
	    mean = getLaplacianMean(i,j,redLap,w);   ////average of the laplacian
	    lapWeightMapR[i+j*w]=mean+std::abs(*redLap);  ///computation of laplacian weight map
	    globWeightMapR[i+j*w]=std::pow((red[i+j*w]-mean),2); ///global weiht map
	    normWeightMapR[i+j*w]= /*globWeightMapR[i+j*w]/*/(lapWeightMapR[i+j*w]+globWeightMapR[i+j*w]);// +  ////computation of normalised weight map
				   // lapWeightMapR[i+j*w]/(lapWeightMapR[i+j*w]+globWeightMapR[i+j*w]);/////mr ancuti didnt reply but i figured it out myself
	 ////for green channel
	    mean = getLaplacianMean(i,j,greenLap,w);
	    lapWeightMapG[i+j*w]=mean+std::abs(*greenLap);
	    globWeightMapG[i+j*w]=std::pow((green[i+j*w]-mean),2);   ///see upwards
	    normWeightMapG[i+j*w]= /*globWeightMapG[i+j*w]/*/(lapWeightMapG[i+j*w]+globWeightMapG[i+j*w]);//+
				   // lapWeightMapG[i+j*w]/(lapWeightMapG[i+j*w]+globWeightMapG[i+j*w]);
	   ////for blue channel
	    mean = getLaplacianMean(i,j,blueLap,w);
	    lapWeightMapB[i+j*w]=mean+std::abs(*blueLap);    ////see upwards
	    globWeightMapB[i+j*w]=std::pow((blue[i+j*w]-mean),2);
	    normWeightMapB[i+j*w]= /*globWeightMapB[i+j*w]/*/(lapWeightMapB[i+j*w]+globWeightMapB[i+j*w]);//+
				  //  lapWeightMapB[i+j*w]/(lapWeightMapB[i+j*w]+globWeightMapB[i+j*w]);
	   
	    max=normWeightMapR[i+j*w]+normWeightMapG[i+j*w]+normWeightMapB[i+j*w]; ///the normalised weight maps must add to 1, 
									/////i must get the max to normalise them
	    redLap++;
	    greenLap++;
	    blueLap++;
	     normWeightMapR[i+j*w]=normWeightMapR[i+j*w]/max;
	      normWeightMapG[i+j*w]=normWeightMapG[i+j*w]/max;
	      normWeightMapB[i+j*w]=normWeightMapB[i+j*w]/max;
	      
	    res=normWeightMapR[i+j*w]*red[i+j*w]+       ////result normalising each weight map and multiply with pixel value from each channel
		  normWeightMapG[i+j*w]*green[i+j*w]+
		  normWeightMapB[i+j*w]*blue[i+j*w];
	    
	    
	 
	   /* *pDestinationData++ =res;
	    *pDestinationData++ = res;
	    *pDestinationData++ =res;*/
	  
	  } 
	}
	
	redLap = redLap - w * h;
	greenLap = greenLap - w * h;  ///reseting the pointers
	blueLap = blueLap - w * h;
	cv::Mat NWMR, NWMG, NWMB,redMat, greenMat, blueMat, endResult, jj,jj2,jj3;
	
	cv::Mat channel[3], dst,tmp, final, tmp2,tmp3;
	NWMR = cv::Mat (h, w, CV_64FC1, normWeightMapR);
	NWMG = cv::Mat (h, w, CV_64FC1, normWeightMapG);
	NWMB = cv::Mat (h, w, CV_64FC1, normWeightMapB);
	
	redMat = cv::Mat (h, w, CV_64FC1, red);
	greenMat = cv::Mat (h, w, CV_64FC1, green);
	blueMat = cv::Mat (h, w, CV_64FC1, blue);
	jj=cv::Mat (h, w, CV_64FC1, redLap);
	std::vector<cv::Mat> finalPyramid;
	
	channel[0] = cv::Mat::zeros(h, w, CV_64FC1);
	channel[2] = cv::Mat::zeros(h, w, CV_64FC1);
	channel[1]=cv::Mat (h, w, CV_64FC1, green);
	cv::Mat A;
	A=cv::Mat::zeros(h, w, CV_64FC3);
	cv::merge(channel,3,A);
	endResult = cv::Mat::zeros(h, w, CV_64FC1);
	 tmp=A;
	 dst=tmp;
	 //cv::imshow(" window" ,greenMat);
	// int key = cv::waitKey(2000);
	 /*for(int i=0; i<3;i++)
	 {
	   //cv::pyrDown(tmp,tmp,cv::Size(tmp.cols/2,tmp.rows/2));
	   //cv::GaussianBlur(tmp,tmp,cv::Size(5,5),0,0);
	   
	 //tmp=dst;
	 }*/
	 cv::GaussianBlur(redMat,tmp,cv::Size(5,5),0,0);
	 cv::GaussianBlur(greenMat,tmp2,cv::Size(5,5),0,0);
	 cv::GaussianBlur(blueMat,tmp3,cv::Size(5,5),0,0);
	 //redMat = redMat - tmp;
	 final=NWMR.mul(jj) +NWMG.mul(greenMat -tmp2 ) +NWMB.mul( blueMat -tmp3);
	      finalPyramid.push_back(final);
	      cv::imshow(" window" ,final);
	 int key = cv::waitKey(2000);
	      
	  for(int i=0; i<7;i++)
	 {
	   cv::pyrDown(NWMR, NWMR, cv::Size(NWMR.cols/2,NWMR.rows/2));
	  cv::pyrDown(redMat, redMat, cv::Size(redMat.cols/2,redMat.rows/2));
	   // cv::GaussianBlur(NWMR,NWMR,cv::Size(5,5),0,0);
	  //  cv::GaussianBlur(redMat,redMat,cv::Size(5,5),0,0);
	   cv::GaussianBlur(redMat,tmp,cv::Size(5,5),0,0);
	   
	   cv::pyrDown(NWMG, NWMG, cv::Size(NWMG.cols/2,NWMG.rows/2));
	  cv::pyrDown(greenMat, greenMat, cv::Size(greenMat.cols/2,greenMat.rows/2));
	   //cv::GaussianBlur(NWMG,NWMG,cv::Size(5,5),0,0);
	  // cv::GaussianBlur(greenMat,greenMat,cv::Size(5,5),0,0);
	   cv::GaussianBlur(greenMat,tmp2,cv::Size(5,5),0,0);
	   
	   cv::pyrDown(NWMB, NWMB, cv::Size(NWMB.cols/2,NWMB.rows/2));
	   cv::pyrDown(blueMat, blueMat, cv::Size(blueMat.cols/2,blueMat.rows/2));
	   //cv::GaussianBlur(NWMB,NWMB,cv::Size(5,5),0,0);
	  // cv::GaussianBlur(blueMat,blueMat,cv::Size(5,5),0,0);
	   cv::GaussianBlur(blueMat,tmp3,cv::Size(5,5),0,0);
 
	   final=NWMR.mul(redMat- tmp)+ NWMG.mul(greenMat - tmp2) + NWMG.mul(blueMat - tmp3);
	// cv::normalize(NWMR,tmp,0.0,255.0,cv::NORM_MINMAX,CV_64FC1);
	   
	  //  
	  finalPyramid.push_back(final);
	 }
	
	 for(int i = 1; i<=7; i++)
	 {
	   tmp = finalPyramid[i];
	   for(int j= 0; j<i;j++)
	   {
	   cv::pyrUp(tmp,tmp,cv::Size(tmp.cols*2,tmp.rows*2));
	
	
 	//cv::imshow(" window" ,tmp);
	//int key = cv::waitKey(3000);
	     
	  //   cv::pyrUp(tmp,tmp,cv::Size(tmp.cols*2,tmp.rows*2));
	  }
	  endResult = endResult + tmp;
	 }
	  
	 // endResult = finalPyramid[0];
	 for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		
	  for (int i = 0; i < pSrc->GetWidth(); i++) ///creating weight maps
	  {
	   *pDestinationData++ =endResult.at<double>(j,i);
	    *pDestinationData++ = endResult.at<double>(j,i);
	    *pDestinationData++ =endResult.at<double>(j,i);
	  }
	}
	  // endResult = endResult + tmp;
	 
	// }
	
	return 0;
}


////probably not needed
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

