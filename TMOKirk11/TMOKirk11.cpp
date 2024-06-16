/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2023-24                                           *
*                                                                              *
*                       Implementation of the TMOKirk11 class                  *
*                       LMSR matrix isnt from verified source - bad results    *
*                                                                              *
*                       Author: Lukas Macejka (xmacej03)                       *
*                       Mail: xmacej03@vutbr.cz                                *
*                                                                              *
*******************************************************************************/
#include "TMOKirk11.h"
#include "./BF/linear_bf.h"
#include <stdlib.h>
#include <math.h> 
#include <cmath>
#include <iostream>
#include <string>
#include "opencv2/opencv.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include "./BF/array.h"
#include <fftw3.h>
#include "./BF/fft_3D.h"
typedef Array_2D<double> image_array;
using namespace cv;
using namespace std;
namespace FFT{
  

  
  unsigned              Support_3D::fftw_flags       = FFTW_EXHAUSTIVE;
  bool                  Support_3D::wisdom_file_set  = false;
  std::string           Support_3D::wisdom_file      = std::string();
  Support_3D::size_type Support_3D::support_number   = 0;
  bool                  Support_3D::auto_save_wisdom = true;
  bool                  Support_3D::wisdom_loaded    = false;
}

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOKirk11::TMOKirk11()
{
	SetName(L"Kirk11");					  
	SetDescription(L"Perceptually Based Tone Mapping for Low-Light Conditions"); 

	dParameter.SetName(L"ParameterName");				
	dParameter.SetDescription(L"ParameterDescription"); 
	dParameter.SetDefault(1);
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); 
	this->Register(dParameter);
}

TMOKirk11::~TMOKirk11()
{
}
/*
 * Convert TMOImage to cv::Mat
 Function taken from TMOParis11 file
 */
Mat TMOImage2Mat(TMOImage* pSrc, Mat *input)
{
	Vec3d* ptrMat;
	int rowsCnt, colsCnt;

	rowsCnt = pSrc->GetHeight();
	colsCnt = pSrc->GetWidth();

	Mat srcConvMat(rowsCnt, colsCnt, CV_32FC3);
	
   Mat tempMat(rowsCnt*colsCnt,1,CV_32FC3);
   int i = 0;
	
		for (int x = 0; x < colsCnt; x++)
   {
   for (int y = 0; y < rowsCnt; y++)
		{
         double R = input->at<float>(0,y*colsCnt + x);
         double G = input->at<float>(1,y*colsCnt + x);
         double B = input->at<float>(2,y*colsCnt + x);
			tempMat.at<Vec3f>(i,0) = Vec3f(B,G,R);
         i++;
		}
	}
      
	/*
	 * If data are continuous, it is possible
	 * to work in one for loop
	 */
	//#pragma omp parallel for collapse(2);

   i=0;
	for (int y = 0; y < rowsCnt ; y++)
	{
		for (int x = 0; x <  colsCnt ; x++)
		{
         srcConvMat.at<Vec3f>(y,x) = tempMat.at<Vec3f>(i,0);
		   i++;
      }
	}
	return srcConvMat;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOKirk11::Transform()
{
    // Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_RGB); // This is format of Y as luminance
	pDst->Convert(TMO_RGB); // x, y as color information
	double *pSourceData = pSrc->GetData();		
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
	double pr, pg, pb;
   int pixels_c = pSrc->GetHeight() * pSrc->GetWidth();
   Mat RGB_M(3,pixels_c,CV_32FC1);
   double *p_source_tmp = pSourceData;

   for(int i=0;i<pSrc->GetWidth();i++){
      for(int j=0;j<pSrc->GetHeight();j++){

         for(int k=0;k<3;k++){
            RGB_M.at<float>(k,(j*pSrc->GetWidth())+i) = (*p_source_tmp);
            p_source_tmp++;
         }

      }
   }
   Mat tester = TMOImage2Mat(pSrc,&RGB_M);

   Mat ConvertMtx(4,3,CV_32FC1);

  /**
   * for ldr photos
  ConvertMtx.at<float>(0,0) = 0.096869562190332;
  ConvertMtx.at<float>(0,1) = 0.318940374720484;
  ConvertMtx.at<float>(0,2) = -0.188428411786113;
  ConvertMtx.at<float>(1,0) = 0.020208210904239;
  ConvertMtx.at<float>(1,1) = 0.291385283197581;
  ConvertMtx.at<float>(1,2) = -0.090918262127325;
  ConvertMtx.at<float>(2,0) = 0.002760510899553;
  ConvertMtx.at<float>(2,1) = -0.008341563564118;
  ConvertMtx.at<float>(2,2) = 0.067213551661950;
  ConvertMtx.at<float>(3,0) = -0.007607045462440;
  ConvertMtx.at<float>(3,1) = 0.122492925567539;
  ConvertMtx.at<float>(3,2) = 0.022445835141881;*/

   /**
	 * HDR photos
	 */
   ConvertMtx.at<float>(0,0) = 0.000167179131;
	ConvertMtx.at<float>(0,1) = 0.021378405283;
	ConvertMtx.at<float>(0,2) = -0.000600300420;
	ConvertMtx.at<float>(1,0) = 0.000037285102;
	ConvertMtx.at<float>(1,1) = 0.017335413184;
	ConvertMtx.at<float>(1,2) = -0.000607142696;
	ConvertMtx.at<float>(2,0) = -0.000384788819;
	ConvertMtx.at<float>(2,1) = 0.011063773501;
	ConvertMtx.at<float>(2,2) = -0.000583754041;
	ConvertMtx.at<float>(3,0) = -0.000211589221;
	ConvertMtx.at<float>(3,1) = 0.018000447150;
	ConvertMtx.at<float>(3,2) = -0.000621909939;

	/**
	 * calculate the LMSR matrix
	 */

   Mat LMSRMtx(4,pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1); 
   LMSRMtx = ConvertMtx*RGB_M;
   
   Mat G_LMS(3,pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1);
   for(int i = 0;i< pSrc->GetWidth()*pSrc->GetHeight();i++){
      //Long 1/(1 + 0.33(𝑞𝐿𝑜𝑛𝑔 + 𝑘1𝑞𝑅𝑜𝑑))^0.5
      G_LMS.at<float>(0,i) = 1/pow(1+0.33*(LMSRMtx.at<float>(0,i) + 0.25 * LMSRMtx.at<float>(3,i)),0.5);
      //Medium 1/(1 + 0.33(𝑞𝑀 𝑒𝑑𝑖𝑢𝑚 + 𝑘1𝑞𝑅𝑜𝑑))0.5
      G_LMS.at<float>(1,i) = 1/pow(1+0.33*(LMSRMtx.at<float>(1,i) + 0.25 * LMSRMtx.at<float>(3,i)),0.5);
      //Short 1/(1 + 0.33(𝑞𝑆ℎ𝑜𝑟𝑡 + 𝑘2𝑞𝑅𝑜𝑑))0.5
      G_LMS.at<float>(2,i) = 1/pow(1+0.33*(LMSRMtx.at<float>(2,i) + 0.4 * LMSRMtx.at<float>(3,i)),0.5); 
   } 
   
   Mat O_RGBYL(3,pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1);
   
   Mat Mat_GStrieska(3,pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1);
   //𝑥 = 𝑦 = 15 a 𝑧 = 5
   //𝜌1 = 1.111, 𝜌2 = 0.939, 𝜌3 = 0.4,𝜌4 = 0.15 a 𝛼 = 0.619
   //𝑙𝑚𝑎𝑥 = 0.637, 𝑚𝑚𝑎𝑥 = 0.392 a 𝑠𝑚𝑎𝑥 = 1.606
   Mat alpha(1,pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1);
   for(int i = 0; i< pSrc->GetWidth()*pSrc->GetHeight();i++){
      //∆o𝑅𝑒𝑑/𝐺𝑟𝑒𝑒𝑛 = 𝑥𝑘1(𝜌1𝑔𝑀𝑒𝑑𝑖𝑢𝑚/𝑚𝑚𝑎𝑥− 𝜌2𝑔𝐿𝑜𝑛𝑔/𝑙𝑚𝑎𝑥)𝑞𝑅𝑜𝑑
      O_RGBYL.at<float>(0,i) = 15 * 0.25 * (1.111 * G_LMS.at<float>(1,i) / 0.392 - 0.939 * G_LMS.at<float>(0,i) / 0.637) * LMSRMtx.at<float>(3,i);
      //∆o𝐵𝑙𝑢𝑒/𝑌𝑒𝑙𝑙𝑜𝑤 = 𝑦(𝜌3𝑔𝑆ℎ𝑜𝑟𝑡/𝑠𝑚𝑎𝑥− 𝜌4(𝛼 𝑔𝐿𝑜𝑛𝑔/𝑙𝑚𝑎𝑥+ (1 − 𝛼)𝑔𝑀𝑒𝑑𝑖𝑢𝑚/𝑚𝑚𝑎𝑥))𝑞𝑅𝑜𝑑
      O_RGBYL.at<float>(1,i) = 15 * (0.4 * (G_LMS.at<float>(2,i) / 1.606) - 0.15 * ( 0.619 * G_LMS.at<float>(0,i) / 0.637 + (1-0.619) * G_LMS.at<float>(1,i) / 0.392)) * LMSRMtx.at<float>(3,i);
      //∆o𝐿𝑢𝑚𝑖𝑛𝑎𝑛𝑐𝑒 = 𝑧(𝛼 𝑔𝐿𝑜𝑛𝑔/𝑙𝑚𝑎𝑥+ (1 − 𝛼) 𝑔𝑀𝑒𝑑𝑖𝑢𝑚/𝑚𝑚𝑎𝑥)𝑞𝑅𝑜𝑑
      O_RGBYL.at<float>(2,i) = 5 * ( 0.619 * G_LMS.at<float>(0,i) / 0.637 + ( 1 - 0.619 ) * G_LMS.at<float>(1,i) / 0.392 ) * LMSRMtx.at<float>(3,i);

      alpha.at<float>(0,i) = 5 * ( 0.619 * G_LMS.at<float>(0,i) / 0.637 + ( 1 - 0.619 ) * G_LMS.at<float>(1,i) / 0.392 );
      if(alpha.at<float>(0,i)<0)alpha.at<float>(0,i)=0;


     float magic = (O_RGBYL.at<float>(2,i) - O_RGBYL.at<float>(0,i))/2;

      Mat_GStrieska.at<float>(0,i) = LMSRMtx.at<float>(0,i) + magic;
      Mat_GStrieska.at<float>(1,i) = LMSRMtx.at<float>(1,i) + O_RGBYL.at<float>(2,i) - magic;
      Mat_GStrieska.at<float>(2,i) = LMSRMtx.at<float>(2,i) + O_RGBYL.at<float>(1,i) + O_RGBYL.at<float>(2,i);

      Mat_GStrieska.at<float>(0,i) += LMSRMtx.at<float>(0,i);
      Mat_GStrieska.at<float>(1,i) += LMSRMtx.at<float>(1,i);
      Mat_GStrieska.at<float>(2,i) += LMSRMtx.at<float>(2,i);
   }
   /**
	 * test monitor settings LMS * RGB
	 */
   
   Mat monitorLMS(3,3,CV_32FC1);
   monitorLMS.at<float>(0,0) = 3.737001334053387; 
  monitorLMS.at<float>(0,1) = 4.474289710481138; 
  monitorLMS.at<float>(0,2) = 0.446763557036099; 
  monitorLMS.at<float>(1,0) = 1.260813115018955; 
  monitorLMS.at<float>(1,1) = 4.394344052644990; 
  monitorLMS.at<float>(1,2) = 0.587779102275466; 
  monitorLMS.at<float>(2,0) = 0.046245589359130; 
  monitorLMS.at<float>(2,1) = 0.244453440088045; 
  monitorLMS.at<float>(2,2) = 1.218879359722785;

   Mat out(3, pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1);
   out = monitorLMS.inv() * Mat_GStrieska;

   /**
	 * Billateral filter aplication
	 */

   Mat bilateral_in = TMOImage2Mat(pSrc,&out);


   cv::Mat outputImage;
   cv::bilateralFilter(bilateral_in, outputImage, 5, 75, 75);  
   
   int i_w = pSrc->GetWidth();
   int i_h = pSrc->GetHeight();

   int rowsCnt = pSrc->GetHeight();
	int colsCnt = pSrc->GetWidth();
   Mat outh(1, pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC3);
   unsigned i=0;
	for (int y = 0; y < rowsCnt ; y++)
	{
		for (int x = 0; x <  colsCnt ; x++)
		{
         outh.at<Vec3f>(0,i) = outputImage.at<Vec3f>(y,x);
		   i++;
      }
	}
   i=0;
   for (int x = 0; x < colsCnt; x++)
   {
   for (int y = 0; y < rowsCnt; y++)
		{
         float r =outh.at<Vec3f>(0,i)[0];
         float g  =outh.at<Vec3f>(0,i)[1];
         float b =outh.at<Vec3f>(0,i)[2];
         out.at<float>(0,y*colsCnt + x) = b;
         out.at<float>(1,y*colsCnt + x) = g;
         out.at<float>(2,y*colsCnt + x) = r;       
         i++;
		}
	}

  ////////////////////////////////////////////////////////////////////////////////////////////
   double x = 0.1f;
   Mat* Idouble = new Mat(3,pSrc->GetWidth()* pSrc->GetHeight(), CV_32FC1);
   Mat* MaxImage = new Mat(1,pSrc->GetWidth()*pSrc->GetHeight(),CV_32FC1);
   Mat* alphaScaling = new Mat(1,pSrc->GetWidth()* pSrc->GetHeight(), CV_32FC1);
   Mat* outputMtx = new Mat(3,pSrc->GetWidth()* pSrc->GetHeight(), CV_32FC1);

   for(int i = 0; i< pSrc->GetWidth()*pSrc->GetHeight();i++){
      Idouble->at<float>(0,i) = ((double) out.at<float>(0,i))/256.0;
      Idouble->at<float>(1,i) = ((double) out.at<float>(1,i))/256.0;
      Idouble->at<float>(2,i) = ((double) out.at<float>(2,i))/256.0;

      MaxImage->at<float>(0,i) = 255.0;
      alphaScaling->at<float>(0,i) = (alpha.at<float>(0,i)*(1-x))+x;
      MaxImage->at<float>(0,i) *= alphaScaling->at<float>(0,i);

      outputMtx->at<float>(0,i) = (Idouble->at<float>(0,i)*MaxImage->at<float>(0,i));
      outputMtx->at<float>(1,i) = (Idouble->at<float>(1,i)*MaxImage->at<float>(0,i));
      outputMtx->at<float>(2,i) = (Idouble->at<float>(2,i)*MaxImage->at<float>(0,i));

      double colorBlendStart = 210.0;
      double intensity = 0.0;
      double bright_factor = 0.0;
      intensity = 0.3 * RGB_M.at<float>(0,i) +  0.6 * RGB_M.at<float>(1,i) + 0.1 * RGB_M.at<float>(2,i);
      if (intensity >= 210.0) 
      {
         bright_factor = 0;
      }
      double finalBlendR = (1-bright_factor)*alpha.at<float>(0,i) + bright_factor;
      double finalBlendG = (1-bright_factor)*alpha.at<float>(1,i) + bright_factor;
      double finalBlendB = (1-bright_factor)*alpha.at<float>(2,i) + bright_factor;

      outputMtx->at<float>(0,i) = (finalBlendR * RGB_M.at<float>(0,i) + ((1- finalBlendR) * outputMtx->at<float>(0,i)));
      outputMtx->at<float>(1,i) = (finalBlendG * RGB_M.at<float>(1,i)+ ((1- finalBlendG) * outputMtx->at<float>(1,i)));
      outputMtx->at<float>(2,i) = (finalBlendB * RGB_M.at<float>(2,i) + ((1- finalBlendB) * outputMtx->at<float>(2,i)));
   }

	for(int i = 0; i < pSrc->GetWidth(); i++)
	{
		for (int j = 0; j < pSrc->GetHeight(); j++)
		{
         if (outputMtx->at<float>(0,i+pSrc->GetWidth()*j) < 0.0f)  outputMtx->at<float>(0,i+pSrc->GetWidth()*j) = 0.0f; 
         if (outputMtx->at<float>(1,i+pSrc->GetWidth()*j) < 0.0f)  outputMtx->at<float>(1,i+pSrc->GetWidth()*j) = 0.0f; 
         if (outputMtx->at<float>(2,i+pSrc->GetWidth()*j) < 0.0f)  outputMtx->at<float>(2,i+pSrc->GetWidth()*j) = 0.0f; 
			
         float w = ( 0.619 * G_LMS.at<float>(0,i) / 0.637 + (1-0.619) * G_LMS.at<float>(1,i) / 0.392);

			*pDestinationData++ = out.at<float>(0,j*pSrc->GetWidth()+i) * w;//r
			*pDestinationData++ = out.at<float>(1,j*pSrc->GetWidth()+i) * w;//g
			*pDestinationData++ = out.at<float>(2,j*pSrc->GetWidth()+i) * w;//b
		}
	}
	return 0;
}
