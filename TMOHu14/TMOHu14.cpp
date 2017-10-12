/* --------------------------------------------------------------------------- *
 * TMOAncuti16.cpp: implementation of the TMOHu14 class.   *
 * Diploma thesis
 * Author : Vladimir Vlkovic, Brno 2017
 * --------------------------------------------------------------------------- */

#include "TMOHu14.h"
#include <boost/concept_check.hpp>
#include <complex>



/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOHu14::TMOHu14()
{
	SetName(L"Hu14");						
	SetDescription(L"Image decolorization using Laplacian operator and multi-scale fusion");

	
}

TMOHu14::~TMOHu14()
{
}
/**
 * returns edge mat using Canny edge detector
 * @param color channel mat
 * @return edge mat of input channel
 */

cv::Mat TMOHu14::getEdgeMat(cv::Mat channel)
{
  double min, max;
  cv::Mat tmp,result;
  int threshold1 = 55;
  int threshold2 = 40;
  cv::minMaxLoc(channel,&min,&max);
	if(min!=max)
	{
	  channel.convertTo(tmp,CV_8U,255*(max-min));
	}
	cv::Canny(tmp,result,threshold1,threshold2);
	
	return result;
}


int TMOHu14::Transform()
{

	
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			 
													
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	
	double min, max;
	cv::Mat redMat, greenMat, blueMat, redEdgeMat, greenEdgeMat, blueEdgeMat, sumEdgeMat, tmpMat;      //////Mat for each color channel
	
	
	
	redMat = cv::Mat::zeros(height, width, CV_32F);
	greenMat = cv::Mat::zeros (height, width, CV_32F);  ////mats for color channels
	blueMat = cv::Mat::zeros (height, width, CV_32F);	
	
	redEdgeMat = cv::Mat::zeros(height, width, CV_8U);
	greenEdgeMat = cv::Mat::zeros (height, width, CV_8U);  ///mats for edges of color channels must be CV_8U because Canny function
	blueEdgeMat = cv::Mat::zeros (height, width, CV_8U);
												
	
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	
 		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			redMat.at<float>(j,i) = *pSourceData++; 
			greenMat.at<float>(j,i) = *pSourceData++;  //getting separate RGB channels
			blueMat.at<float>(j,i) = *pSourceData++;

		}
	}
	
	
	
	redEdgeMat = TMOHu14::getEdgeMat(redMat);
	greenEdgeMat = TMOHu14::getEdgeMat(greenMat); //getting edge mats
	blueEdgeMat = TMOHu14::getEdgeMat(blueMat);
	

	sumEdgeMat = redEdgeMat + greenEdgeMat + blueEdgeMat; ///need to sum all channel edges
	sumEdgeMat.convertTo(tmpMat,CV_32F); ///conversion to float mat
	
	
	redMat=redMat.mul(tmpMat)/255;
	greenMat=greenMat.mul(tmpMat)/255; ///multipling edge map by color channel maps to achieve I_edge (see alg.)
	blueMat=blueMat.mul(tmpMat)/255;
	
	

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
	    pSrc->ProgressBar(j, pSrc->GetHeight());
		
	    for (int i = 0; i < pSrc->GetWidth(); i++) ///result to output, taking only the image correction is discarded
	    {
		  *pDestinationData++ =redMat.at<float>(j,i);
		 *pDestinationData++ = greenMat.at<float>(j,i);
		 *pDestinationData++ =blueMat.at<float>(j,i);
	    }
	}
	  
	
	return 0;
}


