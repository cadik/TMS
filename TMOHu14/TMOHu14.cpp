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
	SetDescription(L"Image and video decolorization in real-time using dominant colors");

	
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
  int threshold1 = 70;
  int threshold2 = 50;
  cv::minMaxLoc(channel,&min,&max);
	if(min!=max)
	{
	  channel.convertTo(tmp,CV_8U,255*(max-min));
	}
	cv::Canny(tmp,result,threshold1,threshold2);
	
	return result;
}


std::map<cv::Vec3f, int, lessVec3b> TMOHu14::getPalette(const cv::Mat& src)
{
    std::map<cv::Vec3f, int, lessVec3b> palette;
    for (int r = 0; r < src.rows; ++r)
    {
        for (int c = 0; c < src.cols; ++c)
        {
	 
	  
	  
            cv::Vec3f color = src.at<cv::Vec3f>(r,c);
	    
           if (palette.count(color) == 0)
            {
                palette[color] = 1;
            }
            else
            {
                palette[color] = palette[color] + 1;
            }
        }
    }
    int pixelCount=src.rows*src.cols;
   /* for (std::map<cv::Vec3f, int, lessVec3b>::iterator it=palette.begin(); it!=palette.end(); ++it)
    {
      float f=(it->second / pixelCount)*100.0;
      if(f <= 0.1)
      {
	palette.erase(it);
      }
      
    }*/
    
    return palette;
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
	
	cv::Mat mergedMat;
	std::vector<cv::Mat> channels;
	channels.push_back(blueMat*255);
	channels.push_back(greenMat*255);
	channels.push_back(redMat*255); /// mergig channels into one mat
	cv::merge(channels,mergedMat);
	

	
	redEdgeMat = TMOHu14::getEdgeMat(redMat);
	greenEdgeMat = TMOHu14::getEdgeMat(greenMat); //getting edge mats
	blueEdgeMat = TMOHu14::getEdgeMat(blueMat);
	

	sumEdgeMat = redEdgeMat + greenEdgeMat + blueEdgeMat; ///need to sum all channel edges
	sumEdgeMat.convertTo(tmpMat,CV_32F); ///conversion to float mat
	
	
	redMat=redMat.mul(tmpMat)/255;
	greenMat=greenMat.mul(tmpMat)/255; ///multipling edge map by color channel maps to achieve I_edge (see alg.)
	blueMat=blueMat.mul(tmpMat)/255;
	
	std::map<cv::Vec3f, int, lessVec3b> palette = getPalette(mergedMat); //gettign the color palette todo: color quantization
	

	
	

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


