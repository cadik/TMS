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

/**
 * quantizes the colors of the I_edge 
 * @param src I_edge Mat
 * @param dst quantized I_edge Mat
 */
void TMOHu14::kmeansColorQuantization(const cv::Mat3b& src, cv::Mat3b& dst)
{
    int K = 256;  ///should be thi but takes too long
    int n = src.rows * src.cols;
    cv::Mat data = src.reshape(1, n);
    data.convertTo(data, CV_32F);

    std::vector<int> labels;
    cv::Mat1f colors;
    cv::kmeans(data, K, labels, cv::TermCriteria(), 1, cv::KMEANS_PP_CENTERS, colors);

    for (int i = 0; i < n; ++i)
    {
        data.at<float>(i, 0) = colors(labels[i], 0);
        data.at<float>(i, 1) = colors(labels[i], 1);
        data.at<float>(i, 2) = colors(labels[i], 2);
    }

    cv::Mat reduced = data.reshape(3, src.rows);
    reduced.convertTo(dst, CV_8U);
    
    cv::imshow("Reduced", dst);
   cv::waitKey();
    
  
  
}
/**
 * Converts BGR vector to XYZ and then to Luv
 * @param bgrVector BGR color vector
 * @return Luv color vector 
 */
cv::Vec3d TMOHu14::rgb2Luv(cv::Vec3b bgrVector)
{
  double x, y, z, L, u, v;
  cv::Vec3d LuvVector;
  
  x = bgrVector[2] * 0.4124 + bgrVector[1] * 0.3576 + bgrVector[0] * 0.1805;
  y = bgrVector[2] * 0.2126 + bgrVector[1] * 0.7152 + bgrVector[0] * 0.0722; //BGR->XYZ converion
  z = bgrVector[2] * 0.0193 + bgrVector[1] * 0.1192 + bgrVector[0] * 0.9505;
  
  
  TMOImage::XyzToLuv(x, y, z, &L, &u, &v);  
  
  LuvVector[0] = L;
  LuvVector[1] = u;
  LuvVector[2] = v;
  
  return LuvVector;
   
}


/**
 * Get feature vector consisting of Luv color and color percentage in image
 * @param src I_edge Mat
 * @return feature vector : color in Luv, color percentage in image
 */
std::map<cv::Vec3d, int, lessVec3b> TMOHu14::getPalette(const cv::Mat& src)
{
    std::map<cv::Vec3b, int, lessVec3b> paletteRGB;
    std::map<cv::Vec3d, int, lessVec3b> paletteLuv;
    float pixelCount=src.rows*src.cols;
    for (int r = 0; r < src.rows; ++r)
    {
        for (int c = 0; c < src.cols; ++c)   ///get every color and pixel count of every color
        {

            cv::Vec3b color = src.at<cv::Vec3b>(r,c);
	    
           if (paletteRGB.count(color) == 0)
            {
                paletteRGB[color] = 1;
            }
            else
            {
                paletteRGB[color] = paletteRGB[color] + 1;
            }
        }
    }
       
    for (std::map<cv::Vec3b, int, lessVec3b>::iterator it=paletteRGB.begin(); it!=paletteRGB.end(); ++it)
    {
      float colorPercentage=0.0;
      colorPercentage=(it->second / pixelCount)*100.0;    
      if(colorPercentage <= 0.1)     /// if color percentage is less then discard
      {
	paletteRGB.erase(it);
      }
      else
      {
	cv::Vec3b tmpBgr = TMOHu14::rgb2Luv(it->first);  ///convert bgr to Luv
	paletteLuv[tmpBgr] = colorPercentage;
      }
      
    }
    
    return paletteLuv;
}



int TMOHu14::Transform()
{

	
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			 
													
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	
	double min, max;
	cv::Mat redMat, greenMat, blueMat, redEdgeMat, greenEdgeMat, blueEdgeMat, sumEdgeMat, tmpMat;      //////Mat for each color channel
	cv::Mat3b reduced;
	
	
	
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
	
	
	redMat=redMat.mul(tmpMat);
	greenMat=greenMat.mul(tmpMat); ///multipling edge map by color channel maps to achieve I_edge (see alg.)
	blueMat=blueMat.mul(tmpMat);
	
	cv::Mat mergedMat;
	std::vector<cv::Mat> channels;
	channels.push_back(blueMat);
	channels.push_back(greenMat);
	channels.push_back(redMat); /// mergig channels into one mat
	cv::merge(channels,mergedMat);

	TMOHu14::kmeansColorQuantization(mergedMat,reduced);  ////decrease number of colors uing color quantization
	
	
	
	std::map<cv::Vec3d, int, lessVec3b> palette = TMOHu14::getPalette(reduced); //gettign the color palette in Luv


	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
	    pSrc->ProgressBar(j, pSrc->GetHeight());
		
	    for (int i = 0; i < pSrc->GetWidth(); i++) ///result to output, taking only the image correction is discarded
	    {
	     
		  *pDestinationData++ =redMat.at<float>(j,i);
		 *pDestinationData++ =greenMat.at<float>(j,i);
		 *pDestinationData++ =blueMat.at<float>(j,i);
	    }
	}
	  
	
	return 0;
}


