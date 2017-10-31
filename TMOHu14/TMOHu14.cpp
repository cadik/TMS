/* --------------------------------------------------------------------------- *
 * TMOHu14.cpp: implementation of the TMOHu14 class.   *
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
  double x, y, z, L, u, v, r,g,b;
  cv::Vec3d LuvVector;
  r = bgrVector[2] / 255.0f;
  g = bgrVector[1] / 255.0f;
  b = bgrVector[0] / 255.0f;
  
  if (r > 0.04045) r = std::pow((r + 0.055)/1.055, 2.4);
  else r = r /12.92;
  
  if (g > 0.04045) g = std::pow((g + 0.055)/1.055, 2.4);
  else g = g /12.92;
  
  if (b > 0.04045) b = std::pow((b + 0.055)/1.055, 2.4);
  else b = b /12.92;
  
  r = r *100;
  g= g*100;
  b=b*100;
  
  x = r * 0.4124 + g * 0.3576 + b * 0.1805;
  y = r * 0.2126 + g * 0.7152 + b * 0.0722; 
  z = r * 0.0193 + g * 0.1192 + b * 0.9505;
  
  
  
  TMOImage::XyzToLuv(x, y, z, &L, &u, &v);  
  
  LuvVector[0] = L;
  LuvVector[1] = u;
  LuvVector[2] = v;
  
  return LuvVector;
   
}


/**
 * Get feature vector consisting of Luv color and color percentage in image
 * @param src I_edge Mat
 * @return feature vector/color pallete/color histogram : color in Luv, color percentage in image
 */
 void TMOHu14::getPalette(std::map<cv::Vec3d, float, lessVec3b>& paletteLuv, cv::Mat& src)
{
   std::map<cv::Vec3b, float, lessVec3b> paletteRGB;
  // std::map<cv::Vec3d, int, lessVec3b> paletteLuv;
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
        int d1e=paletteRGB.size();
	std::map<cv::Vec3b, float, lessVec3b>::iterator it=paletteRGB.begin();
    while( it!=paletteRGB.end())
    {
      float colorPercentage=0.0;
      colorPercentage=(it->second / pixelCount)*100.0f;    
      if(colorPercentage < 0.1f)     /// if color percentage is less then discard
      {
	paletteRGB.erase(it++);
      }
     else
      {
	cv::Vec3d tmpBgr = TMOHu14::rgb2Luv(it->first);  ///convert bgr to Luv
	paletteLuv[tmpBgr] = colorPercentage;
	//paletteRGB[it->first]=colorPercentage;
	//it->second = colorPercentage;
	++it;
      }
      
    }
 
}


/**
 * Get get dominant color descriptr in LUV space
 * merging of perceptual similar colors
 * @param palette Luv color pallete with pixel percentage
 * @return dominat color feature vector : color in Luv, color percentage in image
 */
std::map<cv::Vec3d, float, lessVec3b> TMOHu14::getDominantColorDescriptor(std::map<cv::Vec3d, float, lessVec3b> palette)
{
  int delta = 10; //threshold, empiric value see alg.
  int d=0;
  cv::Vec3d newColor;
  std::map<cv::Vec3d, float, lessVec3b>::iterator it=palette.begin();
  std::map<cv::Vec3d, float, lessVec3b>::iterator it2=palette.begin();
  it++;
  int i=0;
  while( it!=palette.end())
  {
    i++;
    double L,u,v;
    
    L=it->first[0];
    u=it->first[1];
    v=it->first[2];
    
    while ( it2!=palette.end())
    {
      double L2,u2,v2;
      L2=it2->first[0];
      u2=it2->first[1];
      v2=it2->first[2];
      //d=std::abs(r-r2) + std::abs(g-g2) + std::abs(b-b2);
       d = std::sqrt(std::pow(L-L2,2) + std::pow(u-u2,2) + std::pow(v-v2,2));
      
     
	if(d < delta && d > 0)
	{
	  newColor[0] = (L * it->second + L2 * it2->second) / (it->second + it2->second);
	  newColor[1] = (u * it->second + u2 * it2->second) / (it->second + it2->second);
	  newColor[2] = (v * it->second + v2 * it2->second) / (it->second + it2->second);
	  
	  
	  palette[newColor] = it->second + it2->second;
	  palette.erase(it++);
	  palette.erase(it2++);
	  break;
	}
	else ++it2;
      
    
    }
    ++it;
    
    
  }
  
  
  int e = palette.size();
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
	
	
	
	std::map<cv::Vec3d, float, lessVec3b> palette ;
	TMOHu14::getPalette(palette,reduced); //gettign the color palette in rgb
	
	palette =TMOHu14::getDominantColorDescriptor(palette);


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


