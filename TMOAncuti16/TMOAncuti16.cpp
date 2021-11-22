/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                               *
*                                                                              *
*                       Semestral project                                      *
*                       Author: Vladimir Vlkovic                               *
*                       Brno 2017                                              *
*                                                                              *
*                       Implementation of the TMOAncuti16 class                *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOAncuti16.cpp
 * @brief Implementation of the TMOAncuti16 class
 * @author Vladimir Vlkovic
 * @class TMOAncuti16
 */

#include "TMOAncuti16.h"
#include <boost/concept_check.hpp> 
#include <complex>



/**
  *  @brief Constructor
  */
TMOAncuti16::TMOAncuti16()
{
	SetName(L"Ancuti16");						
	//SetDescription(L"Image decolorization using Laplacian operator and multi-scale fusion");
        SetDescription(L"Laplacian-guided image decolorization");
	
}

/**
  *  @brief Destructor
  */
TMOAncuti16::~TMOAncuti16()
{
}

/**
  *  @brief Laplacian-guided image decolorization
  */
int TMOAncuti16::Transform()
{

	
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			 
													
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	
	int correctionWidth, correctionHeight;
	double kdata[]={0,-1,0,-1,4,-1,0,-1,0};
	cv::Mat meanKernel = cv::Mat::ones(3,3,CV_64FC1); /** Used for concolution , to compute the sum of surrounding pixels */
	cv::Mat lapKernel(3,3,CV_64FC1,kdata);/** Laplacian kernel then flipped cause opencv performs a correlation not concolution */
	//cv::flip(lapKernel,lapKernel,-1);
	cv::Mat red, green, blue;      /** Mat for each color channel */
	cv::Mat redLap,greenLap,blueLap; /** Mat fo each laplacian, needed in weight map computation */
	cv::Mat redLapWeightMap,greenLapWeightMap,blueLapWeightMap; /** Laplacian weight map for each channel */
	cv::Mat redGlobalWeightMap, greenGlobalWeightMap,blueGlobalWeightMap; /** Global weight map for each channel */
	cv::Mat redNormalisedWeightMap, greenNormalisedWeightMap, blueNormalisedWeightMap; /** Normalised weigth maps */
	cv::Mat layer; /** Used to store the result for fusion on a certain pyramid layer */
	cv::Mat maxMat; /** Stores the cobined sum of all matrices per element, used for normalisation */
	cv::Mat meanMat; /** Per element lapcian means */
	cv::Mat endResult; /** result */
	cv::Mat tmp, tmp2,tmp3, tmp4; /** Temporary variables */
	
	correctionWidth = std::ceil(log2(width));    /** picture dimensions must be 2^n, */
	correctionHeight = std::ceil(log2(height));
	correctionHeight = std::pow(2 , correctionHeight);
	correctionWidth = std::pow(2 , correctionWidth);
	////////////setting up matrices ////////////////////////////////////////////
	red = cv::Mat::zeros(height, width, CV_64FC1);
	green = cv::Mat::zeros (height, width, CV_64FC1);
	blue = cv::Mat::zeros (height, width, CV_64FC1);
	redLap = cv::Mat::zeros(height, width, CV_64FC1);
	greenLap = cv::Mat::zeros (height, width, CV_64FC1);
	blueLap = cv::Mat::zeros (height, width, CV_64FC1);
	

	endResult = cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1); /** dims must be 2^n for pyramid functions*/
 	//////////////////////////////////////////////////////////////////////////////////////////////
        
	cv::Mat averagingMeanKernel = cv::Mat(height,width,CV_64FC1,cv::Scalar(9)); //used for dividing the sum 
												//sum/9
	
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	
 		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			red.at<double>(j,i) = *pSourceData++; 
			green.at<double>(j,i) = *pSourceData++;  /**getting separate RGB channels */
			blue.at<double>(j,i) = *pSourceData++;

		}
	}
	
	
	cv::filter2D(red,redLap,-1,lapKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);
	cv::filter2D(green,greenLap,-1,lapKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);  /** Laplacian computation */
	cv::filter2D(blue,blueLap,-1,lapKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	
	cv::filter2D(redLap,tmp,-1,meanKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);   /** Getting the sum of neigbor values */
	cv::filter2D(greenLap,tmp2,-1,meanKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);
	cv::filter2D(blueLap,tmp3,-1,meanKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);
	   
	cv::divide(tmp, averagingMeanKernel, meanMat, 1, -1); /** Averaging the sum */
	
	redLapWeightMap = meanMat + cv::abs(redLap);    /** Laplacian weight map computation */
	cv::pow(red - meanMat, 2, redGlobalWeightMap);  /** Global map computation */

	cv::divide(tmp2, averagingMeanKernel, meanMat, 1, -1);   /** See above */
	greenLapWeightMap = meanMat + cv::abs(greenLap);
	cv::pow(green - meanMat, 2, greenGlobalWeightMap);

	cv::divide(tmp3, averagingMeanKernel, meanMat, 1, -1); /** See above */
	blueLapWeightMap = meanMat + cv::abs(blueLap);
	cv::pow(blue - meanMat, 2, blueGlobalWeightMap);
	
	redLap.release();
	greenLap.release();
	blueLap.release();
	meanMat.release();
	
	
	
	cv::normalize(redLapWeightMap, redLapWeightMap, 0.0, 1.0,cv::NORM_MINMAX,CV_64F);   /** Normaling weight maps */
	cv::normalize(greenLapWeightMap, greenLapWeightMap, 0.0, 1.0,cv::NORM_MINMAX,CV_64F);
	cv::normalize(blueLapWeightMap, blueLapWeightMap, 0.0, 1.0,cv::NORM_MINMAX,CV_64F);
	
	
	
	maxMat = redLapWeightMap + redGlobalWeightMap +  greenLapWeightMap + greenGlobalWeightMap + blueLapWeightMap + blueGlobalWeightMap;
	
	cv::divide(redLapWeightMap, maxMat, redLapWeightMap, 1, CV_64FC1);   /** Normalizing in such way taht the sum of the maps is 1 */
	cv::divide(redGlobalWeightMap, maxMat, redGlobalWeightMap, 1, CV_64FC1);
	cv::divide(greenLapWeightMap, maxMat, greenLapWeightMap, 1, CV_64FC1);
	cv::divide(greenGlobalWeightMap, maxMat, greenGlobalWeightMap, 1, CV_64FC1);   /** Normalizing in such way taht the sum of the maps is 1 */
	cv::divide(blueLapWeightMap, maxMat, blueLapWeightMap, 1, CV_64FC1);
	cv::divide(blueGlobalWeightMap, maxMat, blueGlobalWeightMap, 1, CV_64FC1);

	
	
	redNormalisedWeightMap = redLapWeightMap + redGlobalWeightMap;   /** Creating normalised weight maps */
	greenNormalisedWeightMap = greenLapWeightMap + greenGlobalWeightMap;
	blueNormalisedWeightMap = blueLapWeightMap + blueGlobalWeightMap;	
	
	
	
	redGlobalWeightMap.release();
	greenGlobalWeightMap.release();
	blueGlobalWeightMap.release();
	
	redLapWeightMap.release();
	greenLapWeightMap.release();
	blueLapWeightMap.release();

	
	
	
	
	
	maxMat.release();

	layer=redNormalisedWeightMap.mul(red) +greenNormalisedWeightMap.mul(green) +blueNormalisedWeightMap.mul(blue);   /** Creating level 0 of pyramid */
	
	cv::copyMakeBorder(layer,endResult,0,correctionHeight-height,0,correctionWidth-width,cv::BORDER_DEFAULT); /** Resizing matirx to 2^n dimensions  bordes are interpolated */

	for(int i=1; i<=std::floor(log10(width*height));i++) /** Number of levels is the log of the total image size */
	{
	 
	  
	    cv::pyrDown(redNormalisedWeightMap, redNormalisedWeightMap, cv::Size(redNormalisedWeightMap.cols/2,redNormalisedWeightMap.rows/2));
	    cv::pyrDown(red, red, cv::Size(red.cols/2,red.rows/2));
	    cv::GaussianBlur(red,tmp,cv::Size(5,5),0,0); 
	    /** First line gaussian pyramid layer of norm weight maps, second and third preparation for lapalcian layer */
	  
	    cv::pyrDown(greenNormalisedWeightMap, greenNormalisedWeightMap, cv::Size(greenNormalisedWeightMap.cols/2,greenNormalisedWeightMap.rows/2));
	    cv::pyrDown(green, green, cv::Size(green.cols/2,green.rows/2));
	   cv::GaussianBlur(green,tmp2,cv::Size(5,5),0,0);
	  
	  
	    cv::pyrDown(blueNormalisedWeightMap, blueNormalisedWeightMap, cv::Size(blueNormalisedWeightMap.cols/2,blueNormalisedWeightMap.rows/2));
	    cv::pyrDown(blue, blue, cv::Size(blue.cols/2,blue.rows/2));
	   cv::GaussianBlur(blue,tmp3,cv::Size(5,5),0,0);
	    

	    layer=redNormalisedWeightMap.mul(red-tmp)+ greenNormalisedWeightMap.mul(green-tmp2) + greenNormalisedWeightMap.mul(blue-tmp3);
	    
	    
	    /** Gaussian level of weight map multiplied by laplacian level and then summed over 3 channels */
	
	    for(int j= 0; j<i;j++) /** How many times to upsclae current level */
	    {
		cv::pyrUp(layer,layer,cv::Size(layer.cols*2,layer.rows*2));
		/** Upsacling for correct addition to result */
	    }
	    cv::copyMakeBorder(layer,tmp4,0,correctionHeight-layer.rows,0,correctionWidth-layer.cols,cv::BORDER_DEFAULT);/** Enlargement to 2^n for summing, bordes are interpolated */
	    endResult = endResult + tmp4; /** Addition of upsaceled level to result */
	    
	    
	}

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
	    pSrc->ProgressBar(j, pSrc->GetHeight());
		
	    for (int i = 0; i < pSrc->GetWidth(); i++) /** Result to output, taking only the image correction is discarded */
	    {
		  *pDestinationData++ =endResult.at<double>(j,i);
		  *pDestinationData++ = endResult.at<double>(j,i);
		  *pDestinationData++ =endResult.at<double>(j,i);
	    }
	}
	  
	
	return 0;
}
