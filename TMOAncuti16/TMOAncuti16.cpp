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
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	
	int correctionWidth, correctionHeight;
	
	cv::Mat meanKernel = cv::Mat::ones(3,3,CV_64FC1);
	cv::Mat red, green, blue;      //////Mat for each color channel
	cv::Mat redLap,greenLap,blueLap; ///mat fo each laplacian, needed in weight map computation
	cv::Mat redLapWeightMap,greenLapWeightMap,blueLapWeightMap; //laplacian weight map for each channel
	cv::Mat redGlobalWeightMap, greenGlobalWeightMap,blueGlobalWeightMap; //gllobal weight map for each channel
	cv::Mat redNormalisedWeightMap, greenNormalisedWeightMap, blueNormalisedWeightMap; //normalised weigth maps
	cv::Mat layer; //used to store the result for fusion on a certain pyramid layer
	cv::Mat maxMat; //stores the cobined sum of all matrices per element, used for normalisation
	cv::Mat meanMat; //per element lapcian means
	cv::Mat endResult; //result
	cv::Mat tmp, tmp2,tmp3, tmp4; ///temporary variables
	
	correctionWidth = std::ceil(log2(width));    ////picture dimensions must be 2^n, 
	correctionHeight = std::ceil(log2(height));
	correctionHeight = std::pow(2 , correctionHeight);
	correctionWidth = std::pow(2 , correctionWidth);
	////////////setting up matrices for correct dimensions////////////////////////////////////////////
	red = cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	green = cv::Mat::zeros (correctionHeight, correctionWidth, CV_64FC1);
	blue = cv::Mat::zeros (correctionHeight, correctionWidth, CV_64FC1);
	
	redLap=cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	greenLap=cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	blueLap=cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	
	
	redNormalisedWeightMap = cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	greenNormalisedWeightMap = cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	blueNormalisedWeightMap = cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	
	endResult = cv::Mat::zeros(correctionHeight, correctionWidth, CV_64FC1);
	//////////////////////////////////////////////////////////////////////////////////////////////
        
	cv::Mat averagingMeanKernel = cv::Mat::ones(correctionHeight,correctionWidth,CV_64FC1)*9; //used for dividing the sum 
												//sum/9
	
	
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	
 		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			red.at<double>(j,i) = *pSourceData++; 
			green.at<double>(j,i) = *pSourceData++;  //getting separate RGB channels
			blue.at<double>(j,i) = *pSourceData++;

		}
	}
	///computation of laplacian for each chanel lap = img - blurred img////////////////
	cv::GaussianBlur(red,tmp,cv::Size(5,5),0,0);
	cv::GaussianBlur(green,tmp2,cv::Size(5,5),0,0);
	cv::GaussianBlur(blue,tmp3,cv::Size(5,5),0,0);
	
	redLap = redLap - tmp;
	greenLap = greenLap -tmp2;
	blueLap = blueLap -tmp3;
	///////////////////////////////////////////////////////////////////////////////////////
	
	
	cv::filter2D(redLap,redLap,-1,meanKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);   ////getting the sum of neigbor values
	cv::filter2D(greenLap,greenLap,-1,meanKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);
	cv::filter2D(blueLap,blueLap,-1,meanKernel, cv::Point( -1, -1 ), 0, cv::BORDER_DEFAULT);
	   
	cv::divide(redLap, averagingMeanKernel, meanMat, 1, -1); ///averaging the sum
	redLapWeightMap = meanMat + cv::abs(redLap);    ///laplacian weight map computation
	cv::pow(red - meanMat, 2, redGlobalWeightMap);  // global map computation

	cv::divide(greenLap, averagingMeanKernel, meanMat, 1, -1);   ///see above
	greenLapWeightMap = meanMat + cv::abs(greenLap);
	cv::pow(green - meanMat, 2, greenGlobalWeightMap);

	cv::divide(blueLap, averagingMeanKernel, meanMat, 1, -1); ///see above
	blueLapWeightMap = meanMat + cv::abs(blueLap);
	cv::pow(blue - meanMat, 2, blueGlobalWeightMap);

	redNormalisedWeightMap = redLapWeightMap + redGlobalWeightMap;   ///creating normalised weight maps
	greenNormalisedWeightMap = greenLapWeightMap + greenGlobalWeightMap;
	blueNormalisedWeightMap = blueLapWeightMap + blueGlobalWeightMap;

	maxMat = redNormalisedWeightMap + greenNormalisedWeightMap + blueNormalisedWeightMap;

	cv::divide(redNormalisedWeightMap, maxMat, redNormalisedWeightMap, 1, -1);   //normalising normalised weight maps
	cv::divide(greenNormalisedWeightMap, maxMat, greenNormalisedWeightMap, 1, -1);
	cv::divide(blueNormalisedWeightMap, maxMat, blueNormalisedWeightMap, 1, -1);

	layer=redNormalisedWeightMap.mul(red) +greenNormalisedWeightMap.mul(green) +blueNormalisedWeightMap.mul(blue);    
	endResult = layer;  ////computing the 0 level layer and adding it to the resullt
	  
	for(int i=1; i<=std::floor(log10(width*height));i++) //number of levels is the log of the total image size
	{
	    cv::pyrDown(redNormalisedWeightMap, redNormalisedWeightMap, cv::Size(redNormalisedWeightMap.cols/2,redNormalisedWeightMap.rows/2));
	    cv::pyrDown(red, red, cv::Size(red.cols/2,red.rows/2));
	    cv::GaussianBlur(red,tmp,cv::Size(5,5),0,0); 
	    //first line gaussian pyramid layer of norm weight maps, second and third preparation for lapalcian layer
	  
	    cv::pyrDown(greenNormalisedWeightMap, greenNormalisedWeightMap, cv::Size(greenNormalisedWeightMap.cols/2,greenNormalisedWeightMap.rows/2));
	    cv::pyrDown(green, green, cv::Size(green.cols/2,green.rows/2));
	    cv::GaussianBlur(green,tmp2,cv::Size(5,5),0,0);
	  
	    cv::pyrDown(blueNormalisedWeightMap, blueNormalisedWeightMap, cv::Size(blueNormalisedWeightMap.cols/2,blueNormalisedWeightMap.rows/2));
	    cv::pyrDown(blue, blue, cv::Size(blue.cols/2,blue.rows/2));
	    cv::GaussianBlur(blue,tmp3,cv::Size(5,5),0,0);

	    layer=redNormalisedWeightMap.mul(red- tmp)+ greenNormalisedWeightMap.mul(green - tmp2) + greenNormalisedWeightMap.mul(blue - tmp3);
	    ///gaussian level of weight map multiplied by laplacian level and then summed over 3 channels 
	
	    for(int j= 0; j<i;j++) ///how many times to upsclae current level
	    {
		cv::pyrUp(layer,layer,cv::Size(layer.cols*2,layer.rows*2));
		//upsacling for correct addition to result
	    }
	    endResult = endResult + layer; //addition of upsacel level to result
	}

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
	    pSrc->ProgressBar(j, pSrc->GetHeight());
		
	    for (int i = 0; i < pSrc->GetWidth(); i++) ///result to output, taking only the image correction is discarded
	    {
		  *pDestinationData++ =endResult.at<double>(j,i);
		  *pDestinationData++ = endResult.at<double>(j,i);
		  *pDestinationData++ =endResult.at<double>(j,i);
	    }
	}
	  
	
	return 0;
}


