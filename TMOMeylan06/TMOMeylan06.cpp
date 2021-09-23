/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio	                                        *
*                                                                                   *
*                       Author: Jan Kohut [xkohut08 AT stud.fit.vutbr.cz]           *
*                       Brno 2019                                                   *
*                                                                                   *
*                       Implementation of the TMOMeylan06 class                     *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOMeylan06.cpp
 * @brief Implementation of the TMOMeylan06 class
 * @author Jan Kohut
 * @class TMOMeylan06.cpp
 */ 

#include "TMOMeylan06.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMeylan06::TMOMeylan06()
{
	SetName(L"Meylan06");
	SetDescription(L"High Dynamic Range Image Rendering Using a Retinex-Based Adaptive Filter");

	this->kernelRadiusParameter.SetName(L"Kernel Radius");
	this->kernelRadiusParameter.SetDescription(L"Kernel radius of gaussian filter.");
	this->kernelRadiusParameter.SetDefault(49);
	this->kernelRadiusParameter=49;
	this->kernelRadiusParameter.SetRange(3,100);
	this->Register(this->kernelRadiusParameter);

	this->sigmaOrigParameter.SetName(L"Sigma Original");
	this->sigmaOrigParameter.SetDescription(L"Sigma for gaussian filter before edge is detected.");
	this->sigmaOrigParameter.SetDefault(16);
	this->sigmaOrigParameter=16;
	this->sigmaOrigParameter.SetRange(0,20);
	this->Register(this->sigmaOrigParameter);

	this->sigmaEdgeParameter.SetName(L"Sigma Edge");
	this->sigmaEdgeParameter.SetDescription(L"Sigma for gaussian filter after edge is detected.");
	this->sigmaEdgeParameter.SetDefault(5);
	this->sigmaEdgeParameter=5;
	this->sigmaEdgeParameter.SetRange(0,20);
	this->Register(this->sigmaEdgeParameter);

	this->saturationParameter.SetName(L"Saturation");
	this->saturationParameter.SetDescription(L"Color enhancement of two last PCA components.");
	this->saturationParameter.SetDefault(1.6254);
	this->saturationParameter=1.6254;
	this->saturationParameter.SetRange(1.0,10.0);
	this->Register(this->saturationParameter);

	// std::cout << std::endl;
}

TMOMeylan06::~TMOMeylan06()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMeylan06::Transform()
{

	this->numberOfPixels = pSrc->GetHeight() * pSrc->GetWidth();
	this->numberOfPixelsRGB = this->numberOfPixels * 3;

	cv::Mat source(this->numberOfPixels, 3, CV_64F, pSrc->GetData());
	this->Normalize(source);

	// LUMINANCE PROCESSING
	// ***************************************************************************

	// Get luminance as projection to first principal component of PCA analysis
	cv::Mat PCAProjectionForLuminance = this->RGBToPCA(source);
	cv::Mat luminance = this->GetLuminance(PCAProjectionForLuminance);

	this->Normalize(luminance, 0.0, 1.0);

	// Compute global tone mapping of luminance
	std::string type = "exp";
	this->GlobalToneMap(luminance, type);

	// Downsample luminance so the computation is faster
	int maskMaxSize = 200;
	cv::Mat luminanceLowRes = this->ResizeLuminance(luminance, maskMaxSize);

	// Use Canny edge detector to detect edges
	// Lower threshold is set to 0.4 * upperThreshold (MATLAB)
	// Dilatate the edge image so there is always edge in diagonal direction
	double upperThreshold = 0.2;
	cv::Mat edgesLowRes = this->GetEdges(luminanceLowRes, upperThreshold);
	edgesLowRes = this->DilatateEdges(edgesLowRes);

	// Compute mask for local tone mapping
	this->sigmaOrig = this->sigmaOrigParameter;
	this->sigmaEdge = this->sigmaEdgeParameter;
	//this->kernelRadius = (int) 3 * this->sigmaOrig + 1;
	this->kernelRadius = this->kernelRadiusParameter;
	cv::Mat maskLowRes = this->GetMask(luminanceLowRes, edgesLowRes);

	// Upsample mask so it can be substracted from original luminance
	cv::Mat mask(this->iHeight, this->iWidth, CV_64F);
	cv::resize(maskLowRes, mask, cv::Size(this->iWidth, this->iHeight));

	// this->SaveImg("lum.png", luminance);
	// this->SaveImg("mask.png", mask);

	double max = 0.1;
	double scale = 100;
	this->LogMaxScale(luminance, max, scale);
	this->LogMaxScale(mask, max, scale);
	// this->SaveImg("lum_log.png", luminance);
	// this->SaveImg("mask_log.png", mask);

	double sigmaA = 10;
	double sigmaC = 0.5;
	cv::Mat betaFactor = this->GetBetaFactor(luminance, sigmaA, sigmaC);
	mask = this->ElementWiseMul(mask, betaFactor);
	// this->SaveImg("mask_beta.png", mask);

	// Get final luminance by substracting Beta * mask from it
	luminance = this->ElementWiseSub(luminance, mask);
	int numberOfBuckets = 100;
	double minThreshold = 0.01;
	double maxThreshold = 0.99;
	this->HistoClip(luminance, numberOfBuckets, minThreshold, maxThreshold);

	// this->SaveImg("lum_final.png", luminance);

	// std::cout << "LUMINANCE PROCESSING DONE" << std::endl;
	// ***************************************************************************


	// CHROMINANCE PROCESSING
	// ***************************************************************************
	cv::Mat rgb = cv::Mat(source);
	rgb = source.clone();

	// Compute global tone mapping of chrominance
	type = "exp";
	this->GlobalToneMap(rgb, type);

	max = 1;
	scale = 100;
	this->LogMaxScale(rgb, max, scale);

	// Enhance the projection to the second and third principal component representing chrominance
	cv::Mat PCAProjectionForRGB = this->RGBToPCA(rgb);
	double saturationEnhancement = this->saturationParameter;
	this->ScaleSaturation(PCAProjectionForRGB, saturationEnhancement);

	// Normalize luminance so it fits the current PCA analysis
	cv::Mat luminanceRGB = this->GetLuminance(PCAProjectionForRGB);
	double luminanceRGBMin, luminanceRGBMax;
	cv::minMaxLoc(luminanceRGB, &luminanceRGBMin, &luminanceRGBMax);
	this->Normalize(luminance, luminanceRGBMin, luminanceRGBMax);

	// Recompose the projection with processed luminance
	double *PCAProjectionForRGBPtr = PCAProjectionForRGB.ptr<double>(0);
	double* luminancePtr = luminance.ptr<double>(0);
	for (int i = 0; i < this->numberOfPixels; ++i)
	{
		*PCAProjectionForRGBPtr = *luminancePtr;
		PCAProjectionForRGBPtr += 3;
		luminancePtr++;
	}

	// Get the final image as back projection
	cv::Mat finalRGB = this->PCAToRGB(PCAProjectionForRGB);

	numberOfBuckets = 100;
	minThreshold = 0.01;
	maxThreshold = 0.99;
	this->HistoClip(finalRGB, numberOfBuckets, minThreshold, maxThreshold);

	// std::cout << "CHROMINANCE PROCESSING DONE" << std::endl;
	// ***************************************************************************

	double* pDestinationData = pDst->GetData();
	double* finalRGBPtr = finalRGB.ptr<double>(0);
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			*pDestinationData++ = *finalRGBPtr++;
			*pDestinationData++ = *finalRGBPtr++;
			*pDestinationData++ = *finalRGBPtr++;
		}
	}

	// std::cout << "DONE" << std::endl;
	return 0;
}





/* PCA ANALYSIS */
/****************************************************************************/
cv::Mat TMOMeylan06::RGBToPCA(cv::Mat &data)
{
	this->pca = cv::PCA(data, cv::Mat(), cv::PCA::DATA_AS_ROW);
	return this->pca.project(data);
}


cv::Mat TMOMeylan06::PCAToRGB(cv::Mat &PCAProjection)
{
	return this->pca.backProject(PCAProjection);
}


cv::Mat TMOMeylan06::GetLuminance(cv::Mat &PCAProjection)
{
	cv::Mat luminance = cv::Mat(this->iHeight, this->iWidth, CV_64F);
	double* pcaPtr = PCAProjection.ptr<double>(0);
	double* lumPtr = luminance.ptr<double>(0);
	for (int i = 0; i < this->numberOfPixels; ++i)
	{
		*lumPtr++ = *pcaPtr;
		pcaPtr += 3;
	}
	return luminance;
}


void TMOMeylan06::ScaleSaturation(cv::Mat &data, double saturationEnhancement)
{
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.rows; ++i)
	{
		*dataPtr++;
		*dataPtr++ = *dataPtr * saturationEnhancement;
		*dataPtr++ = *dataPtr * saturationEnhancement;
	}
	return;
}
/****************************************************************************/




/* GLOBAL TONE MAPPING */
/****************************************************************************/
void TMOMeylan06::GlobalToneMap(cv::Mat &data, std::string type)
{
	cv::Mat tmpData = cv::Mat(data);
	tmpData = data.clone();
	if (tmpData.channels() == 3)
	{
		double RScale = 0.299;
		double GScale = 0.587;
		double BScale = 0.114;
		this->ScaleRGB(tmpData, RScale, GScale, BScale);
	}
	if (tmpData.channels() == 1 || tmpData.channels() == 3)
	{
		this->Max(tmpData, 0.001);
	}
	else
	{
		std::cerr << "Maylan06::GlobalToneMap: Wrong number of channels." << std::endl;
		return;
	}
	double scale = 100;
	double AL = this->ComputeAL(tmpData, scale);

	double powExp = 1;
	if (type == "exp")
	{
		powExp = exp(AL + 2) / (exp(4) + (1 / 3.0)) + (1 / 3.0);
	}
	else{if (type == "lin")
	{
		powExp = (1 / 6.0) * AL + (2 / 3.0);
	}
	else
	{
		std::cerr << "Maylan06::GlobalToneMap: Wrong type." << std::endl;
		return;
	}}

	powExp = std::min(std::max(powExp, 0.33333), 1.0);
	this->Pow(data, powExp);
}


double TMOMeylan06::ComputeAL(cv::Mat &data, double scale)
{
	double sum = 0;
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		sum += log((*dataPtr) * scale);
		++dataPtr;
	}
	return sum / data.total();
}


void TMOMeylan06::ScaleRGB(cv::Mat &data, double RScale, double GScale, double BScale)
{
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.rows * data.cols; ++i)
	{
		*dataPtr++ = *dataPtr * RScale;
		*dataPtr++ = *dataPtr * GScale;
		*dataPtr++ = *dataPtr * BScale;
	}
	return;
}
/****************************************************************************/




/* EDGE DETECTION */
/****************************************************************************/
cv::Mat TMOMeylan06::GetEdges(cv::Mat &luminance, double upperThresholdRatio)
{
	cv::Mat edges = cv::Mat(luminance);
	edges = luminance.clone();
	this->Normalize(edges, 0.0, 255.0);
	edges.convertTo(edges, CV_8U);
	int upperThreshold = (int) floor(255 * upperThresholdRatio);
	int lowerThreshold = (int) floor(upperThreshold * 0.4);
	// cv::imwrite("input.png", edges);
	cv::Canny(edges, edges, upperThreshold, lowerThreshold);
	// cv::imwrite("canny_edge.png", edges);
	return edges;
}


cv::Mat TMOMeylan06::DilatateEdges(cv::Mat &edges)
{
	cv::Mat dilatatedEdges = cv::Mat(edges);
	dilatatedEdges = edges.clone();
	uint8 data[10] = {0, 1, 0, 1, 0, 1, 0, 0, 0};
	cv::Mat kernel = cv::Mat(3, 3, CV_8U, data);
	dilate(dilatatedEdges, dilatatedEdges, kernel);
	// cv::imwrite("dilatated_edges.png", dilatatedEdges);
	return dilatatedEdges;
}


cv::Mat TMOMeylan06::ResizeLuminance(cv::Mat &luminance, int maskMaxSize)
{
	cv::Mat luminanceLowRes = cv::Mat(luminance);
	luminanceLowRes = luminance.clone();
	if (luminanceLowRes.cols > maskMaxSize | luminanceLowRes.rows > maskMaxSize)
	{
		double resizeFactor = 0;
		if (luminanceLowRes.cols > luminanceLowRes.rows)
		{
			resizeFactor = 1.0 / (luminanceLowRes.cols / maskMaxSize);
		}
		else
		{
			resizeFactor = 1.0 / (luminanceLowRes.rows / maskMaxSize);
		}
		cv::resize(luminanceLowRes, luminanceLowRes, cv::Size(), resizeFactor, resizeFactor);
	}
	return luminanceLowRes;
}
/****************************************************************************/




/* LOCAL TONE MAPPING */
/****************************************************************************/
cv::Mat TMOMeylan06::GetMask(cv::Mat &luminance, cv::Mat &edges)
{

	cv::Mat mask(luminance.rows, luminance.cols, CV_64F, cv::Scalar(0));
	cv::Mat crossCounter(luminance.rows, luminance.cols, CV_32F, cv::Scalar(-1));

	double maskVal;
	int pixelCounter = 0;
	double *maskPtr = mask.ptr<double>(0);
  for (int i = 0; i < mask.rows; i++)
	{
    for (int j = 0; j < mask.cols; j++)
		{
      *maskPtr = this->GetMaskVal(luminance, edges, crossCounter, j, i, pixelCounter);
      ++pixelCounter;
			++maskPtr;
		}
  }

	// std::cout << "MASK DONE" << std::endl;
	return mask;
}


double TMOMeylan06::GetMaskVal(cv::Mat &luminance, cv::Mat &edges, cv::Mat &crossCounter, int x, int y, int counter)
{
	double sigmaCurrent = this->sigmaOrig;

	double weight;
	double sumMask = 0;
	double sumWeights = 0;

	/****************************************************************************/
	/* RIGHT AND LEFT */
	/****************************************************************************/

	int XActual = x;
	int YActual = y;
	int XOffset = 0;
	int YOffset = 0;
	int XOffsetDiag = 0;
	int YOffsetDiag = 0;

	int rightLeftDirection = 0;
	int bottomTopDirection = 0;

  /* go from the center to the RIGHT or LEFT extremity of the surround */
	rightLeftDirection = 0;
	for (int d = 0; d < 2; ++d)
	{
		/* RIGHT */
		if (d == 0)
		{
			rightLeftDirection = +1;
		}
		/* LEFT */
		else
		{
			rightLeftDirection = -1;
		}
		XActual = x;
		YActual = y;
		XOffset = 0;
		YOffset = 0;
		XOffsetDiag = 0;
		YOffsetDiag = 0;
		for (XOffset = 0; XOffset * rightLeftDirection < this->kernelRadius + 1; XOffset+=rightLeftDirection)
		{
	  	XActual = x + XOffset;
			YActual = y + YOffset;
	    if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
			{
	      if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
				{
					sigmaCurrent = this->sigmaEdge;
					if (XActual > 0 & YActual > 0 & (XActual + 1) < crossCounter.cols & (YActual + 1) < crossCounter.rows)
					{
						crossCounter.at<float>(YActual, XActual + rightLeftDirection) = counter;
					}
				}
	      weight = this->GaussDist(sqrt(pow(XOffset, 2) + pow(YOffset, 2)), sigmaCurrent);
	      sumMask += luminance.at<double>(YActual, XActual) * weight;
	      sumWeights += weight;

	      /* go along diagonal direction, toward the BOTTOM or TOP*/
				bottomTopDirection = 0;
				for (int dd = 0; dd < 2; ++dd)
				{
					/* BOTTOM */
					YOffsetDiag = 0;
					if (dd == 0)
					{
						bottomTopDirection = +1;
					}
					/* TOP */
					else
					{
						bottomTopDirection = -1;
					}
					for (XOffsetDiag = XOffset + rightLeftDirection; XOffsetDiag * rightLeftDirection < this->kernelRadius + 1; XOffsetDiag+=rightLeftDirection)
					{
						YOffsetDiag += bottomTopDirection;
						XActual = x + XOffsetDiag;
						YActual = y + YOffsetDiag;
						if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
						{
						  if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
							{
						    sigmaCurrent = this->sigmaEdge;
						    if (XActual > 0 & YActual > 0 & (XActual + 1) < crossCounter.cols & (YActual + 1) < crossCounter.rows)
								{
									crossCounter.at<float>(YActual, XActual + rightLeftDirection) = counter;
						    }
						  }
							weight = this->GaussDist(sqrt(pow(XOffsetDiag, 2) + pow(YOffsetDiag, 2)), sigmaCurrent);
				      sumMask += luminance.at<double>(YActual, XActual) * weight;
				      sumWeights += weight;
						}
					}
					/* reset to initial value when one radial direction is finished */
		      sigmaCurrent = this->sigmaOrig;
				}
	    }
		}
  }


	/****************************************************************************/
	/* BOTTOM AND TOP */
	/****************************************************************************/

  /* go from the center to the BOTTOM or TOP extremity of the surround */
	bottomTopDirection = 0;
	for (int d = 0; d < 2; ++d)
	{
		/* BOTTOM */
		if (d == 0)
		{
			bottomTopDirection = +1;
		}
		/* TOP */
		else
		{
			bottomTopDirection = -1;
		}
		XActual = x;
		YActual = y;
		XOffset = 0;
		YOffset = 0;
		XOffsetDiag = 0;
		YOffsetDiag = 0;
		for (YOffset = bottomTopDirection; YOffset * bottomTopDirection < this->kernelRadius + 1; YOffset+=bottomTopDirection)
		{
	  	XActual = x + XOffset;
			YActual = y + YOffset;
	    if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
			{
	      if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
				{
					sigmaCurrent = this->sigmaEdge;
					if (XActual > 0 & YActual > 0 & (XActual + 1) < crossCounter.cols & (YActual + 1) < crossCounter.rows)
					{
						crossCounter.at<float>(YActual + bottomTopDirection, XActual) = counter;
					}
				}
	      weight = this->GaussDist(sqrt(pow(XOffset, 2) + pow(YOffset, 2)), sigmaCurrent);
	      sumMask += luminance.at<double>(YActual, XActual) * weight;
	      sumWeights += weight;

	      /* go along diagonal direction, toward the BOTTOM or TOP*/
				rightLeftDirection = 0;
				for (int dd = 0; dd < 2; ++dd)
				{
					/* BOTTOM */
					XOffsetDiag = 0;
					if (dd == 0)
					{
						rightLeftDirection = +1;
					}
					/* TOP */
					else
					{
						rightLeftDirection = -1;
					}
					for (YOffsetDiag = YOffset + bottomTopDirection; YOffsetDiag * bottomTopDirection < this->kernelRadius + 1; YOffsetDiag+=bottomTopDirection)
					{
						XOffsetDiag += rightLeftDirection;
						XActual = x + XOffsetDiag;
						YActual = y + YOffsetDiag;
						if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
						{
						  if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
							{
						    sigmaCurrent = this->sigmaEdge;
						    if (XActual > 0 & YActual > 0 & (XActual + 1) < crossCounter.cols & (YActual + 1) < crossCounter.rows)
								{
									crossCounter.at<float>(YActual + bottomTopDirection, XActual) = counter;
						    }
						  }
							weight = this->GaussDist(sqrt(pow(XOffsetDiag, 2) + pow(YOffsetDiag, 2)), sigmaCurrent);
				      sumMask += luminance.at<double>(YActual, XActual) * weight;
				      sumWeights += weight;
						}
					}
					/* reset to initial value when one radial direction is finished */
		      sigmaCurrent = this->sigmaOrig;
				}
	    }
		}
  }
	if (sumWeights == 0)
	{
		sumWeights = 0.00000001;
	}
	return sumMask / sumWeights;
}


bool TMOMeylan06::IsAnEdge(cv::Mat &edges, cv::Mat &crossCounter, int x, int y, int counter)
{
	return (((int)edges.at<uint8>(y, x) > 100) | ((int)crossCounter.at<float>(y, x) == counter));
}


double TMOMeylan06::GaussDist(double d, double s)
{
  return exp(-pow(d, 2) / pow(2 * s, 2));
}
/****************************************************************************/




/* BETA FACTOR */
/****************************************************************************/
cv::Mat TMOMeylan06::GetBetaFactor(cv::Mat &luminance, double c, double a)
{
	cv::Mat betaFactor = cv::Mat(luminance);
  betaFactor = luminance.clone();
	this->Min(betaFactor, 1.0);
	this->Max(betaFactor, 0.0);
	double *betaFactorPtr = betaFactor.ptr<double>(0);
	for (int i = 0; i < this->numberOfPixels; ++i)
	{
		*betaFactorPtr++ = abs(1 - (1 / (1 + exp(-a * (*betaFactorPtr - c)))));
	}
	this->Normalize(betaFactor, 0.0, 1.0);
	return betaFactor;
}
/****************************************************************************/




/* HELPER METHODS */
/****************************************************************************/
void TMOMeylan06::Normalize(cv::Mat &data)
{
	double min, max;
	cv::minMaxLoc(data, &min, &max);
	if (max == 0)
	{
		return;
	}
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		*dataPtr = *dataPtr / max;
		++dataPtr;
	}
	return;
}


void TMOMeylan06::HistoClip(cv::Mat &data, int numberOfBuckets, double minThreshold, double maxThreshold)
{
	double min, max;
	cv::minMaxLoc(data, &min, &max);
	double range = std::abs(max - min);
	double bucketSize = range / numberOfBuckets;

	std::vector<std::vector<double>> histogram(numberOfBuckets);
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		int bucketIndex = (int) (floor(((*dataPtr) - min) / bucketSize));
		bucketIndex = std::min(bucketIndex, numberOfBuckets - 1);
		bucketIndex = std::max(0, bucketIndex);
		histogram[bucketIndex].push_back(*dataPtr);
		++dataPtr;
	}

	std::vector<double> cumulativeHistogram(numberOfBuckets);
	int cumul = 0;
	for (size_t i = 0; i < numberOfBuckets; ++i)
	{
	  cumul += (int) histogram[i].size();
		double fraction = cumul / ((double) data.total());
		cumulativeHistogram[i] = fraction;
	}

	double newMin = min;
	for (int i = 0; i < numberOfBuckets; ++i)
	{
		if (cumulativeHistogram[i] > minThreshold)
		{
			newMin = this->GetMin(histogram[i].data(), (int) histogram[i].size());
			break;
		}
	}

	double newMax = max;
	for (int i = cumulativeHistogram.size() - 1; i >= 0; --i)
	{
		if (cumulativeHistogram[i] < maxThreshold)
		{
			newMax = this->GetMax(histogram[i].data(), (int) histogram[i].size());
			break;
		}
	}

	this->Max(data, newMin);
	this->Min(data, newMax);
	this->Normalize(data, 0.0, 1.0);
}


void TMOMeylan06::LogMaxScale(cv::Mat &data, double max, double scale)
{
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		double tmpScale = (*dataPtr) * scale;
		if (max > tmpScale)
		{
			*dataPtr = log(max) / log(scale);
		}
		else
		{
			*dataPtr = log(tmpScale) / log(scale);
		}
		++dataPtr;
	}
	return;
}


void TMOMeylan06::Normalize(cv::Mat &data, double lowerBound, double upperBound)
{
	cv::normalize(data, data, lowerBound, upperBound, cv::NORM_MINMAX);
	return;
}


cv::Mat TMOMeylan06::ElementWiseMul(cv::Mat &first, cv::Mat &second)
{
	cv::Mat result(first.rows, first.cols, CV_64F, cv::Scalar(0));
	double *firstPtr = first.ptr<double>(0);
	double *secondPtr = second.ptr<double>(0);
	double *resultPtr = result.ptr<double>(0);
	for (int i = 0; i < first.rows * first.cols; ++i)
	{
		*resultPtr = (*firstPtr) * (*secondPtr);
		++firstPtr;
		++secondPtr;
		++resultPtr;
	}
	return result;
}


cv::Mat TMOMeylan06::ElementWiseSub(cv::Mat &first, cv::Mat &second)
{
	cv::Mat result(first.rows, first.cols, CV_64F, cv::Scalar(0));
	double *firstPtr = first.ptr<double>(0);
	double *secondPtr = second.ptr<double>(0);
	double *resultPtr = result.ptr<double>(0);
	for (int i = 0; i < first.rows * first.cols; ++i)
	{
		*resultPtr = (*firstPtr) - (*secondPtr);
		++firstPtr;
		++secondPtr;
		++resultPtr;
	}
	return result;
}


void TMOMeylan06::Pow(cv::Mat &data, double exponent)
{
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		*dataPtr = pow(*dataPtr, exponent);
		++dataPtr;
	}
}


void TMOMeylan06::Max(cv::Mat &data, double max)
{
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		if (*dataPtr < max)
		{
			*dataPtr = max;
		}
		++dataPtr;
	}
	return;
}


void TMOMeylan06::Min(cv::Mat &data, double min)
{
	double *dataPtr = data.ptr<double>(0);
	for (int i = 0; i < data.total(); ++i)
	{
		if (*dataPtr > min)
		{
			*dataPtr = min;
		}
		++dataPtr;
	}
	return;
}


double TMOMeylan06::GetMax(double *data, int dataLength)
{
	double max = data[0];
	for (int i = 1; i < dataLength; ++i)
	{
		if (data[i] > max)
		{
			max = data[i];
		}
	}
	return max;
}


double TMOMeylan06::GetMin(double *data, int dataLength)
{
	double min = data[0];
	for (int i = 1; i < dataLength; ++i)
	{
		if (data[i] < min)
		{
			min = data[i];
		}
	}
	return min;
}


void TMOMeylan06::SaveImg(std::string name, cv::Mat &data)
{
	cv::Mat dataToShow = cv::Mat(data);
	dataToShow = data.clone();
	this->Normalize(dataToShow, 0.0, 255.0);
	if (dataToShow.channels() == 1)
	{
		dataToShow.convertTo(dataToShow, CV_8UC1);
	}
	else
	{
		dataToShow.convertTo(dataToShow, CV_8UC3);
	}
	cv::imwrite(name, dataToShow);
	return;
}
/****************************************************************************/
