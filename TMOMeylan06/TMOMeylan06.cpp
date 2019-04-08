/* --------------------------------------------------------------------------- *
 * TMOMeylan06.cpp: implementation of the TMOMeylan06 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMeylan06.h"
#include "canny.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMeylan06::TMOMeylan06()
{
	SetName(L"Meylan06");
	SetDescription(L"High Dynamic Range Image Rendering Using a Retinex-Based Adaptive Filter");

	this->saturationParameter.SetName(L"Saturation");
	this->saturationParameter.SetDescription(L"Color enhancement of two last PCA components.");
	this->saturationParameter.SetDefault(1.6254);
	this->saturationParameter=1.6254;
	this->saturationParameter.SetRange(1.0,10.0);
	this->Register(this->saturationParameter);

	std::cout << std::endl;
}

TMOMeylan06::~TMOMeylan06()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOMeylan06::Transform()
{

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	this->numberOfPixels = pSrc->GetHeight() * pSrc->GetWidth();
	this->numberOfPixelsRGB = this->numberOfPixels * 3;


	// LUMINANCE PROCESSING
	// ***************************************************************************
	this->Normalize(pSourceData, this->numberOfPixelsRGB);

	cv::Mat PCAProjectionForLuminance = this->RGBToPCA(pSourceData);
	cv::Mat luminance = this->GetLuminance(PCAProjectionForLuminance);

	this->Normalize(luminance.ptr<double>(0), this->numberOfPixels, 0.0, 1.0);

	this->GlobalMapping(luminance.ptr<double>(0), this->numberOfPixels, 1, "exp");

	cv::Mat luminanceLowRes = cv::Mat(luminance);
	luminanceLowRes = luminance.clone();
	luminanceLowRes = this->ReshapeGray(luminanceLowRes);

	int maskMaxSize = 200;

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

	cv::Mat edges = this->GetEdges(luminanceLowRes, 0.2);
	edges = this->DilatateEdges(edges);

	this->sigmaOrig = 16;
	this->sigmaEdge = 1;
	//this->kernelSize = (int) 3 * this->sigmaOrig + 1;
	this->kernelSize = 6;
	cv::Mat mask = this->GetMask(luminanceLowRes, edges);

	cv::Mat maskHighRes(this->iHeight, this->iWidth, CV_64F);
	cv::resize(mask, maskHighRes, cv::Size(this->iWidth, this->iHeight));

	cv::Mat maskToShow = cv::Mat(maskHighRes);
	maskToShow = maskHighRes.clone();
	this->Normalize(maskToShow.ptr<double>(0), maskToShow.cols * maskToShow.rows, 0.0, 255.0);
	cv::imwrite("mask.png", maskToShow);


	this->LogMaxScale(luminance.ptr<double>(0), this->numberOfPixels, 0.1, 100);

	this->HistoClip(luminance.ptr<double>(0), this->numberOfPixels, 100, 0.01, 0.99);


	std::cout << "LUM PROCESS DONE" << std::endl;
	// ***************************************************************************


	// COLOR PROCESSING
	// ***************************************************************************
	std::unique_ptr<double[]> rgbData = std::make_unique<double[]>(this->numberOfPixelsRGB);
	std::memcpy(rgbData.get(), pSourceData, this->numberOfPixels * 3 * sizeof(double));

	this->GlobalMapping(rgbData.get(), this->numberOfPixels * 3, 3, "exp");
	this->LogMaxScale(rgbData.get(), this->numberOfPixels * 3, 1, 100);

	cv::Mat PCAProjectionForRGB = this->RGBToPCA(rgbData.get());

	double saturationEnhancement = this->saturationParameter;
	double* PCAProjectionForRGBPtr = PCAProjectionForRGB.ptr<double>(0);
	// R -> luminance, G and B -> colors
	this->ScaleRGB(PCAProjectionForRGBPtr, 1, saturationEnhancement, saturationEnhancement);

	cv::Mat luminanceRGB = this->GetLuminance(PCAProjectionForRGB);
	double luminanceRGBMax = this->GetMax(luminanceRGB.ptr<double>(0), this->numberOfPixels);
	double luminanceRGBMin = this->GetMin(luminanceRGB.ptr<double>(0), this->numberOfPixels);
	this->Normalize(luminance.ptr<double>(0), this->numberOfPixels, luminanceRGBMin, luminanceRGBMax);

	// Recompose the projection with processed luminance
	PCAProjectionForRGBPtr = PCAProjectionForRGB.ptr<double>(0);
	double* luminancePtr = luminance.ptr<double>(0);
	for (int i = 0; i < this->numberOfPixels; ++i)
	{
		*PCAProjectionForRGBPtr = *luminancePtr;
		PCAProjectionForRGBPtr += 3;
		luminancePtr++;
	}
	cv::Mat finalRGB = this->PCAToRGB(PCAProjectionForRGB);

	this->HistoClip(finalRGB.ptr<double>(0), this->numberOfPixels * 3, 100, 0.01, 0.99);

	std::cout << "COLOR PROCESS DONE" << std::endl;
	// ***************************************************************************

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
	std::cout << "DONE" << std::endl;
	return 0;
}


cv::Mat TMOMeylan06::GetEdges(cv::Mat &luminance, double upperThresholdRatio)
{
	cv::Mat edges = cv::Mat(luminance);
	edges = luminance.clone();
	this->Normalize(edges.ptr<double>(0), edges.rows * edges.cols, 0.0, 255.0);
	edges.convertTo(edges, CV_8U);
	int upperThreshold = (int) floor(255 * upperThresholdRatio);
	int lowerThreshold = (int) floor(upperThreshold * 0.4);
	cv::imwrite("input.png", edges);
	cv::Canny(edges, edges, upperThreshold, lowerThreshold);
	cv::imwrite("canny_edge.png", edges);
	return edges;
}


cv::Mat TMOMeylan06::DilatateEdges(cv::Mat &edges)
{
	cv::Mat dilatatedEdges = cv::Mat(edges);
	dilatatedEdges = edges.clone();
	uint8 data[10] = {0, 1, 0, 1, 0, 1, 0, 0, 0};
	cv::Mat kernel = cv::Mat(3, 3, CV_8U, data);
	dilate(dilatatedEdges, dilatatedEdges, kernel);
	cv::imwrite("dilatated_edges.png", dilatatedEdges);
	return dilatatedEdges;
}


cv::Mat TMOMeylan06::GetMask(cv::Mat &luminance, cv::Mat &edges)
{

	cv::Mat mask(luminance.rows, luminance.cols, CV_64F);
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

	std::cout << "MASK DONE" << std::endl;
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
			//std::cout << "RIGHT" << std::endl;
			rightLeftDirection = +1;
		}
		/* LEFT */
		else
		{
			//std::cout << "LEFT" << std::endl;
			rightLeftDirection = -1;
		}
		XActual = x;
		YActual = y;
		XOffset = 0;
		YOffset = 0;
		XOffsetDiag = 0;
		YOffsetDiag = 0;
		for (XOffset = 0; XOffset * rightLeftDirection < this->kernelSize + 1; XOffset+=rightLeftDirection)
		{
	  	XActual = x + XOffset;
			YActual = y + YOffset;
	    if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
			{
	      if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
				{
					//std::cout << "EDGE" << std::endl;
					sigmaCurrent = this->sigmaEdge;
					if (XActual + rightLeftDirection < crossCounter.cols)
					{
						crossCounter.at<float>(YActual, XActual + rightLeftDirection) = counter;
					}
				}
				//std::cout << XActual << " " << YActual << std::endl;
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
					for (XOffsetDiag = XOffset + rightLeftDirection; XOffsetDiag * rightLeftDirection < this->kernelSize + 1; XOffsetDiag+=rightLeftDirection)
					{
						YOffsetDiag += bottomTopDirection;
						XActual = x + XOffsetDiag;
						YActual = y + YOffsetDiag;
						if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
						{
						  if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
							{
								//std::cout << "EDGE" << std::endl;
						    sigmaCurrent = this->sigmaEdge;
						    if (XActual > 0 & YActual > 0 & (XActual + 1) < crossCounter.cols & (YActual + 1) < crossCounter.rows)
								{
									crossCounter.at<float>(YActual, XActual + rightLeftDirection) = counter;
						    }
						  }
							//std::cout << XActual << " " << YActual << std::endl;
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

	//std::cout << std::endl << std::endl << std::endl;


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
			//std::cout << "BOTTOM" << std::endl;
			bottomTopDirection = +1;
		}
		/* TOP */
		else
		{
			//std::cout << "TOP" << std::endl;
			bottomTopDirection = -1;
		}
		XActual = x;
		YActual = y;
		XOffset = 0;
		YOffset = 0;
		XOffsetDiag = 0;
		YOffsetDiag = 0;
		for (YOffset = bottomTopDirection; YOffset * bottomTopDirection < this->kernelSize + 1; YOffset+=bottomTopDirection)
		{
	  	XActual = x + XOffset;
			YActual = y + YOffset;
	    if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
			{
	      if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
				{
					//std::cout << "EDGE" << std::endl;
					sigmaCurrent = this->sigmaEdge;
					if (YActual + bottomTopDirection < crossCounter.cols)
					{
						crossCounter.at<float>(YActual + bottomTopDirection, XActual) = counter;
					}
				}
				//std::cout << XOffset << " " << YOffset << std::endl;
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
					for (YOffsetDiag = YOffset + bottomTopDirection; YOffsetDiag * bottomTopDirection < this->kernelSize + 1; YOffsetDiag+=bottomTopDirection)
					{
						XOffsetDiag += rightLeftDirection;
						XActual = x + XOffsetDiag;
						YActual = y + YOffsetDiag;
						if (XActual < luminance.cols && YActual < luminance.rows && XActual >= 0 && YActual >= 0)
						{
						  if (this->IsAnEdge(edges, crossCounter, XActual, YActual, counter))
							{
								//std::cout << "EDGE" << std::endl;
						    sigmaCurrent = this->sigmaEdge;
						    if (XActual > 0 & YActual > 0 & (XActual + 1) < crossCounter.cols & (YActual + 1) < crossCounter.rows)
								{
									crossCounter.at<float>(YActual + bottomTopDirection, XActual) = counter;
						    }
						  }
							//std::cout << XOffsetDiag << " " << YOffsetDiag << std::endl;
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

	//std::cout << std::endl << std::endl << std::endl;

	if (sumWeights == 0)
	{
		sumWeights = 0.00000001;
	}
	return sumMask / sumWeights;
}


bool TMOMeylan06::IsAnEdge(cv::Mat &edges, cv::Mat &crossCounter, int x, int y, int counter)
{
	//std::cout << (int)edges.at<uint8>(y, x) << std::endl;
	return (((int)edges.at<uint8>(y, x) > 100) | ((int)crossCounter.at<float>(y, x) == counter));
}


double TMOMeylan06::GaussDist(double d, double s)
{
  return exp(-pow(d, 2) / pow(2 * s, 2));
}


cv::Mat TMOMeylan06::RGBToPCA(double *rgbSourceData)
{
	cv::Mat rgb = cv::Mat(this->numberOfPixels, 3, CV_64F, rgbSourceData);
	this->pca = cv::PCA(rgb, cv::Mat(), cv::PCA::DATA_AS_ROW);
	return this->pca.project(rgb);
}


cv::Mat TMOMeylan06::PCAToRGB(cv::Mat &PCAProjection)
{
	return this->pca.backProject(PCAProjection);
}


cv::Mat TMOMeylan06::GetLuminance(cv::Mat &PCAProjection)
{
	cv::Mat luminance = cv::Mat(this->numberOfPixels, 1, CV_64F);
	double* pcaPtr = PCAProjection.ptr<double>(0);
	double* lumPtr = luminance.ptr<double>(0);
	for (int i = 0; i < this->numberOfPixels; ++i)
	{
		*lumPtr++ = *pcaPtr;
		pcaPtr += 3;
	}
	return luminance;
}


void TMOMeylan06::GlobalMapping(double* data, int dataLength, int numberOfChannels, std::string type)
{
	std::unique_ptr<double[]> tmpData = std::make_unique<double[]>(dataLength);
	std::memcpy(tmpData.get(), data, dataLength * sizeof(double));
	if (numberOfChannels == 3)
	{
		this->ScaleRGB(tmpData.get(), 0.299, 0.587, 0.114);
	}
	if (numberOfChannels == 1 || numberOfChannels == 3)
	{
		this->Max(tmpData.get(), dataLength, 0.001);
	}
	else
	{
		std::cerr << "Maylan06::GlobalMapping: Wrong number of channels." << std::endl;
		return;
	};
	double AL = this->ComputeAL(tmpData.get(), dataLength, 100);

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
		std::cerr << "Maylan06::GlobalMapping: Wrong type." << std::endl;
		return;
	}};

	powExp = std::min(std::max(powExp, 0.33333), 1.0);
	this->Pow(data, dataLength, powExp);
}


double TMOMeylan06::ComputeAL(double *data, int dataLength, double scale)
{
	double sum = 0;
	for (int i = 0; i < dataLength; ++i)
	{
		sum += log(data[i] * scale);
	}
	return sum / dataLength;
}


void TMOMeylan06::HistoClip(double* data, int dataLength, int numberOfBuckets, double minThreshold, double maxThreshold)
{
	double min = this->GetMin(data, dataLength);
	double max = this->GetMax(data, dataLength);;
	double range = std::abs(max - min);
	double bucketSize = range / numberOfBuckets;

	std::vector<std::vector<double>> histogram(numberOfBuckets);
	for (int i = 0; i < dataLength; ++i)
	{
		int bucketIndex = (int) (floor((data[i] - min) / bucketSize));
		bucketIndex = std::min(bucketIndex, numberOfBuckets - 1);
		histogram[bucketIndex].push_back(data[i]);
	}

	std::vector<double> cumulativeHistogram(numberOfBuckets);
	int cumul = 0;
	for (size_t i = 0; i < numberOfBuckets; ++i)
	{
	  cumul += (int) histogram[i].size();
		double fraction = cumul / ((double) dataLength);
		cumulativeHistogram[i] = fraction;
	}

	double newMin = this->GetMin(data, dataLength);
	for (int i = 0; i < numberOfBuckets; ++i)
	{
		if (cumulativeHistogram[i] > minThreshold)
		{
			newMin = this->GetMin(histogram[i].data(), (int) histogram[i].size());
			break;
		}
	}

	double newMax = this->GetMax(data, dataLength);
	for (int i = cumulativeHistogram.size() - 1; i >= 0; --i)
	{
		if (cumulativeHistogram[i] < maxThreshold)
		{
			newMax = this->GetMax(histogram[i].data(), (int) histogram[i].size());
			break;
		}
	}

	this->Max(data, dataLength, newMin);
	this->Min(data, dataLength, newMax);
	this->Normalize(data, dataLength, 0.0, 1.0);
}


void TMOMeylan06::ScaleRGB(double* data, double RScale, double GScale, double BScale)
{
	double *rgbPtr = data;
	for (int i = 0; i < this->numberOfPixels; ++i)
	{
		*data++ = *data * RScale;
		*data++ = *data * GScale;
		*data++ = *data * BScale;
	}
	return;
}


void TMOMeylan06::LogMaxScale(double *data, int dataLength, double max, double scale)
{
	for (int i = 0; i < dataLength; ++i)
	{
		double tmpScale = data[i] * scale;
		if (max > tmpScale)
		{
			data[i] = log(max) / log(scale);
		}
		else
		{
			data[i] = log(tmpScale) / log(scale);
		}
	}
	return;
}


void TMOMeylan06::Pow(double *data, int dataLength, double exponent)
{
	for (int i = 0; i < dataLength; ++i)
	{
		data[i] = pow(data[i], exponent);
	}
}


void TMOMeylan06::Max(double *data, int dataLength, double max)
{
	for (int i = 0; i < dataLength; ++i)
	{
		if (data[i] < max)
		{
			data[i] = max;
		}
	}
	return;
}


void TMOMeylan06::Min(double *data, int dataLength, double min)
{
	for (int i = 0; i < dataLength; ++i)
	{
		if (data[i] > min)
		{
			data[i] = min;
		}
	}
	return;
}


void TMOMeylan06::Normalize(double *data, int dataLength)
{
	double max = this->GetMax(data, dataLength);
	if (max == 0)
	{
		return;
	}
	for (int i = 0; i < dataLength; ++i)
	{
		data[i] = data[i] / max;
	}
	return;
}


void TMOMeylan06::Normalize(double *data, int dataLength, double lowerBound, double upperBound)
{
	double max = this->GetMax(data, dataLength);
	double min = this->GetMin(data, dataLength);
	if (max == 0)
	{
		return;
	}
	for (int i = 0; i < dataLength; ++i)
	{
		data[i] = lowerBound + ((data[i] - min) / (max - min)) * (upperBound - lowerBound);
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

void TMOMeylan06::SaveImg(std::string name, double *data, bool RGB)
{
	int x = this->iWidth;
	int y = this->iHeight;
	if (RGB)
	{
		this->Normalize(data, x * y * 3, 0.0, 255.0);
		int b;
		int g;
		int r;
		cv::Mat imgRGB(y, x, CV_8UC3, cv::Scalar(0, 0, 0));
		uint8 *imgRGBPtr = imgRGB.ptr<uint8>(0);
		for (int i = 0; i < y; ++i)
		{
			for (int j = 0; j < x; ++j)
			{
				*imgRGBPtr++ = (int) *data++;
				*imgRGBPtr++ = (int) *data++;
				*imgRGBPtr++ = (int) *data++;
			}
		}
		cv::imwrite(name, imgRGB);
	}
	else
	{
		this->Normalize(data, x * y, 0.0, 255.0);
		cv::Mat imgGray(y, x, CV_8UC1, cv::Scalar(0));
		uint8 *imgGrayPtr = imgGray.ptr<uint8>(0);
		for (int i = 0; i < y; ++i)
		{
			for (int j = 0; j < x; ++j)
			{
				*imgGrayPtr++ = (int) *data++;
			}
		}
		cv::imwrite(name, imgGray);
	}
	return;
}


cv::Mat TMOMeylan06::ReshapeGray(cv::Mat &source)
{
	int x = this->iWidth;
	int y = this->iHeight;
	double *sourcePtr = source.ptr<double>(0);
	cv::Mat result(y, x, CV_64F, cv::Scalar(0));
	double *resultPtr = result.ptr<double>(0);
	for (int i = 0; i < y; ++i)
	{
		for (int j = 0; j < x; ++j)
		{
			*resultPtr++ = *sourcePtr++;
		}
	}
	return result;
}
