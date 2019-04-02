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

	cv::Mat luminanceForEdgeDetection = cv::Mat(luminance);
	luminanceForEdgeDetection = luminance.clone();
	this->Normalize(luminanceForEdgeDetection.ptr<double>(0), this->numberOfPixels, 0.0, 255.0);
	luminanceForEdgeDetection = this->ResizeGray(luminanceForEdgeDetection);
	luminanceForEdgeDetection.convertTo(luminanceForEdgeDetection, CV_8U);
	int upperThreshold = (int) floor(255 * 0.2);
	int lowerThreshold = (int) floor(upperThreshold * 0.4);
	cv::imwrite("input.jpg", luminanceForEdgeDetection);
	cv::Canny(luminanceForEdgeDetection, luminanceForEdgeDetection, upperThreshold, lowerThreshold);
	cv::imwrite("canny_edge.jpg", luminanceForEdgeDetection);

	/*
	int dilationSize = 1;
	cv::Mat element = cv::getStructuringElement(0, cv::Size(2 * dilationSize + 1, 2 * dilationSize + 1),
																			cv::Point(dilationSize, dilationSize));
	cv::dilate(luminanceForEdgeDetection, luminanceForEdgeDetection, element);
	cv::imwrite("canny_edge_with_dilatation.jpg", luminanceForEdgeDetection);
	*/
	
	this->LogMaxScale(luminance.ptr<double>(0), this->numberOfPixels, 0.1, 100);

	this->HistoClip(luminance.ptr<double>(0), this->numberOfPixels, 100, 0.01, 0.99);



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


cv::Mat TMOMeylan06::ResizeGray(cv::Mat &source)
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
