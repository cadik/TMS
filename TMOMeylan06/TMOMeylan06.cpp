/* --------------------------------------------------------------------------- *
 * TMOMeylan06.cpp: implementation of the TMOMeylan06 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOMeylan06.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOMeylan06::TMOMeylan06()
{
	SetName(L"Meylan06");						// TODO - Insert operator name
	SetDescription(L"Add your TMO description here");	// TODO - Insert description
	/*
	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription");	// TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter=1.;
	dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
	this->Register(dParameter);
	*/
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

	this->Normalize(pSourceData, this->numberOfPixels * 3);

	cv::Mat PCAProjection = this->RGBToPCA(pSourceData);
	cv::Mat luminance = this->GetLuminance(PCAProjection);
	this->Normalize(luminance.ptr<double>(0), this->numberOfPixels, 0.0, 1.0);
	//std::cout << "NORMALIZED LUMINANCE" << std::endl;
	//std::cout << luminance << std::endl;

	this->GlobalMapping(luminance.ptr<double>(0), this->numberOfPixels, 1, "exp");
	this->LogMaxScale(luminance.ptr<double>(0), this->numberOfPixels, 0.1, 100);

	this->HistoClip(luminance.ptr<double>(0), this->numberOfPixels, 100, 0.01, 0.99);

	double* lumPtr = luminance.ptr<double>(0);
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			*pDestinationData++ = *lumPtr;
			*pDestinationData++ = *lumPtr;
			*pDestinationData++ = *lumPtr;
			++lumPtr;
		}
	}
	std::cout << "DONE" << std::endl;
	return 0;
}

cv::Mat TMOMeylan06::RGBToPCA(double *rgbSourceData)
{
	cv::Mat rgb = cv::Mat(this->numberOfPixels, 3, CV_64F);
	double* rgbPtr = rgb.ptr<double>(0);
	for (int i = 0; i < this->numberOfPixels; i++)
	{
		*rgbPtr++ = *rgbSourceData++;
		*rgbPtr++ = *rgbSourceData++;
		*rgbPtr++ = *rgbSourceData++;
	}
	rgbSourceData = pSrc->GetData();
	//std::cout << "INPUT" << std::endl;
	//std::cout << rgb << std::endl;
	this->pca = cv::PCA(rgb, cv::Mat(), cv::PCA::DATA_AS_ROW);
 	cv::Mat meanVector = this->pca.mean;
 	cv::Mat eigenVectors = this->pca.eigenvectors;
	//std::cout << "MEAN" << meanVector << std::endl;
	//std::cout << "EIGEN VECTORS" << eigenVectors << std::endl;
	cv::Mat PCAProjection = this->pca.project(rgb);
	//std::cout << "PROJECTION" << PCAProjection;
	return PCAProjection;
}

cv::Mat TMOMeylan06::PCAToRGB(cv::Mat &PCAProjection)
{
	cv::Mat reconstruction = this->pca.backProject(PCAProjection);
	//std::cout << "BACK PROJECTION" << reconstruction;
	return reconstruction;
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
	//std::cout << "LUMINANCE" << std::endl;
	//std::cout << luminance << std::endl;
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
	/*
	std::cout << min << std::endl;
	std::cout << range << std::endl;
	std::cout << bucketSize << std::endl;
	*/

	//std::cout << "HIST START" << std::endl;
	std::vector<std::vector<double>> histogram(numberOfBuckets);
	for (int i = 0; i < dataLength; ++i)
	{
		int bucketIndex = (int) (floor((data[i] - min) / bucketSize));
		//std::cout << data[i] << std::endl;
		//std::cout << floor(data[i] - min) / bucketSize << std::endl;
		//std::cout << bucketIndex << std::endl;
		//std::cout << std::endl;
		bucketIndex = std::min(bucketIndex, numberOfBuckets - 1);
		histogram[bucketIndex].push_back(data[i]);
	}
	//std::cout << "HIST DONE" << std::endl;

	std::vector<double> cumulativeHistogram(numberOfBuckets);
	int cumul = 0;
	for (size_t i = 0; i < numberOfBuckets; ++i)
	{
		//std::cout << histogram[i].size() << std::endl;
	  cumul += (int) histogram[i].size();
		//std::cout << cumul << std::endl;
		//std::cout << std::endl;
		double fraction = cumul / ((double) dataLength);
		cumulativeHistogram[i] = fraction;
	}

	double newMin = -1000000000.0;
	for (int i = 0; i < numberOfBuckets; ++i)
	{
		if (cumulativeHistogram[i] > minThreshold)
		{
			newMin = this->GetMin(histogram[i].data(), (int) histogram[i].size());
			//std::cout << "NEW MIN: " << newMin << std::endl;
			break;
		}
	}

	double newMax = 1000000000.0;
	for (int i = cumulativeHistogram.size() - 1; i >= 0; --i)
	{
		if (cumulativeHistogram[i] < maxThreshold)
		{
			newMax = this->GetMax(histogram[i].data(), (int) histogram[i].size());
			//std::cout << "NEW MAX: " << newMax << std::endl;
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

//void TMOMeylan06::ClipHist(double *data, int dataLength, )

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
