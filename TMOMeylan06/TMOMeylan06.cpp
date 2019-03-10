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

	cv::Mat PCAProjection = this->RGBToPCA(pSourceData);
	cv::Mat reconstruction = this->PCAToRGB(PCAProjection);

	std::cout << "PCA DONE" << std::endl;

	pSourceData = reconstruction.ptr<double>(0);

	double pY, px, py;

        int j=0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			pY *= dParameter;							// Parameters can be used like
														// simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	std::cout << "DONE" << std::endl;
	pDst->Convert(TMO_RGB);
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
 	//cv::Mat meanVector = this->pca.mean;
 	//cv::Mat eigenVectors = this->pca.eigenvectors;
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
