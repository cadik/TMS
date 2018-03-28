/* --------------------------------------------------------------------------- *
 * TMOEisemann04.cpp: implementation of the TMOEisemann04 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOEisemann04.h"
#include <cmath>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOEisemann04::TMOEisemann04()
{
	SetName(L"Eisemann04");						// TODO - Insert operator name
	SetDescription(L"Flash/No-flash");	// TODO - Insert description

	flashImagePathParameter.SetName(L"Flash");
	flashImagePathParameter.SetDescription(L"Input flash image path");
	flashImagePathParameter.SetDefault("");
	this->Register(flashImagePathParameter);

	
}

TMOEisemann04::~TMOEisemann04()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOEisemann04::Transform()
{

	// Load flash image
	if(flashImagePathParameter.GetString() == ""){
		std::cerr << "No flash image provided" << std::endl;
		return 1;
	}
	flashImage = TMOImage();
	flashImage.Open(flashImagePathParameter.GetString().c_str());

	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_RGB);								// This is format of Y as luminance
	pDst->Convert(TMO_RGB);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components
	double pr, pg, pb;


	double* intensityF = ComputeIntensity(&flashImage, NULL);
	double* intensityNF = ComputeIntensity(&flashImage, pSrc);
	double* intensityNFLarge = ComputeIntensity(pSrc, NULL);

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
		//pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pr = *pSourceData++;
			pg = *pSourceData++;
			pb = *pSourceData++;


			// Here you can use your transform 
			// expressions and techniques...
			//pY *= dParameter;							// Parameters can be used like
														// simple variables

			double intensity = intensityNF[j*pSrc->GetWidth()+i];
			// and store results to the destination image
			*pDestinationData++ = intensity;
			*pDestinationData++ = intensity;//px;
			*pDestinationData++ = intensity;//py;
		}
	//pSrc->ProgressBar(j, pSrc->GetHeight());
	}
	free(intensityF);
	free(intensityNF);
	free(intensityNFLarge);

	pDst->Convert(TMO_RGB);
	return 0;
}

double* TMOEisemann04::ComputeIntensity(TMOImage *inputImage, TMOImage *weightImage){
	double* intensity = (double*)malloc(sizeof(double)*inputImage->GetWidth()*inputImage->GetHeight());
	double* inImageData = inputImage->GetData();
	double* wgImageData = weightImage != NULL ? weightImage->GetData() : NULL;
	if(weightImage == NULL)
		std::cerr << "single "<<std::endl;
	else
		std::cerr << "Two "<<std::endl;

	double pr,pg,pb;
	for (int j = 0; j < inputImage->GetHeight(); j++)
	{
		//pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
		for (int i = 0; i < inputImage->GetWidth(); i++)
		{

			pr = *inImageData++;
			pg = *inImageData++;
			pb = *inImageData++;
			double sum = pr + pg + pb;
			if(sum == 0.)
				std::cerr << pr << " " <<pg<< " "<<pb << " " << "SUM:" <<sum << std::endl;
			if(wgImageData != NULL){

				pr = *wgImageData++;
				pg = *wgImageData++;
				pb = *wgImageData++;	
				
			}

			double value = (pr/sum)*pr + (pg/sum)*pg + (pb/sum)*pb;
			intensity[j*inputImage->GetWidth()+i] = fmax(fmin(value, 1.), 0.);
		}
	//pSrc->ProgressBar(j, pSrc->GetHeight());
	}
	// normalize
	double min = intensity[0]; 
	double max = intensity[0]; 
	for(int i = 1; i < inputImage->GetHeight() * inputImage->GetWidth(); i++){
		if(intensity[i] > max){
			max = intensity[i];
		}
		if(intensity[i] < min){
			min = intensity[i];
		}
	}
	std::cerr << min << " " << max << std::endl;
	for(int i = 0; i < inputImage->GetHeight() * inputImage->GetWidth(); i++){
		intensity[i] = (intensity[i] - min)/ (max-min);// * 255.0;
	}

	return intensity;
}