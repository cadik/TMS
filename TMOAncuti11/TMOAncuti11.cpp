/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       Diploma thesis                                         *
*                       Author: Matus Bicanovsky                               *
*                       Brno 2025                                              *
*                                                                              *
*                       Implementation of the TMOAncuti10 class                *
*                       Enhancing by saliency-guided decolorization               *
*******************************************************************************/
/* --------------------------------------------------------------------------- *
 * TMOAncuti11.cpp: implementation of the TMOAncuti11 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAncuti11.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAncuti11::TMOAncuti11()
{
	SetName(L"Ancuti11");					  // TODO - Insert operator name
	SetDescription(L"Operator for C2G image and video conversion based on paper Enhancing by saliency-guided decolorization"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOAncuti11::~TMOAncuti11()
{
}

cv::Mat TMOAncuti11::decolorization(cv::Mat &input, double eta)
{
	// default values of parameters from the paper
	double mu = 0.1;   //saturation threshold
	double nu = 0.6;   //lightness threshold
	double phi = 300.0; //offset angle in degrees
	double kappa = 2.0; //period
	double gamma = 0.7; //constrast gain
	double chi = 0.9; //global intensity scale
	
	//convert input to HLS
	cv::Mat input8U;
	input.convertTo(input8U, CV_8UC3, 255.0);
	cv::Mat hlsInput;
	cv::cvtColor(input8U, hlsInput, cv::COLOR_BGR2HLS);
	
	//split into channels
	std::vector<cv::Mat> hlsChannels(3);
	cv::split(hlsInput, hlsChannels);
	
	cv::Mat H, L, S;
	hlsChannels[0].convertTo(H, CV_32F);
	hlsChannels[1].convertTo(L, CV_32F);
	hlsChannels[2].convertTo(S, CV_32F);

	//normalize the channels
	H = H * 2.0f;         //H now in range 0-360
	L = L / 255.0f;       //L now in range 0-1
	S = S / 255.0f;       //S now in range 0-1

	//find highlighted regions where S >= mu & L >= nu
	double sumLS = 0.0;
	int countLS = 0;
	cv::Mat mask = cv::Mat::zeros(input.size(), CV_8U);
	for(int i = 0; i < input.rows; i++)
	{
		for(int j = 0; j < input.cols; j++)
		{
			if(S.at<float>(i,j) >= mu && L.at<float>(i,j) >= nu)
			{
				mask.at<uchar>(i,j) = 255;
				sumLS += L.at<float>(i,j) * S.at<float>(i,j);
				countLS++;
			}
		}
	}
	double averageLS = 0.0;
	if(countLS > 0)
	{
		averageLS = sumLS / static_cast<double>(countLS);
	}

	//apply equation (4) to get L_constr
	cv::Mat L_constr = cv::Mat::zeros(input.size(), CV_32F);
	double phiRad = phi * CV_PI / 180.0;   //phi in radians

	for(int i = 0; i < input.rows; i++)
	{
		for(int j = 0; j < input.cols; j++)
		{
			double hueRad = (H.at<float>(i,j) * kappa) * (CV_PI / 180.0); //hue in radians
			double cos = std::cos(hueRad - phiRad);
			if(mask.at<int>(i,j) == 255)
			{
				//highlighted region
				double tmp = L.at<float>(i,j) + gamma * averageLS * cos;
				L_constr.at<float>(i,j) = static_cast<float>(tmp);
			}
			else
			{
				double tmp = L.at<float>(i,j) * (1.0 + gamma * cos * S.at<float>(i,j));
				L_constr.at<float>(i,j) = static_cast<float>(tmp);
			}
		}
	}

	//normalize L_constr -> L_res [0...1], multiply by chi
	double minV, maxV;
	cv::minMaxLoc(L_constr, &minV, &maxV);
	double range = (maxV - minV) > 1e-9 ? (maxV - minV) : 1e-9;

	cv::Mat L_res = cv::Mat::zeros(input.size(), CV_32F);
	for(int i = 0; i < input.rows; i++)
	{
		for(int j = 0; j < input.cols; j++)
		{
			double normalized = (L_constr.at<float>(i,j) - minV) / range;
			normalized *= chi;
			L_res.at<float>(i,j) = static_cast<float>(normalized);
		}
	}

	//restrict intesity to [min(R,G,B), max(R,G,B)] for each pixel
	cv::Mat output = cv::Mat::zeros(input.size(), CV_32F);
	for(int i = 0; i < input.rows; i++)
	{
		for(int j = 0; j < input.cols; j++)
		{
			//how do i select specific channel value at i,j?
			unsigned char Bval = input8U.at<cv::Vec3b>(i,j)[0];
			unsigned char Gval = input8U.at<cv::Vec3b>(i,j)[1];
			unsigned char Rval = input8U.at<cv::Vec3b>(i,j)[2];

			float minRGB = std::min({Bval, Gval, Rval}) / 255.0f;
			float maxRGB = std::max({Bval, Gval, Rval}) / 255.0f;

			float val = L_res.at<float>(i,j);
			//clamp
			if(val < minRGB) val = minRGB;
			if(val > maxRGB) val = maxRGB;

			output.at<float>(i,j) = val;
		}
	}

	//blend with original, equation (6)
	cv::Mat finalGray = cv::Mat::zeros(input.size(), CV_32F);
	for(int i = 0; i < input.rows; i++)
	{
		for(int j = 0; j < input.cols; j++)
		{
			double blended = (output.at<float>(i,j) + eta * L.at<float>(i,j)) / (1.0 + eta);
			finalGray.at<float>(i,j) = static_cast<float>(blended);
		}
	}
	return finalGray;
}
/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOAncuti11::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	cv::Mat input(height, width, CV_64FC3);
	for(int i = 0; i < height; i++)
	{
		for(int j = 0; j < width; j++)
		{
			double R = *pSourceData++;
			double G = *pSourceData++;
			double B = *pSourceData++;
			input.at<cv::Vec3d>(i, j)[0] = B; // to create BGR input matrix for opencv operations
			input.at<cv::Vec3d>(i, j)[1] = G;
			input.at<cv::Vec3d>(i, j)[2] = R;
		}
	}
	double eta = 0.2; //blend factor for image from paper
	cv::Mat result = decolorization(input, eta);

	double pY, px, py;

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			float val = result.at<float>(j, i);
			*pDestinationData++ = val;
			*pDestinationData++ = val;
			*pDestinationData++ = val;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	return 0;
}

int TMOAncuti11::TransformVideo()
{
	cv::VideoCapture vid = vSrc->getVideoCaptureObject();
	cv::Mat currentFrame, output;
	std::vector<cv::Mat> channels(3);
	double eta = 1.1; //blend factor for video from paper
	for(int i = 0; i < vSrc->GetTotalNumberOfFrames(); i++)
	{
		vSrc->GetMatVideoFrame(vid, i, currentFrame);
		cv::Mat result = decolorization(currentFrame, eta);
		
		channels[0] = result;
		channels[1] = result;
		channels[2] = result;
		
		cv::merge(channels, output);
		vDst->setMatFrame(vDst->getVideoWriterObject(), output);
		fprintf(stderr, "\rFrame %d/%d processed", i, vSrc->GetTotalNumberOfFrames());
		fflush(stdout);
	}
	fprintf(stderr, "\n");
	return 0;
}
