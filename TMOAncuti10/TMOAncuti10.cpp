/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       Diploma thesis                                         *
*                       Author: Matus Bicanovsky                               *
*                       Brno 2024                                              *
*                                                                              *
*                       Implementation of the TMOAncuti10 class                *
*                       Image and Video Decolorization by Fusion               *
*******************************************************************************/
/* --------------------------------------------------------------------------- *
 * TMOAncuti10.cpp: implementation of the TMOAncuti10 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOAncuti10.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOAncuti10::TMOAncuti10()
{
	SetName(L"Ancuti10");					 
	SetDescription(L"Color to grayscale operator for images and video, method from paper: Image and Video Decolorization by Fusion"); 
}

TMOAncuti10::~TMOAncuti10()
{
}


//function to apply separable binomial kernel
void TMOAncuti10::applySeparableBinomialKernel(double* input, double* output, int width, int height)
{
	double kernel[5] = {1.0, 4.0, 6.0, 4.0, 1.0};    // 1D kernel
	double factor = 1.0 / 16.0;                      // normalization factor

	double* tmp = new double[width * height];
	//apply horizontal kernel
	for(int y = 0; y < height; y++){
		for(int x = 0; x < width; x++)
		{
			double sum = 0.0;
			for(int k = -2; k <= 2; k++){
				int idx = std::min(std::max(x + k, 0), width - 1);
				sum += input[y * width + idx] * kernel[k + 2];
			}
			tmp[y * width + x] = sum * factor;
		}
	}
	//apply vertical kernel
	for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double sum = 0.0;
            for (int k = -2; k <= 2; k++) {
                int idy = std::min(std::max(y + k, 0), height - 1);
                sum += tmp[idy * width + x] * kernel[k + 2];
            }
            output[y * width + x] = sum * factor;
        }
    }
	delete[] tmp;
}

//function to compute saliency map, inputs are channel data, width and height and output is saliency map
void TMOAncuti10::computeSaliencyMap(double* channel, double* saliencyMap, int width, int height)
{
	double meanValue = calculateAPV(channel, width, height);             // mean value of the channel
	double* blurredChannel = new double[width * height];                  
	applySeparableBinomialKernel(channel, blurredChannel, width, height);
	//compute map
	for (int i = 0; i < width * height; i++)
	{
		saliencyMap[i] = fabs(meanValue - blurredChannel[i]);			//paper equation (2)
	}
	delete[] blurredChannel;
}

//function to compute exposedness map, inputs are channel data, width and height and output is exposedness map
void TMOAncuti10::computeExposednessMap(double* channel, double* exposednessMap, int width, int height)
{
	double sigma = 0.25;
	double sigma2 = 2 * sigma * sigma;

	for(int i = 0; i < width * height; i++)
	{
		exposednessMap[i] = exp(- ((channel[i] - 0.5) * (channel[i] - 0.5)/(sigma2)));  			//formula from the paper, equation (3)
	}
}

//function to compute chromatic map, inputs are channel data, chromatic map, saturation values, width and height
void TMOAncuti10::computeChromaticMap(double* channel, double* chromaticMap, double* saturation, int width, int height)
{
	int totalCount = width * height;
	//compute chromatic map
	for (int i = 0; i < totalCount; i++)
	{
		chromaticMap[i] = fabs(channel[i] - saturation[i]) / sqrt(2.0);         
	}
}

//function to compute Laplacian pyramid, inputs are input image, number of levels and output is Laplacian pyramid
void TMOAncuti10::computeLaplacianPyramid(const cv::Mat& input, std::vector<cv::Mat>& pyramid, int levels)
{
	pyramid.clear();
	cv::Mat current = input.clone();
	for(int i = 0; i < levels; i++)
	{
		cv::Mat down, up, laplacian;                 //temporary images
		cv::pyrDown(current, down);                  //downsample the image
		cv::pyrUp(down, up, current.size());         //upsample the downsampled image
		laplacian = current - up;                    //compute the Laplacian
		pyramid.push_back(laplacian);                
		current = down;
	}
	pyramid.push_back(current); 		//add the smallest level
}

//computation of Gaussian pyramid, inputs are input frame, number of levels and output is Gaussian pyramid
void TMOAncuti10::computeGaussianPyramid(const cv::Mat& input, std::vector<cv::Mat>& pyramid, int levels)
{
	pyramid.clear();
	cv::Mat current = input.clone();
	for(int i = 0; i < levels; i++)
	{
		pyramid.push_back(current);
		cv::pyrDown(current, current);
	}
	pyramid.push_back(current); 		//add the smallest level
}

//function to normalize weight maps, input is vector of weight maps
void TMOAncuti10::normalizeWeightMaps(std::vector<cv::Mat>& weightMaps)
{
	cv::Mat sum = cv::Mat::zeros(weightMaps[0].size(), weightMaps[0].type());
	for(const auto& weightMap : weightMaps){
		sum += weightMap;
	}
	for(auto& weightMap : weightMaps){
		weightMap /= sum;
	}
}

//function to fuse Laplacian and Gaussian pyramids to form final pyramid
void TMOAncuti10::fusePyramids(const std::vector<std::vector<cv::Mat>>& laplacianPyramids, const std::vector<std::vector<cv::Mat>>& gaussianPyramids, std::vector<cv::Mat>& fusedPyramid)
{
	int levels = laplacianPyramids[0].size();					//number of levels
	fusedPyramid.resize(levels);								//resize the output pyramid
	for(int i = 0; i < levels; i++){
		fusedPyramid[i] = cv::Mat::zeros(laplacianPyramids[0][i].size(), laplacianPyramids[0][i].type());			//initialize the output pyramid
		for(int k = 0; k < laplacianPyramids.size(); k++){
			fusedPyramid[i] += gaussianPyramids[k][i].mul(laplacianPyramids[k][i]);					//fuse the pyramids according to the formula in paper, equation (5)
		}
	}
}

//function for reconstruction of image from pyramid
cv::Mat TMOAncuti10::reconstructFromPyramid(const std::vector<cv::Mat>& pyramid)
{
	cv::Mat current = pyramid.back();					//initialize the current image
	for(int i = pyramid.size() - 2; i >= 0; i--){			
		cv::Mat up;
		cv::pyrUp(current, up, pyramid[i].size());			//upsample the current image
		current = up + pyramid[i];							//add the upsampled image to the current image
	}
	return current;
}

//main function of the c2g operator
int TMOAncuti10::Transform()
{

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												
	int width = pSrc->GetWidth();
    int height = pSrc->GetHeight();
	double *R_Channel = new double[width * height];
    double *G_Channel = new double[width * height];
    double *B_Channel = new double[width * height];
    double *HKL_Channel = new double[width * height];

	double *pSourceDataIter = pSourceData;
	//extract the RGB channels
	for(int i = 0; i < width * height; i++)
	{
		R_Channel[i] = *pSourceDataIter++;
		G_Channel[i] = *pSourceDataIter++;
		B_Channel[i] = *pSourceDataIter++;
	}
	//compute HKL channel
	for (int i = 0; i < width * height; i++)
	{
		//firstly we need to convert RGB intput into CIELch color space, that is done as RGB -> XYZ -> CIELab -> CIELch
		double X, Y, Z, L, a, b, C, h;
		convertRGBtoXYZ(R_Channel[i], G_Channel[i], B_Channel[i], X, Y, Z);
		convertXYZtoCIELAB(X, Y, Z, L, a, b);
		convertCIELABtoCIELCh(L, a, b, C, h);
		//then we compute the L_hk channel frame
		HKL_Channel[i] = computeL_HK(L, C, h);
	}
	double *saliencyMapR = new double[width * height];
    double *saliencyMapG = new double[width * height];
    double *saliencyMapB = new double[width * height];
    double *saliencyMapHKL = new double[width * height];
	//compute saliency maps
	computeSaliencyMap(R_Channel, saliencyMapR, width, height);
    computeSaliencyMap(G_Channel, saliencyMapG, width, height);
    computeSaliencyMap(B_Channel, saliencyMapB, width, height);
    computeSaliencyMap(HKL_Channel, saliencyMapHKL, width, height);

	double *exposednessMapR = new double[width * height];
    double *exposednessMapG = new double[width * height];
    double *exposednessMapB = new double[width * height];
    double *exposednessMapHKL = new double[width * height];
	//compute exposedness maps
	computeExposednessMap(R_Channel, exposednessMapR, width, height);
    computeExposednessMap(G_Channel, exposednessMapG, width, height);
    computeExposednessMap(B_Channel, exposednessMapB, width, height);
    computeExposednessMap(HKL_Channel, exposednessMapHKL, width, height);

	double *chromaticMapR = new double[width * height];
    double *chromaticMapG = new double[width * height];
    double *chromaticMapB = new double[width * height];
    double *chromaticMapHKL = new double[width * height];

	//compute saturation values
	double *saturation = new double[width * height];
    for (int i = 0; i < width * height; i++)
    {
        double H, S, L;
        convertRGBtoHSL(R_Channel[i], G_Channel[i], B_Channel[i], H, S, L);
        saturation[i] = S;
    }
	//compute chromatic maps
	computeChromaticMap(R_Channel, chromaticMapR, saturation, width, height);
    computeChromaticMap(G_Channel, chromaticMapG, saturation, width, height);
    computeChromaticMap(B_Channel, chromaticMapB, saturation, width, height);
    computeChromaticMap(HKL_Channel, chromaticMapHKL, saturation, width, height);

	//normalize the weight maps
    std::vector<cv::Mat> weightMapsR = {
        cv::Mat(height, width, CV_64F, saliencyMapR),
        cv::Mat(height, width, CV_64F, exposednessMapR),
        cv::Mat(height, width, CV_64F, chromaticMapR)
    };
    std::vector<cv::Mat> weightMapsG = {
        cv::Mat(height, width, CV_64F, saliencyMapG),
        cv::Mat(height, width, CV_64F, exposednessMapG),
        cv::Mat(height, width, CV_64F, chromaticMapG)
    };
    std::vector<cv::Mat> weightMapsB = {
        cv::Mat(height, width, CV_64F, saliencyMapB),
        cv::Mat(height, width, CV_64F, exposednessMapB),
        cv::Mat(height, width, CV_64F, chromaticMapB)
    };
    std::vector<cv::Mat> weightMapsHKL = {
        cv::Mat(height, width, CV_64F, saliencyMapHKL),
        cv::Mat(height, width, CV_64F, exposednessMapHKL),
        cv::Mat(height, width, CV_64F, chromaticMapHKL)
    };
    normalizeWeightMaps(weightMapsR);
	normalizeWeightMaps(weightMapsG);
    normalizeWeightMaps(weightMapsB);
    normalizeWeightMaps(weightMapsHKL);

	//combine the normalized maps to form the final weight map for each channel
    cv::Mat finalWeightMapR = weightMapsR[0] + weightMapsR[1] + weightMapsR[2];
    cv::Mat finalWeightMapG = weightMapsG[0] + weightMapsG[1] + weightMapsG[2];
    cv::Mat finalWeightMapB = weightMapsB[0] + weightMapsB[1] + weightMapsB[2];
    cv::Mat finalWeightMapHKL = weightMapsHKL[0] + weightMapsHKL[1] + weightMapsHKL[2];

	//compute Laplacian pyramids for each input channel
    std::vector<cv::Mat> inputs = {
        cv::Mat(height, width, CV_64F, R_Channel),
        cv::Mat(height, width, CV_64F, G_Channel),
        cv::Mat(height, width, CV_64F, B_Channel),
        cv::Mat(height, width, CV_64F, HKL_Channel)
    };
	int levels = 5; // number of pyramid levels
	std::vector<std::vector<cv::Mat>> laplacianPyramids(inputs.size());
	for (int i = 0; i < inputs.size(); i++) {
        computeLaplacianPyramid(inputs[i], laplacianPyramids[i], levels);
    }
	fprintf(stderr, "Laplacian pyramids computed\n");
	
    std::vector<cv::Mat> finalWeightMaps = {
        finalWeightMapR,
        finalWeightMapG,
        finalWeightMapB,
        finalWeightMapHKL
    };
	//compute Gaussian pyramids for each normalized weight map
    std::vector<std::vector<cv::Mat>> gaussianPyramids(finalWeightMaps.size());
    for (int i = 0; i < finalWeightMaps.size(); i++) {
        computeGaussianPyramid(finalWeightMaps[i], gaussianPyramids[i], levels);
    }
	fprintf(stderr, "Gaussian pyramids computed\n");
	//fuse the pyramids
    std::vector<cv::Mat> fusedPyramid;
    fusePyramids(laplacianPyramids, gaussianPyramids, fusedPyramid);
	//reconstruct the final image from the fused pyramid
	fprintf(stderr, "Pyramids fused\n");
    cv::Mat fusedImage = reconstructFromPyramid(fusedPyramid);
	fprintf(stderr, "Image reconstructed\n");
	cv::normalize(fusedImage, fusedImage, 0, 1, cv::NORM_MINMAX);
	fprintf(stderr, "Reconstruction done.\n");

	
	double pY, px, py;
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			
			double pixel = fusedImage.at<double>(j, i);
			*pDestinationData++ = pixel;
			*pDestinationData++ = pixel;
			*pDestinationData++ = pixel;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	
	// Clean up
    delete[] R_Channel;
    delete[] G_Channel;
    delete[] B_Channel;
    delete[] HKL_Channel;
    delete[] saliencyMapR;
    delete[] saliencyMapG;
    delete[] saliencyMapB;
    delete[] saliencyMapHKL;
    delete[] exposednessMapR;
    delete[] exposednessMapG;
    delete[] exposednessMapB;
    delete[] exposednessMapHKL;
    delete[] chromaticMapR;
    delete[] chromaticMapG;
    delete[] chromaticMapB;
    delete[] chromaticMapHKL;
    delete[] saturation;
	return 0;
}

//function to conver RGB to HSL color space
void TMOAncuti10::convertRGBtoHSL(double R, double G, double B, double &H, double &S, double &L)
{
	//find the maximum and minimum values of R, G, B
	double max = std::max({R, G, B});
	double min = std::min({R, G, B});
	//calculate luminance
	L = (max + min) / 2.0;
	if(max == min){
		H = S = 0.0; 			//achromatic
	}
	else{
		double tmp = max - min;
		S = (L > 0.5) ? tmp / (2.0 - max - min) : tmp / (max + min);
		double del_R = (((max - R) / 6.0) + (tmp / 2.0)) / tmp;
		double del_G = (((max - G) / 6.0) + (tmp / 2.0)) / tmp;
		double del_B = (((max - B) / 6.0) + (tmp / 2.0)) / tmp;
		if(R == max)
			H = del_B - del_G;
		else if(G == max)
			H = (1.0 / 3.0) + del_R - del_B;
		else if(B == max)
			H = (2.0 / 3.0) + del_G - del_R;
		if(H < 0.0)
			H += 1.0;
		if(H > 1.0)
			H -= 1.0;
	}
}
//function to convert RGB to XYZ color space
void TMOAncuti10::convertRGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z)
{

	// //apply gamma correction
	R = (R/255.0 > 0.04045) ? pow((R + 0.055) / 1.055, 2.4) : R / 12.92;
    G = (G/255.0 > 0.04045) ? pow((G + 0.055) / 1.055, 2.4) : G / 12.92;
    B = (B/255.0 > 0.04045) ? pow((B + 0.055) / 1.055, 2.4) : B / 12.92;
	
	//convert RGB to XYZ using D65 illuminant
	X = R * 0.4124564 + G * 0.3575761 + B * 0.1804375;
    Y = R * 0.2126729 + G * 0.7151522 + B * 0.0721750;
    Z = R * 0.0193339 + G * 0.1191920 + B * 0.9503041;
}
//function to convert XYZ to CIELAB color space
void TMOAncuti10::convertXYZtoCIELAB(double X, double Y, double Z, double &L, double &a, double &b)
{
	//reference white points
	const double ref_X = 95.047;
	const double ref_Y = 100.000;
	const double ref_Z = 108.883;
	//normalize XYZ values
	X /= ref_X;
	Y /= ref_Y;
	Z /= ref_Z;
	
	X = (X > 0.008856) ? pow(X, 1.0 / 3.0) : (7.787 * X) + (16.0 / 116.0);
    Y = (Y > 0.008856) ? pow(Y, 1.0 / 3.0) : (7.787 * Y) + (16.0 / 116.0);
    Z = (Z > 0.008856) ? pow(Z, 1.0 / 3.0) : (7.787 * Z) + (16.0 / 116.0);

	L = (116.0 * Y) - 16.0;
	a = 500.0 * (X - Y);
	b = 200.0 * (Y - Z);
}
//function to convert CIELAB to CIELCh color space
void TMOAncuti10::convertCIELABtoCIELCh(double L, double a, double b, double &C, double &h)
{
	C = sqrt(a * a + b * b);
	double var_h = atan2(b, a);
	if (var_h > 0)
		var_h = (var_h / M_PI) * 180.0;
	else
		var_h = 360 - (abs(var_h) / M_PI) * 180.0;
	h = var_h;
}
//function to compute L_HK value for given pixel
double TMOAncuti10::computeL_HK(double L, double C, double H)
{
	return L + (2.5 - 0.025 * L)*(0.116 * fabs(sin((H-90)/2)) + 0.085) * C;               //computation of L_hk value for one pixel based on equation (1)
}

//function to calculate average pixel value
double TMOAncuti10::calculateAPV(double* channel, int width, int height)
{
	double sum = 0.0;
	int totalCount = width * height;
	for (int i = 0; i < totalCount; i++)
	{
		sum += channel[i];
	}
	return sum / totalCount;
}
