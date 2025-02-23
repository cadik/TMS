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
////////////////////////////////////////////////////////////////////////////////////////////////////////
//                         Itti-koch saliency model part                                           /////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// function for computing feature maps using rules: c in {2,3,4}, s = c+3 or c+4 as in the paper, c is center scale, s is surround scale
// the surround is upsampled to the size of center scale and then absolute difference is computed
// input is vector of gaussian pyramid levels and output is vector of scaleDifference structures, each with difference map, c and s values
std::vector<scaleDifference> TMOAncuti11::mapsDifference(std::vector<cv::Mat>& input)
{
	//according to paper c is [2,3,4], s = c+3 or c+4 so delta is [3,4]
	std::vector<scaleDifference> result;
	for(int c = 2; c <= 4; c++)
	{
		for(int delta = 3; delta <= 4; delta++)
		{
			int s = c + delta;
			if(s >= input.size()) break;
			cv::Mat center = input[c];
			cv::Mat surround = input[s];
			//upsample surround to match center size and compute absolute difference
			cv::Mat upSurround;
			cv::resize(surround, upSurround, center.size(), 0, 0, cv::INTER_LINEAR);
			//abs diff
			cv::Mat diff;
			cv::absdiff(center, upSurround, diff);
			scaleDifference tmp;
			tmp.diff = diff;               //storing difference map
			tmp.c = c;
			tmp.s = s;
			result.push_back(tmp);
		}
	}
	return result;
}

// calculating intensity channel according to formula from paper, I = (R + G + B) / 3
// input is BGR image and output is single channel intensity image
cv::Mat TMOAncuti11::intensityChannel(cv::Mat& input)
{
	cv::Mat tmp;
	std::vector<cv::Mat> channels(3);
	cv::split(input, channels);
	tmp = (channels[0] + channels[1] + channels[2]) / 3.0f;
	return tmp;
}

// building 9-level Gaussian pyramid for input image using openCV pyrDwon function
// input is the base image and levels which is this case is 9 according to the paper, output is vector of pyramid levels
std::vector<cv::Mat> TMOAncuti11::buildPyramid(cv::Mat& input, int levels=9)
{
	std::vector<cv::Mat> pyramid(levels);
	pyramid[0] = input.clone();               //scale 0 = full sized image
	for(int i = 1; i < levels; i++)
	{
		cv::Mat down;
		cv::pyrDown(pyramid[i-1], down);
		pyramid[i] = down;
	}
	return pyramid;
}

//implementation of normalization operator N(.) from the paper, steps:
// 1. rescale the map so global max is equal to N_scale variable
// 2. find local maxima in scaled map, then compute the average amplitude of those local maxima
// 3. normalize whole map by formula from paper -> (N_scale - avg)^2
// input is single channel map, output is normalized map
cv::Mat TMOAncuti11::ittiNormalize(cv::Mat& input)
{
	float N_scale = 8.f;  //scaling factor from paper
	double minVal, maxVal;
	cv::minMaxLoc(input, &minVal, &maxVal);
	if(maxVal < 1e-6) return input.clone();   // if the map is empty
	// step 1. scale so global max is N_scale
	float alpha = N_scale / static_cast<float>(maxVal);
	cv::Mat scaled;
	input.convertTo(scaled, CV_32F, alpha, 0.f);

	//find local maxima
	int w = input.cols;
	int h = input.rows;
	// kernel size for dilitation
	int kernelW = std::max(w/32, 1);
	int kernelH = std::max(h/32, 1);

	cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(kernelW, kernelH));
	cv::Mat dilated;
	cv::dilate(scaled, dilated, kernel);
	cv::Mat maxima = (scaled == dilated);  // maxima[i,j] = 255 if local maxima, 0 otherwise, as in paper only local maxima are considered

	// find amount of local maxima and their sum
	double maxSum = 0.0;
	int maxCount = 0;
	for(int i = 0; i < h; i++)
	{
		for(int j = 0; j < w; j++)
		{
			if(maxima.at<uchar>(i,j) == 255)
			{
				maxSum += scaled.at<float>(i,j);
				maxCount++;
			}
		}
	}
	if(maxCount < 1) return scaled;     //no local maxima found

	float avg = static_cast<float>(maxSum / static_cast<float>(maxCount));
	float tmp = (N_scale - avg);
	tmp *= tmp;   // (N_scale - avg)^2
	
	cv::Mat output;
	scaled.convertTo(output, CV_32F, tmp, 0.f);    //multiplying entire map by calculated factor
	return output;
}

//function sums multiple maps (results of mapsDifference function) for given channel and create resulting conspicuity map
cv::Mat TMOAncuti11::sumNormalized(std::vector<scaleDifference>& input, cv::Size size)
{
	//choosing the smallest center scale among diffs
	int minC = 99;
	for(auto &tmp : input)
	{
		if(tmp.c < minC) minC = tmp.c;
	}
	// final size should be scale(minC)
	cv::Size finalSize = size;
	for(int i = 0; i<minC; i++)
	{
		finalSize.width = (finalSize.width + 1) / 2;
		finalSize.height = (finalSize.height + 1) / 2;
	}
	cv::Mat sumMaps = cv::Mat::zeros(finalSize, CV_32F);
	for(auto &tmp : input)
	{
		//resize to finalSize
		cv::Mat resized;
		cv::resize(tmp.diff, resized, finalSize, 0, 0, cv::INTER_LINEAR);
		//normalize
		cv::Mat normalized = ittiNormalize(resized);
		//add to sumMaps
		sumMaps += normalized;
	}
	//final normalization of sum, creating conspicuity map for given channel
	sumMaps = ittiNormalize(sumMaps);
	return sumMaps;
}

//function builds four color opponency channgels (red, green, blue, yellow) from input image
// with formulas from paper: R' = R - (G + B) / 2, 
// 							 G' = G - (R + B) / 2, 
// 							 B' = B - (R + G) / 2, 
// 							 Y' = (R + G) / 2 - |R - G| / 2 - B
// also thresholding low intensities as in paper I < 1/10 of max(I)
void TMOAncuti11::computeColorChannels(cv::Mat& input, cv::Mat& Rp, cv::Mat& Gp, cv::Mat& Bp, cv::Mat& Yp)
{
	//input to float [0..1]
	cv::Mat inputF;
	input.convertTo(inputF, CV_32F, 1.0/255.0);
	//split into channels
	std::vector<cv::Mat> channels(3);
	cv::split(inputF, channels);
	cv::Mat R = channels[2];
	cv::Mat G = channels[1];
	cv::Mat B = channels[0];

	//in paper threshold low intensities => 0 if I < 1/10 of max(I)
	cv::Mat I = (R + G + B) / 3.0f;
	double maxVal;
	cv::minMaxLoc(I, nullptr, &maxVal);
	float threshold = static_cast<float>(maxVal / 10.f);
	//clamp the R,G,B matrixes to 0 if I < threshold
	cv::threshold(R, R, threshold, 0.f, cv::THRESH_TOZERO);
	cv::threshold(G, G, threshold, 0.f, cv::THRESH_TOZERO);
	cv::threshold(B, B, threshold, 0.f, cv::THRESH_TOZERO);
	// calculate color opponency channels based on formulas from paper
	Rp = R - (G + B) / 2.0f;
	Gp = G - (R + B) / 2.0f;
	Bp = B - (R + G) / 2.0f;
	Yp = (R + G) / 2.0f - cv::abs(R - G) / 2.0f - B;

	//clamp negative values to 0 as in paper
	cv::threshold(Rp, Rp, 0.f, 0.f, cv::THRESH_TOZERO);
	cv::threshold(Gp, Gp, 0.f, 0.f, cv::THRESH_TOZERO);
	cv::threshold(Bp, Bp, 0.f, 0.f, cv::THRESH_TOZERO);
	cv::threshold(Yp, Yp, 0.f, 0.f, cv::THRESH_TOZERO);
}

// function for applying gabor filter to input image with given angle to detect local orientation
// input angle is in degrees, output is filtered image
cv::Mat TMOAncuti11::gaborFilter(cv::Mat& input, float angle)
{
	int kernelSize = 11;   //kernel size
	float sigma = 4.0f;    //standard deviation
	float lambda = 10.0f; //wavelength
	float gamma = 0.5f;  //aspect ratio
	float psi = 0.0f;   //phase offset

	cv::Mat kernel(kernelSize, kernelSize, CV_32F);
	int halfSize = kernelSize / 2;
	for(int i = 0; i < kernelSize; i++)
	{
		for(int j = 0; j < kernelSize; j++)
		{
			float x = j - halfSize;
			float y = halfSize - i;
			float xprime = x * std::cos(angle) + y * std::sin(angle);
			float yprime = -x * std::sin(angle) + y * std::cos(angle);
			float val = std::exp(-(xprime*xprime + gamma*gamma*yprime*yprime)/(2.f*sigma*sigma))
						* std::cos(2.f*(float)CV_PI*xprime/lambda + psi);
			kernel.at<float>(i,j) = val;
		}
	}
	cv::Mat result;
	cv::filter2D(input, result, CV_32F, kernel);
	return result;
}

// main function for computing saliency map using Itti-Koch model
// if input is grayscale we discard color channels computation (this is specific for usage in decolorization algorithm)
// this function realizes main steps of model:
// create intensity channel -> create its difference maps -> sum and normalize them -> resulting intensityNorm
// compute color channels -> create RG and BY difference maps -> sum and normalize them -> resulting Cmap
// compute orientation maps for angles {0, 45, 90, 135} -> difference maps -> sum and normalize them -> resulting orientationMap
// create final saliency map as average of intensityNorm, Cmap and orientationMap and lastly upsample it to original size
cv::Mat TMOAncuti11::computeSaliencyMap(cv::Mat &input, bool color)
{	
	//get intesity
	cv::Mat intensity;
	if(color)
	{
		intensity = intensityChannel(input);
	}
	else{
		input.convertTo(intensity, CV_32F, 1.0/255.0);
	}
	auto intensityPyramid = buildPyramid(intensity, 9);       //build pyramid
	auto intensityDiff = mapsDifference(intensityPyramid);    //compute difference maps
	cv::Mat intensityNorm = sumNormalized(intensityDiff, input.size());  //sum and normalize maps

	cv::Mat Cmap;
	if(color)
	{
		//color channels
		cv::Mat Rp, Gp, Bp, Yp;
		computeColorChannels(input, Rp, Gp, Bp, Yp);

		//building pyramids for RG channel (R' - G') and BY channel (B' - Y'), computing difference maps, summing and normalizing them
		cv::Mat RG, BY;
		cv::absdiff(Rp, Gp, RG);
		cv::absdiff(Bp, Yp, BY);
		auto RGPyramid = buildPyramid(RG, 9);
		auto RGDiff = mapsDifference(RGPyramid);
		cv::Mat RGnorm = sumNormalized(RGDiff, input.size());
	
		auto BYPyramid = buildPyramid(BY, 9);
		auto BYDiff = mapsDifference(BYPyramid);
		cv::Mat BYnorm = sumNormalized(BYDiff, input.size());
	
		Cmap = RGnorm + BYnorm;
	}
	//calculate orientation maps for each angle in {0, 45, 90, 135}
	std::vector<float> angles = {0.0, 45.0, 90.0, 135.0};
	cv::Mat orientationSum;
	bool first = true;
	for(float angle : angles)
	{
		cv::Mat filtered = gaborFilter(intensity, angle);
		
		auto orientationPyramid = buildPyramid(filtered, 9);
		auto orientationDiff = mapsDifference(orientationPyramid);
		cv::Mat orientationNorm = sumNormalized(orientationDiff, input.size());
		if(first)
		{
			first = false;
			orientationSum = orientationNorm.clone();
		}
		else{
			orientationSum += orientationNorm;
		}
	}
	//final saliency is average of all normalized maps intensityNorm, Cmap, orientationMap
	cv::Mat saliencyMap;
	if(color)
	{
		saliencyMap = (intensityNorm + Cmap + orientationSum) / 3.0f;
	}
	else
	{
		saliencyMap = (intensityNorm + orientationSum) / 2.0f;
	}
	
	cv::Mat result;
	cv::resize(saliencyMap, result, input.size(), 0, 0, cv::INTER_LINEAR); // upsample to original size
	return result;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/////                        end of Itti-koch saliency model part                              /////
////////////////////////////////////////////////////////////////////////////////////////////////////


//function for finding most salient points in saliency map, returns topN points with highest saliency value
std::vector<cv::Point> TMOAncuti11::findSalientPoint(cv::Mat &saliencyMap, int topN, int radius)
{
	std::vector<std::pair<float, cv::Point>> allPoints;
	for(int i = 0; i < saliencyMap.rows; i++)
	{
		for(int j = 0; j < saliencyMap.cols; j++)
		{
			float val = saliencyMap.at<float>(i,j);
			allPoints.push_back(std::make_pair(val, cv::Point(j,i)));
		}
	}
	//sort by saliency value
	std::sort(allPoints.begin(), allPoints.end(), [](auto &a, auto &b){
		return a.first > b.first; 
	});

	std::vector<cv::Point> result;
	for(int i = 0; i < topN; i++)
	{
		result.push_back(allPoints[i].second);
	}
	return result;
}
//function for computing average saliency in region around center point with given radius
double TMOAncuti11::averageSaliencyInRegion(cv::Mat& saliencyMap, cv::Point center, int radius)
{
	double sum = 0.0;
	int cnt = 0;
	for(int dy = -radius; dy <= radius; dy++)
	{
		for(int dx = -radius; dx <= radius; dx++)
		{
			int xx = center.x + dx;
			int yy = center.y + dy;
			//check for out of bounds
			if(xx >= 0 && xx < saliencyMap.cols && yy >= 0 && yy < saliencyMap.rows)
			{
				//check if current (dx,dy) is inside circle
				if(dx*dx + dy*dy <= radius*radius)
				{
					float val = saliencyMap.at<float>(yy, xx);
					sum += val;
					cnt++;
				}
			}
		}
	}
	if(cnt > 0)
	{
		return sum / static_cast<double>(cnt);    //return average saliency in region
	}
	else
	{
		return 0.0;
	}
}

double TMOAncuti11::offsetAngleSelection(cv::Mat &input)
{
	cv::Mat colorInput;
	input.convertTo(colorInput, CV_8UC3, 255.0);
	cv::Mat colorSal = computeSaliencyMap(colorInput, true);    //compute color saliency map
	cv::Mat baselineGray;
	cv::cvtColor(colorInput, baselineGray, cv::COLOR_BGR2GRAY);
	cv::Mat graySal = computeSaliencyMap(baselineGray, false);  //compute saliency map for naive grayscale image
	//save saliency maps
	double minVal, maxVal, grayMin, grayMax;
	cv::minMaxLoc(colorSal, &minVal, &maxVal);
	cv::minMaxLoc(graySal, &grayMin, &grayMax);
	cv::Mat storeColorSal, storeGraySal;
	colorSal.convertTo(storeColorSal, CV_8UC1, 255.0/(maxVal - minVal), -255.0*minVal/(maxVal - minVal));
	graySal.convertTo(storeGraySal, CV_8UC1, 255.0/(grayMax - grayMin), -255.0*grayMin/(grayMax - grayMin));
	cv::imwrite("colorSal.jpg", storeColorSal);
	cv::imwrite("graySal.jpg", storeGraySal);
	int topN = 5;        //amount of top salient regions to consider
	int radius = 25;     //radius for non-maximum suppression
	std::vector<cv::Point> mostSalientPoints = findSalientPoint(colorSal, topN, radius);

	int regionRadius = 50;
	float lostRatio = 0.5;
	std::vector<cv::Point> lostRegions;

	for(int i = 0; i < mostSalientPoints.size(); i++)
	{
		cv::Point p = mostSalientPoints[i];
		double colorAvg = averageSaliencyInRegion(colorSal, p, regionRadius);
		double grayAvg = averageSaliencyInRegion(graySal, p, regionRadius);
		float tmpRatio = (float)((grayAvg + 1e-6)/(colorAvg + 1e-6));
		fprintf(stderr, "Color avg: %f, Gray avg: %f, Ratio: %f Position: x:%d y:%d\n", colorAvg, grayAvg, tmpRatio, p.x, p.y);
		if(tmpRatio < lostRatio)
		{
			lostRegions.push_back(p);
		}
	}
	if(lostRegions.empty()) return 250.0; //default value

	std::vector<double> angleOptions = {200.0, 250.0, 300.0, 320.0, 350.0};
	double bestAngle = 200.0;
	double bestScore = 0.0;
	for(int i = 0; i < angleOptions.size(); i++)
	{
		double angle = angleOptions[i];
		cv::Mat tmp = decolorization(input, 0.2, angle);
		cv::Mat tmpSal = computeSaliencyMap(tmp, false);
		//calculate total saliency in lost regions
		double totalSal = 0.0;
		for(int j = 0; j < lostRegions.size(); j++)
		{
			cv::Point p = lostRegions[j];
			totalSal += averageSaliencyInRegion(tmpSal, p, regionRadius);
		}
		fprintf(stderr, "Angle: %f, Total saliency: %f\n", angle, totalSal);
		if(totalSal > bestScore)
		{
			bestScore = totalSal;
			bestAngle = angle;
		}
	}
	return bestAngle;
}


cv::Mat TMOAncuti11::decolorization(cv::Mat &input, double eta, double phi)
{
	// default values of parameters from the paper
	double mu = 0.1;   //saturation threshold
	double nu = 0.6;   //lightness threshold
	//double phi = 300.0; //offset angle in degrees
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
			double cos = std::cos(hueRad + phiRad);
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
	double phi = offsetAngleSelection(input);
	cv::Mat result = decolorization(input, eta, phi);


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
		cv::Mat result = decolorization(currentFrame, eta, 300);
		
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
