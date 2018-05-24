/* --------------------------------------------------------------------------- *
 * TMOEisemann04.cpp: implementation of the TMOEisemann04 class.   *
 * ---------------------------------------------------------------------------*/
/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*           Semestral project for VYF - Computation Photography course         *
*                   Implementation of method described in article              *
*              Flash Photography Enhancement via Intrinsic Relighting          *
*                      by Elmar Eisemann and Frédo Durand                      *
*                                                                              *
*   Loosely based on previous Java implementation by Anna Ostroukh from 2017   *
*             Author: Roman Jaška [xjaska00 AT stud.fit.vutbr.cz]              *
*                                  Brno 2018                                   *
*                                                                              *
*******************************************************************************/

#include "TMOEisemann04.h"
//#define SAVE_ALL // Saves all intermediary image outputs

 /* ---------------------------------------------------------------------------*
  * Constructor serves for describing a technique and input parameters         *
  * --------------------------------------------------------------------------*/
TMOEisemann04::TMOEisemann04()
{
	SetName(L"Eisemann04");
	SetDescription(L"Flash/No-flash");

	flashImagePathParameter.SetName(L"Flash");
	flashImagePathParameter.SetDescription(L"Input flash image path");
	flashImagePathParameter.SetDefault("");
	this->Register(flashImagePathParameter);

	shadowCorrectionParameter.SetName(L"ShadowCorr");
	shadowCorrectionParameter.SetDescription(L"Enables shadow correction.");
	shadowCorrectionParameter.SetDefault(false);

	this->Register(shadowCorrectionParameter);

	useOptimized();
}

TMOEisemann04::~TMOEisemann04(){}

/**Helper functions **********************************************************/
// Creates a new 3-channel double precision Mat from a TMOImage
Mat toMat64C3(TMOImage* image) {
	return
		Mat(image->GetHeight(),
		image->GetWidth(),
		CV_64FC3,
		image->GetData());
}

// Creates a new 3-channel single precision Mat from an input RGB fp Mat
Mat toLuv(Mat& image) {
	Mat imageLuv;
	image.convertTo(imageLuv, CV_32FC3);
	cvtColor(imageLuv, imageLuv, COLOR_RGB2Luv);
	return imageLuv;
}

// Creates a new 3-channel double precision Mat from an input LUV fp Mat
Mat fromLuv(Mat image) {
	Mat imageRgb;
	cvtColor(image, imageRgb, COLOR_Luv2RGB);
	imageRgb.convertTo(imageRgb, CV_64FC3);
	return imageRgb;
}

// Scales a floating point image to unsigned byte range and saves it as jpg.
void saveImage(string name, Mat image) {
	Mat output;
	image.convertTo(output, CV_8UC3, 255., 0);
	cvtColor(output, output, COLOR_RGB2BGR);
	imwrite(name + ".jpg", output);
}

// Turns single channel mat into 3-channel, duplicating the channels.
Mat to3C(Mat&gray) {
	Mat color;
	merge(vector<Mat>{gray, gray, gray}, color);
	return color;
}

// Finds the max value in an image and sums its occurences
double sumMaxValues(Mat&channel) {
	double min, max, sum;
	minMaxLoc(channel, &min, &max);

	for (int j = 0; j < channel.size().width; j++)
		for (int i = 0; i < channel.size().height; i++)
			if (channel.at<double>(i, j) == max)
				sum += max;
	return sum;
}

// Creates a 2D gaussian kernel of given sigma value
Mat getGaussianKernel2D(double sigma) {
	Mat kernel = getGaussianKernel(sigma * 6 - 1, sigma);
	Mat kernel2d = kernel * kernel.t();

#ifdef SAVE_ALL
	Mat kernelNormalized;
	normalize(kernel2d, kernelNormalized, 0, 1, NORM_MINMAX);
	saveImage("kernel", to3C(kernelNormalized));
#endif
	return kernel2d;
}

// Locates the minimum and maximum values of a masked image
void minMaxLocMask(Mat&image, Mat&mask, double *min, double *max) {
	bool isMinSet = false, isMaxSet = false;
	for (int j = 0; j < image.size().height; j++) {
		for (int i = 0; i < image.size().width; i++) {
			bool maskValue = mask.at<double>(j, i) > 0.888;
			if (maskValue == false) {
				double imagePix = image.at<double>(j, i);
				if (!isMinSet) { *min = imagePix; isMinSet = true; }
				if (!isMaxSet) { *max = imagePix; isMaxSet = true; }

				if (imagePix > *max) {
					*max = imagePix;
				}
				if (imagePix < *min) {
					*min = imagePix;
				}
			}
		}
	}
}
/******************************************************************************/

// Computes the intensity component of an image. Second image may be provided
// to be used for weight values.
Mat computeIntensity(Mat&img, Mat* img2 = NULL) {
	Mat channels[3]; split(img, channels);
	Mat sums = channels[0];
	add(sums, channels[1], sums);
	add(sums, channels[2], sums);

	Mat temp;
	Mat weights = img2 == NULL ? img.clone() : img2->clone();

	divide(img, to3C(sums), temp);
	multiply(temp, weights, temp);

	split(temp, channels);

	Mat intensity = channels[0];
	add(intensity, channels[1], intensity);
	add(intensity, channels[2], intensity);

	intensity.convertTo(intensity, CV_64FC1);
	return intensity;
}

// Computes the large scale image component of the given grayscale image.
// Requires the size of the image's diagonal
// May be provided with a mask to exclude certain areas from filtering
Mat computeLargeScale(Mat&intensity, double diagonal, Mat *mask = NULL) {
	Mat output;
	Mat intensity32;
	intensity.convertTo(intensity32, CV_32F);
	for (int i = 0; i < 1; i++) {
		bilateralFilter(intensity32, output, -1, 1./255.*7, 0.015* diagonal,
			BORDER_DEFAULT);
		intensity32 = output.clone();
	}
	output.convertTo(output, CV_64F);
	if (mask != NULL) {
		Mat maskInv = 1 - (*mask > 0.0);
		maskInv.convertTo(maskInv, CV_64FC1);
		divide(output, maskInv, output);
	}
	return output;
}

// Computes the detail component from the intensity and larg scale components
Mat computeDetail(Mat&intensity, Mat&largeScale) {
	Mat detail;
	divide(intensity, largeScale, detail, CV_64FC1);
	return detail;
}

// Computes the color component from the src. image and its intensity component
Mat computeColor(Mat&image, Mat&intensity) {
	Mat color;
	divide(image, to3C(intensity), color);
	return color;
}

// Calculates the channel weight for white ballance correction
// Needs to be provided with corresponding flash and no flash channels
double getWhiteBallanceChannelWeight(Mat&channelF, Mat&channelNF) {
	return pow(0.4 * sumMaxValues(channelF) + 0.6*sumMaxValues(channelNF), 0.2);
}

// Corrects the white ballance of the output image, using both input images
Mat whiteBalanceCorr(Mat&image, Mat&imageF, Mat&imageNF) {
	Mat channels[3]; split(image, channels);
	Mat channelsF[3]; split(imageF, channelsF);
	Mat channelsNF[3]; split(imageNF, channelsNF);
	for (int i = 0; i < 3; i++) {
		double weight =
			getWhiteBallanceChannelWeight(channelsF[i], channelsNF[i]);
		multiply(channels[i], weight, channels[i]);
	}
	Mat corrected = Mat();
	merge(vector<Mat>{channels[0], channels[1], channels[2]}, corrected);

	normalize(corrected, corrected, 0., 1., NORM_MINMAX, CV_64FC3);
	return corrected;
}

// Calculates the color similarity if a pixel's neighborhood
Mat getColorSimilarity(Vec3f pixel, Mat&image2) {
	Mat output;
	Mat image1 = Mat(image2.size(), CV_32FC3);
	image1 = Scalar(pixel[0], pixel[1], pixel[2]);
	subtract(image2, image1, output);
	pow(output, 2, output);
	Mat channels[3];
	split(output, channels);

	add(channels[0], channels[1], output);
	add(output, channels[2], output);
	Mat delta;
	sqrt(output, delta);

	normalize(delta, delta, 0, 1, NORM_MINMAX);
	subtract(1, delta, delta);
	return delta;
}

// Advanced shadow correction kernel applied to each shadow pixel
// i, j - pixel coordinates
// kernelHalf - kernel radius
// colorF_Luv - Luv version of the color component of flash image
// imageF_Luv - Luv version of the flash image
// imageNF_Luv - Luv version of the non-flash image
// invMask - inverted version of the shadow mask
// kernel - a gaussian kernel
// corrected - the image beign corrected
void localShadowCorrection(int i, int j, int kernelHalf, Mat*colorF_Luv,
	Mat*imageF_Luv, Mat*imageNF_Luv, Mat*invMask, Mat*kernel, Mat*corrected) {

	// Find and apply boundaries
	int miny = max(j - kernelHalf, 0);
	int maxy = min(j + kernelHalf + 1, colorF_Luv->size().height);

	int minx = max(i - kernelHalf, 0);
	int maxx = min(i + kernelHalf + 1, colorF_Luv->size().width);

	int currentWidth = maxx - minx;
	if (currentWidth % 2 == 1) currentWidth--;
	int currentHeight = maxy - miny;
	if (currentHeight % 2 == 1) currentHeight--;
	int kernelLeft = (kernelHalf - currentWidth / 2);
	int kernelTop = (kernelHalf - currentHeight / 2);
	Rect roiK = Rect(kernelLeft, kernelTop, currentWidth, currentHeight);
	if (currentHeight == 0 || currentWidth == 0) return;

	// Get regions of interest of each input component
	Rect roi = Rect(minx, miny, currentWidth, currentHeight);
	Mat colorFroi = (*colorF_Luv)(roi);
	Mat imNFroi = (*imageNF_Luv)(roi);
	Mat maskRoi = (*invMask)(roi);
	Mat kernelRoi = (*kernel)(roiK);
	Mat colorSimilarity =
		getColorSimilarity(imageF_Luv->at<Vec3f>(j, i), imNFroi);

	GaussianBlur(colorSimilarity, colorSimilarity, Size(0, 0), 0.01);

	// Calculate the final value
	Mat weightsRoiWithKernel;
	multiply(colorSimilarity, kernelRoi, weightsRoiWithKernel);
	multiply(weightsRoiWithKernel, maskRoi, weightsRoiWithKernel);

	Scalar weightSum = sum(weightsRoiWithKernel);

	Mat colorFWeighted;
	multiply(colorFroi, to3C(weightsRoiWithKernel), colorFWeighted);

	Scalar valueSum = sum(colorFWeighted);

	corrected->at<Vec3f>(j, i)[1] = valueSum[1] / weightSum[0];
	corrected->at<Vec3f>(j, i)[2] = valueSum[2] / weightSum[0];
}

// Applies the advanced local shadow correction
// colorF_Luv - Luv version of the color component of flash image
// imageF_Luv - Luv version of the flash image
// imageNF_Luv - Luv version of the non-flash image
// shadowMask - mask of the shadows
// diagonal - image's diagonal
Mat shadowCorrection(Mat colorF_Luv, Mat imageF_Luv, Mat imageNF_Luv,
	Mat& shadowMask, double diagonal) {

	// Set up output Mat
	Mat corrected = colorF_Luv.clone();
	// Invert the shadow mask
	Mat invMask = 1 - (shadowMask > 0);
	invMask.convertTo(invMask, CV_32F);
	// Prepare gaussian kernel
	Mat kernel = getGaussianKernel2D(diagonal*0.025);
	int kernelHalf = kernel.size().height / 2;
	kernel.convertTo(kernel, CV_32F);
	// Run the local correction
	int width = colorF_Luv.size().width;
	int height = colorF_Luv.size().height;
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			if (shadowMask.at<double>(j, i) > 0.)
				localShadowCorrection(i, j, kernelHalf, &colorF_Luv,
					&imageF_Luv, &imageNF_Luv, &invMask, &kernel, &corrected);
		}
	}
	return corrected;
}

// Computes the gradient magnitude map of an image
Mat computeGradientImage(Mat&image) {
	Mat temp = Mat(), gradX, gradY, grad;
	image.convertTo(temp, CV_32F);
	Sobel(temp, gradX, CV_32F, 1, 0);
	Sobel(temp, gradY, CV_32F, 0, 1);
	gradX = abs(gradX);
	gradY = abs(gradY);
	magnitude(gradX, gradY, grad);
	threshold(grad, grad, 0.1, 0., THRESH_TOZERO);
	temp.release();
	grad.convertTo(grad, CV_64FC1);
	return grad;
}

// Calculates the histogram for shadow detection
Mat computeHistogram(Mat&image) {
	Mat image32; image.convertTo(image32, CV_32FC1);
	Mat histogram;
	int channels[] = { 0 };
	int histSize = 128;
	float ranges[] = { 0.,1. };
	const float * rang[] = { ranges };
	calcHist(&image32, 1, 0, Mat(), histogram, 1, &histSize, rang, true, false);
	GaussianBlur(histogram, histogram, Size(0, 0), 2);
#ifdef SAVE_ALL
	Mat histOut; normalize(histogram, histOut, 0, 1, NORM_MINMAX, CV_64FC1);
	saveImage("hist", to3C(histOut));
#endif
	return histogram;
}

// Computes the umbra component of shadow
Mat computeUmbra(Mat&intensityF, Mat&intensityNF) {
	// We expect the difference image ∆I between flash and no-flash 
	// to tell how much additional light was received from the flash.
	Mat deltaI;
	subtract(intensityF, intensityNF, deltaI);	
	deltaI = abs(deltaI);

#ifdef SAVE_ALL
	saveImage("deltaI", to3C(deltaI));
#endif

	// We compute the histogram of pixels ∆I. We use 128 bins and smooth it with
	// a Gaussian blur of variance two bins. //We start with a coarse threshold
	// of 0.2 and discard all pixels where ∆I is above this value.
	Mat histogram = computeHistogram(deltaI);				 
	float thresholdDeltaI = 0.2;								
	float thresholdHistIndex = (int)(128 * thresholdDeltaI) - 1;
	for (int i = thresholdHistIndex - 1; i >= 0; i--) {
		if (histogram.at<float>(0, i) >
			histogram.at<float>(0, thresholdHistIndex)) {
			break;
		}
		else
		{
			thresholdHistIndex = i;
		}
	}
	thresholdDeltaI = (double)thresholdHistIndex / 128.;
	Mat umbra = deltaI <= thresholdDeltaI;
	umbra.convertTo(umbra, CV_64FC1, 1. / 255.);
#ifdef SAVE_ALL
	saveImage("umbra", to3C(umbra));
#endif
	return umbra;
}

// Computes the penumbra component of an image based on the F and NF intensities
// and the detected umbra component
Mat computePenumbra(Mat& umbra, Mat& intensityF, Mat& intensityNF,
	double diagonal, bool advanced) {
	// We compute the magnitude of the gradient ∇I^f and ∇I^nf...
	Mat gradF = computeGradientImage(intensityF);
	Mat gradNF = computeGradientImage(intensityNF);
	// ...and smooth it with a Gaussian of variance 2 pixels to remove noise.
	GaussianBlur(gradF, gradF, Size(0, 0), 2);
	GaussianBlur(gradNF, gradNF, Size(0, 0), 2);
	// We identify candidate penumbra pixels as pixels where the gradient is
	// stronger in the flash image.
	Mat strongerFlashGradients = Mat::zeros(umbra.size(), CV_64FC1);
	for (int j = 0; j < umbra.size().height; j++) {
		for (int i = 0; i < umbra.size().width; i++) {
			double value = gradF.at<double>(j, i);
			if (value > gradNF.at<double>(j, i)) {
				strongerFlashGradients.at<double>(j, i) = 1;
			}
		}
	}
	// We then keep only pixels that are “close” to umbra pixels, that is,
	// such that at least one of their neighbors is in umbra. In practice, 
	// we use a square neighborhood of size 1% of the photo’s diagonal. This
	// computation can be performed efficiently by convolving the
	// binary umbra map with a box filter.

	int boxWidth = (int)(0.01*diagonal);
	Mat umbraNeighborhoods;
	blur(umbra, umbraNeighborhoods, Size(boxWidth, boxWidth));
	threshold(umbraNeighborhoods, umbraNeighborhoods, 0.01, 1., THRESH_BINARY);
	Mat penumbraTouchingUmbra;
	multiply(strongerFlashGradients, umbraNeighborhoods, penumbraTouchingUmbra);

if(advanced) return penumbraTouchingUmbra;

#ifdef SAVE_ALL
	saveImage("gradF", to3C(gradF));
	saveImage("gradNF", to3C(gradNF));
	saveImage("umbraNeighborhoods", to3C(umbraNeighborhoods));
	saveImage("strongerFlashGradients", to3C(strongerFlashGradients));
	saveImage("penumbraTouchingUmbra", to3C(penumbraTouchingUmbra));
#endif	

	// We also must account for shadows cast by tiny objects such
	// as pieces of fur, since these might have a pure penumbra without umbra.
	// We use a similar strategy and consider as
	// shadow pixels that have a large number of neighbors with
	// higher gradient in the flash image. We use a threshold of 80%
	// on a square neighborhood of size 0.7% of the photo’s diagonal.
	boxWidth = (int)(0.007*diagonal);
	Mat higherGradientNeighborhoods;
	blur(strongerFlashGradients, higherGradientNeighborhoods,
		Size(boxWidth, boxWidth));
	higherGradientNeighborhoods = higherGradientNeighborhoods > 0.8;
	Mat purePenumbraPixels;
	higherGradientNeighborhoods.convertTo(purePenumbraPixels,CV_64FC1,
		1. / 255.);

	Mat penumbra; add(penumbraTouchingUmbra, purePenumbraPixels, penumbra);
	threshold(penumbra, penumbra, 1, 1, THRESH_TRUNC);

#ifdef SAVE_ALL
	saveImage("higherGradientNeighborhoods", to3C(higherGradientNeighborhoods));
	saveImage("purePenumbraPixels", to3C(purePenumbraPixels));
	saveImage("penumbra", to3C(penumbra));
#endif
	return penumbra;
}


// The body of the implemented method
int TMOEisemann04::Transform()
{
	if (flashImagePathParameter.GetString() == "") {
		cerr << "No flash image provided" << endl;
		return 1;
	}

	flashImage = TMOImage();
	flashImage.Open(flashImagePathParameter.GetString().c_str());

	Mat imageF = toMat64C3(&flashImage);
	Mat imageNF = toMat64C3(pSrc);

	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();
	double diagonal = sqrt(width*width + height * height);

	Mat intensityF = computeIntensity(imageF);
	Mat intensityNF = computeIntensity(imageNF);
#ifdef SAVE_ALL
	saveImage("intensityF", to3C(intensityF));
	saveImage("intensityNF", to3C(intensityNF));
#endif
	Mat umbra = computeUmbra(intensityF, intensityNF);
	Mat penumbra = computePenumbra(umbra, intensityF, intensityNF, diagonal,
		shadowCorrectionParameter.GetBool());
	Mat shadowMask; add(umbra, penumbra, shadowMask);
	threshold(shadowMask, shadowMask, 1, 1, THRESH_TRUNC);

#ifdef SAVE_ALL	
	saveImage("shadowMask", to3C(shadowMask));
#endif

	Mat largeScaleF = computeLargeScale(intensityF, diagonal, &shadowMask);
	Mat largeScaleNF = computeLargeScale(intensityNF, diagonal);

#ifdef SAVE_ALL
	saveImage("largeScaleF", to3C(largeScaleF));
	saveImage("largeScaleNF", to3C(largeScaleNF));
#endif

	Mat colorF = computeColor(imageF, intensityF);
	Mat colorNF = computeColor(imageNF, intensityNF);

#ifdef SAVE_ALL
	saveImage("colorF", colorF);
	saveImage("colorNF", colorNF);
#endif

	Mat detailF = computeDetail(intensityF, largeScaleF);
	Mat detailNF = computeDetail(intensityNF, largeScaleNF);

#ifdef SAVE_ALL
	saveImage("detailF", to3C(detailF));
	saveImage("detailNF", to3C(detailNF));
#endif
	double minFlashOutsideShadow, maxFlashOutsideShadow;
	minMaxLocMask(detailF, shadowMask, &minFlashOutsideShadow,
		&maxFlashOutsideShadow);
	Mat scaledNonFlashDetail;
	normalize(detailNF, scaledNonFlashDetail, minFlashOutsideShadow,
		maxFlashOutsideShadow, NORM_MINMAX, CV_64F, shadowMask > 0.0);
#ifdef SAVE_ALL
	saveImage("detailNFscaled", to3C(scaledNonFlashDetail));
#endif

	Mat detailFinal = detailF.clone();
	detailNF.copyTo(detailFinal, shadowMask > 0.0);

#ifdef SAVE_ALL	
	saveImage("detailFinal", to3C(detailFinal));
#endif

	Mat result, basicReconstruction;
	multiply(detailFinal, largeScaleNF, basicReconstruction);
	multiply(to3C(basicReconstruction), colorF, basicReconstruction);
	basicReconstruction =
		whiteBalanceCorr(basicReconstruction, imageF, imageNF);
#ifdef SAVE_ALL
	saveImage("basicReconstruction", basicReconstruction);
#endif
	if(shadowCorrectionParameter.GetBool()){
		result = fromLuv(shadowCorrection(
					toLuv(basicReconstruction),
					toLuv(imageF),
					toLuv(imageNF),
					shadowMask, 
					diagonal));
		result.convertTo(result, CV_64FC3);
#ifdef SAVE_ALL
		saveImage("shadowCorrection", result);
#endif
	}
	else{
		result = basicReconstruction;
	}

	// Write the output data
	double* pDestinationData = pDst->GetData();
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			Vec3d resultPixel = result.at<Vec3d>(j, i);
			for (int c = 0; c < 3; c++) {
				*pDestinationData++ = resultPixel[c];
			}
		}
	}
	pDst->Convert(TMO_RGB);
	return 0;
}
