/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                  *
*                                                                              *
*                       Semestral project                                      *
*                       Author: Lucie Svobodova                                *
*                       Brno 2025                                              *
*                                                                              *
*                       Implementation of the TMOAncuti19 class                *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOAncuti19.cpp
 * @brief Implementation of the TMOAncuti19 class
 * @author Lucie Svobodova
 * @class TMOAncuti19
 */

#include "TMOAncuti19.h"

/**
 * @brief Constructor.
 */
TMOAncuti19::TMOAncuti19()
{
	SetName(L"Ancuti19");
	SetDescription(L"Image decolorization based on information theory.");

	alphaParameter.SetName(L"alpha");
	alphaParameter.SetDescription(L"Scaling factor for local entropy threshold; controls influence of local detail.");
	alphaParameter.SetDefault(0.85);
	alphaParameter = 0.85;
	alphaParameter.SetRange(0.0, 1.0);

	this->Register(alphaParameter);
}

/**
 * @brief Destructor.
 */
TMOAncuti19::~TMOAncuti19()
{
}

/**
 * @brief Main method implementing the tone mapping operator.
 *
 * Decomposes the input image, computes global and local weight maps,
 * builds multi-scale pyramids, fuses them, reconstructs the final image,
 * normalizes it, and writes the output.
 *
 * @return 0 on success
 * @todo ica_model.yml shall be copied to final dir (install)
 */
int TMOAncuti19::Transform()
{
	// Load ICA model used for global weight map computation
	if (!LoadICAModel("ica_model.yml")) {
		std::cerr << "Error: ICA model not loaded." << std::endl;
		return 1;
	}

	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();

	// Decompose the input color image into R, G, and B channels
	cv::Mat red, green, blue;
	DecomposeInput(pSrc, red, green, blue);

	// Compute global (saliency) and local weight maps for each channel
	cv::Mat normW_R, normW_G, normW_B;
	ComputeWeightMaps(red, green, blue, normW_R, normW_G, normW_B);

	// Build Laplacian pyramids for each channel
	std::vector<cv::Mat> lapPyr_R, lapPyr_G, lapPyr_B;
	int pyramidLevels = 6;
	BuildPyramids(red, lapPyr_R, pyramidLevels);
	BuildPyramids(green, lapPyr_G, pyramidLevels);
	BuildPyramids(blue, lapPyr_B, pyramidLevels);

	// Build Gaussian pyramids for each normalized weight map
	std::vector<cv::Mat> gaussPyr_W_R, gaussPyr_W_G, gaussPyr_W_B;
	BuildWeightPyramids(normW_R, gaussPyr_W_R, pyramidLevels);
	BuildWeightPyramids(normW_G, gaussPyr_W_G, pyramidLevels);
	BuildWeightPyramids(normW_B, gaussPyr_W_B, pyramidLevels);

	// Fuse the pyramids across channels
	int numLevels = static_cast<int>(lapPyr_R.size());
	std::vector<cv::Mat> fusedPyr;
	for (int lvl = 0; lvl < numLevels; lvl++) {
		cv::Mat fused = gaussPyr_W_R[lvl].mul(lapPyr_R[lvl]) +
						gaussPyr_W_G[lvl].mul(lapPyr_G[lvl]) +
						gaussPyr_W_B[lvl].mul(lapPyr_B[lvl]);
		fusedPyr.push_back(fused);
	}

	// Reconstruct the final image from the fused pyramid
	cv::Mat result = ReconstructFromPyramid(fusedPyr);
	
	// Normalize the final result
	cv::Mat normalizedResult = NormalizeImage(result);

	// Write the final result to the destination image
	WriteOutput(pDst, normalizedResult);

	return 0;
}

/**
 * @brief Decomposes the input image into red, green, and blue channels.
 *
 * @param src Pointer to the input image
 * @param red Output red channel
 * @param green Output green channel
 * @param blue Output blue channel
 */
void TMOAncuti19::DecomposeInput(TMOImage* src, cv::Mat &red, cv::Mat &green, cv::Mat &blue)
{
	int width = src->GetWidth();
	int height = src->GetHeight();
	red = cv::Mat(height, width, CV_64F);
	green = cv::Mat(height, width, CV_64F);
	blue = cv::Mat(height, width, CV_64F);

	double *pSourceData = src->GetData();
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			int idx = (j * width + i) * 3;
			red.at<double>(j, i)   = pSourceData[idx];
			green.at<double>(j, i) = pSourceData[idx + 1];
			blue.at<double>(j, i)  = pSourceData[idx + 2];
		}
	}
}

/**
 * @brief Loads the pre-trained ICA model from a YAML file.
 *
 * @param modelFile Path to the model YAML file
 * @return true if the model is loaded successfully
 */
bool TMOAncuti19::LoadICAModel(const std::string &modelFile)
{
	cv::FileStorage fs(modelFile, cv::FileStorage::READ);
	if (!fs.isOpened()) {
		std::cerr << "Error: Could not open ICA model file: " << modelFile << std::endl;
		return false;
	}
	fs["icaBasis"] >> mIcaBasis;
	fs["mu"] >> mMu;
	fs["sigma"] >> mSigma;
	fs.release();
	return true;
}
 
/**
 * @brief Creates global and local weight maps for each channel, combines and normalizes them.
 *
 * @param red Red channel
 * @param green Green channel
 * @param blue Blue channel
 * @param normW_R Output normalized weight map for red
 * @param normW_G Output normalized weight map for green
 * @param normW_B Output normalized weight map for blue
 */
void TMOAncuti19::ComputeWeightMaps(const cv::Mat &red, const cv::Mat &green, const cv::Mat &blue,
									cv::Mat &normW_R, cv::Mat &normW_G, cv::Mat &normW_B)
{
	// Compute global weight maps
	cv::Mat globalW_R = ComputeGlobalWeightMap(red, 31);
	cv::Mat globalW_G = ComputeGlobalWeightMap(green, 31);
	cv::Mat globalW_B = ComputeGlobalWeightMap(blue, 31);

	// Compute local weight maps
	cv::Mat localW_R = ComputeLocalWeightMap(red, 5);
	cv::Mat localW_G = ComputeLocalWeightMap(green, 5);
	cv::Mat localW_B = ComputeLocalWeightMap(blue, 5);

	// Combine the weight maps for each channel
	cv::Mat norm_R = globalW_R + localW_R;
	cv::Mat norm_G = globalW_G + localW_G;
	cv::Mat norm_B = globalW_B + localW_B;

	// Normalize the weight maps so that for every pixel, the R+G+B=1
	cv::Mat sumWeights = norm_R + norm_G + norm_B;
	sumWeights.setTo(1, sumWeights < 1e-6);
	normW_R = norm_R / sumWeights;
	normW_G = norm_G / sumWeights;
	normW_B = norm_B / sumWeights;
}

/**
 * @brief Computes the entropy over a given patch.
 *
 * @param patch Input patch
 * @return Entropy value
 */
double TMOAncuti19::ComputeEntropy(const cv::Mat &patch)
{
	int numBins = 256;
	std::vector<int> hist(numBins, 0);
	int total = patch.rows * patch.cols;
	for (int r = 0; r < patch.rows; r++) {
		for (int c = 0; c < patch.cols; c++) {
			double val = patch.at<double>(r, c);
			int bin = static_cast<int>(val * (numBins - 1));
			hist[bin]++;
		}
	}
	double entropy = 0.0;
	for (int i = 0; i < numBins; i++) {
		if (hist[i] > 0) {
			double p = static_cast<double>(hist[i]) / total;
			entropy += p * log(1.0 / p);
		}
	}
	return entropy;
}
 
/**
 * @brief Computes a local weight map using entropy over a sliding patch.
 *
 * @param channel Input channel
 * @param patchSize Size of the patch
 * @return Local weight map
 */
cv::Mat TMOAncuti19::ComputeLocalWeightMap(const cv::Mat &channel, int patchSize)
{
	int width = channel.cols;
	int height = channel.rows;
	cv::Mat localEntropy = cv::Mat::zeros(height, width, CV_64F);
	int halfPatch = patchSize / 2;
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++) {
			int r0 = std::max(0, r - halfPatch);
			int r1 = std::min(height - 1, r + halfPatch);
			int c0 = std::max(0, c - halfPatch);
			int c1 = std::min(width - 1, c + halfPatch);
			cv::Rect roi(c0, r0, c1 - c0 + 1, r1 - r0 + 1);
			cv::Mat patch = channel(roi);
			double entropy = ComputeEntropy(patch);
			localEntropy.at<double>(r, c) = entropy;
		}
	}
	double maxEntropy;
	cv::minMaxLoc(localEntropy, nullptr, &maxEntropy);
	double theta = alphaParameter * maxEntropy;
	cv::Mat weight = localEntropy / theta;
	cv::threshold(weight, weight, 1.0, 1.0, cv::THRESH_TRUNC);
	return weight;
}

/**
 * @brief Computes the global weight (self-information) at a given pixel.
 *
 * @param channel Input channel
 * @param r Row coordinate
 * @param c Column coordinate
 * @param patchSize Size of the patch
 * @return Global weight at the pixel
 */
double TMOAncuti19::ComputeGlobalWeightAtPixel(const cv::Mat &channel, int r, int c, int patchSize)
{
	int halfPatch = patchSize / 2;
	int r0 = std::max(0, r - halfPatch);
	int r1 = std::min(channel.rows - 1, r + halfPatch);
	int c0 = std::max(0, c - halfPatch);
	int c1 = std::min(channel.cols - 1, c + halfPatch);
	cv::Rect roi(c0, r0, c1 - c0 + 1, r1 - r0 + 1);
	cv::Mat patch = channel(roi);

	if (!patch.isContinuous())
		patch = patch.clone();

	if (patch.rows != patchSize || patch.cols != patchSize)
		cv::resize(patch, patch, cv::Size(patchSize, patchSize));

	cv::Mat patchVec = patch.reshape(1, patchSize * patchSize);
	cv::Mat coeffs = mIcaBasis * patchVec;

	double logLikelihoodSum = 0.0;
	int numComponents = coeffs.rows;
	for (int i = 0; i < numComponents; i++) {
		double c_val = coeffs.at<double>(i, 0);
		double m = mMu.at<double>(i, 0);
		double s = mSigma.at<double>(i, 0);
		double logL = -std::log(s * std::sqrt(2 * M_PI))
					- ((c_val - m) * (c_val - m)) / (2 * s * s);
		logLikelihoodSum += logL;
	}
	return -logLikelihoodSum;
}

/**
 * @brief Computes a global weight map for one channel.
 *
 * @param channel Input channel
 * @param patchSize Size of the patch
 * @return Global weight map
 */
cv::Mat TMOAncuti19::ComputeGlobalWeightMap(const cv::Mat &channel, int patchSize)
{
	cv::Mat weightMap(channel.size(), CV_64F);
	for (int r = 0; r < channel.rows; r++) {
		for (int c = 0; c < channel.cols; c++) {
			double weight = ComputeGlobalWeightAtPixel(channel, r, c, patchSize);
			weightMap.at<double>(r, c) = weight;
		}
	}
	return weightMap;
}

/**
 * @brief Constructs a Laplacian pyramid for the given channel.
 *
 * @param channel Input channel
 * @param lapPyr Output Laplacian pyramid
 * @param numLevels Number of pyramid levels
 */
void TMOAncuti19::BuildPyramids(const cv::Mat &channel, std::vector<cv::Mat> &lapPyr, int numLevels)
{
	std::vector<cv::Mat> gaussPyr;
	gaussPyr.push_back(channel);
	for (int lvl = 1; lvl < numLevels; lvl++) {
		cv::Mat down;
		cv::pyrDown(gaussPyr[lvl - 1], down);
		if (down.cols < 2 || down.rows < 2)
			break;
		gaussPyr.push_back(down);
	}
	for (size_t lvl = 0; lvl < gaussPyr.size() - 1; lvl++) {
		cv::Mat up;
		cv::pyrUp(gaussPyr[lvl + 1], up, gaussPyr[lvl].size());
		cv::Mat lap = gaussPyr[lvl] - up;
		lapPyr.push_back(lap);
	}
	lapPyr.push_back(gaussPyr.back());
}

/**
 * @brief Constructs a Gaussian pyramid for the given weight map.
 *
 * @param weight Input weight map
 * @param gaussPyr Output Gaussian pyramid
 * @param numLevels Number of pyramid levels
 */
void TMOAncuti19::BuildWeightPyramids(const cv::Mat &weight, std::vector<cv::Mat> &gaussPyr, int numLevels)
{
	gaussPyr.push_back(weight);
	for (int lvl = 1; lvl < numLevels; lvl++) {
		cv::Mat down;
		cv::pyrDown(gaussPyr[lvl - 1], down);
		if (down.cols < 2 || down.rows < 2)
			break;
		gaussPyr.push_back(down);
	}
}

/**
 * @brief Reconstructs the image from its fused Laplacian pyramid.
 *
 * @param fusedPyr Input fused pyramid
 * @return Reconstructed image
 */
cv::Mat TMOAncuti19::ReconstructFromPyramid(const std::vector<cv::Mat> &fusedPyr)
{
	cv::Mat reconstructed = fusedPyr.back();
	for (int lvl = static_cast<int>(fusedPyr.size()) - 2; lvl >= 0; lvl--) {
		cv::Mat up;
		cv::pyrUp(reconstructed, up, fusedPyr[lvl].size());
		reconstructed = up + fusedPyr[lvl];
	}
	return reconstructed;
}

/**
 * @brief Normalizes a grayscale image to the range [0, 1].
 *
 * @param input Input image
 * @return Normalized image
 */
cv::Mat TMOAncuti19::NormalizeImage(const cv::Mat &input)
{
    double minVal, maxVal;
    cv::minMaxLoc(input, &minVal, &maxVal);
    cv::Mat output;
    if (maxVal > minVal)
        output = (input - minVal) / (maxVal - minVal);
    else
        output = cv::Mat::zeros(input.size(), input.type());
    return output;
}

/**
 * @brief Writes the final grayscale result to the destination image.
 *
 * @param dst Destination image
 * @param result Final result image
 */
void TMOAncuti19::WriteOutput(TMOImage* dst, const cv::Mat &result)
{
	int width = dst->GetWidth();
	int height = dst->GetHeight();
	double *pDestinationData = dst->GetData();
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			double val = result.at<double>(j, i);
			*pDestinationData++ = val;
			*pDestinationData++ = val;
			*pDestinationData++ = val;
		}
	}
}
