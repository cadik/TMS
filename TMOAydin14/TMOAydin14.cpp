/*******************************************************************************
*                                                                              *
*                        Brno University of Technology                         *
*                      Faculty of Information Technology                       *
*                                                                              *
*                              Tone Mapping Methods                            *
*                                                                              *
*            Author: Hugo Bohácsek [xbohach00 AT stud.fit.vutbr.cz]            *
*                                   Brno 2026                                  *
*                                                                              *
*                    Implementation of the TMOAydin14 class                    *
*            Temporally Coherent Local Tone Mapping of HDR Video               *
*                                                                              *
*******************************************************************************/

#include "TMOAydin14.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <iostream>

// Static definitions
constexpr double TMOAydin14::EPSILON;

TMOAydin14::TMOAydin14() {
	SetName(L"Aydin14");
	SetDescription(L"Temporally Coherent Local Tone Mapping of HDR Video (Aydin, Stefanoski, Croci, Gross, Smolic 2014)");

	// Spatial filtering parameters
	dSigmaS.SetName(L"SigmaS");
	dSigmaS.SetDescription(L"Spatial permeability sigma: controls edge sensitivity. Lower values preserve more edges. Paper range: [0.1, 1.0].");
	dSigmaS.SetDefault(0.4);
	dSigmaS.SetRange(0.01, 2.0);
	this->Register(dSigmaS);
	dSigmaS = 0.4;

	dAlpha.SetName(L"Alpha");
	dAlpha.SetDescription(L"Edge-stopping exponent for Lorentzian function. Paper uses alpha=2.");
	dAlpha.SetDefault(2.0);
	dAlpha.SetRange(1.0, 4.0);
	this->Register(dAlpha);
	dAlpha = 2.0;

	iNumIter.SetName(L"NumIter");
	iNumIter.SetDescription(L"Number of spatial filtering iterations. Paper uses 20. More iterations capture coarser detail in the base layer.");
	iNumIter.SetDefault(20);
	iNumIter.SetRange(1, 100);
	this->Register(iNumIter);
	iNumIter = 20;

	// Temporal filtering parameters
	iTemporalRadius.SetName(L"TemporalRadius");
	iTemporalRadius.SetDescription(L"Temporal neighborhood radius (frames in each direction). Paper uses 10 (total window = 21). Set to 0 for single-image mode.");
	iTemporalRadius.SetDefault(10);
	iTemporalRadius.SetRange(0, 30);
	this->Register(iTemporalRadius);
	iTemporalRadius = 10;

	dSigmaPC.SetName(L"SigmaPC");
	dSigmaPC.SetDescription(L"Photo-constancy sigma for temporal permeability. Paper: 0.1. Controls sensitivity to warping errors.");
	dSigmaPC.SetDefault(0.1);
	dSigmaPC.SetRange(0.01, 1.0);
	this->Register(dSigmaPC);
	dSigmaPC = 0.1;

	dSigmaFG.SetName(L"SigmaFG");
	dSigmaFG.SetDescription(L"Flow gradient sigma for temporal permeability. Paper: 2.0. Penalizes complex motion with high flow gradients.");
	dSigmaFG.SetDefault(2.0);
	dSigmaFG.SetRange(0.1, 10.0);
	this->Register(dSigmaFG);
	dSigmaFG = 2.0;

	// Tone mapping parameters
	dCompression.SetName(L"Compression");
	dCompression.SetDescription(L"Base layer compression factor (log-domain). Paper range: [0.2, 0.4]. Lower = stronger compression.");
	dCompression.SetDefault(0.3);
	dCompression.SetRange(0.05, 1.0);
	this->Register(dCompression);
	dCompression = 0.3;

	dDragoBias.SetName(L"DragoBias");
	dDragoBias.SetDescription(L"Drago et al. 2003 bias parameter. Paper range: [0.2, 0.5]. Controls contrast distribution between darks and brights.");
	dDragoBias.SetDefault(0.35);
	dDragoBias.SetRange(0.05, 1.0);
	this->Register(dDragoBias);
	dDragoBias = 0.35;

	bUseDrago.SetName(L"UseDrago");
	bUseDrago.SetDescription(L"Use Drago adaptive logarithmic tone curve instead of a simple compression factor.");
	bUseDrago.SetDefault(false);
	this->Register(bUseDrago);
	bUseDrago = false;

	dSaturation.SetName(L"Saturation");
	dSaturation.SetDescription(L"Color saturation control. 0.0 = grayscale, 1.0 = full color.");
	dSaturation.SetDefault(1.0);
	dSaturation.SetRange(0.0, 2.0);
	this->Register(dSaturation);
	dSaturation = 1.0;

	dGamma.SetName(L"Gamma");
	dGamma.SetDescription(L"Display gamma for final output. Typical: 2.2.");
	dGamma.SetDefault(2.2);
	dGamma.SetRange(1.0, 3.0);
	this->Register(dGamma);
	dGamma = 2.2;
}

TMOAydin14::~TMOAydin14() {
}

int TMOAydin14::Transform() {
	// Work in linear RGB space throughout (avoid clipping)
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	const double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();
	int numPixels = width * height;

	// Retrieve parameters
	double sigmaS = dSigmaS.GetDouble();
	double alphaEdge = dAlpha.GetDouble();
	int numIter = iNumIter.GetInt();
	int temporalRadius = iTemporalRadius.GetInt();
	double sigmaPC = dSigmaPC.GetDouble();
	double sigmaFG = dSigmaFG.GetDouble();
	double compression = dCompression.GetDouble();
	double dragoBias = dDragoBias.GetDouble();
	bool useDrago = bUseDrago.GetBool();
	double saturation = dSaturation.GetDouble();
	double gamma = dGamma.GetDouble();

	// Step 1: Extract current frame data
	FrameData currentFrame;
	currentFrame.width = width;
	currentFrame.height = height;
	currentFrame.logLum.resize(numPixels);
	currentFrame.logR.resize(numPixels);
	currentFrame.logG.resize(numPixels);
	currentFrame.logB.resize(numPixels);

	// Extract linear RGB and compute log-domain values
	for (int i = 0; i < numPixels; i++) {
		double R = std::max(pSourceData[i * 3], EPSILON);
		double G = std::max(pSourceData[i * 3 + 1], EPSILON);
		double B = std::max(pSourceData[i * 3 + 2], EPSILON);

		double Y = 0.2126 * R + 0.7152 * G + 0.0722 * B;
		Y = std::max(Y, EPSILON);

		currentFrame.logLum[i] = std::log10(Y);
		currentFrame.logR[i] = std::log10(R);
		currentFrame.logG[i] = std::log10(G);
		currentFrame.logB[i] = std::log10(B);
	}

	// Prepare normalized grayscale image for optical flow computation
	double minLogLum = *std::min_element(currentFrame.logLum.begin(), currentFrame.logLum.end());
	double maxLogLum = *std::max_element(currentFrame.logLum.begin(), currentFrame.logLum.end());
	double rangeLogLum = std::max(maxLogLum - minLogLum, EPSILON);

	currentFrame.grayFloat = cv::Mat(height, width, CV_32F);
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			float val = (float)((currentFrame.logLum[y * width + x] - minLogLum) / rangeLogLum);
			currentFrame.grayFloat.at<float>(y, x) = val;
		}
	}

	pSrc->ProgressBar(1, 10);

	// Step 2: Spatial filtering of current frame's log luminance
	spatialFilter(currentFrame.logLum, currentFrame.spatialFiltered, width, height, sigmaS, alphaEdge, numIter);

	pSrc->ProgressBar(3, 10);

	// Step 3: Temporal cache management
	if (width != lastWidth || height != lastHeight) { // Reset cache if frame dimensions changed
		cache.clear();
		lastWidth = width;
		lastHeight = height;
	}

	// Set max cache size = 2 * temporalRadius + 1
	cache.maxFrames = 2 * temporalRadius + 1;

	// Compute optical flow with previous frame (if available)
	if (!cache.frames.empty()) {
		cv::Mat flow;
		computeOpticalFlow(cache.frames.back().grayFloat, currentFrame.grayFloat, flow);
		cache.forwardFlows.push_back(flow);
	}

	// Add current frame to cache
	cache.frames.push_back(currentFrame);

	// Trim cache to max size
	while ((int)cache.frames.size() > cache.maxFrames) {
		cache.frames.pop_front();
		if (!cache.forwardFlows.empty()) {
			cache.forwardFlows.pop_front();
		}
	}

	pSrc->ProgressBar(4, 10);

	// Step 4: Base layer computation
	std::vector<double> baseLayer;

	// Center frame index in cache (current frame = last in cache)
	int centerIdx = (int)cache.frames.size() - 1;
	int numFramesInCache = (int)cache.frames.size();

	if (temporalRadius == 0 || numFramesInCache <= 1) { // Single-image mode: base layer = spatially filtered log luminance
		baseLayer = currentFrame.spatialFiltered;
	} else { // Video mode: warp spatially filtered frames to center, temporal filter

		// Warp all spatially filtered frames to center frame
		std::vector<std::vector<double>> warpedSpatial(numFramesInCache);
		warpedSpatial[centerIdx] = cache.frames[centerIdx].spatialFiltered;

		for (int i = 0; i < numFramesInCache; i++) {
			if (i == centerIdx) continue;
			cv::Mat cumulFlow = accumulateFlow(i, centerIdx);
			if (cumulFlow.empty()) { // Fallback: use unwarped
				warpedSpatial[i] = cache.frames[i].spatialFiltered;
			} else {
				warpFrame(cache.frames[i].spatialFiltered, warpedSpatial[i], cumulFlow, width, height);
			}
		}

		// Calculate temporal permeabilities between adjacent warped frames
		std::vector<std::vector<double>> temporalPerms(numFramesInCache - 1);
		for (int i = 0; i < numFramesInCache - 1; i++) {
			// Get the flow between frames i and i+1
			cv::Mat flowBetween;
			if (i < (int)cache.forwardFlows.size()) {
				flowBetween = cache.forwardFlows[i];
			}
			computeTemporalPermeability(warpedSpatial[i], warpedSpatial[i + 1], flowBetween, temporalPerms[i], width, height, sigmaPC, sigmaFG, alphaEdge);
		}

		// Temporal filter -> base layer
		temporalFilter1D(warpedSpatial, temporalPerms, centerIdx, baseLayer, numPixels);
	}

	pSrc->ProgressBar(6, 10);

	// Step 5: Detail layer computation
	std::vector<double> detailR(numPixels), detailG(numPixels), detailB(numPixels);

	if (temporalRadius == 0 || numFramesInCache <= 1) { // Single-image mode: detail = logChannel - base
		for (int i = 0; i < numPixels; i++) {
			detailR[i] = currentFrame.logR[i] - baseLayer[i];
			detailG[i] = currentFrame.logG[i] - baseLayer[i];
			detailB[i] = currentFrame.logB[i] - baseLayer[i];
		}
	} else { // Video mode - temporal filter each color channel, then subtract base

		// Warp original (non-spatially-filtered) color channels to center
		auto temporalFilterChannel = [&](
			const std::function<const std::vector<double>&(int)>& getChannel,
			std::vector<double>& filteredChannel)
		{
			std::vector<std::vector<double>> warpedChannels(numFramesInCache);
			for (int i = 0; i < numFramesInCache; i++) {
				if (i == centerIdx) {
					warpedChannels[i] = getChannel(i);
				} else {
					cv::Mat cumulFlow = accumulateFlow(i, centerIdx);
					if (cumulFlow.empty()) {
						warpedChannels[i] = getChannel(i);
					} else {
						warpFrame(getChannel(i), warpedChannels[i], cumulFlow, width, height);
					}
				}
			}

			// Compute temporal permeabilities (same confidence as base layer)
			std::vector<std::vector<double>> temporalPerms(numFramesInCache - 1);
			for (int i = 0; i < numFramesInCache - 1; i++) {
				cv::Mat flowBetween;
				if (i < (int)cache.forwardFlows.size()) {
					flowBetween = cache.forwardFlows[i];
				}
				computeTemporalPermeability(warpedChannels[i], warpedChannels[i + 1], flowBetween, temporalPerms[i], width, height, sigmaPC, sigmaFG, alphaEdge);
			}

			temporalFilter1D(warpedChannels, temporalPerms, centerIdx, filteredChannel, numPixels);
		};

		std::vector<double> filteredR, filteredG, filteredB;

		temporalFilterChannel([&](int i) -> const std::vector<double>& {
			return cache.frames[i].logR;
		}, filteredR);

		temporalFilterChannel([&](int i) -> const std::vector<double>& {
			return cache.frames[i].logG;
		}, filteredG);

		temporalFilterChannel([&](int i) -> const std::vector<double>& {
			return cache.frames[i].logB;
		}, filteredB);

		// Detail = temporally filtered color - base
		for (int i = 0; i < numPixels; i++) {
			detailR[i] = filteredR[i] - baseLayer[i];
			detailG[i] = filteredG[i] - baseLayer[i];
			detailB[i] = filteredB[i] - baseLayer[i];
		}
	}

	pSrc->ProgressBar(8, 10);

	// Step 6: Tone curve application to base layer
	if (useDrago) {
		applyDragoToneCurve(baseLayer, numPixels, dragoBias, 100.0);
	} else {
		applyCompressionFactor(baseLayer, numPixels, compression);
	}

	// Step 7: Recombine base + detail, convert to display values
	for (int i = 0; i < numPixels; i++) {
		// Recombine in log domain
		double outLogR = baseLayer[i] + saturation * detailR[i];
		double outLogG = baseLayer[i] + saturation * detailG[i];
		double outLogB = baseLayer[i] + saturation * detailB[i];

		// Convert from log10 to linear
		double R = std::pow(10.0, outLogR);
		double G = std::pow(10.0, outLogG);
		double B = std::pow(10.0, outLogB);

		// Clamp negatives
		R = std::max(0.0, R);
		G = std::max(0.0, G);
		B = std::max(0.0, B);

		pDestinationData[i * 3]     = R;
		pDestinationData[i * 3 + 1] = G;
		pDestinationData[i * 3 + 2] = B;
	}

	// Normalize output to [0, 1] using percentile-based approach
	{
		// Collect all values
		std::vector<double> allValues(numPixels * 3);
		for (int i = 0; i < numPixels * 3; i++) {
			allValues[i] = pDestinationData[i];
		}
		std::sort(allValues.begin(), allValues.end());

		// Use 0.5th and 99.5th percentiles for better normalization
		int lowIdx = (int)(0.005 * allValues.size());
		int highIdx = (int)(0.995 * allValues.size());
		highIdx = std::min(highIdx, (int)allValues.size() - 1);
		double lowVal = allValues[lowIdx];
		double highVal = allValues[highIdx];
		double range = std::max(highVal - lowVal, EPSILON);

		for (int i = 0; i < numPixels; i++) {
			double R = (pDestinationData[i * 3]     - lowVal) / range;
			double G = (pDestinationData[i * 3 + 1] - lowVal) / range;
			double B = (pDestinationData[i * 3 + 2] - lowVal) / range;

			// Clamp
			R = std::max(0.0, std::min(1.0, R));
			G = std::max(0.0, std::min(1.0, G));
			B = std::max(0.0, std::min(1.0, B));

			// Apply display gamma
			double invGamma = 1.0 / gamma;
			R = std::pow(R, invGamma);
			G = std::pow(G, invGamma);
			B = std::pow(B, invGamma);

			pDestinationData[i * 3]     = R;
			pDestinationData[i * 3 + 1] = G;
			pDestinationData[i * 3 + 2] = B;
		}
	}

	pSrc->ProgressBar(10, 10);
	return 0;
}

// Compute horizontal permeability map
void TMOAydin14::computePermeabilityH(const std::vector<double>& input, std::vector<double>& perm, int width, int height, double sigma, double alpha) {
	perm.resize(height * (width - 1));

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width - 1; x++) {
			double diff = input[y * width + x] - input[y * width + (x + 1)];
			perm[y * (width - 1) + x] = lorentzian(diff, sigma, alpha);
		}
	}
}

// Compute vertical permeability map
void TMOAydin14::computePermeabilityV(const std::vector<double>& input, std::vector<double>& perm, int width, int height, double sigma, double alpha) {
	perm.resize((height - 1) * width);

	for (int y = 0; y < height - 1; y++) {
		for (int x = 0; x < width; x++) {
			double diff = input[y * width + x] - input[(y + 1) * width + x];
			perm[y * width + x] = lorentzian(diff, sigma, alpha);
		}
	}
}

// Single iteration of permeability-based diffusion filter
void TMOAydin14::spatialFilterIteration(const std::vector<double>& original, std::vector<double>& current, const std::vector<double>& perm, int width, int height, bool horizontal) {
	std::vector<double> result(width * height);

	if (horizontal) { // Process each row independently
		for (int y = 0; y < height; y++) {
			int rowOffset = y * width;
			int permRowOffset = y * (width - 1);

			// Forward pass (left to right) - accumulate left contributions
			std::vector<double> leftSum(width, 0.0);
			std::vector<double> leftWeight(width, 0.0);

			for (int x = 1; x < width; x++) {
				double pi = perm[permRowOffset + (x - 1)];
				leftSum[x] = pi * (leftSum[x - 1] + current[rowOffset + (x - 1)]);
				leftWeight[x] = pi * (leftWeight[x - 1] + 1.0);
			}

			// Backward pass (right to left) - accumulate right contributions
			std::vector<double> rightSum(width, 0.0);
			std::vector<double> rightWeight(width, 0.0);

			for (int x = width - 2; x >= 0; x--) {
				double pi = perm[permRowOffset + x];
				rightSum[x] = pi * (rightSum[x + 1] + current[rowOffset + (x + 1)]);
				rightWeight[x] = pi * (rightWeight[x + 1] + 1.0);
			}

			// Combine
			for (int x = 0; x < width; x++) {
				double S = leftWeight[x] + rightWeight[x] + 1.0;
				result[rowOffset + x] = (leftSum[x] + rightSum[x] + original[rowOffset + x]) / S;
			}
		}
	} else { // Process each column independently (vertical filtering)
		for (int x = 0; x < width; x++) {
			// Forward pass (top to bottom)
			std::vector<double> topSum(height, 0.0);
			std::vector<double> topWeight(height, 0.0);

			for (int y = 1; y < height; y++) {
				double pi = perm[(y - 1) * width + x];
				topSum[y] = pi * (topSum[y - 1] + current[(y - 1) * width + x]);
				topWeight[y] = pi * (topWeight[y - 1] + 1.0);
			}

			// Backward pass (bottom to top)
			std::vector<double> botSum(height, 0.0);
			std::vector<double> botWeight(height, 0.0);

			for (int y = height - 2; y >= 0; y--) {
				double pi = perm[y * width + x];
				botSum[y] = pi * (botSum[y + 1] + current[(y + 1) * width + x]);
				botWeight[y] = pi * (botWeight[y + 1] + 1.0);
			}

			// Combine
			for (int y = 0; y < height; y++) {
				double S = topWeight[y] + botWeight[y] + 1.0;
				result[y * width + x] = (topSum[y] + botSum[y] + original[y * width + x]) / S;
			}
		}
	}

	current = result;
}

// Iterative spatial filtering
void TMOAydin14::spatialFilter(const std::vector<double>& input, std::vector<double>& output, int width, int height, double sigma, double alphaEdge, int numIter) {
	// Compute permeability maps from the original input
	std::vector<double> permH, permV;
	computePermeabilityH(input, permH, width, height, sigma, alphaEdge);
	computePermeabilityV(input, permV, width, height, sigma, alphaEdge);

	// Initialize J^(0) = I
	output = input;

	// Iterate, alternating between horizontal and vertical diffusion
	for (int k = 0; k < numIter; k++) {
		bool horizontal = (k % 2 == 0);
		if (horizontal) {
			spatialFilterIteration(input, output, permH, width, height, true);
		} else {
			spatialFilterIteration(input, output, permV, width, height, false);
		}
	}
}

// Compute dense optical flow using OpenCV Farneback
void TMOAydin14::computeOpticalFlow(const cv::Mat& prev, const cv::Mat& next, cv::Mat& flow) {
	// Convert to 8-bit for Farneback (better stability)
	cv::Mat prev8u, next8u;
	prev.convertTo(prev8u, CV_8U, 255.0);
	next.convertTo(next8u, CV_8U, 255.0);

	cv::calcOpticalFlowFarneback(prev8u, next8u, flow,
								 0.5,   // pyr_scale
								 5,     // levels
								 15,    // winsize
								 3,     // iterations
								 7,     // poly_n
								 1.5,   // poly_sigma
								 0);    // flags
}

// Accumulate flow vectors from frame "fromIdx" to frame "toIdx"
cv::Mat TMOAydin14::accumulateFlow(int fromIdx, int toIdx) {
	if (fromIdx == toIdx) {
		return cv::Mat();
	}

	int width = cache.frames[0].width;
	int height = cache.frames[0].height;

	cv::Mat cumulFlow = cv::Mat::zeros(height, width, CV_32FC2);

	if (fromIdx < toIdx) { // Forward direction - accumulate forward flows
		for (int i = fromIdx; i < toIdx; i++) {
			if (i >= (int)cache.forwardFlows.size()) break;
			cv::add(cumulFlow, cache.forwardFlows[i], cumulFlow);
		}
	} else { // Backward direction - accumulate negated forward flows
		for (int i = toIdx; i < fromIdx; i++) {
			if (i >= (int)cache.forwardFlows.size()) break;
			cv::subtract(cumulFlow, cache.forwardFlows[i], cumulFlow);
		}
	}

	return cumulFlow;
}

// Warp frame using flow field with bilinear interpolation
void TMOAydin14::warpFrame(const std::vector<double>& src, std::vector<double>& dst, const cv::Mat& flow, int width, int height) {
	int numPixels = width * height;
	dst.resize(numPixels);

	// Build remap coordinates
	// 		for each pixel (x,y) in dst, sample src at (x+dx, y+dy)
	// 		where (dx,dy) = -flow(x,y) (warp source to align with destination)
	cv::Mat mapX(height, width, CV_32F);
	cv::Mat mapY(height, width, CV_32F);

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			cv::Vec2f f = flow.at<cv::Vec2f>(y, x);
			// Flow goes from source to center; warp source by shifting with flow
			mapX.at<float>(y, x) = (float)x + f[0];
			mapY.at<float>(y, x) = (float)y + f[1];
		}
	}

	// Convert source for remapping
	cv::Mat srcMat(height, width, CV_64F);
	for (int i = 0; i < numPixels; i++) {
		srcMat.at<double>(i / width, i % width) = src[i];
	}

	cv::Mat dstMat;
	cv::remap(srcMat, dstMat, mapX, mapY, cv::INTER_LINEAR, cv::BORDER_REPLICATE);

	for (int i = 0; i < numPixels; i++) {
		dst[i] = dstMat.at<double>(i / width, i % width);
	}
}

// Compute temporal permeability between two adjacent warped frames
void TMOAydin14::computeTemporalPermeability(const std::vector<double>& frame1, const std::vector<double>& frame2, const cv::Mat& flow, std::vector<double>& perm, int width, int height, double sigmaPC, double sigmaFG, double alpha) {
	int numPixels = width * height;
	perm.resize(numPixels);

	for (int i = 0; i < numPixels; i++) {
		// Photo-constancy: Lorentzian of color difference between warped frames
		double colorDiff = frame1[i] - frame2[i];
		double pcPerm = lorentzian(colorDiff, sigmaPC, alpha);

		// Flow gradient: Lorentzian of flow gradient magnitude
		double fgPerm = 1.0; // Default: no flow gradient penalty

		if (!flow.empty()) {
			int y = i / width;
			int x = i % width;

			// Compute gradient of flow field using finite differences
			double gradMag = 0.0;

			if (x > 0 && x < width - 1 && y > 0 && y < height - 1) {
					cv::Vec2f fRight = flow.at<cv::Vec2f>(y, x + 1);
					cv::Vec2f fLeft = flow.at<cv::Vec2f>(y, x - 1);
					cv::Vec2f fDown = flow.at<cv::Vec2f>(y + 1, x);
					cv::Vec2f fUp = flow.at<cv::Vec2f>(y - 1, x);

					// Central differences for flow gradient
					double dFx_dx = (fRight[0] - fLeft[0]) * 0.5;
					double dFy_dx = (fRight[1] - fLeft[1]) * 0.5;
					double dFx_dy = (fDown[0] - fUp[0]) * 0.5;
					double dFy_dy = (fDown[1] - fUp[1]) * 0.5;

					gradMag = std::sqrt(dFx_dx * dFx_dx + dFy_dx * dFy_dx + dFx_dy * dFx_dy + dFy_dy * dFy_dy);
			}

			fgPerm = lorentzian(gradMag, sigmaFG, alpha);
		}

		// Final temporal permeability: product of both measures
		perm[i] = pcPerm * fgPerm;
	}
}

// 1D temporal filtering along the temporal dimension
void TMOAydin14::temporalFilter1D(const std::vector<std::vector<double>>& warpedFrames, const std::vector<std::vector<double>>& temporalPerms, int centerIdx, std::vector<double>& output, int numPixels) {
	int numFrames = (int)warpedFrames.size();
	output.resize(numPixels);

	// For each pixel independently, apply 1D permeability-based filter along time
	for (int p = 0; p < numPixels; p++) {
		// Left pass: accumulate contributions from past frames
		double lSum = 0.0;
		double lWeight = 0.0;

		for (int t = 0; t < centerIdx; t++) {
			double pi = (t < (int)temporalPerms.size() && !temporalPerms[t].empty()) ? temporalPerms[t][p] : EPSILON;
			lSum = pi * (lSum + warpedFrames[t][p]);
			lWeight = pi * (lWeight + 1.0);
		}

		// Right pass: accumulate contributions from future frames
		double rSum = 0.0;
		double rWeight = 0.0;

		for (int t = numFrames - 2; t >= centerIdx; t--) {
			double pi = (t >= 0 && t < (int)temporalPerms.size() && !temporalPerms[t].empty()) ? temporalPerms[t][p] : EPSILON;
			rSum = pi * (rSum + warpedFrames[t + 1][p]);
			rWeight = pi * (rWeight + 1.0);
		}

		// Combine with fidelity term
		double S = lWeight + rWeight + 1.0;
		output[p] = (lSum + rSum + warpedFrames[centerIdx][p]) / S;
	}
}

// Log-domain compression (Durand & Dorsey 2002)
void TMOAydin14::applyCompressionFactor(std::vector<double>& base, int numPixels, double c) {
	// Find max of base layer
	double maxBase = *std::max_element(base.begin(), base.begin() + numPixels);

	for (int i = 0; i < numPixels; i++) {
		base[i] = c * (base[i] - maxBase);
	}
}

// Drago et al. 2003 Adaptive Logarithmic Mapping
void TMOAydin14::applyDragoToneCurve(std::vector<double>& base, int numPixels, double bias, double Ldmax) {
	// Convert base from log10 to linear luminance
	std::vector<double> Lw(numPixels);
	double Lwmax = 0.0;
	for (int i = 0; i < numPixels; i++) {
		Lw[i] = std::pow(10.0, base[i]);
		Lw[i] = std::max(Lw[i], EPSILON);
		Lwmax = std::max(Lwmax, Lw[i]);
	}
	Lwmax = std::max(Lwmax, EPSILON);

	double biasExponent = std::log(bias) / std::log(0.5);
	double logLwmax1 = std::log10(Lwmax + 1.0);
	logLwmax1 = std::max(logLwmax1, EPSILON);

	for (int i = 0; i < numPixels; i++) {
		double lwNorm = Lw[i] / Lwmax;
		double logBase = 2.0 + std::pow(lwNorm, biasExponent) * 8.0;
		double Ld = Ldmax * 0.01 / logLwmax1 * std::log(Lw[i] + 1.0) / std::log(logBase);
		Ld = std::max(Ld, EPSILON);

		// Convert back to log10 domain
		base[i] = std::log10(Ld);
	}
}
