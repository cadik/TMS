/*******************************************************************************
*                                                                              *
*                        Brno University of Technology                         *
*                      Faculty of Information Technology                       *
*                                                                              *
*                              Tone Mapping Methods                            *
*                                                                              *
*            Author: Hugo Bohácsek [xbohach00 AT stud.fit.vutbr.cz]            *
*                                   Brno 2025                                  *
*                                                                              *
*                   Implementation of the TMOEilertsen15 class                 *
*                        Real-time noise-aware tone mapping                    *
*                                                                              *
*******************************************************************************/

#include "TMOEilertsen15.h"
#include <cstring>
#include <cassert>
#include <cfloat>
#include <iostream>

// Minimum luminance to avoid log(0)
static const double LUM_EPS = 1e-6;

TMOEilertsen15::TMOEilertsen15() {
	SetName(L"Eilertsen15");
	SetDescription(L"Real-time noise-aware tone mapping (Eilertsen, Mantiuk, Unger 2015)");

	dPeakLuminance.SetName(L"Peak Luminance");
	dPeakLuminance.SetDescription(L"Lmax: Peak display luminance.");
	dPeakLuminance.SetDefault(200.0);
	dPeakLuminance.SetRange(50.0, 10000.0);
	this->Register(dPeakLuminance);
	dPeakLuminance = 200.0;

	dBlackLevel.SetName(L"Black Level");
	dBlackLevel.SetDescription(L"Lblack: Display black level.");
	dBlackLevel.SetDefault(0.5);
	dBlackLevel.SetRange(0.001, 5.0);
	this->Register(dBlackLevel);
	dBlackLevel = 0.5;

	dGamma.SetName(L"Display Gamma");
	dGamma.SetDescription(L"Gamma: Display gamma, typically 2.2.");
	dGamma.SetDefault(2.2);
	dGamma.SetRange(1.0, 3.0);
	this->Register(dGamma);
	dGamma = 2.2;

	dAmbientLight.SetName(L"Ambient Light");
	dAmbientLight.SetDescription(L"Eamb: Ambient illuminance.");
	dAmbientLight.SetDefault(50.0);
	dAmbientLight.SetRange(0.0, 50000.0);
	this->Register(dAmbientLight);
	dAmbientLight = 50.0;

	dReflectivity.SetName(L"Reflectivity");
	dReflectivity.SetDescription(L"k: Display panel reflectivity, 0.005-0.01 for LCD.");
	dReflectivity.SetDefault(0.01);
	dReflectivity.SetRange(0.001, 0.05);
	this->Register(dReflectivity);
	dReflectivity = 0.01;

	dDetailScaling.SetName(L"Detail Scaling");
	dDetailScaling.SetDescription(L"e: Detail enhancement factor. 1.0=preserve, >1=enhance");
	dDetailScaling.SetDefault(1.0);
	dDetailScaling.SetRange(0.0, 10.0);
	this->Register(dDetailScaling);
	dDetailScaling = 1.0;

	dTonePriority.SetName(L"Tone Priority");
	dTonePriority.SetDescription(L"Balance bright/dark tones. -1=prioritize bright, 0=neutral, 1=prioritize dark");
	dTonePriority.SetDefault(0.0);
	dTonePriority.SetRange(-1.0, 1.0);
	this->Register(dTonePriority);
	dTonePriority = 0.0;

	bLocalToneCurves.SetName(L"Local Tone Curves");
	bLocalToneCurves.SetDescription(L"Use locally adaptive tone curves.");
	bLocalToneCurves.SetDefault(true);
	this->Register(bLocalToneCurves);
	bLocalToneCurves = true;

	iFilterIterations.SetName(L"Filter Iterations");
	iFilterIterations.SetDescription(L"N: Number of diffusion iterations. Paper uses N=12.");
	iFilterIterations.SetDefault(12);
	iFilterIterations.SetRange(1, 30);
	this->Register(iFilterIterations);
	iFilterIterations = 12;

	dSigma.SetName(L"Filter Sigma");
	dSigma.SetDescription(L"Sigma: Starting kernel size for diffusion filter. Paper uses 3.0.");
	dSigma.SetDefault(3.0);
	dSigma.SetRange(0.5, 10.0);
	this->Register(dSigma);
	dSigma = 3.0;

	dEdgeStop.SetName(L"Edge Stop Lambda");
	dEdgeStop.SetDescription(L"Lambda: Edge-stop threshold for Tukey biweight. Paper uses 0.5.");
	dEdgeStop.SetDefault(0.5);
	dEdgeStop.SetRange(0.05, 2.0);
	this->Register(dEdgeStop);
	dEdgeStop = 0.5;

	dNoiseA.SetName(L"Noise A");
	dNoiseA.SetDescription(L"a: Signal-dependent noise parameter (photon noise).");
	dNoiseA.SetDefault(0.001);
	dNoiseA.SetRange(0.0, 1.0);
	this->Register(dNoiseA);
	dNoiseA = 0.001;

	dNoiseB.SetName(L"Noise B");
	dNoiseB.SetDescription(L"b: Signal-independent noise parameter (read-out noise).");
	dNoiseB.SetDefault(0.0001);
	dNoiseB.SetRange(0.0, 0.1);
	this->Register(dNoiseB);
	dNoiseB = 0.0001;

	bEstimateNoise.SetName(L"Estimate Noise");
	bEstimateNoise.SetDescription(L"Automatically estimate noise parameters from image.");
	bEstimateNoise.SetDefault(true);
	this->Register(bEstimateNoise);
	bEstimateNoise = true;

	// --- Temporal parameter ---
	dFrameRate.SetName(L"Frame Rate");
	dFrameRate.SetDescription(L"Video frame rate (fps) for temporal IIR filter. Set 0 to disable temporal filtering.");
	dFrameRate.SetDefault(24.0);
	dFrameRate.SetRange(0.0, 120.0);
	this->Register(dFrameRate);
	dFrameRate = 24.0;

	// --- Initialize persistent state ---
	prevWidth = 0;
	prevHeight = 0;
	mDelta = 0.2;
	mMinLog = 0.0;
	mMaxLog = 6.0;
}

TMOEilertsen15::~TMOEilertsen15() {
}


// Calculate noise variance
double TMOEilertsen15::noiseVariance(double I, const NoiseModel &nm) const {
	return nm.a * std::max(I, 0.0) + nm.b;
}

// Calculate noise
double TMOEilertsen15::noiseLogLevel(double I, const NoiseModel &nm) const {
	if (I <= LUM_EPS) { // very noisy at near-zero luminance
		return 1.0;
	}
	double sigma = std::sqrt(std::max(noiseVariance(I, nm), 0.0));
	return std::log10((I + sigma) / I);
}

// Noise estimation using Foi et al. 2008
TMOEilertsen15::NoiseModel TMOEilertsen15::estimateNoiseModel(const double *luminance, int width, int height) {
	const int blockSize = 8;
	const double uniformThreshold = 0.15;

	// Collect (mean, variance) pairs from uniform blocks
	std::vector<double> means, variances;

	for (int by = 0; by + blockSize <= height; by += blockSize) {
		for (int bx = 0; bx + blockSize <= width; bx += blockSize) {
			// Compute block mean
			double sum = 0.0;
			for (int dy = 0; dy < blockSize; dy++) {
					for (int dx = 0; dx < blockSize; dx++) {
						sum += luminance[(by + dy) * width + (bx + dx)];
					}
			}
			double mean = sum / (blockSize * blockSize);

			if (mean < LUM_EPS) {
				continue;
			}

			// Compute block variance
			double varSum = 0.0;
			for (int dy = 0; dy < blockSize; dy++) {
				for (int dx = 0; dx < blockSize; dx++) {
					double diff = luminance[(by + dy) * width + (bx + dx)] - mean;
					varSum += diff * diff;
				}
			}
			double var = varSum / (blockSize * blockSize - 1);

			// Uniformity check: coefficient of variation
			double cv = std::sqrt(var) / std::max(mean, LUM_EPS);
			if (cv < uniformThreshold) {
				means.push_back(mean);
				variances.push_back(var);
			}
		}
	}

	NoiseModel nm;
	nm.a = dNoiseA;
	nm.b = dNoiseB;

	if (means.size() < 10) { // Not enough data, fall back to defaults
		return nm;
	}

	// Linear regression: variance = a * mean + b
	double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
	int n = (int)means.size();
	for (int i = 0; i < n; i++) {
		sumX += means[i];
		sumY += variances[i];
		sumXX += means[i] * means[i];
		sumXY += means[i] * variances[i];
	}
	double denom = n * sumXX - sumX * sumX;
	if (std::abs(denom) > 1e-12) {
		nm.a = (n * sumXY - sumX * sumY) / denom;
		nm.b = (sumY - nm.a * sumX) / n;
		nm.a = std::max(nm.a, 0.0);
		nm.b = std::max(nm.b, 0.0);
	}

	return nm;
}

// Returns peak contrast sensitivity
double TMOEilertsen15::csfPeak(double Ld) const {
	if (Ld <= 0) {
		Ld = 1e-6;
	}
	double logL = std::log10(Ld);

	// Piecewise linear model
	double logCSF;
	if (logL < -2.0) {
		logCSF = 0.7 + 0.25 * (logL + 2.0);
	} else if (logL < 2.0) {
		logCSF = 0.7 + 0.4 * (logL + 2.0);
	} else {
		logCSF = 2.3;
	}

	logCSF = std::max(0.3, std::min(2.5, logCSF));
	return std::pow(10.0, logCSF);
}

// Calculate visibility threshold
double TMOEilertsen15::visibilityThreshold(double Ld) const {
	double Ct = 1.0 / csfPeak(Ld);
	Ct = std::min(Ct, 0.999);
	return 0.5 * std::log10((Ct + 1.0) / (1.0 - Ct));
}

// Get reflected light
double TMOEilertsen15::reflectedLight() {
	return (double)dReflectivity / M_PI * (double)dAmbientLight;
}

// Calculate display model
double TMOEilertsen15::displayModel(double Lprime) {
	Lprime = std::max(0.0, std::min(1.0, Lprime));
	double Lrefl = reflectedLight();
	double gamma = dGamma.GetDouble();
	double Lmax = dPeakLuminance.GetDouble();
	double Lblack = dBlackLevel.GetDouble();
	return std::pow(Lprime, gamma) * (Lmax - Lblack) + Lblack + Lrefl;
}

// Get the display's dynamic range
double TMOEilertsen15::displayDynamicRange() {
	double Ld1 = displayModel(1.0);
	double Ld0 = displayModel(0.0);
	if (Ld0 <= 0) {
		Ld0 = 1e-6;
	}
	return std::log10(Ld1 / Ld0);
}

// Separable Gaussian blur
void TMOEilertsen15::gaussianBlur(const double *input, double *output, int width, int height, double sigma) const {
	if (sigma < 0.1) {
		std::memcpy(output, input, width * height * sizeof(double));
		return;
	}

	int radius = (int)std::ceil(3.0 * sigma);
	int kernelSize = 2 * radius + 1;

	// Build 1D kernel
	std::vector<double> kernel(kernelSize);
	double sum = 0.0;
	for (int i = 0; i < kernelSize; i++) {
		double x = i - radius;
		kernel[i] = std::exp(-x * x / (2.0 * sigma * sigma));
		sum += kernel[i];
	}
	for (int i = 0; i < kernelSize; i++) {
		kernel[i] /= sum;
	}

	// Temp buffer for horizontal pass
	std::vector<double> temp(width * height);

	// Horizontal pass
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double val = 0.0;
			double wsum = 0.0;
			for (int k = -radius; k <= radius; k++) {
				int xx = std::max(0, std::min(width - 1, x + k));
				double w = kernel[k + radius];
				val += input[y * width + xx] * w;
				wsum += w;
			}
			temp[y * width + x] = val / wsum;
		}
	}

	// Vertical pass
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double val = 0.0;
			double wsum = 0.0;
			for (int k = -radius; k <= radius; k++) {
				int yy = std::max(0, std::min(height - 1, y + k));
				double w = kernel[k + radius];
				val += temp[yy * width + x] * w;
				wsum += w;
			}
			output[y * width + x] = val / wsum;
		}
	}
}

// Tukey's biweight edge-stop function
double TMOEilertsen15::tukeyBiweight(double x, double lambda) const {
	if (std::abs(x) >= lambda) {
		return 0.0;
	}
	double t = x / lambda;
	double u = 1.0 - t * t;
	return u * u;
}

// Gradient magnitude using linear ramp formulation
void TMOEilertsen15::computeGradientMagnitude(const double *input, double *gradMag, int width, int height, int radius) const {
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double gx = 0.0, gy = 0.0;
			for (int d = -radius; d <= radius; d++) {
				int xx = std::max(0, std::min(width - 1, x + d));
				int yy = std::max(0, std::min(height - 1, y + d));
				gx += d * input[y * width + xx];
				gy += d * input[yy * width + x];
			}
			gradMag[y * width + x] = std::sqrt(gx * gx + gy * gy);
		}
	}
}

// Fast detail extraction diffusion
void TMOEilertsen15::detailExtractionDiffusion(const double *input, double *baseLayer, int width, int height, int N, double sigma, double lambda) const {
	int size = width * height;

	std::memcpy(baseLayer, input, size * sizeof(double));

	std::vector<double> ln(size);
	std::vector<double> gradMag(size);

	for (int k = 1; k <= N; k++) {
		double sigma_k = sigma * std::sqrt(2.0 * k - 1.0);

		// Gaussian blur of currently filtered
		gaussianBlur(baseLayer, ln.data(), width, height, sigma_k);

		// Gradient magnitude of lf
		int radius = (int)std::ceil(3.0 * sigma_k);
		radius = std::max(1, radius);
		computeGradientMagnitude(baseLayer, gradMag.data(), width, height, radius);

		// Constrain gradient magnitude
		for (int i = 0; i < size; i++) {
			double constraint = k * std::abs(ln[i] - input[i]);
			gradMag[i] = std::max(gradMag[i], constraint);
		}

		// Edge-stop and iterative update
		for (int i = 0; i < size; i++) {
			double wr = tukeyBiweight(gradMag[i], lambda);
			baseLayer[i] = (1.0 - wr) * baseLayer[i] + wr * ln[i];
		}
	}
}

// Calculate local contrast
void TMOEilertsen15::computeLocalContrast(const double *logLum, int width, int height, double sigma, double *contrast) const {
	int size = width * height;

	std::vector<double> filtered(size);
	std::vector<double> logLumSq(size);
	std::vector<double> filteredSq(size);

	for (int i = 0; i < size; i++) {
		logLumSq[i] = logLum[i] * logLum[i];
	}

	gaussianBlur(logLum, filtered.data(), width, height, sigma);
	gaussianBlur(logLumSq.data(), filteredSq.data(), width, height, sigma);

	for (int i = 0; i < size; i++) {
		double var = filteredSq[i] - filtered[i] * filtered[i];
		contrast[i] = std::sqrt(std::max(0.0, var));
	}
}

// Create a content/noise aware histogram
void TMOEilertsen15::computeContentHistogram(const double *logLum, const double *linLum, int width, int height, const NoiseModel &nm, double minLog, double maxLog, std::vector<double> &histogram) const {
	int N = NUM_SEGMENTS;
	double delta = mDelta;
	int size = width * height;

	histogram.assign(N, 0.0);

	// Local contrast
	std::vector<double> contrast(size);
	computeLocalContrast(logLum, width, height, 3.0, contrast.data());

	// Noise-aware weighted histogram
	double totalWeight = 0.0;

	for (int i = 0; i < size; i++) {
		double l = logLum[i];
		double I = linLum[i];

		//Noise level at this pixel
		double n = noiseLogLevel(I, nm);

		if (contrast[i] <= n) {
			continue;
		}

		// Find histogram bin k
		int k = (int)((l - minLog) / delta);
		if (k < 0 || k >= N) {
			continue;
		}

		// Weight by contrast value
		histogram[k] += contrast[i];
		totalWeight += contrast[i];
	}

	// Normalize to probability distribution
	if (totalWeight > 0) {
		for (int k = 0; k < N; k++) {
			histogram[k] /= totalWeight;
		}
	} else {
		// Fallback to uniform histogram if no pixels pass noise threshold
		for (int k = 0; k < N; k++) {
			histogram[k] = 1.0 / N;
		}
	}
}

// Apply tone priority
void TMOEilertsen15::applyTonePriority(std::vector<double> &histogram, double priority) const {
	if (std::abs(priority) < 1e-6) {
		return;
	}

	int N = (int)histogram.size();
	for (int k = 0; k < N; k++) {
		// Position: 0 = darkest, 1 = brightest
		double position = (double)k / (N - 1);
		double weight;
		if (priority > 0.0) {
			weight = 1.0 + priority * (1.0 - position); // emphasize dark
		} else {
			weight = 1.0 - priority * position; // emphasize bright
		}
		histogram[k] *= weight;
	}

	// Renormalize
	double sum = std::accumulate(histogram.begin(), histogram.end(), 0.0);
	if (sum > 0) {
		for (auto &h : histogram) {
			h /= sum;
		}
	}
}

// Minimum contrast distortion tone curve
std::vector<double> TMOEilertsen15::computeSlopes(const std::vector<double> &histogram, double r) const {
	int N = NUM_SEGMENTS;
	double delta = mDelta;

	std::vector<double> slopes(N, 0.0);

	// Iterative threshold to enforce sk >= 0
	// Start with a small threshold p0 = 0.0001
	double pt = 0.0001;
	const int maxIter = 20;

	for (int iter = 0; iter < maxIter; iter++) {
		std::vector<int> activeSet;
		for (int k = 0; k < N; k++) {
			if (histogram[k] > pt) {
				activeSet.push_back(k);
			}
		}

		if (activeSet.empty()) {
			break;
		}

		double invSum = 0.0;
		for (int k : activeSet) {
			invSum += 1.0 / histogram[k];
		}

		double newPt = ((double)activeSet.size() - r / delta) / invSum;

		if (std::abs(newPt - pt) < 1e-10) {
			break;
		}
		pt = newPt;
	}

	// Final active set with converged threshold
	std::vector<int> activeSet;
	for (int k = 0; k < N; k++) {
		if (histogram[k] > pt) {
			activeSet.push_back(k);
		}
	}

	if (activeSet.empty()) { // Fallback: all segments active with uniform slope
		double uniformSlope = r / (delta * N);
		std::fill(slopes.begin(), slopes.end(), uniformSlope);
		return slopes;
	}

	// Compute sum of reciprocals for active set
	double invSum = 0.0;
	for (int k : activeSet) {
		invSum += 1.0 / histogram[k];
	}

	int M = (int)activeSet.size();

	for (int k = 0; k < N; k++) {
		slopes[k] = 0.0;
	}

	for (int k : activeSet) {
		double sk = 1.0 + (r / delta - M) / (histogram[k] * invSum);
		slopes[k] = std::max(0.0, sk);
	}

	// Verify total slope budget: sum(sk * delta) <= r
	double totalSlope = 0.0;
	for (int k = 0; k < N; k++) {
		totalSlope += slopes[k] * delta;
	}

	// Rescale if exceeding budget (for safety)
	if (totalSlope > r + 1e-6) {
		double scale = r / totalSlope;
		for (int k = 0; k < N; k++) {
			slopes[k] *= scale;
		}
	}

	return slopes;
}

// Build tone curve node values from slopes
std::vector<double> TMOEilertsen15::buildToneCurveNodes(const std::vector<double> &slopes, double r) const {
	int N = NUM_SEGMENTS;
	double delta = mDelta;

	std::vector<double> nodes(N + 1);
	nodes[0] = -r; // minimum display log-luminance
	for (int k = 0; k < N; k++) {
		nodes[k + 1] = nodes[k] + slopes[k] * delta;
	}

	return nodes;
}

// Apply piecewise linear tone curve to a single log-luminance value
double TMOEilertsen15::applyToneCurve(double logLum, double minLog, const std::vector<double> &nodeOutputs) const {
	int N = NUM_SEGMENTS;
	double delta = mDelta;

	// Find segment
	double pos = (logLum - minLog) / delta;
	int segment = (int)std::floor(pos);
	segment = std::max(0, std::min(N - 1, segment));

	// Linear interpolation within segment
	double t = pos - segment;
	t = std::max(0.0, std::min(1.0, t));

	return nodeOutputs[segment] + t * (nodeOutputs[segment + 1] - nodeOutputs[segment]);
}

// Calculate and apply local tone curves if enabled
void TMOEilertsen15::computeAndApplyLocalToneCurves(double *baseLayer, const double *logLum, const double *linLum, int width, int height, const NoiseModel &nm, double minLog, double maxLog, double r, double priority) {
	int N = NUM_SEGMENTS;
	int tileSize = DEFAULT_TILE_SIZE;
	int tilesX = (width + tileSize - 1) / tileSize;
	int tilesY = (height + tileSize - 1) / tileSize;

	// Compute IIR alpha for temporal filtering
	double fps = dFrameRate;
	double alpha = 1.0; // no filtering by default
	if (fps > 0) {
		alpha = 1.0 - std::exp(-2.0 * M_PI * TEMPORAL_CUTOFF_HZ / fps);
   }

	// Compute global histogram
	std::vector<double> globalHist;
	computeContentHistogram(logLum, linLum, width, height, nm, minLog, maxLog, globalHist);
	applyTonePriority(globalHist, priority);

	// Reset temporal cache if dimensions changed
	if (width != prevWidth || height != prevHeight) {
		prevLocalNodes.clear();
		prevGlobalNodes.clear();
		prevWidth = width;
		prevHeight = height;
	}

	// Compute per-tile tone curves
	int totalTiles = tilesX * tilesY;
	std::vector<std::vector<double>> tileNodes(totalTiles);

	for (int ty = 0; ty < tilesY; ty++) {
		for (int tx = 0; tx < tilesX; tx++) {
			int tileIdx = ty * tilesX + tx;

			// Tile bounds
			int x0 = tx * tileSize;
			int y0 = ty * tileSize;
			int x1 = std::min(x0 + tileSize, width);
			int y1 = std::min(y0 + tileSize, height);
			int tileW = x1 - x0;
			int tileH = y1 - y0;
			int tilePixels = tileW * tileH;

			// Extract tile data
			std::vector<double> tileLog(tilePixels);
			std::vector<double> tileLin(tilePixels);
			for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
               int ti = (y - y0) * tileW + (x - x0);
               int si = y * width + x;
               tileLog[ti] = logLum[si];
               tileLin[ti] = linLum[si];
            }
			}

			// Compute tile histogram
			std::vector<double> tileHist;
			computeContentHistogram(tileLog.data(), tileLin.data(), tileW, tileH, nm, minLog, maxLog, tileHist);

			// Blend: 10% global + 90% local
			std::vector<double> blendedHist(N);
			for (int k = 0; k < N; k++) {
				blendedHist[k] = LOCAL_GLOBAL_RATIO * globalHist[k] + (1.0 - LOCAL_GLOBAL_RATIO) * tileHist[k];
			}

			applyTonePriority(blendedHist, priority);

			// Compute tone curve slopes and nodes
			std::vector<double> slopes = computeSlopes(blendedHist, r);
			tileNodes[tileIdx] = buildToneCurveNodes(slopes, r);

			// Temporal IIR filtering on node positions
			if ((int)prevLocalNodes.size() == totalTiles && alpha < 1.0) {
            for (int k = 0; k <= N; k++) {
               tileNodes[tileIdx][k] = alpha * tileNodes[tileIdx][k] + (1.0 - alpha) * prevLocalNodes[tileIdx][k];
            }
			}
		}
	}

	// Store for next frame's temporal filtering
	prevLocalNodes = tileNodes;

	// Apply local tone curves with bilinear interpolation between tiles
	// Each tile's tone curve is centered at the tile's center
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int idx = y * width + x;
			double logVal = baseLayer[idx];

			// Find position relative to tile centers
			double tileXf = ((double)x / tileSize) - 0.5;
			double tileYf = ((double)y / tileSize) - 0.5;

			int tx0 = (int)std::floor(tileXf);
			int ty0 = (int)std::floor(tileYf);
			int tx1 = tx0 + 1;
			int ty1 = ty0 + 1;

			double fx = tileXf - tx0;
			double fy = tileYf - ty0;

			// Clamp tile indices
			tx0 = std::max(0, std::min(tilesX - 1, tx0));
			tx1 = std::max(0, std::min(tilesX - 1, tx1));
			ty0 = std::max(0, std::min(tilesY - 1, ty0));
			ty1 = std::max(0, std::min(tilesY - 1, ty1));

			// Bilinear interpolation of tone-mapped values from 4 nearest tiles
			double v00 = applyToneCurve(logVal, minLog, tileNodes[ty0 * tilesX + tx0]);
			double v10 = applyToneCurve(logVal, minLog, tileNodes[ty0 * tilesX + tx1]);
			double v01 = applyToneCurve(logVal, minLog, tileNodes[ty1 * tilesX + tx0]);
			double v11 = applyToneCurve(logVal, minLog, tileNodes[ty1 * tilesX + tx1]);

			baseLayer[idx] = (1.0 - fx) * (1.0 - fy) * v00 + fx * (1.0 - fy) * v10 + (1.0 - fx) * fy * v01 + fx * fy * v11;
		}
	}
}

// Noise aware detail control
void TMOEilertsen15::noiseAwareDetailControl(double *detailLayer, const double *baseOrig, const double *baseToneMapped, const double *linLuminance, int width, int height, const NoiseModel &nm, double enhancement) const {
	int size = width * height;

	double peakLd = 200.0;

	for (int i = 0; i < size; i++) {
		// Noise level from original base (log domain, depends on input luminance)
		double I = linLuminance[i];
		double n = noiseLogLevel(I, nm);

		// Visibility threshold at displayed luminance
		double displayLum = peakLd * std::pow(10.0, baseToneMapped[i]);
		displayLum = std::max(displayLum, LUM_EPS);
		double V = visibilityThreshold(displayLum);

		// Modulation factor
		double modulation = 1.0;
		if (n > 0) {
			modulation = std::min(1.0, V / n);
      }

		detailLayer[i] = enhancement * modulation * detailLayer[i];
	}
}


// Main Transform function (Figure 4.1 pipeline)
int TMOEilertsen15::Transform() {
	// Work in RGB, extract luminance manually for full control
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double *pSourceData = pSrc->GetData();
	double *pDestData = pDst->GetData();

	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();
	int size = width * height;

	// Read parameters
	int N = iFilterIterations;
	double sigma = dSigma;
	double lambda = dEdgeStop;
	double enhancement = dDetailScaling;
	double priority = dTonePriority;
	bool useLocal = bLocalToneCurves;

	pSrc->ProgressBar(0, 100);

	// Step 1: Extract luminance from RGB
	std::vector<double> srcR(size), srcG(size), srcB(size);
	std::vector<double> luminance(size);

	for (int i = 0; i < size; i++) {
		srcR[i] = pSourceData[i * 3 + 0];
		srcG[i] = pSourceData[i * 3 + 1];
		srcB[i] = pSourceData[i * 3 + 2];

		// CIE Y luminance
		luminance[i] = 0.2126 * srcR[i] + 0.7152 * srcG[i] + 0.0722 * srcB[i];
		luminance[i] = std::max(luminance[i], LUM_EPS);
	}

	pSrc->ProgressBar(5, 100);

	// Step 2: Noise model
	NoiseModel nm;
	if (bEstimateNoise) {
		nm = estimateNoiseModel(luminance.data(), width, height);
   } else {
		nm.a = dNoiseA;
		nm.b = dNoiseB;
	}

	pSrc->ProgressBar(10, 100);

	// Step 3: Transform to log domain
	std::vector<double> logLum(size);
	for (int i = 0; i < size; i++) {
		logLum[i] = std::log10(luminance[i]);
   }

	// Compute log-luminance range for tone curve binning
	double minLog = *std::min_element(logLum.begin(), logLum.end());
	double maxLog = *std::max_element(logLum.begin(), logLum.end());

	// Ensure a minimum range to avoid degenerate cases
	double dataRange = maxLog - minLog;
	if (dataRange < 0.1) {
		double mid = (maxLog + minLog) / 2.0;
		minLog = mid - 0.5;
		maxLog = mid + 0.5;
		dataRange = 1.0;
	}

	// Adapt bin width
	mDelta = dataRange / NUM_SEGMENTS;
	mMinLog = minLog;
	mMaxLog = maxLog;

	pSrc->ProgressBar(15, 100);

	// Step 4: Edge-stopping spatial filter
	std::vector<double> baseLayer(size);
	std::vector<double> detailLayer(size);

	detailExtractionDiffusion(logLum.data(), baseLayer.data(), width, height, N, sigma, lambda);

	// Detail layer
	for (int i = 0; i < size; i++) {
		detailLayer[i] = logLum[i] - baseLayer[i];
   }

	pSrc->ProgressBar(50, 100);


	// Step 5: Local tone curves
	double r = displayDynamicRange();

	// Keep a copy of original base for noise-aware detail control
	std::vector<double> baseOrig(baseLayer.begin(), baseLayer.end());

	if (useLocal) { // Use local tone curves with tile-based processing
		computeAndApplyLocalToneCurves(baseLayer.data(), logLum.data(), luminance.data(), width, height, nm, minLog, maxLog, r, priority);
	} else { // Use global tone curve
		std::vector<double> histogram;
		computeContentHistogram(logLum.data(), luminance.data(), width, height, nm, minLog, maxLog, histogram);
		applyTonePriority(histogram, priority);

		std::vector<double> slopes = computeSlopes(histogram, r);
		std::vector<double> nodes = buildToneCurveNodes(slopes, r);

		// Temporal IIR filtering
		double fps = dFrameRate;
		if (fps > 0 && !prevGlobalNodes.empty()) {
			double alpha = 1.0 - std::exp(-2.0 * M_PI * TEMPORAL_CUTOFF_HZ / fps);
			for (int k = 0; k <= NUM_SEGMENTS; k++) {
				nodes[k] = alpha * nodes[k] + (1.0 - alpha) * prevGlobalNodes[k];
         }
		}
		prevGlobalNodes = nodes;

		// Apply tone curve to base layer
		for (int i = 0; i < size; i++) {
			baseLayer[i] = applyToneCurve(baseLayer[i], minLog, nodes);
      }
	}

	pSrc->ProgressBar(70, 100);

	// Step 6: Noise-aware detail control
	noiseAwareDetailControl(detailLayer.data(), baseOrig.data(), baseLayer.data(), luminance.data(), width, height, nm, enhancement);

	pSrc->ProgressBar(80, 100);


	// Step 6: Combine base and detail, restore color
	double LdMax = displayModel(1.0);  // peak display luminance
	double invLrefl = reflectedLight();
	double invLblack = dBlackLevel.GetDouble();
	double invLmax = dPeakLuminance.GetDouble();
	double invGamma = dGamma.GetDouble();
	double invDenom = invLmax - invLblack;
	if (invDenom <= 0) {
      invDenom = 1.0;
   }
	double invGammaExp = 1.0 / invGamma;

	for (int i = 0; i < size; i++) {
		double tmLogLum = baseLayer[i] + detailLayer[i];

		// Convert from normalized display log-luminance to actual Ld
		double Ld = LdMax * std::pow(10.0, tmLogLum);

		// Inverse display model
		double numerator = std::max(0.0, Ld - invLblack - invLrefl);
		double Lprime = std::pow(numerator / invDenom, invGammaExp);
		Lprime = std::max(0.0, std::min(1.0, Lprime));

		// Color restoration as per Mantiuk et al. 2009
		double Y = luminance[i];
		double s = COLOR_SATURATION;

		double rRatio = srcR[i] / Y;
		double gRatio = srcG[i] / Y;
		double bRatio = srcB[i] / Y;

		pDestData[i * 3 + 0] = std::max(0.0, std::min(1.0, std::pow(rRatio, s) * Lprime));
		pDestData[i * 3 + 1] = std::max(0.0, std::min(1.0, std::pow(gRatio, s) * Lprime));
		pDestData[i * 3 + 2] = std::max(0.0, std::min(1.0, std::pow(bRatio, s) * Lprime));
	}

	pSrc->ProgressBar(100, 100);
	return 0;
}
