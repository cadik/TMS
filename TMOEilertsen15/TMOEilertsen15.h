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
*                  Implementation of the TMOEilertsen15 header                 *
*                        Real-time noise-aware tone mapping                    *
*                                                                              *
*******************************************************************************/

#ifndef TMO_EILERTSEN15_H
#define TMO_EILERTSEN15_H

#include "TMO.h"
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <numeric>

class TMOEilertsen15 : public TMO {
public:
	TMOEilertsen15();
	virtual ~TMOEilertsen15();
	virtual int Transform();

protected:
	struct NoiseModel {
		double a;
		double b;
	};

	double noiseVariance(double I, const NoiseModel &nm) const;

	double noiseLogLevel(double I, const NoiseModel &nm) const;

	// Estimate noise parameters from image
	NoiseModel estimateNoiseModel(const double *luminance, int width, int height);

	double visibilityThreshold(double Ld) const;

	// Peak CSF at adaptation luminance
	double csfPeak(double Ld) const;

	double displayModel(double Lprime);

	double reflectedLight();

	double displayDynamicRange();

	// Compute optimal slopes sk for given histogram p(lk)
	std::vector<double> computeSlopes(const std::vector<double> &histogram, double r) const;

	// Build piecewise linear tone curve node values from slopes
	std::vector<double> buildToneCurveNodes(const std::vector<double> &slopes, double r) const;

	// Apply piecewise linear tone curve to a single value
	double applyToneCurve(double logLum, double minLog, const std::vector<double> &nodeOutputs) const;

	// Compute content-weighted, noise-aware histogram
	void computeContentHistogram(const double *logLum, const double *linLum, int width, int height, const NoiseModel &nm, double minLog, double maxLog, std::vector<double> &histogram) const;

	// Local contrast as Gaussian-windowed standard deviation
	void computeLocalContrast(const double *logLum, int width, int height, double sigma, double *contrast) const;

	// Tone priority weighting
	void applyTonePriority(std::vector<double> &histogram, double priority) const;

	// Compute per-tile histograms, blend with global, compute local tone curves + apply with bilinear interpolation between tiles
	void computeAndApplyLocalToneCurves(double *baseLayer, const double *logLum, const double *linLum, int width, int height, const NoiseModel &nm, double minLog, double maxLog, double r, double priority);

	// Fast detail extraction diffusion
	void detailExtractionDiffusion(const double *input, double *baseLayer, int width, int height, int N, double sigma, double lambda) const;

	// Separable Gaussian blur
	void gaussianBlur(const double *input, double *output, int width, int height, double sigma) const;

	// Gradient magnitude with linear ramp formulation
	void computeGradientMagnitude(const double *input, double *gradMag, int width, int height, int radius) const;

	// Tukey's biweight edge-stop function
	double tukeyBiweight(double x, double lambda) const;

	void noiseAwareDetailControl(double *detailLayer, const double *baseOrig, const double *baseToneMapped, const double *linLuminance, int width, int height, const NoiseModel &nm, double enhancement) const;

private:
	// Display
	TMODouble dPeakLuminance;  // Lmax (default 200)
	TMODouble dBlackLevel;     // Lblack (default 0.5)
	TMODouble dGamma;          // Display gamma (default 2.2)
	TMODouble dAmbientLight;   // Eamb (default 50)
	TMODouble dReflectivity;   // Panel reflectivity (default 0.01)

	// Processing 
	TMODouble dDetailScaling;   // Detail enhancement (default 1.0)
	TMODouble dTonePriority;    // Tone priority: -1=bright, 0=neutral, 1=dark
	TMOBool bLocalToneCurves;   // Use local tone curves (default true)
	TMOInt iFilterIterations;   // N: diffusion iterations (default 12)
	TMODouble dSigma;           // sigma: starting kernel size (default 3.0)
	TMODouble dEdgeStop;        // lambda: edge-stop threshold (default 0.5)

	// Noise
	TMODouble dNoiseA;          // a: signal-dependent noise
	TMODouble dNoiseB;          // b: signal-independent noise
	TMOBool bEstimateNoise;     // auto-estimate noise

	// Temporal
	TMODouble dFrameRate;       // fps for IIR filter (default 24)

	// Constants
	static const int NUM_SEGMENTS = 30;               // N segments in tone curve
	static constexpr double TEMPORAL_CUTOFF_HZ = 0.5;  // IIR cutoff frequency
	static const int DEFAULT_TILE_SIZE = 230;           // approx. 5 visual degrees
	static constexpr double LOCAL_GLOBAL_RATIO = 0.1;   // 10% global, 90% local

	double mDelta;
	double mMinLog, mMaxLog;

	// Persistent state across frames
	std::vector<double> prevGlobalNodes;
	// For local tone curves (per tile):
	std::vector<std::vector<double>> prevLocalNodes;
	int prevWidth, prevHeight;

	static constexpr double COLOR_SATURATION = 0.6;
};

#endif // TMO_EILERTSEN15_H
