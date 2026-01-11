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
*                    Implementation of the TMOEilertsen15 header               *
*                         Real-time noise-aware tone mapping                   *
*                     https://computergraphics.on.liu.se/rntm/                 *
*                                                                              *
*******************************************************************************/

#ifndef TMO_EILERTSEN15_H
#define TMO_EILERTSEN15_H

#include "TMO.h"
#include <vector>
#include <cmath>
#include <deque>

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
    
    // Compute noise variance
    double computeNoiseVariance(double luminance, const NoiseModel& noise);
    
    // Compute noise in log domain
    double computeNoiseLevel(double luminance, const NoiseModel& noise);
    
    // Estimate noise using Foi et al. 2008 method
    NoiseModel estimateNoise();
    void fitLinearModel(const std::vector<std::pair<double, double>>& samples, double& a, double& b);

    // HDR-VDP-2 contrast sensitivity function
    double computeContrastSensitivity(double luminance, double frequency);
    
    // Detection threshold
    double computeDetectionThreshold(double displayLuminance);
    
    // Visibility threshold
    double computeVisibilityThreshold(double displayLuminance);

    // Forward: Ld(L') = (L')^γ (Lmax - Lblack) + Lblack + Lrefl
    double displayModel(double pixelValue);
    
    // Inverse display model
    double inverseDisplayModel(double displayLuminance);
    
    // Reflected ambient light
    double computeReflectedLight();
    
    // Dynamic range
    double computeDisplayDynamicRange();

    // Compute tone curve slopes with iterative threshold
    void computeToneCurveSlopes(const std::vector<double>& histogram, const NoiseModel& noise, std::vector<double>& slopes);
    
    // Build and apply piecewise linear tone curve
    void applyToneCurve(double* baseLayer, int width, int height, const std::vector<double>& slopes, double minLog, double maxLog);

    // Local contrast
    void computeLocalContrast(const double* logLuminance, int width, int height, double* contrast);
    
    // Weighted histogram based on contrast above noise
    void computeNoiseAwareHistogram(const double* logLuminance, const double* contrast, int width, int height, const NoiseModel& noise, std::vector<double>& histogram, double minLog, double maxLog);

    struct IIRFilter {
        std::deque<std::vector<double>> history;
        double alpha;
        int maxHistory;
        IIRFilter() : alpha(0.105), maxHistory(10) {}
        void reset() { history.clear(); }
    };
    
    std::vector<IIRFilter> temporalFilters;
    int lastTilesX = 0;
    int lastTilesY = 0;
    
    // Apply temporal filtering to tone curve node positions
    void applyTemporalFilter(std::vector<double>& nodeValues, IIRFilter& filter);

    // Compute local tone curves with 90% local + 10% global blending
    void computeLocalToneCurves(const double* logLuminance, const double* contrast, int width, int height, const NoiseModel& noise, std::vector<std::vector<double>>& localSlopes, int& tilesX, int& tilesY, double& minLog, double& maxLog);
    
    // Apply local tone curves with bilinear interpolation
    void applyLocalToneCurves(double* baseLayer, int width, int height, const std::vector<std::vector<double>>& localSlopes, int tilesX, int tilesY, double minLog, double maxLog);

    // Main diffusion algorithm
    void fastDetailExtractionDiffusion(const double* input, double* output, int width, int height);
    
    // Separable Gaussian blur
    void gaussianBlur(const double* input, double* output, int width, int height, double sigma);
    
    // Gradient magnitude with linear ramp
    void computeGradientMagnitude(const double* input, double* gradMag, int width, int height, int radius);
    
    // Tukey's biweight edge-stop function
    double tukeyBiweight(double x, double lambda);

    // dout = e · min(1, V(btm)/n(b)) · d
    void noiseAwareDetailControl(double* detailLayer, const double* baseOrig, const double* baseTM, int width, int height, const NoiseModel& noise);

    // Block-based operations for noise estimation
    double computeBlockMean(const double* luminance, int width, int height, int bx, int by, int blockSize);
    double computeBlockVariance(const double* luminance, int width, int height, int bx, int by, int blockSize, double mean);
    bool isUniformBlock(const double* luminance, int width, int height, int bx, int by, int blockSize);

private:

    // Display parameters
    TMODouble dPeakLuminance;      // Lmax: Peak display luminance [cd/m2]
    TMODouble dBlackLevel;         // Lblack: Display black level [cd/m2]
    TMODouble dGamma;              // γ: Display gamma (typically 2.2)
    TMODouble dAmbientLight;       // Eamb: Ambient illuminance [lux]
    TMODouble dReflectivity;       // k: Display reflectivity [0-1%]
    
    // Processing parameters
    TMODouble dDetailScaling;      // e: Detail enhancement factor (Eq. 28)
    TMODouble dNoiseControl;       // Additional noise visibility control
    TMODouble dTonePriority;       // Tone priority: -1=bright, 0=neutral, 1=dark
    TMOBool bLocalToneCurves;      // Use local tone curves (Section 4.3)
    TMOInt iFilterIterations;      // N: Number of diffusion iterations (default 12)
    TMODouble dSigma;              // σ: Starting kernel size for diffusion (default 3.0)
    
    // Noise model parameters
    TMODouble dNoiseA;             // a: Signal-dependent noise parameter
    TMODouble dNoiseB;             // b: Signal-independent noise parameter
    TMOBool bEstimateNoise;        // Automatically estimate noise (Foi et al. 2008)

    static const int NUM_TONE_CURVE_NODES = 30;
    static constexpr double TONE_CURVE_BIN_WIDTH = 0.2;
    static const int DEFAULT_TILE_SIZE = 230;
    static constexpr double EDGE_STOP_LAMBDA = 0.5;
    static constexpr double TEMPORAL_CUTOFF_HZ = 0.5;
    
    void applyTonePriority(std::vector<double>& histogram);
};

#endif // TMO_EILERTSEN15_H