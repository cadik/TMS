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
*                    Implementation of the TMOAydin14 header                   *
*            Temporally Coherent Local Tone Mapping of HDR Video               *
*                                                                              *
*******************************************************************************/

#ifndef TMO_Aydin14_H
#define TMO_Aydin14_H

#include "TMO.h"
#include <vector>
#include <deque>
#include <cmath>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include <opencv2/opencv.hpp>
#include <opencv2/video/tracking.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#undef EPS

class TMOAydin14 : public TMO {
public:
    TMOAydin14();
    virtual ~TMOAydin14();
    virtual int Transform();

protected:
    // Per-frame cached data
    struct FrameData {
        std::vector<double> logLum;          // log10(luminance) per pixel
        std::vector<double> logR, logG, logB;// log10(R), log10(G), log10(B)
        std::vector<double> spatialFiltered; // Spatially filtered log luminance
        cv::Mat grayFloat;                   // Normalized luminance [0,1] CV_32F for optical flow
        int width, height;
    };

    // Temporal frame cache for video processing
    struct TemporalCache {
        std::deque<FrameData> frames;
        // Forward optical flow: forwardFlows[i] = flow from frames[i] to frames[i+1]
        std::deque<cv::Mat> forwardFlows;
        int maxFrames;
        int cachedWidth, cachedHeight;

        TemporalCache() : maxFrames(21), cachedWidth(0), cachedHeight(0) {}

        void clear() {
            frames.clear();
            forwardFlows.clear();
            cachedWidth = 0;
            cachedHeight = 0;
        }
    };

    TemporalCache cache;
    int lastWidth = 0;
    int lastHeight = 0;


    // Compute horizontal permeability - Lorentzian edge-stopping
    void computePermeabilityH(const std::vector<double>& input, std::vector<double>& perm, int width, int height, double sigma, double alpha);

    // Compute vertical permeability (vertically)
    void computePermeabilityV(const std::vector<double>& input, std::vector<double>& perm, int width, int height, double sigma, double alpha);

    // Single spatial filter iteration using permeability-based diffusion
    void spatialFilterIteration(const std::vector<double>& original, std::vector<double>& current, const std::vector<double>& perm, int width, int height, bool horizontal);

    // Full iterative spatial filtering - alternates H/V for numIter iterations
    void spatialFilter(const std::vector<double>& input, std::vector<double>& output, int width, int height, double sigma, double alphaEdge, int numIter);

    // Compute dense optical flow between two frames (OpenCV Farneback)
    void computeOpticalFlow(const cv::Mat& prev, const cv::Mat& next, cv::Mat& flow);

    // Accumulate flow vectors from frame at index "fromIdx" to "toIdx"
    cv::Mat accumulateFlow(int fromIdx, int toIdx);

    // Warp a single-channel frame using flow field (bilinear interpolation)
    void warpFrame(const std::vector<double>& src, std::vector<double>& dst, const cv::Mat& flow, int width, int height);

    // Compute temporal permeability between two adjacent warped frames
    void computeTemporalPermeability(const std::vector<double>& frame1, const std::vector<double>& frame2, const cv::Mat& flow, std::vector<double>& perm, int width, int height, double sigmaPC, double sigmaFG, double alpha);

    // 1D temporal filtering along warped temporal path
    void temporalFilter1D(const std::vector<std::vector<double>>& warpedFrames, const std::vector<std::vector<double>>& temporalPerms, int centerIdx, std::vector<double>& output, int numPixels);

    // Log-domain compression factor (Durand & Dorsey 2002)
    void applyCompressionFactor(std::vector<double>& base, int numPixels, double c);

    // Drago et al. 2003 adaptive logarithmic tone curve
    void applyDragoToneCurve(std::vector<double>& base, int numPixels, double bias, double Ldmax);

    // Lorentzian edge-stopping function
    inline double lorentzian(double diff, double sigma, double alpha) const {
        double ratio = std::fabs(diff) / sigma;
        return 1.0 / (1.0 + std::pow(ratio, alpha));
    }

private:

    // Spatial permeability sigma: controls edge sensitivity
    TMODouble dSigmaS;

    // Edge-stopping exponent alpha
    TMODouble dAlpha;

    // Number of spatial filtering iterations
    TMOInt iNumIter;

    // Temporal radius: frames in each direction
    TMOInt iTemporalRadius;

    // Photo-constancy sigma for temporal permeability
    TMODouble dSigmaPC;

    // Flow gradient sigma for temporal permeability
    TMODouble dSigmaFG;

    // Compression factor for log-domain base compression
    TMODouble dCompression;

    // Drago bias parameter
    TMODouble dDragoBias;

    // Toggle between compression factor and Drago tone curve
    TMOBool bUseDrago;

    // Saturation control (0.0 = grayscale, 1.0 = full color)
    TMODouble dSaturation;

    // Display gamma
    TMODouble dGamma;

    static constexpr double EPSILON = 1e-6;
};

#endif // TMO_Aydin14_H