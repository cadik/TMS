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
*                      Implementation of the TMOShan12 header                  *
*               Tone Mapping High Dynamic Range Videos using Wavelets          *
*                                                                              *
*******************************************************************************/

#ifndef TMO_SHAN12_H
#define TMO_SHAN12_H

#include "TMO.h"
#include <vector>
#include <deque>
#include <cmath>

class TMOShan12 : public TMO {
public:
    TMOShan12();
    virtual ~TMOShan12();
    virtual int Transform();

protected:
    
    // Wavelet sub-band data for a single frame
    struct WaveletPyramid {
        std::vector<std::vector<double>> bands;    // High-pass detail bands
        std::vector<double> lowpass;               // Coarsest low-pass signal
        int width;
        int height;
        int numBands;
    };

    // Temporal cache: stores log-luminance of recent frames
    struct TemporalFrameCache {
        std::deque<std::vector<double>> frameHistory;
        int maxFrames;
        int cachedWidth;
        int cachedHeight;

        TemporalFrameCache() : maxFrames(5), cachedWidth(0), cachedHeight(0) {}

        void clear() {
            frameHistory.clear();
            cachedWidth = 0;
            cachedHeight = 0;
        }
    };

    TemporalFrameCache temporalCache;
    int lastWidth = 0;
    int lastHeight = 0;

    // 2D redundant wavelet analysis with B-spline wavelets
    void waveletAnalysis2D(const std::vector<double> &input, int width, int height, int numBands, WaveletPyramid &pyramid);
    
    // 2D redundant wavelet synthesis
    static void waveletSynthesis2D(const WaveletPyramid &pyramid, std::vector<double> &output);
    
    // 1D separable low-pass filtering with B-spline filter [0.5, 1, 0.5]
    void lowPassFilter1D(const std::vector<double> &input, std::vector<double> &output, int length, int stride, int numLines, int lineStride, int level);
    
    // 2D low-pass filter
    void lowPassFilter2D(const std::vector<double> &input, std::vector<double> &output, int width, int height, int level);
    
    // Temporal edge-avoiding low-pass filter
    void temporalEdgeAvoidingLowPass(const std::deque<std::vector<double>> &frames, int centerIdx, std::vector<double> &output, int width, int height, double alpha);
    
    // Temporal edge-avoiding decomposition for temporally coherent wavelet bands
    void temporalEdgeAvoidingDecomposition(int width, int height, int numBands, double alpha, WaveletPyramid &pyramid);
    
    // Compute band activity from Gaussian-blurred absolute values
    void computeBandActivity(const std::vector<double> &band, std::vector<double> &activity, int width, int height, int level);
    
    // Compute aggregate activity
    void computeAggregateActivity(const std::vector<std::vector<double>> &activities, std::vector<double> &aggregate, int width, int height);
    
    // Compute Li et al. 2005 gain map
    void computeGainMapLi(const std::vector<double> &aggregate, int width, int height, int bandIdx, int numBands, double gamma, double delta, std::vector<double> &gainMap);
    
    // Compute improved gain map
    void computeGainMapImproved(const std::vector<double> &aggregate, int width, int height, int bandIdx, int numBands, double gamma, double delta, std::vector<double> &gainMap);
    
    // Level-dependent multipliers
    static double getLevelMultiplier(int bandIdx, int numBands);
    
    // Binary edge function
    inline bool isEdge(double val1, double val2, double alpha) const {
        return std::fabs(val1 - val2) >= alpha;
    }
    
    // Gaussian blur with sigma proportional to wavelet level
    void gaussianBlur(const std::vector<double> &input, std::vector<double> &output, int width, int height, double sigma);

private:

    // gamma: Controls linearity of gain function
    TMODouble dGamma;
    
    // delta: Normalization constant for gain function
    TMODouble dDelta;
    
    // N: Number of wavelet decomposition bands
    TMOInt iNumBands;
    
    // alpha: Edge detection threshold for temporal edge avoidance
    TMODouble dAlpha;
    
    // Enable improved gain control
    TMOBool bImprovedGain;
    
    // Number of temporal frames for 3D wavelet (window size)
    TMOInt iTemporalWindow;

    static constexpr double EPSILON = 1e-6;
    static constexpr double LP_FILTER[3] = {0.5, 1.0, 0.5};
    static constexpr double HP_FILTER[5] = {1.0 / 3.0, -6.0 / 3.0, 10.0 / 3.0, -6.0 / 3.0, 1.0 / 3.0};
    static constexpr double RECON_LP[3] = {0.5, 1.0, 0.5};
    static constexpr double RECON_HP[5] = {1.0 / 3.0, -6.0 / 3.0, 10.0 / 3.0, -6.0 / 3.0, 1.0 / 3.0};
};

#endif // TMO_SHAN12_H