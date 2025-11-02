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

class TMOEilertsen15 : public TMO {
    public:
        TMOEilertsen15();
        virtual ~TMOEilertsen15();
        virtual int Transform();

    protected:
        // Parameters
        TMODouble dPeakLuminance;      // Peak display luminance [cd/m2]
        TMODouble dBlackLevel;         // Display black level [cd/m2]
        TMODouble dAmbientLight;       // Ambient light [lux]
        TMODouble dReflectivity;       // Display reflectivity [0-1%]
        TMODouble dDetailScaling;      // Detail enhancement factor
        TMODouble dNoiseControl;       // Noise visibility control
        TMODouble dTonePriority;       // Tone priority <-1; 1>
        TMOBool bLocalToneCurves;      // Use local tone curves
        TMOInt iFilterIterations;      // Number of diffusion iterations

    struct NoiseModel {
        double a;
        double b;
    };

    double computeNoiseLevel(double luminance, const NoiseModel& noise);
    double computeVisibilityThreshold(double displayLuminance);
    double computeDisplayDynamicRange();
    
    double displayModel(double pixelValue);
    double inverseDisplayModel(double luminance);
    
    void computeToneCurve(const std::vector<double>& histogram, const NoiseModel& noise, std::vector<double>& toneCurveSlopes);
    void computeLocalToneCurves(double* luminance, int width, int height, std::vector<std::vector<double>>& localCurves);
    
    void applyToneCurve(double* base, int width, int height, const std::vector<double>& slopes);
    void applyLocalToneCurves(double* base, int width, int height, const std::vector<std::vector<double>>& localCurves);
    
    void detailExtractionDiffusion(double* input, double* output, int width, int height, int iterations, double sigma);
    void gaussianBlur(double* input, double* output, int width, int height, double sigma);
    void computeGradientMagnitude(double* input, double* gradMag, int width, int height, int radius);
    
    void noiseAwareDetailControl(double* detail, double* baseOrig, double* baseTM, int width, int height, const NoiseModel& noise, double detailScale);

    void computeHistogram(double* luminance, int width, int height, std::vector<double>& histogram, int numBins);
    void computeWeightedHistogram(double* luminance, int width, int height, const NoiseModel& noise, std::vector<double>& histogram, int numBins);

    double computeBlockMean(int bx, int by, int blockSize);
    double computeBlockVariance(int bx, int by, int blockSize);
    bool isUniform(int bx, int by, int blockSize);
    void fitLinearModel(const std::vector<std::pair<double, double>>& samples, double& a, double& b);
    
    private:
        static const int NUM_TONE_CURVE_NODES = 30;
        static const double TONE_CURVE_BIN_WIDTH;
        static const int TILE_SIZE = 230;
        
        NoiseModel estimateNoise();

};

#endif // TMO_EILERTSEN15_H