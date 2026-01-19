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
*                     Implementation of the TMOEilertsen15 class               *
*                         Real-time noise-aware tone mapping                   *
*                     https://computergraphics.on.liu.se/rntm/                 *
*                                                                              *
*******************************************************************************/

#include "TMOEilertsen15.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <numeric>

const int TMOEilertsen15::NUM_TONE_CURVE_NODES;
constexpr double TMOEilertsen15::TONE_CURVE_BIN_WIDTH;
const int TMOEilertsen15::DEFAULT_TILE_SIZE;
constexpr double TMOEilertsen15::EDGE_STOP_LAMBDA;
constexpr double TMOEilertsen15::TEMPORAL_CUTOFF_HZ;


TMOEilertsen15::TMOEilertsen15() {
   SetName(L"Eilertsen15");
   SetDescription(L"Real-time noise-aware tone mapping (Eilertsen et al. 2015)");

   dPeakLuminance.SetName(L"Peak Luminance");
   dPeakLuminance.SetDescription(L"Lmax: Peak display luminance [cd/m2]");
   dPeakLuminance.SetDefault(200.0);
   dPeakLuminance.SetRange(50.0, 10000.0);
   this->Register(dPeakLuminance);
   dPeakLuminance = 200.;

   dBlackLevel.SetName(L"Black Level");
   dBlackLevel.SetDescription(L"Lblack: Display black level [cd/m2]");
   dBlackLevel.SetDefault(0.5);
   dBlackLevel.SetRange(0.01, 10.0);
   this->Register(dBlackLevel);
   dBlackLevel = .5;

   dGamma.SetName(L"Gamma");
   dGamma.SetDescription(L"y: Display gamma");
   dGamma.SetDefault(2.2);
   dGamma.SetRange(1.0, 3.0);
   this->Register(dGamma);
   dGamma = 2.2;

   dAmbientLight.SetName(L"Ambient Light");
   dAmbientLight.SetDescription(L"Eamb: Ambient illuminance [lux]");
   dAmbientLight.SetDefault(100.0);
   dAmbientLight.SetRange(0.0, 5000.0);
   this->Register(dAmbientLight);
   dAmbientLight = 100.;

   dReflectivity.SetName(L"Reflectivity");
   dReflectivity.SetDescription(L"k: Display reflectivity [%]");
   dReflectivity.SetDefault(0.8);
   dReflectivity.SetRange(0.1, 5.0);
   this->Register(dReflectivity);
   dReflectivity = .8;

   dDetailScaling.SetName(L"Detail Scaling");
   dDetailScaling.SetDescription(L"e: Detail enhancement factor");
   dDetailScaling.SetDefault(1.0);
   dDetailScaling.SetRange(0.0, 4.0);
   this->Register(dDetailScaling);
   dDetailScaling = 1.;

   dNoiseControl.SetName(L"Noise Control");
   dNoiseControl.SetDescription(L"Additional noise visibility control");
   dNoiseControl.SetDefault(1.0);
   dNoiseControl.SetRange(0.0, 2.0);
   this->Register(dNoiseControl);
   dNoiseControl = 1.;

   dTonePriority.SetName(L"Tone Priority");
   dTonePriority.SetDescription(L"Priority: -1=bright, 0=neutral, 1=dark");
   dTonePriority.SetDefault(0.0);
   dTonePriority.SetRange(-1.0, 1.0);
   this->Register(dTonePriority);
   dTonePriority = 0.;

   bLocalToneCurves.SetName(L"Local Tone Curves");
   bLocalToneCurves.SetDescription(L"Use local tone curves");
   bLocalToneCurves.SetDefault(true);
   this->Register(bLocalToneCurves);
   bLocalToneCurves = true;

   iFilterIterations.SetName(L"Filter Iterations");
   iFilterIterations.SetDescription(L"N: Number of diffusion iterations");
   iFilterIterations.SetDefault(12);
   iFilterIterations.SetRange(1, 20);
   this->Register(iFilterIterations);
   iFilterIterations = 12;

   dSigma.SetName(L"Sigma");
   dSigma.SetDescription(L"o: Starting kernel size for diffusion");
   dSigma.SetDefault(3.0);
   dSigma.SetRange(1.0, 10.0);
   this->Register(dSigma);
   dSigma = 3.;

   dNoiseA.SetName(L"Noise A");
   dNoiseA.SetDescription(L"a: Signal-dependent noise parameter");
   dNoiseA.SetDefault(0.001);
   dNoiseA.SetRange(0.0, 1.0);
   this->Register(dNoiseA);
   dNoiseA = .001;

   dNoiseB.SetName(L"Noise B");
   dNoiseB.SetDescription(L"b: Signal-independent noise parameter");
   dNoiseB.SetDefault(0.0001);
   dNoiseB.SetRange(0.0, 0.1);
   this->Register(dNoiseB);
   dNoiseB = .0001;

   bEstimateNoise.SetName(L"Estimate Noise");
   bEstimateNoise.SetDescription(L"Automatically estimate noise (Foi et al. 2008)");
   bEstimateNoise.SetDefault(true);
   this->Register(bEstimateNoise);
   bEstimateNoise = true;
}

TMOEilertsen15::~TMOEilertsen15() {
}

int TMOEilertsen15::Transform() {
   // Convert to Yxy color space for luminance processing
   pSrc->Convert(TMO_Yxy);
   pDst->Convert(TMO_Yxy, true);
   
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int size = width * height;
   
   double* luminance = new double[size];
   double* logLuminance = new double[size];
   double* baseLayer = new double[size];
   double* detailLayer = new double[size];
   double* contrast = new double[size];
   
   // Extract luminance and convert to log domain
   const double* pSourceData = pSrc->GetData();
   for (int i = 0; i < size; i++) {
      luminance[i] = std::max(pSourceData[i * 3], 1e-6);
      logLuminance[i] = std::log10(luminance[i]);
   }
   
   // Get or estimate noise model
   NoiseModel noise;
   if (bEstimateNoise) {
      noise = estimateNoise();
   } else {
      noise.a = dNoiseA;
      noise.b = dNoiseB;
   }
   
   // Step 1: Detail extraction using fast diffusion
   fastDetailExtractionDiffusion(logLuminance, baseLayer, width, height);
   
   // Compute detail layer
   for (int i = 0; i < size; i++) {
      detailLayer[i] = logLuminance[i] - baseLayer[i];
   }
   
   // Compute local contrast for histogram weighting
   computeLocalContrast(logLuminance, width, height, contrast);
   
   // Step 2: Compute and apply tone curves
   if (bLocalToneCurves) {
      std::vector<std::vector<double>> localSlopes;
      int tilesX, tilesY;
      double minLog, maxLog;
      
      computeLocalToneCurves(logLuminance, contrast, width, height, noise, localSlopes, tilesX, tilesY, minLog, maxLog);
      applyLocalToneCurves(baseLayer, width, height, localSlopes, tilesX, tilesY, minLog, maxLog);
   } else {
      double minLog = baseLayer[0], maxLog = baseLayer[0];
      for (int i = 1; i < size; i++) {
         minLog = std::min(minLog, baseLayer[i]);
         maxLog = std::max(maxLog, baseLayer[i]);
      }
      
      std::vector<double> histogram(NUM_TONE_CURVE_NODES, 0.0);
      computeNoiseAwareHistogram(logLuminance, contrast, width, height, noise, histogram, minLog, maxLog);
      
      std::vector<double> slopes(NUM_TONE_CURVE_NODES);
      computeToneCurveSlopes(histogram, noise, slopes);
      applyToneCurve(baseLayer, width, height, slopes, minLog, maxLog);
   }
   
   // Step 3: Noise-aware detail control
   noiseAwareDetailControl(detailLayer, logLuminance, baseLayer, width, height, noise);
   
   // Combine base and detail layers
   for (int i = 0; i < size; i++) {
      double toneMappedLog = baseLayer[i] + detailLayer[i];
      double toneMappedLin = std::pow(10.0, toneMappedLog);
      double pixelValue = std::max(0.0, std::min(1.0, toneMappedLin));
      
      pDst->GetData()[i * 3 + 0] = pixelValue;  // Y
      pDst->GetData()[i * 3 + 1] = pSourceData[i * 3 + 1];  // x
      pDst->GetData()[i * 3 + 2] = pSourceData[i * 3 + 2];  // y
   }
   
   delete[] luminance;
   delete[] logLuminance;
   delete[] baseLayer;
   delete[] detailLayer;
   delete[] contrast;
   
   pDst->Convert(TMO_RGB);
   return 0;
}

// Equation 1
double TMOEilertsen15::computeNoiseVariance(double luminance, const NoiseModel& noise) const{
   return noise.a * luminance + noise.b;
}

// Equation 2
double TMOEilertsen15::computeNoiseLevel(double luminance, const NoiseModel& noise) {
   double variance = computeNoiseVariance(luminance, noise);
   double sigma = std::sqrt(variance);
   return std::log10((luminance + sigma) / luminance);
}

// Foi et al. 2008 noise estimation
TMOEilertsen15::NoiseModel TMOEilertsen15::estimateNoise() {
   NoiseModel noise;
   std::vector<std::pair<double, double>> samples;
   
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   const double* pData = pSrc->GetData();
   
   // Extract luminance
   int size = width * height;
   double* lum = new double[size];
   for (int i = 0; i < size; i++) {
      double r = pData[i * 3];
      double g = pData[i * 3 + 1];
      double b = pData[i * 3 + 2];
      lum[i] = 0.299 * r + 0.587 * g + 0.114 * b;
   }
   
   // Analyze blocks
   int blockSize = 32;
   for (int by = 0; by < height; by += blockSize) {
      for (int bx = 0; bx < width; bx += blockSize) {
         if (isUniformBlock(lum, width, height, bx, by, blockSize)) {
               double mean = computeBlockMean(lum, width, height, bx, by, blockSize);
               double variance = computeBlockVariance(lum, width, height, bx, by, blockSize, mean);
               if (mean > 0 && variance > 0) {
                  samples.push_back({mean, variance});
               }
         }
      }
   }
   
   delete[] lum;
   
   // Fit linear model: variance = a * mean + b
   if (samples.size() >= 10) {
      fitLinearModel(samples, noise.a, noise.b);
   } else {
      noise.a = dNoiseA;
      noise.b = dNoiseB;
   }
   
   return noise;
}

void TMOEilertsen15::fitLinearModel(const std::vector<std::pair<double, double>>& samples, double& a, double& b) {
   int n = samples.size();
   double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
   
   for (const auto& sample : samples) {
      double x = sample.first;
      double y = sample.second;
      sumX += x;
      sumY += y;
      sumXY += x * y;
      sumX2 += x * x;
   }
   
   // Linear regression: y = ax + b
   double denominator = n * sumX2 - sumX * sumX;
   if (std::abs(denominator) > 1e-10) {
      a = (n * sumXY - sumX * sumY) / denominator;
      b = (sumY - a * sumX) / n;
      
      a = std::max(0.0, a);
      b = std::max(0.0, b);
   } else {
      a = 0.001;
      b = 0.0001;
   }
}

double TMOEilertsen15::computeBlockMean(const double* luminance, int width, int height, int bx, int by, int blockSize) {
   double sum = 0.0;
   int count = 0;
   
   int endX = std::min(bx + blockSize, width);
   int endY = std::min(by + blockSize, height);
   
   for (int y = by; y < endY; y++) {
      for (int x = bx; x < endX; x++) {
         sum += luminance[y * width + x];
         count++;
      }
   }
   
   return (count > 0) ? sum / count : 0.0;
}

double TMOEilertsen15::computeBlockVariance(const double* luminance, int width, int height, int bx, int by, int blockSize, double mean) {
   double sumSquaredDiff = 0.0;
   int count = 0;
   
   int endX = std::min(bx + blockSize, width);
   int endY = std::min(by + blockSize, height);
   
   for (int y = by; y < endY; y++) {
      for (int x = bx; x < endX; x++) {
         double diff = luminance[y * width + x] - mean;
         sumSquaredDiff += diff * diff;
         count++;
      }
   }
   
   return (count > 1) ? sumSquaredDiff / (count - 1) : 0.0;
}

bool TMOEilertsen15::isUniformBlock(const double* luminance, int width, int height, int bx, int by, int blockSize) {
   double mean = computeBlockMean(luminance, width, height, bx, by, blockSize);
   double variance = computeBlockVariance(luminance, width, height, bx, by, blockSize, mean);
   
   // Block is uniform if coefficient of variation is low
   double stdDev = std::sqrt(variance);
   double cv = (mean > 0) ? stdDev / mean : 0;
   
   return cv < 0.1;  // Threshold for uniformity
}

// HDR-VDP-2 contrast sensitivity function - Mantiuk et al. 2011
double TMOEilertsen15::computeContrastSensitivity(double luminance, double frequency) {
   luminance = std::max(luminance, 0.01);
   
   double A = 75.0;
   double B = 0.2;
   
   double S_max = A * std::pow(luminance / 100.0, B);
   
   double f_max = 4.0;
   double f_norm = frequency / f_max;
   double spatial_term = std::exp(-0.5 * std::pow(f_norm - 1.0, 2));
   
   double sensitivity = S_max * spatial_term;
   
   return std::max(sensitivity, 1.0);
}

// Equation 3: Detection threshold Ct(Ld) = 1/CSF(Ld)
double TMOEilertsen15::computeDetectionThreshold(double displayLuminance) {
   // Use peak CSF (frequency at peak around 4 cpd)
   double csf = computeContrastSensitivity(displayLuminance, 4.0);
   return 1.0 / csf;
}

// Equation 3: V(Ld)
double TMOEilertsen15::computeVisibilityThreshold(double displayLuminance) {
   double Ct = computeDetectionThreshold(displayLuminance);
   
   Ct = std::min(Ct, 0.99);
   
   return 0.5 * std::log10((Ct + 1.0) / (1.0 - Ct));
}

// Equation 4: Ld(L')
double TMOEilertsen15::displayModel(double pixelValue) {
   double gamma = dGamma;
   double Lmax = dPeakLuminance;
   double Lblack = dBlackLevel;
   double Lrefl = computeReflectedLight();
   
   pixelValue = std::max(0.0, std::min(1.0, pixelValue));
   
   return std::pow(pixelValue, gamma) * (Lmax - Lblack) + Lblack + Lrefl;
}

// Equation 5
double TMOEilertsen15::computeReflectedLight() {
   double k = dReflectivity / 100.0;
   double Eamb = dAmbientLight;
   return (k / M_PI) * Eamb;
}

// Inverse display model
/*
* Note: Not used in main tone mapping pipeline - the tone curve already
* produces pixel values in range [0,1]. This function is provided for
* completeness of the display model implementation
*/
// cppcheck-suppress unusedFunction
double TMOEilertsen15::inverseDisplayModel(double displayLuminance) {
   double gamma = dGamma;
   double Lmax = dPeakLuminance;
   double Lblack = dBlackLevel;
   double Lrefl = computeReflectedLight();
   
   double numerator = displayLuminance - Lblack - Lrefl;
   double denominator = Lmax - Lblack;
   
   if (denominator <= 0) return 0.0;
   
   double ratio = numerator / denominator;
   ratio = std::max(0.0, ratio);
   
   return std::pow(ratio, 1.0 / gamma);
}

// Equation 6
double TMOEilertsen15::computeDisplayDynamicRange() {
   double Ld1 = displayModel(1.0);
   double Ld0 = displayModel(0.0);
   
   if (Ld0 <= 0) Ld0 = 1e-6;
   
   double ratio = Ld1 / Ld0;
   double r = std::log10(ratio);

   return r;
}

// Equations 14-17
void TMOEilertsen15::computeToneCurveSlopes(const std::vector<double>& histogram, const NoiseModel& noise, std::vector<double>& slopes) {
   int N = NUM_TONE_CURVE_NODES;
   double r = computeDisplayDynamicRange();
   double delta = TONE_CURVE_BIN_WIDTH;
   
   slopes.resize(N);
   
   std::vector<double> weightedHist = histogram;
   applyTonePriority(weightedHist);
   
   // 15-17
   double pt = 0.0001;
   std::vector<bool> validSegments(N, false);
   
   // Iteratively find segments with positive slopes
   for (int iteration = 0; iteration < 10; iteration++) {
      int numValid = 0;
      double sumInvProb = 0.0;
      
      for (int k = 0; k < N; k++) {
         if (weightedHist[k] > pt) {
               validSegments[k] = true;
               numValid++;
               sumInvProb += 1.0 / weightedHist[k];
         } else {
               validSegments[k] = false;
         }
      }
      
      if (numValid == 0) break;
      
      // 17
      double newPt = (numValid - r / delta) / sumInvProb;
      
      // Check convergence
      if (std::abs(newPt - pt) < 1e-6) break;
      pt = std::max(0.0, newPt);
   }
   
   // 14
   double sumInvProb = 0.0;
   for (int k = 0; k < N; k++) {
      if (validSegments[k]) {
         sumInvProb += 1.0 / weightedHist[k];
      }
   }
   
   for (int k = 0; k < N; k++) {
      if (validSegments[k] && weightedHist[k] > 0) {
         slopes[k] = 1.0 + (r / delta - N) / (weightedHist[k] * sumInvProb);
         slopes[k] = std::max(0.0, slopes[k]);
      } else {
         slopes[k] = 0.0;
      }
   }
}

// Apply tone priority weighting
void TMOEilertsen15::applyTonePriority(std::vector<double>& histogram) {
   if (std::abs(dTonePriority) < 1e-6) return;
   
   int N = histogram.size();
   for (int k = 0; k < N; k++) {
      double position = (double)k / (N - 1);
      
      double weight = 1.0;
      if (dTonePriority > 0.0) {
         weight = 1.0 + dTonePriority * (1.0 - position);
      } else {
         weight = 1.0 - dTonePriority * position;
      }
      
      histogram[k] *= weight;
   }
   
   double sum = std::accumulate(histogram.begin(), histogram.end(), 0.0);
   if (sum > 0) {
      std::for_each(histogram.begin(), histogram.end(), [sum](double& h) { h /= sum; });
   }
}

// Build and apply piecewise linear tone curve
void TMOEilertsen15::applyToneCurve(double* baseLayer, int width, int height, const std::vector<double>& slopes, double minLog, double maxLog) {
   int N = NUM_TONE_CURVE_NODES;
   double delta = TONE_CURVE_BIN_WIDTH;
   double r = computeDisplayDynamicRange();
   int size = width * height;
   
   std::vector<double> nodeInputs(N + 1);
   std::vector<double> nodeOutputs(N + 1);
   
   nodeInputs[0] = minLog;
   nodeOutputs[0] = -r;
   
   for (int k = 0; k < N; k++) {
      nodeInputs[k + 1] = nodeInputs[k] + delta;
      nodeOutputs[k + 1] = nodeOutputs[k] + slopes[k] * delta;
   }
   
   // Apply tone curve via linear interpolation
   for (int i = 0; i < size; i++) {
      double val = baseLayer[i];
      
      int segment = (int)((val - minLog) / delta);
      segment = std::max(0, std::min(N - 1, segment));
      
      double t = (val - nodeInputs[segment]) / delta;
      t = std::max(0.0, std::min(1.0, t));
      
      baseLayer[i] = nodeOutputs[segment] + t * (nodeOutputs[segment + 1] - nodeOutputs[segment]);
   }
}

// Equation 18
void TMOEilertsen15::computeLocalContrast(const double* logLuminance, int width, int height, double* contrast) {
   int size = width * height;
   double sigma = 3.0;
   
   double* filtered = new double[size];
   double* filteredSquared = new double[size];
   
   double* logLumSquared = new double[size];
   for (int i = 0; i < size; i++) {
      logLumSquared[i] = logLuminance[i] * logLuminance[i];
   }
   
   gaussianBlur(logLuminance, filtered, width, height, sigma);
   gaussianBlur(logLumSquared, filteredSquared, width, height, sigma);
   
   for (int i = 0; i < size; i++) {
      double variance = filteredSquared[i] - filtered[i] * filtered[i];
      contrast[i] = std::sqrt(std::max(0.0, variance));
   }
   
   delete[] filtered;
   delete[] filteredSquared;
   delete[] logLumSquared;
}

// Equations 19-21
void TMOEilertsen15::computeNoiseAwareHistogram(const double* logLuminance, const double* contrast, int width, int height, const NoiseModel& noise, std::vector<double>& histogram, double minLog, double maxLog) {
   int size = width * height;
   int numBins = histogram.size();
   double binWidth = (maxLog - minLog) / numBins;
   
   std::fill(histogram.begin(), histogram.end(), 0.0);
   
   // 20-21
   double totalWeight = 0.0;
   
   for (int i = 0; i < size; i++) {
      double lumLin = std::pow(10.0, logLuminance[i]);
      
      double noiseLevel = computeNoiseLevel(lumLin, noise);
      
      // 20
      if (contrast[i] > noiseLevel) {
         // 21
         int bin = (int)((logLuminance[i] - minLog) / binWidth);
         bin = std::max(0, std::min(numBins - 1, bin));
         
         // 19
         double weight = contrast[i] * dNoiseControl;
         histogram[bin] += weight;
         totalWeight += weight;
      }
   }
   
   // 19
   if (totalWeight > 0) {
      for (int i = 0; i < numBins; i++) {
         histogram[i] /= totalWeight;
      }
   } else {
      for (int i = 0; i < numBins; i++) {
         histogram[i] = 1.0 / numBins;
      }
   }
}

// Apply IIR filter for flicker reduction
void TMOEilertsen15::applyTemporalFilter(std::vector<double>& nodeValues, IIRFilter& filter) {
   if (filter.history.empty()) {
      filter.alpha = 0.105;
      filter.history.push_back(nodeValues);
      return;
   }
   
   std::vector<double> filtered = nodeValues;
   const std::vector<double>& prev = filter.history.back();
   
   for (size_t i = 0; i < nodeValues.size(); i++) {
      filtered[i] = filter.alpha * nodeValues[i] + (1.0 - filter.alpha) * prev[i];
   }
   
   filter.history.push_back(filtered);
   if (filter.history.size() > filter.maxHistory) {
      filter.history.pop_front();
   }
   
   nodeValues = filtered;
}

// Compute local tone curves with 90% local + 10% global blending
void TMOEilertsen15::computeLocalToneCurves(const double* logLuminance, const double* contrast, int width, int height, const NoiseModel& noise, std::vector<std::vector<double>>& localSlopes, int& tilesX, int& tilesY, double& minLog, double& maxLog) {
   int size = width * height;
   
   minLog = logLuminance[0];
   maxLog = logLuminance[0];
   for (int i = 1; i < size; i++) {
      minLog = std::min(minLog, logLuminance[i]);
      maxLog = std::max(maxLog, logLuminance[i]);
   }
   
   int tileSize = DEFAULT_TILE_SIZE;
   tilesX = (width + tileSize - 1) / tileSize;
   tilesY = (height + tileSize - 1) / tileSize;

   if (temporalFilters.empty() || tilesX != lastTilesX || tilesY != lastTilesY) {
      temporalFilters.clear();
      temporalFilters.resize(tilesX * tilesY);
      lastTilesX = tilesX;
      lastTilesY = tilesY;
   }
   
   std::vector<double> globalHist(NUM_TONE_CURVE_NODES, 0.0);
   computeNoiseAwareHistogram(logLuminance, contrast, width, height, noise, globalHist, minLog, maxLog);
   
   localSlopes.resize(tilesY);
   for (int ty = 0; ty < tilesY; ty++) {
      localSlopes[ty].resize(tilesX * NUM_TONE_CURVE_NODES);
   }
   
   for (int ty = 0; ty < tilesY; ty++) {
      for (int tx = 0; tx < tilesX; tx++) {
         int x0 = tx * tileSize;
         int y0 = ty * tileSize;
         int x1 = std::min(x0 + tileSize, width);
         int y1 = std::min(y0 + tileSize, height);
         
         int tileWidth = x1 - x0;
         int tileHeight = y1 - y0;
         int tilePixels = tileWidth * tileHeight;
         
         double* tileLog = new double[tilePixels];
         double* tileContrast = new double[tilePixels];
         
         int idx = 0;
         for (int y = y0; y < y1; y++) {
               for (int x = x0; x < x1; x++) {
                  tileLog[idx] = logLuminance[y * width + x];
                  tileContrast[idx] = contrast[y * width + x];
                  idx++;
               }
         }
         
         std::vector<double> localHist(NUM_TONE_CURVE_NODES, 0.0);
         computeNoiseAwareHistogram(tileLog, tileContrast, tileWidth, tileHeight, noise, localHist, minLog, maxLog);
         
         // Blend with global
         for (int k = 0; k < NUM_TONE_CURVE_NODES; k++) {
               localHist[k] = 0.9 * localHist[k] + 0.1 * globalHist[k];
         }
         
         std::vector<double> slopes(NUM_TONE_CURVE_NODES);
         computeToneCurveSlopes(localHist, noise, slopes);
         
         int tileIdx = ty * tilesX + tx;
         applyTemporalFilter(slopes, temporalFilters[tileIdx]);
         
         for (int k = 0; k < NUM_TONE_CURVE_NODES; k++) {
               localSlopes[ty][tx * NUM_TONE_CURVE_NODES + k] = slopes[k];
         }
         
         delete[] tileLog;
         delete[] tileContrast;
      }
   }
}

// Apply local tone curves with bilinear interpolation
void TMOEilertsen15::applyLocalToneCurves(double* baseLayer, int width, int height, const std::vector<std::vector<double>>& localSlopes, int tilesX, int tilesY, double minLog, double maxLog) {
   int N = NUM_TONE_CURVE_NODES;
   double delta = TONE_CURVE_BIN_WIDTH;
   double r = computeDisplayDynamicRange();
   int tileSize = DEFAULT_TILE_SIZE;
   
   // For each pixel, interpolate tone curve from surrounding tiles
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double val = baseLayer[y * width + x];
         
         double tileX = (double)x / tileSize;
         double tileY = (double)y / tileSize;
         
         int tx0 = (int)std::floor(tileX);
         int ty0 = (int)std::floor(tileY);
         int tx1 = std::min(tx0 + 1, tilesX - 1);
         int ty1 = std::min(ty0 + 1, tilesY - 1);
         
         double wx = tileX - tx0;
         double wy = tileY - ty0;
         
         // Find segment in tone curve
         int segment = (int)((val - minLog) / delta);
         segment = std::max(0, std::min(N - 1, segment));
         double t = (val - (minLog + segment * delta)) / delta;
         t = std::max(0.0, std::min(1.0, t));
         
         auto getSlope = [&](int ty, int tx, int seg) -> double {
               return localSlopes[ty][tx * N + seg];
         };
         
         double slope00 = getSlope(ty0, tx0, segment);
         double slope01 = getSlope(ty0, tx1, segment);
         double slope10 = getSlope(ty1, tx0, segment);
         double slope11 = getSlope(ty1, tx1, segment);
         
         double slope0 = slope00 * (1 - wx) + slope01 * wx;
         double slope1 = slope10 * (1 - wx) + slope11 * wx;
         double slope = slope0 * (1 - wy) + slope1 * wy;
         
         double output = -r;
         
         for (int s = 0; s < segment; s++) {
               double s00 = getSlope(ty0, tx0, s);
               double s01 = getSlope(ty0, tx1, s);
               double s10 = getSlope(ty1, tx0, s);
               double s11 = getSlope(ty1, tx1, s);
               
               double s0 = s00 * (1 - wx) + s01 * wx;
               double s1 = s10 * (1 - wx) + s11 * wx;
               double sInterp = s0 * (1 - wy) + s1 * wy;
               
               output += sInterp * delta;
         }
         
         output += slope * delta * t;
         
         baseLayer[y * width + x] = output;
      }
   }
}

// Fast detail extraction diffusion
void TMOEilertsen15::fastDetailExtractionDiffusion(const double* input, double* output, int width, int height) {
   int size = width * height;
   int N = iFilterIterations;
   double sigma = dSigma;
   double lambda = EDGE_STOP_LAMBDA;
   double* temp = new double[size];
   double* gradMag = new double[size];
   
   std::memcpy(output, input, size * sizeof(double));

   for (int iter = 0; iter < N; iter++) {
      double sigma_k = sigma * std::sqrt((double)(iter + 1));
      
      gaussianBlur(output, temp, width, height, sigma_k);
      
      int radius = (int)std::ceil(3.0 * sigma_k);
      computeGradientMagnitude(output, gradMag, width, height, radius);
      
      for (int i = 0; i < size; i++) {
         double dist = std::abs(temp[i] - input[i]);
         gradMag[i] = std::max(gradMag[i], (double)(iter + 1) * dist);
      }
      
      std::vector<double> wr(size);
      for (int i = 0; i < size; i++) {
         wr[i] = tukeyBiweight(gradMag[i], lambda);
      }
      
      for (int i = 0; i < size; i++) {
         output[i] = (1.0 - wr[i]) * output[i] + wr[i] * temp[i];
      }
   }
   
   delete[] temp;
   delete[] gradMag;
}

// Separable Gaussian
void TMOEilertsen15::gaussianBlur(const double* input, double* output, int width, int height, double sigma) {
   int radius = (int)std::ceil(3.0 * sigma);
   int kernelSize = 2 * radius + 1;
   std::vector<double> kernel(kernelSize);
   
   double sum = 0.0;
   for (int i = 0; i < kernelSize; i++) {
      int x = i - radius;
      kernel[i] = std::exp(-0.5 * (x * x) / (sigma * sigma));
      sum += kernel[i];
   }
   
   for (int i = 0; i < kernelSize; i++) {
      kernel[i] /= sum;
   }
   
   double* temp = new double[width * height];
   
   // Horizontal pass
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double pixelSum = 0.0;
         double weightSum = 0.0;
         
         for (int k = -radius; k <= radius; k++) {
               int xx = x + k;
               if (xx >= 0 && xx < width) {
                  pixelSum += input[y * width + xx] * kernel[k + radius];
                  weightSum += kernel[k + radius];
               }
         }
         
         temp[y * width + x] = pixelSum / weightSum;
      }
   }
   
   // Vertical pass
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double pixelSum = 0.0;
         double weightSum = 0.0;
         
         for (int k = -radius; k <= radius; k++) {
               int yy = y + k;
               if (yy >= 0 && yy < height) {
                  pixelSum += temp[yy * width + x] * kernel[k + radius];
                  weightSum += kernel[k + radius];
               }
         }
         
         output[y * width + x] = pixelSum / weightSum;
      }
   }
   
   delete[] temp;
}

// Gradient magnitude with linear ramp
void TMOEilertsen15::computeGradientMagnitude(const double* input, double* gradMag, int width, int height, int radius) {
   // 27
   
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double gradX = 0.0;
         double weightSumX = 0.0;
         
         for (int delta = -radius; delta <= radius; delta++) {
               int xx = x + delta;
               if (xx >= 0 && xx < width) {
                  gradX += delta * input[y * width + xx];
                  weightSumX += std::abs(delta);
               }
         }
         
         if (weightSumX > 0) {
               gradX /= weightSumX;
         }
         
         double gradY = 0.0;
         double weightSumY = 0.0;
         
         for (int delta = -radius; delta <= radius; delta++) {
               int yy = y + delta;
               if (yy >= 0 && yy < height) {
                  gradY += delta * input[yy * width + x];
                  weightSumY += std::abs(delta);
               }
         }
         
         if (weightSumY > 0) {
               gradY /= weightSumY;
         }
         
         gradMag[y * width + x] = std::sqrt(gradX * gradX + gradY * gradY);
      }
   }
}

// Tukey's biweight edge-stop function
double TMOEilertsen15::tukeyBiweight(double x, double lambda) {
   if (x <= lambda) {
      double ratio = x / lambda;
      double term = 1.0 - ratio * ratio;
      return term * term;
   } else {
      return 0.0;
   }
}

// Equation 28
void TMOEilertsen15::noiseAwareDetailControl(double* detailLayer, const double* baseOrig, const double* baseTM, int width, int height, const NoiseModel& noise) {
   int size = width * height;
   double e = dDetailScaling;
   
   for (int i = 0; i < size; i++) {
      double lumOrig = std::pow(10.0, baseOrig[i]);
      double lumTM = std::pow(10.0, baseTM[i]);
      
      double noiseLevelOrig = computeNoiseLevel(lumOrig, noise);

      double displayLum = lumTM;
      double visibilityThreshold = computeVisibilityThreshold(displayLum);
      
      double scaleFactor = 1.0;
      if (noiseLevelOrig > 0) {
         scaleFactor = std::min(1.0, visibilityThreshold / noiseLevelOrig);
      }
      
      detailLayer[i] *= e * scaleFactor;
   }
}
