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

TMOEilertsen15::TMOEilertsen15(){
   
   SetName(L"Eilertsen15");
   SetDescription(L"Real-time noise-aware tone mapping");

   // Display parameters
   dPeakLuminance.SetName(L"Peak Luminance");
   dPeakLuminance.SetDescription(L"Peak display luminance in cd/m2");
   dPeakLuminance.SetDefault(200.0);
   dPeakLuminance.SetRange(50.0, 10000.0);
   this->Register(dPeakLuminance);

   dBlackLevel.SetName(L"Black Level");
   dBlackLevel.SetDescription(L"Display black level in cd/m2");
   dBlackLevel.SetDefault(0.5);
   dBlackLevel.SetRange(0.01, 10.0);
   this->Register(dBlackLevel);

   dAmbientLight.SetName(L"Ambient Light");
   dAmbientLight.SetDescription(L"Ambient illuminance in lux");
   dAmbientLight.SetDefault(100.0);
   dAmbientLight.SetRange(0.0, 5000.0);
   this->Register(dAmbientLight);

   dReflectivity.SetName(L"Reflectivity");
   dReflectivity.SetDescription(L"Display reflectivity in %");
   dReflectivity.SetDefault(0.8);
   dReflectivity.SetRange(0.1, 5.0);
   this->Register(dReflectivity);

   // Processing parameters
   dDetailScaling.SetName(L"Detail Scaling");
   dDetailScaling.SetDescription(L"Detail enhancement factor");
   dDetailScaling.SetDefault(1.0);
   dDetailScaling.SetRange(0.0, 4.0);
   this->Register(dDetailScaling);

   dNoiseControl.SetName(L"Noise Control");
   dNoiseControl.SetDescription(L"Noise visibility control");
   dNoiseControl.SetDefault(1.0);
   dNoiseControl.SetRange(0.0, 2.0);
   this->Register(dNoiseControl);

   dTonePriority.SetName(L"Tone Priority");
   dTonePriority.SetDescription(L"Priority: -1=bright, 0=neutral, 1=dark");
   dTonePriority.SetDefault(0.0);
   dTonePriority.SetRange(-1.0, 1.0);
   this->Register(dTonePriority);

   bLocalToneCurves.SetName(L"Local Tone Curves");
   bLocalToneCurves.SetDescription(L"Use local tone curves");
   bLocalToneCurves.SetDefault(true);
   this->Register(bLocalToneCurves);

   iFilterIterations.SetName(L"Filter Iterations");
   iFilterIterations.SetDescription(L"Number of diffusion iterations");
   iFilterIterations.SetDefault(12);
   iFilterIterations.SetRange(1, 20);
   this->Register(iFilterIterations);
}

TMOEilertsen15::~TMOEilertsen15(){
}

int TMOEilertsen15::Transform(){
   // Convert to luminance
   pSrc->Convert(TMO_Yxy);
   pDst->Convert(TMO_Yxy);
   
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int size = width * height;
   
   double* luminance = new double[size];
   double* logLuminance = new double[size];
   double* baseLayer = new double[size];
   double* detailLayer = new double[size];
   
   // Extract luminance and convert to log domain
   double* pSourceData = pSrc->GetData();
   for (int i = 0; i < size; i++) {
      luminance[i] = pSourceData[i * 3];
      logLuminance[i] = std::log10(std::max(luminance[i], 1e-6));
   }
   
   NoiseModel noise = estimateNoise();
   
   // Step 1: Detail extraction using fast diffusion
   detailExtractionDiffusion(logLuminance, baseLayer, width, height, iFilterIterations, 3.0);
   
   // Compute detail layer
   for (int i = 0; i < size; i++) {
      detailLayer[i] = logLuminance[i] - baseLayer[i];
   }
   
   // Step 2: Compute and apply tone curves
   if (bLocalToneCurves) {
      std::vector<std::vector<double>> localCurves;
      computeLocalToneCurves(baseLayer, width, height, localCurves);
      applyLocalToneCurves(baseLayer, width, height, localCurves);
   } else {
      std::vector<double> histogram(NUM_TONE_CURVE_NODES);
      computeWeightedHistogram(baseLayer, width, height, noise, 
                              histogram, NUM_TONE_CURVE_NODES);
      std::vector<double> slopes(NUM_TONE_CURVE_NODES);
      computeToneCurve(histogram, noise, slopes);
      applyToneCurve(baseLayer, width, height, slopes);
   }
   
   // Step 3: Noise-aware detail control
   noiseAwareDetailControl(detailLayer, logLuminance, baseLayer, width, height, noise, dDetailScaling);
   
   // Combine base and detail layers
   for (int i = 0; i < size; i++) {
      double toneMappedLog = baseLayer[i] + detailLayer[i];
      double toneMappedLum = std::pow(10.0, toneMappedLog);
      
      // Convert back to linear and apply to output
      double scale = toneMappedLum / std::max(luminance[i], 1e-6);
      
      double* pDestData = pDst->GetData();
      pDestData[i * 3] = toneMappedLum;
      pDestData[i * 3 + 1] = pSourceData[i * 3 + 1];  // x chromaticity
      pDestData[i * 3 + 2] = pSourceData[i * 3 + 2];  // y chromaticity
   }
   
   delete[] luminance;
   delete[] logLuminance;
   delete[] baseLayer;
   delete[] detailLayer;
   
   pDst->Convert(TMO_RGB);
   return 0;
}

// Display model
double TMOEilertsen15::displayModel(double pixelValue){
   double gamma = 2.2;
   double Lrefl = (dReflectivity / 100.0) * dAmbientLight / M_PI;
   return std::pow(pixelValue, gamma) * (dPeakLuminance - dBlackLevel) + dBlackLevel + Lrefl;
}

double TMOEilertsen15::inverseDisplayModel(double luminance){
   double gamma = 2.2;
   double Lrefl = (dReflectivity / 100.0) * dAmbientLight / M_PI;
   double val = (luminance - dBlackLevel - Lrefl) / (dPeakLuminance - dBlackLevel);
   return std::pow(std::max(val, 0.0), 1.0 / gamma);
}

// Effective display dynamic range
double TMOEilertsen15::computeDisplayDynamicRange(){
   double Ld_max = displayModel(1.0);
   double Ld_min = displayModel(0.0);
   return std::log10(Ld_max / Ld_min);
}

// Noise model
double TMOEilertsen15::computeNoiseLevel(double luminance, const NoiseModel& noise){
   double variance = noise.a * luminance + noise.b;
   double sigma = std::sqrt(std::max(variance, 0.0));
   return std::log10((luminance + sigma) / std::max(luminance, 1e-6));
}

// Visibility threshold
double computeVisibilityThreshold(double displayLuminance, double frequency) //TODO: full HDR-VDP-2 CSF - separate file?
{
   double L = std::max(displayLuminance, 0.01);
   
   double log_La = std::log10(L);
   
   // CSF parameters
   double A = 75.0;
   double B = 0.2;
   double f_max = 4.0;
   
   double S_max = A * std::pow(L / 100.0, B);
   
   // Spatial frequency component
   double f_norm = frequency / f_max;
   double spatial_term = std::exp(-std::pow(f_norm - 1.0, 2) / 2.0);
   
   double sensitivity = S_max * spatial_term;
   double contrastThreshold = 1.0 / sensitivity;
   
   return 0.5 * std::log10((contrastThreshold + 1.0) / (1.0 - contrastThreshold));
}

// Tone curve slopes
void TMOEilertsen15::computeToneCurve(const std::vector<double>& histogram, const NoiseModel& noise, std::vector<double>& slopes){
   int N = NUM_TONE_CURVE_NODES;
   double r = computeDisplayDynamicRange();
   double delta = TONE_CURVE_BIN_WIDTH;
   
   // Sum of inverse probabilities
   double sumInvProb = 0.0;
   std::vector<double> validHist;
   double threshold = 0.0001;
   
   for (int k = 0; k < N; k++) {
      if (histogram[k] > threshold) {
         sumInvProb += 1.0 / histogram[k];
         validHist.push_back(histogram[k]);
      }
   }
   
   // Slopes
   for (int k = 0; k < N; k++) {
      if (histogram[k] > threshold) {
         slopes[k] = 1.0 + (r / delta - N) / (histogram[k] * sumInvProb);
         slopes[k] = std::max(0.0, slopes[k]);
      } else {
         slopes[k] = 0.0;
      }
   }
}

// Apply tone curve to base layer
void TMOEilertsen15::applyToneCurve(double* base, int width, int height, const std::vector<double>& slopes){
   double delta = TONE_CURVE_BIN_WIDTH;
   double r = computeDisplayDynamicRange();
   int size = width * height;
   
   // Find min/max of base layer
   double minLog = base[0], maxLog = base[0];
   for (int i = 1; i < size; i++) {
      minLog = std::min(minLog, base[i]);
      maxLog = std::max(maxLog, base[i]);
   }
   
   // Build piecewise linear tone curve
   std::vector<double> nodePositions(NUM_TONE_CURVE_NODES + 1);
   std::vector<double> nodeValues(NUM_TONE_CURVE_NODES + 1);
   
   nodePositions[0] = minLog;
   nodeValues[0] = -r;
   
   for (int k = 0; k < NUM_TONE_CURVE_NODES; k++) {
      nodePositions[k + 1] = nodePositions[k] + delta;
      nodeValues[k + 1] = nodeValues[k] + slopes[k] * delta;
   }
   
   // Apply tone curve
   for (int i = 0; i < size; i++) {
      double val = base[i];
      
      int segment = (int)((val - minLog) / delta);
      segment = std::max(0, std::min(NUM_TONE_CURVE_NODES - 1, segment));
      
      double t = (val - nodePositions[segment]) / delta;
      base[i] = nodeValues[segment] + t * (nodeValues[segment + 1] - nodeValues[segment]);
   }
}

// Fast detail extraction diffusion
void TMOEilertsen15::detailExtractionDiffusion(double* input, double* output, int width, int height, int iterations, double sigma){
   int size = width * height;
   double* temp = new double[size];
   double* gradMag = new double[size];
   
   std::memcpy(output, input, size * sizeof(double));
   
   // Iterative diffusion
   for (int iter = 0; iter < iterations; iter++) {
      double sigma_k = sigma * std::sqrt((double)(iter + 1) / iterations);
      
      gaussianBlur(output, temp, width, height, sigma_k);
      
      int radius = (int)(3 * sigma_k);
      computeGradientMagnitude(output, gradMag, width, height, radius);
      
      for (int i = 0; i < size; i++) {
         double dist = std::abs(temp[i] - input[i]);
         gradMag[i] = std::max(gradMag[i], (double)iter * dist);
      }
      
      // Edge-stop function (Tukey's biweight)
      double lambda = 0.5;
      for (int i = 0; i < size; i++) {
         double x = gradMag[i] / lambda;
         double wr;
         if (x <= 1.0) {
               wr = (1.0 - x * x) * (1.0 - x * x);
         } else {
               wr = 0.0;
         }
         
         output[i] = (1.0 - wr) * output[i] + wr * temp[i];
      }
   }
   
   delete[] temp;
   delete[] gradMag;
}

// Gaussian blur (separable)
void TMOEilertsen15::gaussianBlur(double* input, double* output, int width, int height, double sigma){ //TODO: rewise
   int radius = (int)(3 * sigma);
   int kernelSize = 2 * radius + 1;
   
   // Build Gaussian kernel
   std::vector<double> kernel(kernelSize);
   double sum = 0.0;
   for (int i = 0; i < kernelSize; i++) {
      int x = i - radius;
      kernel[i] = std::exp(-(x * x) / (2.0 * sigma * sigma));
      sum += kernel[i];
   }
   for (int i = 0; i < kernelSize; i++) {
      kernel[i] /= sum;
   }
   
   // Temporary buffer
   double* temp = new double[width * height];
   
   // Horizontal pass
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double val = 0.0;
         for (int k = -radius; k <= radius; k++) {
               int xx = std::max(0, std::min(width - 1, x + k));
               val += input[y * width + xx] * kernel[k + radius];
         }
         temp[y * width + x] = val;
      }
   }
   
   // Vertical pass
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double val = 0.0;
         for (int k = -radius; k <= radius; k++) {
               int yy = std::max(0, std::min(height - 1, y + k));
               val += temp[yy * width + x] * kernel[k + radius];
         }
         output[y * width + x] = val;
      }
   }
   
   delete[] temp;
}

// Gradient magnitude
void TMOEilertsen15::computeGradientMagnitude(double* input, double* gradMag, int width, int height, int radius){
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double gradX = 0.0, gradY = 0.0;
         
         for (int dx = -radius; dx <= radius; dx++) {
               int xx = std::max(0, std::min(width - 1, x + dx));
               gradX += dx * input[y * width + xx];
         }
         
         for (int dy = -radius; dy <= radius; dy++) {
               int yy = std::max(0, std::min(height - 1, y + dy));
               gradY += dy * input[yy * width + x];
         }
         
         gradMag[y * width + x] = std::sqrt(gradX * gradX + gradY * gradY);
      }
   }
}

// Noise-aware detail control
void TMOEilertsen15::noiseAwareDetailControl(double* detail, double* baseOrig, double* baseTM, int width, int height, const NoiseModel& noise, double detailScale){
   for (int i = 0; i < width * height; i++) {
      double lumOrig = std::pow(10.0, baseOrig[i]);
      double noiseLevel = computeNoiseLevel(lumOrig, noise);
      
      double lumTM = std::pow(10.0, baseTM[i]);
      double lumDisplay = displayModel(inverseDisplayModel(lumTM));
      double visThreshold = computeVisibilityThreshold(lumDisplay);
      
      double ratio = std::min(1.0, visThreshold / (noiseLevel * dNoiseControl));
      detail[i] *= detailScale * ratio;
   }
}

// Weighted histogram
void TMOEilertsen15::computeWeightedHistogram(double* luminance, int width, int height, const NoiseModel& noise, std::vector<double>& histogram, int numBins){
   std::fill(histogram.begin(), histogram.end(), 0.0);
   
   double minLog = luminance[0], maxLog = luminance[0];
   for (int i = 1; i < width * height; i++) {
      minLog = std::min(minLog, luminance[i]);
      maxLog = std::max(maxLog, luminance[i]);
   }
   
   double range = maxLog - minLog;
   double binWidth = range / numBins;
   
   // Local contrast
   std::vector<double> contrast(width * height);
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double sum = 0.0, sumSq = 0.0;
         int count = 0;
         
         for (int dy = -3; dy <= 3; dy++) {
               for (int dx = -3; dx <= 3; dx++) {
                  int xx = x + dx, yy = y + dy;
                  if (xx >= 0 && xx < width && yy >= 0 && yy < height) {
                     double val = luminance[yy * width + xx];
                     sum += val;
                     sumSq += val * val;
                     count++;
                  }
               }
         }
         
         double mean = sum / count;
         double variance = sumSq / count - mean * mean;
         contrast[y * width + x] = std::sqrt(std::max(0.0, variance));
      }
   }
   
   // Build weighted histogram
   double totalWeight = 0.0;
   for (int i = 0; i < width * height; i++) {
      double lumLin = std::pow(10.0, luminance[i]);
      double noiseLevel = computeNoiseLevel(lumLin, noise);
      
      if (contrast[i] > noiseLevel) {
         int bin = (int)((luminance[i] - minLog) / binWidth);
         bin = std::max(0, std::min(numBins - 1, bin));
         
         double weight = contrast[i];
         histogram[bin] += weight;
         totalWeight += weight;
      }
   }
   
   // Normalize
   if (totalWeight > 0) {
      for (int i = 0; i < numBins; i++) {
         histogram[i] /= totalWeight;
      }
   }
}

// Local tone curves
void TMOEilertsen15::computeLocalToneCurves(double* luminance, int width, int height, std::vector<std::vector<double>>& localCurves){
   // Divide image into tiles
   int tilesX = (width + TILE_SIZE - 1) / TILE_SIZE;
   int tilesY = (height + TILE_SIZE - 1) / TILE_SIZE;
   
   NoiseModel noise = estimateNoise();
   localCurves.resize(tilesY);
   for (int ty = 0; ty < tilesY; ty++) {
      localCurves[ty].resize(tilesX * NUM_TONE_CURVE_NODES);
   }
   
   std::vector<double> globalHist(NUM_TONE_CURVE_NODES);
   computeWeightedHistogram(luminance, width, height, noise, globalHist, NUM_TONE_CURVE_NODES);
   
   // Compute tone curve for each tile
   for (int ty = 0; ty < tilesY; ty++) {
      for (int tx = 0; tx < tilesX; tx++) {
         int x0 = tx * TILE_SIZE;
         int y0 = ty * TILE_SIZE;
         int x1 = std::min(x0 + TILE_SIZE, width);
         int y1 = std::min(y0 + TILE_SIZE, height);
         
         std::vector<double> tileData;
         for (int y = y0; y < y1; y++) {
               for (int x = x0; x < x1; x++) {
                  tileData.push_back(luminance[y * width + x]);
               }
         }
         
         std::vector<double> localHist(NUM_TONE_CURVE_NODES);
         computeWeightedHistogram(tileData.data(), x1 - x0, y1 - y0, noise, localHist, NUM_TONE_CURVE_NODES);
         
         // Blend with global (10% global, 90% local)
         for (int k = 0; k < NUM_TONE_CURVE_NODES; k++) {
               localHist[k] = 0.1 * globalHist[k] + 0.9 * localHist[k];
         }
         
         std::vector<double> slopes(NUM_TONE_CURVE_NODES);
         computeToneCurve(localHist, noise, slopes);
         
         for (int k = 0; k < NUM_TONE_CURVE_NODES; k++) {
               localCurves[ty][tx * NUM_TONE_CURVE_NODES + k] = slopes[k];
         }
      }
   }
}

// Apply local tone curves with interpolation
void TMOEilertsen15::applyLocalToneCurves(double* base, int width, int height, const std::vector<std::vector<double>>& localCurves){ //TODO: Split image into squares as the paper suggests
   std::vector<double> globalSlopes(NUM_TONE_CURVE_NODES);
   for (int k = 0; k < NUM_TONE_CURVE_NODES; k++) {
      globalSlopes[k] = 1.0;
   }
   applyToneCurve(base, width, height, globalSlopes);
}

// Estimate noise
TMOEilertsen15::NoiseModel TMOEilertsen15::estimateNoise(){
   NoiseModel noise;
   std::vector<std::pair<double, double>> samples; // (mean, variance)
   
   int blockSize = 32;
   for (int by = 0; by < iHeight; by += blockSize) {
      for (int bx = 0; bx < iWidth; bx += blockSize) {
         double mean = computeBlockMean(bx, by, blockSize);
         double variance = computeBlockVariance(bx, by, blockSize);
         
         if (isUniform(bx, by, blockSize)) {
               samples.push_back({mean, variance});
         }
      }
   }
   
   fitLinearModel(samples, noise.a, noise.b);
   
   return noise;
}

double TMOEilertsen15::computeBlockMean(int bx, int by, int blockSize){
   double sum = 0.0;
   int count = 0;
   
   int endX = std::min(bx + blockSize, iWidth);
   int endY = std::min(by + blockSize, iHeight);
   
   for (int y = by; y < endY; y++) {
      for (int x = bx; x < endX; x++) {
         double r = pSrc->GetPixel(x, y)[0];
         double g = pSrc->GetPixel(x, y)[1];
         double b = pSrc->GetPixel(x, y)[2];
         double luminance = 0.299 * r + 0.587 * g + 0.114 * b;
         
         sum += luminance;
         count++;
      }
   }
   
   return (count > 0) ? sum / count : 0.0;
}

double TMOEilertsen15::computeBlockVariance(int bx, int by, int blockSize){
   double mean = computeBlockMean(bx, by, blockSize);
   double sumSquaredDiff = 0.0;
   int count = 0;
   
   int endX = std::min(bx + blockSize, iWidth);
   int endY = std::min(by + blockSize, iHeight);
   
   for (int y = by; y < endY; y++) {
      for (int x = bx; x < endX; x++) {
         double r = pSrc->GetPixel(x, y)[0];
         double g = pSrc->GetPixel(x, y)[1];
         double b = pSrc->GetPixel(x, y)[2];
         double luminance = 0.299 * r + 0.587 * g + 0.114 * b;
         
         double diff = luminance - mean;
         sumSquaredDiff += diff * diff;
         count++;
      }
   }
   
   return (count > 1) ? sumSquaredDiff / (count - 1) : 0.0;
}

bool TMOEilertsen15::isUniform(int bx, int by, int blockSize){
   double variance = computeBlockVariance(bx, by, blockSize);
   double mean = computeBlockMean(bx, by, blockSize);
   
   if (mean < 1e-6) return false;
   
   double stdDev = std::sqrt(variance);
   double coeffVariation = stdDev / mean;
   
   return coeffVariation < 0.1;
}

void TMOEilertsen15::fitLinearModel(const std::vector<std::pair<double, double>>& samples, double& a, double& b){
   if (samples.empty()) {
      a = 0.0;
      b = 0.0;
      return;
   }
   
   int n = samples.size();
   
   double meanX = 0.0;
   double meanY = 0.0;
   
   for (const auto& sample : samples) {
      meanX += sample.first;
      meanY += sample.second;
   }
   meanX /= n;
   meanY /= n;
   
   double covariance = 0.0;
   double varianceX = 0.0;
   
   for (const auto& sample : samples) {
      double dx = sample.first - meanX;
      double dy = sample.second - meanY;
      
      covariance += dx * dy;
      varianceX += dx * dx;
   }
   
   // Slope: a = Cov(X,Y) / Var(X)
   if (varianceX < 1e-10) {
      a = 0.0;
      b = meanY;
   } else {
      a = covariance / varianceX;
      b = meanY - a * meanX;
   }
   
   a = std::max(0.0, a);
   b = std::max(0.0, b);
}
