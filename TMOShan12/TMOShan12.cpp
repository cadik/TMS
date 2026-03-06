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
*                      Implementation of the TMOShan12 class                   *
*               Tone Mapping High Dynamic Range Videos using Wavelets          *
*                                                                              *
*******************************************************************************/

#include "TMOShan12.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <numeric>

// Static definitions
constexpr double TMOShan12::EPSILON;
constexpr double TMOShan12::LP_FILTER[3];
constexpr double TMOShan12::HP_FILTER[5];
constexpr double TMOShan12::RECON_LP[3];
constexpr double TMOShan12::RECON_HP[5];

TMOShan12::TMOShan12() {
   SetName(L"Shan12");
   SetDescription(L"Tone Mapping HDR Videos using Wavelets (Shan, Meyer, DeRose, Anderson 2012");


   dGamma.SetName(L"Gamma");
   dGamma.SetDescription(L"Gamma: Controls gain function linearity. As gamma->1, gain approaches constant.");
   dGamma.SetDefault(0.6);
   dGamma.SetRange(0.1, 0.95);
   this->Register(dGamma);
   dGamma = 0.6;

   dDelta.SetName(L"Delta");
   dDelta.SetDescription(L"Delta: Normalization constant for gain function.");
   dDelta.SetDefault(1.0);
   dDelta.SetRange(0.01, 100.0);
   this->Register(dDelta);
   dDelta = 1.0;

   iNumBands.SetName(L"NumBands");
   iNumBands.SetDescription(L"N: Number of wavelet decomposition bands. Paper uses 4-9, typically 4 or 5.");
   iNumBands.SetDefault(5);
   iNumBands.SetRange(2, 9);
   this->Register(iNumBands);
   iNumBands = 5;

   dAlpha.SetName(L"Alpha");
   dAlpha.SetDescription(L"Alpha: Edge detection threshold for temporal edge avoidance.");
   dAlpha.SetDefault(0.5);
   dAlpha.SetRange(0.01, 2.0);
   this->Register(dAlpha);
   dAlpha = 0.5;

   bImprovedGain.SetName(L"Improved Gain");
   bImprovedGain.SetDescription(L"Use improved gain control. Reduces haloing compared to Li et al. 2005.");
   bImprovedGain.SetDefault(true);
   this->Register(bImprovedGain);
   bImprovedGain = true;

   iTemporalWindow.SetName(L"Temporal Window");
   iTemporalWindow.SetDescription(L"Number of frames in temporal cache. LP filter needs >=3, HP filter needs >=5. Set to 1 for single-image mode (no temporal).");
   iTemporalWindow.SetDefault(5);
   iTemporalWindow.SetRange(1, 15);
   this->Register(iTemporalWindow);
   iTemporalWindow = 5;
}

TMOShan12::~TMOShan12() {
}

int TMOShan12::Transform() {
   // Work in RGB space: extract luminance, tone-map, scale RGB by ratio
   pSrc->Convert(TMO_RGB);
   pDst->Convert(TMO_RGB);

   const double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int numPixels = width * height;
   int N = iNumBands.GetInt();
   double gamma = dGamma.GetDouble();
   double delta = dDelta.GetDouble();
   double alpha = dAlpha.GetDouble();
   bool useImproved = bImprovedGain.GetBool();
   int temporalWindow = iTemporalWindow.GetInt();


   // Step 1: Extract RGB and compute luminance
   std::vector<double> srcR(numPixels), srcG(numPixels), srcB(numPixels);
   std::vector<double> luminance(numPixels);

   for (int i = 0; i < numPixels; i++)
   {
      srcR[i] = pSourceData[i * 3];
      srcG[i] = pSourceData[i * 3 + 1];
      srcB[i] = pSourceData[i * 3 + 2];

      // CIE luminance from RGB
      luminance[i] = 0.2126 * srcR[i] + 0.7152 * srcG[i] + 0.0722 * srcB[i];

      // Ensure positive luminance for log transform
      if (luminance[i] < EPSILON)
         luminance[i] = EPSILON;
   }

   // Work in log domain as is standard for wavelet-based TMOs
   std::vector<double> logLum(numPixels);
   for (int i = 0; i < numPixels; i++)
   {
      logLum[i] = std::log(luminance[i]);
   }

   // Reset cache if frame dimensions changed
   if (width != lastWidth || height != lastHeight)
   {
      temporalCache.clear();
      lastWidth = width;
      lastHeight = height;
   }

   // Update cache parameters
   temporalCache.maxFrames = temporalWindow;
   temporalCache.cachedWidth = width;
   temporalCache.cachedHeight = height;

   temporalCache.frameHistory.push_back(logLum);

   while ((int)temporalCache.frameHistory.size() > temporalCache.maxFrames)
   {
      temporalCache.frameHistory.pop_front();
   }

   pSrc->ProgressBar(0, 10);

   // Step 2: Wavelet decomposition
   WaveletPyramid pyramid;
   bool useTemporalFiltering = (temporalWindow > 1) && ((int)temporalCache.frameHistory.size() >= 3);

   if (useTemporalFiltering) {
      // 3D edge-avoiding wavelet decomposition
      temporalEdgeAvoidingDecomposition(width, height, N, alpha, pyramid);
   } else {
      // Pure 2D spatial wavelet decomposition
      waveletAnalysis2D(logLum, width, height, N, pyramid);
   }

   pSrc->ProgressBar(3, 10);

   // Step 3: Gain control

   // Compute activity maps for each band
   std::vector<std::vector<double>> activities(N);
   for (int i = 0; i < N; i++) {
      activities[i].resize(numPixels);
      computeBandActivity(pyramid.bands[i], activities[i], width, height, N - 1 - i);
   }

   // Compute aggregate activity
   std::vector<double> aggregate(numPixels, 0.0);
   computeAggregateActivity(activities, aggregate, width, height);

   pSrc->ProgressBar(5, 10);

   // Apply gain to each band
   for (int i = 0; i < N; i++) {
      std::vector<double> gainMap(numPixels);

      if (useImproved) {
         // Improved gain control
         computeGainMapImproved(aggregate, width, height, i, N, gamma, delta, gainMap);
      } else {
         // Li et al. 2005 gain control
         computeGainMapLi(aggregate, width, height, i, N, gamma, delta, gainMap);
      }

      // Apply gain with level multiplier
      double mi = getLevelMultiplier(i, N);
      for (int p = 0; p < numPixels; p++) {
         pyramid.bands[i][p] *= gainMap[p] * mi;
      }
   }

   pSrc->ProgressBar(7, 10);

   // Compress lowpass to target display range using linear scaling in log domain
   {
      double lpMin = *std::min_element(pyramid.lowpass.begin(), pyramid.lowpass.end());
      double lpMax = *std::max_element(pyramid.lowpass.begin(), pyramid.lowpass.end());
      double lpRange = lpMax - lpMin;

      // Compute the mean of the lowpass as the anchor point
      double lpMean = 0.0;
      for (int i = 0; i < numPixels; i++)
      {
         lpMean += pyramid.lowpass[i];
      }
      lpMean /= numPixels;

      if (lpRange > EPSILON)
      {
         // Compress to target range in log domain
         double targetRange = 3.5;
         double scale = targetRange / lpRange;

         for (int i = 0; i < numPixels; i++)
         {
               // Compress around the mean to preserve average scene brightness
               pyramid.lowpass[i] = (pyramid.lowpass[i] - lpMean) * scale + lpMean;
         }
      }
   }


   // Step 4: Synthesis - Reconstruct tone-mapped luminance
   std::vector<double> logLumTM(numPixels);
   waveletSynthesis2D(pyramid, logLumTM);

   pSrc->ProgressBar(8, 10);

   // Step 5: Convert from log domain, normalize, and apply to RGB
   std::vector<double> linLumTM(numPixels);
   for (int i = 0; i < numPixels; i++) {
      linLumTM[i] = std::exp(logLumTM[i]);
   }

   // Percentile-based normalization to avoid outlier sensitivity
   std::vector<double> sorted = linLumTM;
   std::sort(sorted.begin(), sorted.end());

   int idxLow = (int)(0.005 * numPixels);
   int idxHigh = (int)(0.995 * numPixels);
   idxLow = std::max(0, std::min(numPixels - 1, idxLow));
   idxHigh = std::max(0, std::min(numPixels - 1, idxHigh));

   double pLow = sorted[idxLow];
   double pHigh = sorted[idxHigh];
   double normRange = pHigh - pLow;
   if (normRange < EPSILON) {
      normRange = 1.0;
   }

   // Display gamma for perceptual uniformity
   double displayGamma = 1.0 / 2.2;

   pSrc->ProgressBar(9, 10);

   for (int i = 0; i < numPixels; i++) {
      // Normalize tone-mapped luminance to [0, 1] in linear space
      double YtmLinear = (linLumTM[i] - pLow) / normRange;
      YtmLinear = std::max(0.0, std::min(1.0, YtmLinear));

      // Scale original RGB channels by the linear luminance ratio
      double Yorig = luminance[i];
      double ratio = (Yorig > EPSILON) ? (YtmLinear / Yorig) : 0.0;

      double R = srcR[i] * ratio;
      double G = srcG[i] * ratio;
      double B = srcB[i] * ratio;

      // Clamp negatives before gamma (negative^fractional = NaN)
      R = std::max(0.0, R);
      G = std::max(0.0, G);
      B = std::max(0.0, B);

      // Apply display gamma to the final RGB channels
      R = std::pow(R, displayGamma);
      G = std::pow(G, displayGamma);
      B = std::pow(B, displayGamma);

      double maxCh = std::max(R, std::max(G, B));
      if (maxCh > 1.0)
      {
         double YtmGamma = std::pow(YtmLinear, displayGamma);
         if (maxCh > YtmGamma + EPSILON)
         {
               // Blend toward achromatic luminance so max channel == 1.0
               double t = (1.0 - YtmGamma) / (maxCh - YtmGamma);
               t = std::max(0.0, std::min(1.0, t));
               R = YtmGamma + t * (R - YtmGamma);
               G = YtmGamma + t * (G - YtmGamma);
               B = YtmGamma + t * (B - YtmGamma);
         }
         else
         {
               R = std::min(1.0, R);
               G = std::min(1.0, G);
               B = std::min(1.0, B);
         }
      }

      pDestinationData[i * 3]     = std::max(0.0, std::min(1.0, R));
      pDestinationData[i * 3 + 1] = std::max(0.0, std::min(1.0, G));
      pDestinationData[i * 3 + 2] = std::max(0.0, std::min(1.0, B));
   }

   pSrc->ProgressBar(10, 10);
   return 0;
}

// Temporal edge-avoiding low-pass filter
void TMOShan12::temporalEdgeAvoidingLowPass(const std::deque<std::vector<double>> &frames, int centerIdx, std::vector<double> &output, int width, int height, double alpha) {
   int numPixels = width * height;
   int numFrames = (int)frames.size();
   output.resize(numPixels);

   const double f[3] = {0.25, 0.50, 0.25};
   const int filterHalfWidth = 1;

   const std::vector<double> &centerFrame = frames[centerIdx];

   for (int p = 0; p < numPixels; p++) {
      double sum = 0.0;
      double weightSum = 0.0;

      for (int t = -filterHalfWidth; t <= filterHalfWidth; t++) {
         int frameIdx = centerIdx + t;

         // Mirror boundary in temporal direction
         if (frameIdx < 0) {
            frameIdx = -frameIdx;
         }
         if (frameIdx >= numFrames){
            frameIdx = 2 * (numFrames - 1) - frameIdx;
         }
         frameIdx = std::max(0, std::min(numFrames - 1, frameIdx));

         double filterCoeff = f[t + filterHalfWidth];

         // Edge-avoiding: ghost particle if temporal discontinuity
         double value;
         if (isEdge(centerFrame[p], frames[frameIdx][p], alpha)) {
            value = centerFrame[p];
         }
         else {
            value = frames[frameIdx][p];
         }

         sum += value * filterCoeff;
         weightSum += filterCoeff;
      }

      output[p] = sum / weightSum;
   }
}

// Temporal edge-avoiding band decomposition
void TMOShan12::temporalEdgeAvoidingDecomposition(int width, int height, int numBands, double alpha, WaveletPyramid &pyramid) {
   int numPixels = width * height;
   int numFrames = (int)temporalCache.frameHistory.size();
   int centerIdx = numFrames - 1;

   // Step 1: Perform 2D spatial wavelet analysis on each cached frame
   std::vector<WaveletPyramid> framePyramids(numFrames);
   for (int f = 0; f < numFrames; f++) {
      waveletAnalysis2D(temporalCache.frameHistory[f], width, height, numBands, framePyramids[f]);
   }

   // Step 2: For each spatial band, apply temporal edge-avoiding filtering
   pyramid.width = width;
   pyramid.height = height;
   pyramid.numBands = numBands;
   pyramid.bands.resize(numBands);

   // For each band, collect the band data across frames and filter temporally
   for (int b = 0; b < numBands; b++) {
      pyramid.bands[b].resize(numPixels);

      std::deque<std::vector<double>> bandAcrossFrames;
      for (int f = 0; f < numFrames; f++) {
         bandAcrossFrames.push_back(framePyramids[f].bands[b]);
      }

      std::vector<double> temporalLP(numPixels);
      temporalEdgeAvoidingLowPass(bandAcrossFrames, centerIdx, temporalLP, width, height, alpha);

      pyramid.bands[b] = temporalLP;
   }

   // Also temporally filter the low-pass residual
   {
      std::deque<std::vector<double>> lpAcrossFrames;
      for (int f = 0; f < numFrames; f++) {
         lpAcrossFrames.push_back(framePyramids[f].lowpass);
      }

      temporalEdgeAvoidingLowPass(lpAcrossFrames, centerIdx, pyramid.lowpass, width, height, alpha);
   }
}

// Wavelet Analysis 2D
void TMOShan12::waveletAnalysis2D(const std::vector<double> &input, int width, int height, int numBands, WaveletPyramid &pyramid) {
   int numPixels = width * height;
   pyramid.width = width;
   pyramid.height = height;
   pyramid.numBands = numBands;
   pyramid.bands.resize(numBands);

   // Start from the finest level
   std::vector<double> current = input;

   for (int level = numBands - 1; level >= 0; level--) {
      // Compute low-pass at this level
      std::vector<double> lowpass(numPixels);
      lowPassFilter2D(current, lowpass, width, height, numBands - 1 - level);

      // Detail band = current - lowpass
      pyramid.bands[level].resize(numPixels);
      for (int i = 0; i < numPixels; i++) {
         pyramid.bands[level][i] = current[i] - lowpass[i];
      }

      current = lowpass;
   }

   pyramid.lowpass = current;
}

// Wavelet Synthesis 2D
void TMOShan12::waveletSynthesis2D(const WaveletPyramid &pyramid, std::vector<double> &output) {
   int numPixels = pyramid.width * pyramid.height;
   output = pyramid.lowpass;

   // Add all detail bands (from coarsest to finest)
   for (int level = 0; level < pyramid.numBands; level++) {
      for (int i = 0; i < numPixels; i++) {
         output[i] += pyramid.bands[level][i];
      }
   }
}

// 2D Low-pass filter (separable)
void TMOShan12::lowPassFilter2D(const std::vector<double> &input, std::vector<double> &output, int width, int height, int level) {
   int numPixels = width * height;
   std::vector<double> temp(numPixels);

   // Horizontal pass
   lowPassFilter1D(input, temp, width, 1, height, width, level);

   // Vertical pass
   lowPassFilter1D(temp, output, height, width, width, 1, level);
}

// 1D Low-pass filter
void TMOShan12::lowPassFilter1D(const std::vector<double> &input, std::vector<double> &output, int length, int stride, int numLines, int lineStride, int level) {
   // Dilation factor for stationary wavelets
   int dilation = 1 << level;

   const double f0 = 0.25;
   const double f1 = 0.50;
   const double f2 = 0.25;

   for (int line = 0; line < numLines; line++) {
      for (int i = 0; i < length; i++) {
         // Mirror boundary handling
         int idx_m1 = i - dilation;
         int idx_p1 = i + dilation;

         // Reflect at boundaries
         if (idx_m1 < 0) {
            idx_m1 = -idx_m1;
         }
         if (idx_p1 >= length) {
            idx_p1 = 2 * (length - 1) - idx_p1;
         }

         // Clamp to valid range
         idx_m1 = std::max(0, std::min(length - 1, idx_m1));
         idx_p1 = std::max(0, std::min(length - 1, idx_p1));

         double val = f0 * input[line * lineStride + idx_m1 * stride] + f1 * input[line * lineStride + i * stride] + f2 * input[line * lineStride + idx_p1 * stride];

         output[line * lineStride + i * stride] = val;
      }
   }
}

// Compute band activity
void TMOShan12::computeBandActivity(const std::vector<double> &band, std::vector<double> &activity, int width, int height, int level) {
   int numPixels = width * height;

   // Compute absolute value of band
   std::vector<double> absBand(numPixels);
   for (int i = 0; i < numPixels; i++) {
      absBand[i] = std::fabs(band[i]);
   }

   double sigma = std::pow(2.0, level + 1);
   gaussianBlur(absBand, activity, width, height, sigma);
}

// Compute aggregate activity
void TMOShan12::computeAggregateActivity(const std::vector<std::vector<double>> &activities, std::vector<double> &aggregate, int width, int height) {
   int numPixels = width * height;
   std::fill(aggregate.begin(), aggregate.end(), 0.0);

   for (size_t i = 0; i < activities.size(); i++) {
      for (int p = 0; p < numPixels; p++) {
         aggregate[p] += activities[i][p];
      }
   }
}

// Gain map: Li et al. 2005
void TMOShan12::computeGainMapLi(const std::vector<double> &aggregate, int width, int height, int bandIdx, int numBands, double gamma, double delta, std::vector<double> &gainMap) {
   int numPixels = width * height;
   double exponent = gamma - 1.0;

   for (int i = 0; i < numPixels; i++) {
      double A = aggregate[i];
      gainMap[i] = std::pow((A + EPSILON) / delta, exponent);

      // Clamp gain to prevent extreme amplification
      if (!std::isfinite(gainMap[i]) || gainMap[i] < 0.0) {
         gainMap[i] = 0.0;
      }
      if (gainMap[i] > 10.0) {
         gainMap[i] = 10.0;
      }
   }
}

// Gain map: Improved
void TMOShan12::computeGainMapImproved(const std::vector<double> &aggregate, int width, int height, int bandIdx, int numBands, double gamma, double delta, std::vector<double> &gainMap) {
   int numPixels = width * height;

   double gammaP = std::min(gamma + 0.05 * (numBands - 1 - bandIdx), 0.9);
   double exponentP = gammaP - 1.0;
   double exponentOrig = gamma - 1.0;

   double sumAG_orig = 0.0;
   double sumAG_new_numer = 0.0;

   for (int i = 0; i < numPixels; i++) {
      double A = aggregate[i];
      double Aeps = A + EPSILON;

      sumAG_orig += A * std::pow(Aeps / delta, exponentOrig);
      sumAG_new_numer += A * std::pow(Aeps, exponentP);
   }

   double deltaP;
   if (std::fabs(exponentP) > 1e-10 && sumAG_orig > 1e-20) {
      double ratio = sumAG_new_numer / sumAG_orig;
      deltaP = std::pow(ratio, 1.0 / exponentP);
   }
   else {
      deltaP = delta;
   }

   if (deltaP < EPSILON || !std::isfinite(deltaP)) {
      deltaP = delta;
   }

   // Compute the gain map with adjusted parameters
   for (int i = 0; i < numPixels; i++) {
      double A = aggregate[i];
      gainMap[i] = std::pow((A + EPSILON) / deltaP, exponentP);

      if (!std::isfinite(gainMap[i]) || gainMap[i] < 0.0) {
         gainMap[i] = 0.0;
      }
      if (gainMap[i] > 10.0) {
         gainMap[i] = 10.0;
      }
   }
}

// Level-dependent multipliers (finest=1.0, next=0.7, rest=0.4)
double TMOShan12::getLevelMultiplier(int bandIdx, int numBands) {
   int fromFinest = numBands - 1 - bandIdx;

   if (fromFinest == 0) {
      return 1.0;   // Finest level
   } else if (fromFinest == 1) {
      return 0.7;   // Second finest
   } else {
      return 0.4;   // All coarser levels
   }
}

// Gaussian blur
void TMOShan12::gaussianBlur(const std::vector<double> &input, std::vector<double> &output, int width, int height, double sigma) {
   int numPixels = width * height;
   output.resize(numPixels);

   int kernelRadius = std::max(1, (int)std::ceil(3.0 * sigma));
   kernelRadius = std::min(kernelRadius, std::min(width / 2, height / 2));

   // Precompute 1D Gaussian kernel
   std::vector<double> kernel(2 * kernelRadius + 1);
   double kernelSum = 0.0;
   double twoSigmaSq = 2.0 * sigma * sigma;

   for (int k = -kernelRadius; k <= kernelRadius; k++) {
      kernel[k + kernelRadius] = std::exp(-(double)(k * k) / twoSigmaSq);
      kernelSum += kernel[k + kernelRadius];
   }
   // Normalize
   for (int k = 0; k < (int)kernel.size(); k++) {
      kernel[k] /= kernelSum;
   }

   // Horizontal pass
   std::vector<double> temp(numPixels);
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         double val = 0.0;
         for (int k = -kernelRadius; k <= kernelRadius; k++) {
            int nx = x + k;
            // Mirror boundary
            if (nx < 0) nx = -nx;
            if (nx >= width) nx = 2 * (width - 1) - nx;
            nx = std::max(0, std::min(width - 1, nx));

            val += input[y * width + nx] * kernel[k + kernelRadius];
         }
         temp[y * width + x] = val;
      }
   }

   // Vertical pass
   for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
         double val = 0.0;
         for (int k = -kernelRadius; k <= kernelRadius; k++) {
            int ny = y + k;
            // Mirror boundary
            if (ny < 0) ny = -ny;
            if (ny >= height) ny = 2 * (height - 1) - ny;
            ny = std::max(0, std::min(height - 1, ny));

            val += temp[ny * width + x] * kernel[k + kernelRadius];
         }
         output[y * width + x] = val;
      }
   }
}
