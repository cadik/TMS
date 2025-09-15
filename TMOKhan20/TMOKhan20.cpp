/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                  *
*                                                                              *
*                       Semestral project                                      *
*                       Author: Milan Tichavsky                                *
*                       Brno 2025                                              *
*                                                                              *
*                       Implementation of the TMOKhan20 class                  *
*                                                                              *
*******************************************************************************/

#include "TMOKhan20.h"

#include <iostream>
#include <vector>
#include <algorithm>

// Constants from the paper, see equation #1
#define M (1305.0/8192.0)
#define N (2523.0/32.0)
#define UNDEF_INT_PARAM 0
#define EPSILON 0.01


/**
 * Construct a luminance lookup table (LUT) based on histogram binning.
 *
 * @param bins Number of bins in the histogram, in the paper N
 * @param truncation Truncation factor to limit bin counts, in the paper k
 * @param luminance Input image luminance values
 * @param lum_len Length of the luminance array
 * @param min_luminance Minimum observed luminance in the input image
 * @param max_luminance Maximum observed luminance in the input image
 */
LuminanceLUT::LuminanceLUT(int bins, int truncation, const std::vector<double> &luminance, int lum_len,
                           double min_luminance, double max_luminance) {
    // At first, the LUT with N + 1 rows is implemented as two vectors
    std::vector<double> edges(bins + 1); // T(1) calculation
    std::vector<int> bin_count(bins + 1, 0); // T(2) calculation

    // Equation #3 in the paper
    for (int i = 0; i <= bins; i++) {
        edges[i] = min_luminance + i * (max_luminance - min_luminance) / bins;
    }
    // Otherwise for max luminance, we would be getting out of bounds index in upper_bound function
    edges[bins] += EPSILON;

    // Equation #4 in the paper
    for (int i = 0; i < lum_len; i++) {
        long bin_index = std::upper_bound(edges.begin(), edges.end(), luminance[i]) - edges.begin();
        if (bin_index >= 0 && bin_index < bins + 1) {
            bin_count[bin_index]++;
        }
    }
    int prev = 0;
    int max_bin_count = static_cast<int>(lum_len * static_cast<double>(truncation) / bins);
    for (int i = 0; i < bins + 1; i++) {
        if (bin_count[i] > max_bin_count) { // truncate to limit max bin count
            bin_count[i] = max_bin_count;
        }
        bin_count[i] += prev;
        prev = bin_count[i];
    }
    // Equation #4, case for i=0
    assert(bin_count[0] == 0 && "LUT creation poorly implemented, see equation #4 case i=0.");

    // Equation #5 in the paper
    for (int i = 0; i < bins + 1; i++) {
        bin_count[i] = 255 * bin_count[i] / bin_count[bins];
    }

    for (int i = 0; i < bins + 1; i++) {
        lut.emplace_back(edges[i], bin_count[i]);
    }
    assert(lut.size() == bins + 1 && "Implementation error (unexpected LUT size).");
}

/**
 * Return the bin edges for the given key luminance value. Keys represent HDR values, and the corresponding values
 * represent mapped LDR values.
 *
 * @param key Input luminance value to search for.
 * @return Tuple containing (lower key, lower value, upper key, upper value).
 */
std::tuple<double, int, double, int> LuminanceLUT::getValue(double key) const {
    auto it = std::upper_bound(lut.begin(), lut.end(), key,
                               [](double val, const std::pair<double, int> &elem) {
                                   return val < elem.first;
                               });

    if (it == lut.begin()) {
        throw std::runtime_error("Invalid lookup key: " + std::to_string(key));
    } else if (it == lut.end()) {
        return std::make_tuple(lut.back().first, lut.back().second, lut.back().first, lut.back().second);
    }

    auto upper = it;
    auto lower = std::prev(it);
    return std::make_tuple(lower->first, lower->second, upper->first, upper->second);
}

/** Print LUT to stderr for debug purposes */
[[maybe_unused]] void LuminanceLUT::printLUT() const {
    for (const auto &[edge, value]: lut) {
        std::cerr << "Edge: " << edge << " -> Value: " << value << "\n";
    }
}

TMOKhan20::TMOKhan20() {
    SetName(L"Khan20");
    SetDescription(L"Tone-Mapping Using Perceptual-Quantizer and Image Histogram");

    binParameter.SetName(L"Bins");
    binParameter.SetDescription(L"Labeled as N; the number of bins in the histogram");
    binParameter.SetDefault(256);
    binParameter = UNDEF_INT_PARAM;
    binParameter.SetRange(1, 100000);  // max chosen arbitrarily
    this->Register(binParameter);
    truncationParameter.SetName(L"Truncation");
    truncationParameter.SetDescription(
        L"Labeled as k, where `(k/N)*nof_pixels` is the maximum (truncated) count of pixels in one histogram bin"
    );
    truncationParameter.SetDefault(5);
    truncationParameter = UNDEF_INT_PARAM;
    truncationParameter.SetRange(1, 100000);  // max chosen arbitrarily
    this->Register(truncationParameter);
}

TMOKhan20::~TMOKhan20()
= default;

/**
 * Apply tone mapping to the luminance values using LUT-based interpolation.
 *
 * @param luminance Input/output luminance values. Note that the input (HDR) luminance is in range [0, 1], after calling
 *                  the function it is transformed (LDR) to [0, 255]
 * @param lum_len Length of the luminance array.
 * @param lut Precomputed luminance LUT.
 */
void TMOKhan20::ToneMap(std::vector<double> &luminance, int lum_len, LuminanceLUT &lut) {
    for (int i = 0; i < lum_len; i++){
        // Equation #6 in the paper
        auto [k1, v1, k2, v2] = lut.getValue(luminance[i]);
        if (v2 - v1 == 0) {
            luminance[i] = v1;
        } else if (k2 - k1 == 0) {
            throw std::runtime_error("Trying to divide by zero in ToneMap method");
        } else {
            luminance[i] = v1 + (v2 - v1) * (luminance[i] - k1) / (k2 - k1);
        }
    }
}

/**
 * Transform luminance using a PQ equation from the paper.
 *
 * @param pSrc Input image.
 * @param luminance Output luminance values in range [0, 1]
 * @return Tuple containing (min_luminance, max_luminance) in the output luminance vector.
 */
std::tuple<double, double> TMOKhan20::ApplyPerceptualQuantizer(std::vector<double> &luminance) {
    double *pSourceData = pSrc->GetData();
    int k = 0;
    double max_luminance = 0;
    double min_luminance = 1;
    for (int j = 0; j < pSrc->GetHeight(); j++) {
        for (int i = 0; i < pSrc->GetWidth(); i++) {
            double red = *pSourceData++;
            double green = *pSourceData++;
            double blue = *pSourceData++;

            // Equation #2 from the paper
            double l_in = 0.2126 * red + 0.7152 * green + 0.0722 * blue;
            if (l_in > 10000) {
                l_in = 10000;
            } else if (l_in < 0) {
                l_in = 0;
            }

            // Equation #1 from the paper
            double luminance_out = ((107.0 + 2413.0 * pow(l_in / 10000.0, M)) /
                                    (128.0 + 2392.0 * pow(l_in / 10000.0, M)));
            luminance_out = pow(luminance_out, N);

            if (luminance_out > max_luminance) {
                max_luminance = luminance_out;
            }
            if (luminance_out < min_luminance) {
                min_luminance = luminance_out;
            }

            luminance[k++] = luminance_out;
        }
    }
    return std::make_tuple(min_luminance, max_luminance);
}

/**
 * Keeping the original color, just changing the luminance in Yxy model. This step wasn't properly discussed
 * in the paper, so doing the final conversion to RGB what I consider to be the standard way.
 *
 * @param luminance the tone mapped luminance in [0, 255] range
 */
void TMOKhan20::ApplyTransformationToDstImage(std::vector<double> luminance) {
    pSrc->Convert(TMO_Yxy);
    pDst->Convert(TMO_Yxy);
    double *pSourceData = pSrc->GetData();
    double *destination = pDst->GetData();
    for (int i = 0; i < pSrc->GetWidth() * pDst->GetHeight(); i++) {
        pSourceData++;
        // Transforming to Yxy [0,1] range while copying to the output image
        *destination++ = luminance[i] / 255;
        *destination++ = *pSourceData++;
        *destination++ = *pSourceData++;
    }
}

int TMOKhan20::Transform() {
    if (binParameter == 0) binParameter = binParameter.GetDefault();
    if (truncationParameter == 0) truncationParameter = truncationParameter.GetDefault();

    pSrc->Convert(TMO_RGB);
    int luminance_len = pSrc->GetHeight() * pSrc->GetWidth();
    std::vector<double> luminance(luminance_len);
    auto [min_luminance, max_luminance] = TMOKhan20::ApplyPerceptualQuantizer(luminance);

    LuminanceLUT lut(binParameter, truncationParameter, luminance, luminance_len, min_luminance, max_luminance);
    ToneMap(luminance, luminance_len, lut);
    ApplyTransformationToDstImage(luminance);
    return 0;
}

