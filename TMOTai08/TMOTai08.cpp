/**
 * A Two-Stage Contrast Enhancement Algorithm for Digital Images
 * https://doi.org/10.1109/CISP.2008.400
 *
 * Author: Michal Vlnas
 */

#include "TMOTai08.h"
#include <memory>
#include <opencv2/imgproc/imgproc.hpp>

TMOTai08::TMOTai08()
{
    SetName(L"Tai08");
    SetDescription(L"A Two-Stage Contrast Enhancement Algorithm for Digital Images");

    tc.SetName(L"tc");
    tc.SetDescription(L"lumimance of the threshold level");
    tc.SetDefault(0.35);
    tc = 0.35;
    this->Register(tc);

    kc.SetName(L"kc");
    kc.SetDescription(L"denotes the luminance level at the knee point of the conventional knee curve");
    kc.SetDefault(0.9);
    kc = 0.9;
    this->Register(kc);

    gamma.SetName(L"gamma");
    gamma.SetDescription(L"denotes camera-gamma");
    gamma.SetDefault(2.2);
    gamma = 2.2;
    this->Register(gamma);

    th.SetName(L"th");
    th.SetDescription(L"edge threshold");
    th.SetDefault(0.5);
    th = 0.5;
    this->Register(th);
}

int TMOTai08::Transform()
{
    double* source = pSrc->GetData();
    double* destination = pDst->GetData();

    // set these internal variables as in tmogui they are not set, but required in this algorithm
    iWidth = pSrc->GetWidth();
    iHeight = pSrc->GetHeight();

    double max = -1;

    pSrc->Convert(TMO_Yxy);
    pDst->Convert(TMO_Yxy);
    double pY, px, py;

    t = tc;
    k = kc;
    // set fixed value, according to the article
    m = 2.0;

    // get manually min and max luminance from Yxy
    for (int j = 0; j < pSrc->GetHeight(); j++)
    {
        for (int i = 0; i < pSrc->GetWidth(); i++)
        {
            double lum = pSrc->GetLuminanceYxy(i, j);
            if (lum > max)
                max = lum;
        }
    }

    // get max in log space
    gmax = std::log10(max + 1.0);

    // compute cubic curve coefficients
    double s = (1.0 - k) / (t - k);
    double tm = std::pow(t - m, 3.0);

    a = ((s - 1.0) * t + (2.0 - m - m * s)) / tm;
    b = (2.0 * (1.0 - s) * t * t + (m * s + 2.0 * m - 3.0) * (t + m)) / tm;
    c = (s * t * t * t + (s - 4.0) * m * t * t + (6.0 - m - 2.0 * m * s) * m * t - m * m * m) / tm;
    d = ((1.0 - m * s) * t * t * t + (m * s + 2.0 * m - 3.0) * m * t * t) / tm;

    // compute Sobel filter on the whole image
    cv::Mat sobel_src(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC1);
    lum_avg = cv::Mat(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC1);

    for (int j = 0; j < pSrc->GetHeight(); j++)
    {
        for (int i = 0; i < pSrc->GetWidth(); i++)
        {
            sobel_src.at<double>(j, i) = NL(pSrc->GetLuminanceYxy(i, j));
            lum_avg.at<double>(j, i) = NL(pSrc->GetLuminanceYxy(i, j, 1));
        }
    }

    // move 't' to luminance domain
    t = std::pow(tc, 1.0 / gamma);

    cv::Mat grad_x, grad_y;
    cv::Sobel(sobel_src, grad_x, CV_64F, 1, 0, 3);
    cv::Sobel(sobel_src, grad_y, CV_64F, 0, 1, 3);
    cv::addWeighted(grad_x, 0.5, grad_y, 0.5, 0, grad);
    BuildBetaMap();

    for (int j = 0; j < pSrc->GetHeight(); j++)
    {
        pSrc->ProgressBar(j, pSrc->GetHeight());
        for (int i = 0; i < pSrc->GetWidth(); i++)
        {
            pY = *source++;
            px = *source++;
            py = *source++;

            double w = NL(pY);

            double w_avg = lum_avg.at<double>(j, i);
            double beta_value = beta.at<double>(j, i);

            // get alpha coefficient according to mean around the pixel
            double alpha = GetAlpha(j, i, w, w_avg);
            // finish exponent with curve differential
            double exponent = alpha * (1.0 - (w / _ApproxKneeCurve(w)) * _ApproxKneeCurveDifferential(w));
            // alpha mapping
            double E = std::pow((w / w_avg), exponent);
            // contrast-enhancement with log curve and beta mapping
            double frac = std::log(_ApproxKneeCurve(w) * (beta_value - 1.0) + 1.0) / std::log(beta_value);
            double v = 2 * frac * E;

            // and store results to the destination image
            *destination++ = 2.0 * v;
            *destination++ = px;
            *destination++ = py;
        }
    }

    // convert to RGB and correct gamma
    pDst->Convert(TMO_RGB);
    pDst->CorrectGamma(gamma);
    return 0;
}

double TMOTai08::_ApproxKneeCurve(double value)
{
    if (value < t)
        return value;
    else
        return std::min(std::pow((a * std::pow(value, 3.0) + b * value * value + c * value + d), 1.0 / gamma), 1.0);
}

double TMOTai08::_ApproxKneeCurveDifferential(double value)
{
    if (value < t)
        return 1.0;
    else
        return std::min(1.0 / gamma * std::pow((a * std::pow(value, 3.0) + b * value * value + c * value + d), 1.0 / gamma - 1.0)
            * (3 * a * value * value + 2 * b * value + c), 1.0);
}

double TMOTai08::GetAlpha(int x, int y, double value, double avg)
{
    double sobel = std::abs(grad.at<double>(x, y));
    double kweight = 1.0 / (std::exp(sobel));

    return value > avg ? 1.0 + (1.0 - kweight) * 0.75 : kweight;
}

void TMOTai08::BuildBetaMap()
{
    beta = cv::Mat(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC1);

    // use uint8 matrix for pixel predictions
    predictor = cv::Mat(pSrc->GetHeight(), pSrc->GetWidth(), CV_8U);

    // 1. setup beta map for non-edge pixels
    for (int j = 0; j < pSrc->GetHeight(); j++)
    {
        for (int i = 0; i < pSrc->GetWidth(); i++)
        {
            double sobel = std::abs(grad.at<double>(j, i));
            if (sobel < th)
            {
                // compute mean around the selected pixel
                double lum = NL(pSrc->GetLuminanceYxy(i, j));
                double mean = GetMeanAroundPixel(i, j, lum);

                if (lum >= mean)
                {
                    // isnan is required sanity check as for very low lum here could be possible nan value
                    double c_beta = 10.0 * (lum - mean) / lum + 2.0;
                    beta.at<double>(j, i) = std::isnan(c_beta) ? 2.0 : c_beta;
                }
                else
                    beta.at<double>(j, i) = 2;

                predictor.at<uint8_t>(j, i) = BETA_DEFINED;
            }
            else // undefined
            {
                predictor.at<uint8_t>(j, i) = BETA_UNDEFINED;
                beta.at<double>(j, i) = -2.0;
            }
        }
    }

    // 2. finish beta map with edge pixels
    uint64 undefined_count = 0;
    do
    {
        // reset count
        undefined_count = 0;

        for (int j = 0; j < pSrc->GetHeight(); j++) // y
        {
            for (int i = 0; i < pSrc->GetWidth(); i++) // x
            {
                if (predictor.at<uint8_t>(j, i) == BETA_UNDEFINED)
                {
                    double lum = NL(pSrc->GetLuminanceYxy(i, j));
                    double beta_value = GetBetaAroundPixel(i, j, lum);
                    if (beta_value >= 0.0)
                    {
                        // update beta map and predictors
                        beta.at<double>(j, i) = beta_value;
                        predictor.at<uint8_t>(j, i) = BETA_DEFINED;
                    }
                    else
                        undefined_count++;
                }
            }
        }
    }
    while (undefined_count > 0);
}

double TMOTai08::GetMeanAroundPixel(int x, int y, double value)
{
    int xp, yp, count = 0, offset;
    double retval = .0;

    int r = 5;

    for (int i = -r; i <= r; i++)
    {
        for (int j = -r; j <= r; j++)
        {
            // skip itself
            if (i == 0 && j == 0)
                continue;

            xp = x + j;
            yp = y + i;
            offset = yp * iWidth + xp;
            if ((xp < 0) || (xp >= iWidth))
                continue;
            if ((yp < 0) || (yp >= iHeight))
                continue;

            double lum = NL(pSrc->GetOffset(offset)[0]);

            // also check, if luminance difference is not greater than th
            if (std::abs(value - lum) > th)
                continue;

            retval += lum;
            count++;
        }
    }

    // if not valid pixel around, return itself as mean to prevent errors
    return count ? retval / count : value;
}

double TMOTai08::GetBetaAroundPixel(int x, int y, double value)
{
    int xp, yp, offset;
    double diff = 10e7;

    int r = 7;

    int bx = -1;
    int by = -1;

    for (int i = -r; i <= r; i++)
    {
        for (int j = -r; j <= r; j++)
        {
            // skip itself
            if (i == 0 && j == 0)
                continue;

            xp = x + j;
            yp = y + i;
            offset = yp * iWidth + xp;
            if ((xp < 0) || (xp >= iWidth))
                continue;
            if ((yp < 0) || (yp >= iHeight))
                continue;

            if (predictor.at<uint8_t>(yp, xp) == BETA_UNDEFINED)
                continue;

            double lum = NL(pSrc->GetOffset(offset)[0]);
            if (std::abs(value - lum) < diff)
            {
                bx = xp;
                by = yp;
                diff = std::abs(value - lum);
            }
        }
    }

    if (bx >= 0 && by >= 0)
        return beta.at<double>(by, bx);

    return -1.0;
}

double TMOTai08::NL(double value)
{
    // dynamic range normalization, this part is guessed, as it's not in the article
    return 2.0 * std::log10(value + 1.0) / gmax;
}
