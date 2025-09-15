/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       VYF -- term project                                    *
*                       Author: Branislav Dubec                                *
*                       Brno 2024                                              *
*                                                                              *
*                       Image Decolorization Algorithm Based                   *
*                       on Cumulative  Distribution Function                   *
*                       https://www.mdpi.com/2313-433X/10/3/51                 *                                       
********************************************************************************/


/**
 * @file TMOWu24_2.cpp
 * @brief Image Decolorization Algorithm Based on Cumulative  Distribution Function
 * @author Branislav Dubec
 * @class TMOWu24_2
 */

#include "TMOWu24_2.h"

/**
  *  @brief Constructor
  */
TMOWu24_2::TMOWu24_2()
{
	SetName(L"Wu24_2");
	SetDescription(L"Image Decolorization Algorithm" 
    L"Based on Cumulative Distribution Function"); 

}

/**
  *  @brief Destructor
  */

TMOWu24_2::~TMOWu24_2()
{
}


/**
 * @brief transformation function
 * @return exit code
 */
int TMOWu24_2::Transform()
{


// Source image is stored in local parameter pSrc
// Destination image is in pDst

double *pDestinationData = pDst->GetData();
double *pSourceData = pSrc->GetData();




int histSize = 256;

// declare variables to store original value for each color channel 
cv::Mat srcB(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3);
cv::Mat srcG(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3);
cv::Mat srcR(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3);

// declare variables to store histograms for each color channel
std::vector<int> histB(256, 0), histG(256, 0), histR(256, 0);

// calculate histograms for each channel and store original values of 
// image for each channel
for (int i = 0; i < pSrc->GetHeight(); ++i) {
    for (int j = 0; j < pSrc->GetWidth(); ++j) {
        float r,g,b;
        b = *pSourceData++;
        g = *pSourceData++;
        r = *pSourceData++;
        srcB.at<double>(i,j) = b;
        srcG.at<double>(i,j) = g;
        srcR.at<double>(i,j) = r;
        // set minimal and maximal value, 0 and 255
        histB[std::clamp(static_cast<int>(b*255), 0, 256 - 1)]++;
        histG[std::clamp(static_cast<int>(g*255), 0, 256 - 1)]++;
        histR[std::clamp(static_cast<int>(r*255), 0, 256 - 1)]++;
    }
}

// normalize histograms into <0,1> by dividing every value 
// with sum of its values

// calculate sum for each histogram
float sumB = std::accumulate(histB.begin(), histB.end(), 0);
float sumG = std::accumulate(histG.begin(), histG.end(), 0);
float sumR = std::accumulate(histR.begin(), histR.end(), 0);

// declare variables to store normalized histograms
std::vector<float> histBNorm(256, 0), histGNorm(256, 0), histRNorm(256, 0);

// compute normalized histograms
for (int i = 0; i < histSize; i++) {
    histBNorm[i] = histB[i] / sumB;
    histGNorm[i] = histG[i] / sumG;
    histRNorm[i] = histR[i] /sumR;
}



// calculate cumulative distribution function for each channel
std::vector<float> cdfB(histSize, 0.0f), cdfG(histSize, 0.0f), cdfR(histSize, 0.0f);
cdfB[0] = histBNorm[0];
cdfG[0] = histGNorm[0];
cdfR[0] = histRNorm[0];
for (int i = 1; i < histSize; i++) {
    cdfB[i] = cdfB[i-1] + histBNorm[i];
    cdfG[i] = cdfG[i-1] + histGNorm[i];
    cdfR[i] = cdfR[i-1] + histRNorm[i];
}

// calculate optimal cumulative distribution function
std::vector<float> cdfOpt(histSize, 0.0f);
for (int i = 0; i < histSize; i++) {
    cdfOpt[i] = (cdfB[i] + cdfG[i] + cdfR[i]) / 3;
}


// calculate weights for each color channel for its intensity
std::vector<float> weightB(histSize, 0.0f), weightG(histSize, 0.0f), weightR(histSize, 0.0f);
for (int i = 0; i < histSize; i++) {
    weightB[i] = exp(-fabs(cdfB[i] - cdfOpt[i]));
    weightG[i] = exp(-fabs(cdfG[i] - cdfOpt[i]));
    weightR[i] = exp(-fabs(cdfR[i] - cdfOpt[i]));
}

// normalize weights
sumB = std::accumulate(weightB.begin(), weightB.end(), 0.0f);
sumG = std::accumulate(weightG.begin(), weightG.end(), 0.0f);
sumR = std::accumulate(weightR.begin(), weightR.end(), 0.0f);
for (int i = 0; i < histSize; i++) {
    weightB[i] /= sumB;
    weightG[i] /= sumG;
    weightR[i] /= sumR;
}


// calculate output picture
// get original value for each channel
// get pixel intensity for every channel
// multiply original channel value with the weight of the particular intensity
// output is sum of all three channel values normalized by dividing the sum of weights
for (int i = 0; i < pSrc->GetHeight(); ++i) {
    for (int j = 0; j < pSrc->GetWidth(); ++j) {
        

        // get original value for each channel
        float blueValue = srcB.at<double>(i,j);
        float greenValue = srcG.at<double>(i,j);
        float redValue = srcR.at<double>(i,j);

        // get pixel intensity for each channel
        int blueIdx = std::clamp(static_cast<int>(blueValue*255), 0, histSize - 1);
        int greenIdx = std::clamp(static_cast<int>(greenValue*255), 0, histSize - 1);
        int redIdx = std::clamp(static_cast<int>(redValue*255), 0, histSize - 1);


        float weightSum = weightB[blueIdx] + weightG[greenIdx] + weightR[redIdx];
        // each channel value is multiplied with its weight on its intensity
        float blueWeighted = weightB[blueIdx] * blueValue;
        float greenWeighted = weightG[greenIdx] * greenValue;
        float redWeighted = weightR[redIdx] * redValue;

        // normalize output
        float output = (blueWeighted + greenWeighted + redWeighted) / weightSum;
        // for each channel is the same output
        *pDestinationData++ = output;
        *pDestinationData++ = output;
        *pDestinationData++ = output;
    }
}

pSrc->ProgressBar(pSrc->GetHeight(), pSrc->GetWidth());

return 0;

}
