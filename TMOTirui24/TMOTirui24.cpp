/* --------------------------------------------------------------------------- *
 * TMOTirui24.cpp: implementation of the TMOTirui24 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOTirui24.h"
/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOTirui24::TMOTirui24()
{
	SetName(L"Tirui24");					  // TODO - Insert operator name
	SetDescription(L"Image Decolorization Algorithm Based on Cumulative Distribution Function"); // TODO - Insert description


	//this->Register(dParameter);
}

TMOTirui24::~TMOTirui24()
{
}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOTirui24::Transform()
{


// Source image is stored in local parameter pSrc
// Destination image is in pDst

//pSrc->Convert(TMO_RGB);
//pDst->Convert(TMO_RGB);
double *pDestinationData = pDst->GetData();
double *pSourceData = pSrc->GetData();

int histSize = 256;
float range[] = {0.f, 1.f};
const float* histRange[] = {range};
bool uniform = true, accumulate = false;
cv::Mat srcB(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3);
cv::Mat srcG(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3);
cv::Mat srcR(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3);

std::vector<int> histB(256, 0), histG(256, 0), histR(256, 0);
for (int i = 0; i < pSrc->GetHeight(); ++i) {
    for (int j = 0; j < pSrc->GetWidth(); ++j) {
        float r,g,b;
        b = *pSourceData++;
        g = *pSourceData++;
        r = *pSourceData++;
        srcB.at<double>(i,j) = b;
        srcG.at<double>(i,j) = g;
        srcR.at<double>(i,j) = r;
        histB[std::clamp(static_cast<int>(b*255), 0, 256 - 1)]++;
        histG[std::clamp(static_cast<int>(g*255), 0, 256 - 1)]++;
        histR[std::clamp(static_cast<int>(r*255), 0, 256 - 1)]++;
    }
}

float sumB = std::accumulate(histB.begin(), histB.end(), 0);
float sumG = std::accumulate(histG.begin(), histG.end(), 0);
float sumR = std::accumulate(histR.begin(), histR.end(), 0);
std::vector<float> histBNorm(256, 0), histGNorm(256, 0), histRNorm(256, 0);
for (int i = 0; i < histSize; i++) {
    histBNorm[i] = histB[i] / sumB;
    histGNorm[i] = histG[i] / sumG;
    histRNorm[i] = histR[i] /sumR;
}



std::vector<float> cdfB(histSize, 0.0f), cdfG(histSize, 0.0f), cdfR(histSize, 0.0f);
cdfB[0] = histBNorm[0];
cdfG[0] = histGNorm[0];
cdfR[0] = histRNorm[0];

for (int i = 1; i < histSize; i++) {
    cdfB[i] = cdfB[i-1] + histBNorm[i];
    cdfG[i] = cdfG[i-1] + histGNorm[i];
    cdfR[i] = cdfR[i-1] + histRNorm[i];
}

std::vector<float> cdfOpt(histSize, 0.0f);
for (int i = 0; i < histSize; i++) {
    cdfOpt[i] = (cdfB[i] + cdfG[i] + cdfR[i]) / 3;
}

std::vector<float> weightB(histSize, 0.0f), weightG(histSize, 0.0f), weightR(histSize, 0.0f);
for (int i = 0; i < histSize; i++) {
    weightB[i] = exp(-fabs(cdfB[i] - cdfOpt[i]));
    weightG[i] = exp(-fabs(cdfG[i] - cdfOpt[i]));
    weightR[i] = exp(-fabs(cdfR[i] - cdfOpt[i]));
}

 sumB = std::accumulate(weightB.begin(), weightB.end(), 0.0f);
 sumG = std::accumulate(weightG.begin(), weightG.end(), 0.0f);
sumR = std::accumulate(weightR.begin(), weightR.end(), 0.0f);
for (int i = 0; i < histSize; i++) {
    weightB[i] /= sumB;
    weightG[i] /= sumG;
    weightR[i] /= sumR;
}

for (int i = 0; i < pSrc->GetHeight(); ++i) {
    for (int j = 0; j < pSrc->GetWidth(); ++j) {

        float blueValue = srcB.at<double>(i,j);
        float greenValue = srcG.at<double>(i,j);
        float redValue = srcR.at<double>(i,j);

        int blueIdx = std::clamp(static_cast<int>(blueValue*255), 0, histSize - 1);
        int greenIdx = std::clamp(static_cast<int>(greenValue*255), 0, histSize - 1);
        int redIdx = std::clamp(static_cast<int>(redValue*255), 0, histSize - 1);

        float weightSum = weightB[blueIdx] + weightG[greenIdx] + weightR[redIdx];

        float blueWeighted = weightB[blueIdx] * blueValue;
        float greenWeighted = weightG[greenIdx] * greenValue;
        float redWeighted = weightR[redIdx] * redValue;

        float output = (blueWeighted + greenWeighted + redWeighted) / weightSum;

        *pDestinationData++ = output;
        *pDestinationData++ = output;
        *pDestinationData++ = output;
    }
}

pSrc->ProgressBar(pSrc->GetHeight(), pSrc->GetWidth());

return 0;

}
