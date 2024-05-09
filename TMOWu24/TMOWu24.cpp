/* --------------------------------------------------------------------------- *
 * TMOWu24.cpp: implementation of the TMOWu24 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOWu24.h"
#include <numeric>
#include <opencv2/opencv.hpp>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOWu24::TMOWu24()
{
	SetName(L"Wu24");					  // TODO - Insert operator name
	SetDescription(L"Efficient and Effective Image Decolorization Algorithm Based on Cumulative Distribution Function"); // TODO - Insert description

}

TMOWu24::~TMOWu24()
{
}

void normalizeWeights(std::vector<double>& weights) {

   double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
   for (double& w : weights) {
      w /= sum;
   }
}
   

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOWu24::Transform()
{

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();

   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();

   cv::Mat src (height, width, CV_64FC3, pSourceData);
   cv::Mat R,G,B;
   cv::Mat channels[3];
   split(src, channels);
   R = channels[0];
   G = channels[1];
   B = channels[2];

   std::vector<int> histR(256, 0), histG(256, 0), histB(256, 0);
   std::vector<double> cdfR(256), cdfG(256), cdfB(256);

   //compute histograms
   for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
         histR[static_cast<int>(R.at<double>(j, i))]++;
         histG[static_cast<int>(G.at<double>(j, i))]++;
         histB[static_cast<int>(B.at<double>(j, i))]++;
      }
   }

   //compute CDFs
   std::partial_sum(histR.begin(), histR.end(), cdfR.begin());
   std::partial_sum(histG.begin(), histG.end(), cdfG.begin());
   std::partial_sum(histB.begin(), histB.end(), cdfB.begin());

   //normalize CDFs
   for (int i = 0; i < 256; ++i) {
      cdfR[i] /= cdfR.back();
      cdfG[i] /= cdfG.back();
      cdfB[i] /= cdfB.back();
   }

   std::vector<double> cdfOpt(256,0);

   for (int i = 0; i < 256; ++i) {
      cdfOpt[i] = (cdfR[i] + cdfG[i] + cdfB[i]) / 3.0;
   }

   std::vector<double> weightR(256, 0), weightG(256, 0), weightB(256, 0);
   for (int i = 0; i < 256; i++) {
      weightR[i] = exp(-abs(cdfR[i] - cdfOpt[i]));
      weightG[i] = exp(-abs(cdfG[i] - cdfOpt[i]));
      weightB[i] = exp(-abs(cdfB[i] - cdfOpt[i]));
   }

   normalizeWeights(weightR);
   normalizeWeights(weightG);
   normalizeWeights(weightB);

   for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
            
            //sum the weights
            double sumWeights = weightR[static_cast<int>(R.at<double>(j, i))] + weightG[static_cast<int>(G.at<double>(j, i))] + weightB[static_cast<int>(B.at<double>(j, i))];

            double red = (R.at<double>(j, i) * weightR[static_cast<int>(R.at<double>(j, i))]);
            double green = (G.at<double>(j, i) * weightG[static_cast<int>(G.at<double>(j, i))]);
            double blue = (B.at<double>(j, i) * weightB[static_cast<int>(B.at<double>(j, i))]);

            double res = (red + green + blue) / sumWeights; 

            *pDestinationData++ = res; //red channel
            *pDestinationData++ = res; //green channel
            *pDestinationData++ = res; //blue channel
      }
   }

   pSrc->ProgressBar(height, height);
   return 0;
}


