/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                                                                                   *
*                       Author: Filip Sapak [xsapak05@stud.fit.vutbr.cz]            *
*                       Brno 2024                                                   *
*                                                                                   *
*                       Implementation of the TMOWu24 class                         *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOWu24.cpp
 * @brief Implementation of the TMOWu24 class
 * @author Filip Sapak
 * @class TMOWu24.cpp
 */

/* --------------------------------------------------------------------------- *
 * TMOWu24.cpp: implementation of the TMOWu24 class.                           *
 * --------------------------------------------------------------------------- */

#include "TMOWu24.h"
#include <numeric>
#include <opencv2/opencv.hpp>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOWu24::TMOWu24()
{
	SetName(L"Wu24");					  
	SetDescription(L"Efficient and Effective Image Decolorization Algorithm Based on Cumulative Distribution Function"); 

}

TMOWu24::~TMOWu24()
{
}


/* --------------------------------------------------------------------------- *
 * @brief This function normalizes weights                                            *
 * --------------------------------------------------------------------------- */
void normalizeWeights(std::vector<double>& weights) {

   double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
   for (double& w : weights) {
      w /= sum;
   }
}
   

/**
 * @brief Converts image
 * 
 * Source image is stored in local parameter pSrc
 * Destination image is in pDst
 * 
 */
int TMOWu24::Transform()
{

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();

   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();

   //convert source data to cv::Mat
   cv::Mat src (height, width, CV_64FC3, pSourceData);
   cv::Mat R,G,B;
   cv::Mat channels[3];

   //split into channels
   split(src, channels);
   R = channels[0];
   G = channels[1];
   B = channels[2];

   //histograms and CDFs
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

   //compute optimal CDF
   for (int i = 0; i < 256; ++i) {
      cdfOpt[i] = (cdfR[i] + cdfG[i] + cdfB[i]) / 3.0;
   }

   std::vector<double> weightR(256, 0), weightG(256, 0), weightB(256, 0);

   //compute weights
   for (int i = 0; i < 256; i++) {
      weightR[i] = exp(-abs(cdfR[i] - cdfOpt[i]));
      weightG[i] = exp(-abs(cdfG[i] - cdfOpt[i]));
      weightB[i] = exp(-abs(cdfB[i] - cdfOpt[i]));
   }

   //normalize weights
   normalizeWeights(weightR);
   normalizeWeights(weightG);
   normalizeWeights(weightB);

   //compute resulting image
   for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
            
            
            //compute sum of weights on pixel
            double sumWeights = weightR[static_cast<int>(R.at<double>(j, i))] + weightG[static_cast<int>(G.at<double>(j, i))] + weightB[static_cast<int>(B.at<double>(j, i))];

            //compute resulting channel value
            double red = (R.at<double>(j, i) * weightR[static_cast<int>(R.at<double>(j, i))]);
            double green = (G.at<double>(j, i) * weightG[static_cast<int>(G.at<double>(j, i))]);
            double blue = (B.at<double>(j, i) * weightB[static_cast<int>(B.at<double>(j, i))]);

            //sum and normalize
            double res = (red + green + blue) / sumWeights; 

            *pDestinationData++ = res; //red channel
            *pDestinationData++ = res; //green channel
            *pDestinationData++ = res; //blue channel
      }
   }

   pSrc->ProgressBar(height, height);
   return 0;
}


