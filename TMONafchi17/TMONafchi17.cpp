/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*             Author: Tomas Krivanek [xkriva29 AT stud.fit.vutbr.cz]           *
*                                    Brno 2024                                 *
*                                                                              *
*******************************************************************************/

#include "TMONafchi17.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMONafchi17::TMONafchi17()
{
	SetName(L"Nafchi17");					  // TODO - Insert operator name
	SetDescription(L"Global mapping approach that estimates the three linear weighting parameters for RGB based on correlation."); // TODO - Insert description

	r.SetName(L"r");				// TODO - Insert parameters names
	r.SetDescription(L"downsampling factor"); // TODO - Insert parameter descriptions
	r.SetDefault(256);							// TODO - Add default values
	r = 256;
	r.SetRange(1.0, 1024.0); // TODO - Add acceptable range if needed
	this->Register(r);

   imageType.SetName(L"invertContrastMap");
   imageType.SetDescription(L"Using inverted contrast map.");
   imageType.SetDefault(true);
   imageType = true;
   this->Register(imageType);
}

TMONafchi17::~TMONafchi17()
{
}


cv::Mat TMONafchi17::calculateMeanImage(cv::Mat &in01, int width, int height) {

   cv::Mat singleChannelImage;
   cv::transform(in01, singleChannelImage, cv::Matx13f(1,1,1)); // Summing 3 channels into one
   singleChannelImage /= in01.channels(); // Mean image
   cv::Mat meanImage(height, width, CV_64FC3);
   
   cv::merge(
      std::vector<cv::Mat>{
         singleChannelImage,
         singleChannelImage,
         singleChannelImage
      },
      meanImage); // Mean image
   return meanImage;
}

cv::Mat TMONafchi17::calculateStdDevImage(cv::Mat &in01, cv::Mat &meanImage, int width, int height){
   // Calculation of stddev image
   cv::Mat difference = cv::abs(in01 - meanImage);
   cv::Mat squaredDif;
   cv::pow(difference, 2, squaredDif);
   cv::Mat singleChannelImage;
   cv::transform(squaredDif, singleChannelImage, cv::Matx13f(1,1,1)); // Summing 3 channels into one
   singleChannelImage /= 2;
   cv::Mat sqrtImage(height, width, CV_64FC1);
   cv::sqrt(singleChannelImage,sqrtImage);

   cv::Mat result(height, width, CV_64FC3);
   cv::merge(
      std::vector<cv::Mat>{
         sqrtImage,
         sqrtImage,
         sqrtImage
      },
      result); // Mean image
   result /= 0.5774; // Constant is taken from original paper's source code
   
   return result;
}

std::tuple<double,double,double> TMONafchi17::calculatePearsonCoeff(cv::Mat &in01,cv::Mat &contrastMap, int width, int height){
   
   std::vector<cv::Mat> contrastMapChannels;
   cv::split(contrastMap, contrastMapChannels);
   cv::Mat contrastMapR = contrastMapChannels[0];
   auto m0 = cv::mean(contrastMapR)[0];
   
   /* Splitting the input image to channels*/
   std::vector<cv::Mat> channels;
   cv::split(in01, channels);

   auto mr = cv::mean(channels[0])[0];
   auto mg = cv::mean(channels[1])[0];
   auto mb = cv::mean(channels[2])[0];


   cv::Mat d1 = contrastMapR.reshape(1) - m0;
   cv::Mat d12;
   cv::pow(d1, 2, d12);

   cv::Mat dr = channels[0].reshape(1) - mr;
   cv::Mat dr2;
   cv::pow(dr, 2, dr2);   

   cv::Mat dg = channels[1].reshape(1) - mg;
   cv::Mat dg2;
   cv::pow(dg, 2, dg2);

   cv::Mat db = channels[2].reshape(1) - mb;
   cv::Mat db2;
   cv::pow(db, 2, db2);

   auto sumd1 = cv::sum(d12)[0];
   auto sumdr = cv::sum(dr2)[0];
   auto sumdg = cv::sum(dg2)[0];
   auto sumdb = cv::sum(db2)[0];

   auto Rho1 = sum(d1.mul(dr))[0] / std::pow(sumd1 * sumdr,0.5);
   auto Rho2 = sum(d1.mul(dg))[0] / std::pow(sumd1 * sumdg,0.5);
   auto Rho3 = sum(d1.mul(db))[0] / std::pow(sumd1 * sumdb,0.5); 

   return std::make_tuple(Rho1, Rho2, Rho3);
}

std::tuple<double,double,double> TMONafchi17::calculateLambda(cv::Mat &in01, std::tuple<double,double,double> rhos, int width, int height){
   auto rho1 = std::get<0>(rhos);
   auto rho2 = std::get<1>(rhos);
   auto rho3 = std::get<2>(rhos);
   cv::Mat matFromTuple(1,3,CV_64F);
   matFromTuple.at<double>(0,0) = rho1;
   matFromTuple.at<double>(0,1) = rho2;
   matFromTuple.at<double>(0,2) = rho3;
   auto minRho = std::min(rho1, std::min(rho2, rho3));
   auto maxRho = std::max(rho1, std::max(rho2, rho3));

   
   cv::Mat Gamma = ((matFromTuple - minRho) / (maxRho - minRho)) - 0.5;
   

   cv::Mat beta = cv::abs(matFromTuple);
   beta = beta / sum(beta)[0]; 
   cv::Mat lambda = beta + min(beta,Gamma);
   lambda = cv::abs(lambda);
   lambda = lambda / sum(lambda)[0];

   return std::make_tuple(lambda.at<double>(0,0), lambda.at<double>(0,1), lambda.at<double>(0,2));

}
bool TMONafchi17::haveImageGreaterValuesThen(cv::Mat &in01, double value){
   for(int i = 0; i < in01.rows; i++){
      for(int j = 0; j < in01.cols; j++){
         if(in01.at<cv::Vec3d>(i,j)[0] > value || in01.at<cv::Vec3d>(i,j)[1] > value || in01.at<cv::Vec3d>(i,j)[2] > value){
            return true;
         }
      }
   }
   return false;
}
cv::Mat TMONafchi17::convertToDouble(cv::Mat &in01){
   cv::Mat outputImage;
   in01.convertTo(outputImage, CV_64FC3, 1.0 / 255.0); // Convert to double
   return outputImage;

}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMONafchi17::Transform()
{
   pSrc->Convert(TMO_RGB,false);
   pDst->Convert(TMO_RGB,false); 
   
   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();
   
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   
   cv::Mat in01(height, width, CV_64FC3, pSourceData);
   
   if(haveImageGreaterValuesThen(in01, 1.0)){
      in01 = convertToDouble(in01);
   }

   cv::Mat Mu = calculateMeanImage(in01, width, height);
   cv::Mat Sigma = calculateStdDevImage(in01, Mu, width, height);

   cv::Mat Q = Mu.mul(Sigma);
   auto koefs1 = calculatePearsonCoeff(in01, Q, width, height);

   cv::Mat Q_minus1= Mu.mul(1-Sigma);
   auto koefs2 = calculatePearsonCoeff(in01, Q_minus1, width, height);

   double Lambda1a, Lambda2a, Lambda3a,Lambda1b, Lambda2b, Lambda3b;
   std::tie(Lambda1a, Lambda2a, Lambda3a) = calculateLambda(in01, koefs1, width, height);
   std::tie(Lambda1b, Lambda2b, Lambda3b) = calculateLambda(in01, koefs2, width, height);


   std::vector<cv::Mat> channels;
   cv::split(in01, channels);
   cv::Mat R = channels[0];
   cv::Mat G = channels[1];
   cv::Mat B = channels[2];

   cv::Mat result1 = (R.mul(Lambda1a) + G.mul(Lambda2a) + B.mul(Lambda3a));
   cv::merge(
      std::vector<cv::Mat>{
         result1,
         result1,
         result1,
      },
      result1
   );
   cv::Mat result2 = (R.mul(Lambda1b) + G.mul(Lambda2b) + B.mul(Lambda3b));
   cv::merge(
      std::vector<cv::Mat>{
         result2,
         result2,
         result2,
      },
      result2
   );

   if(imageType){
      memcpy(pDestinationData, result1.data, width*height*3*sizeof(double));
   } else {
      memcpy(pDestinationData, result2.data, width*height*3*sizeof(double));
   }
   
	return 0;
}
