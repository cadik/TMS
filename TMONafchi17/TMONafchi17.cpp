/* --------------------------------------------------------------------------- *
 * TMONafchi17.cpp: implementation of the TMONafchi17 class.   *
 * --------------------------------------------------------------------------- */
/***********************************************************************************
* Author: Bc. Tomas Krivanek (xkriva29)   xkriva29@stud.fit.vutbr.cz               *
************************************************************************************/
#include "TMONafchi17.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMONafchi17::TMONafchi17()
{
	SetName(L"Nafchi17");					  // TODO - Insert operator name
	SetDescription(L"Global mapping approach that estimates the three linear weighting parameters for RGB based on correlation."); // TODO - Insert description

	dParameter.SetName(L"r");				// TODO - Insert parameters names
	dParameter.SetDescription(L"downsampling factor"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(256);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(1.0, 1024.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMONafchi17::~TMONafchi17()
{
}


cv::Mat TMONafchi17::calculateMeanImage(cv::Mat &in01, int width, int height) {

   cv::Mat singleChannelImage;
   cv::transform(in01, singleChannelImage, cv::Matx13f(1,1,1));
   singleChannelImage /= in01.channels();
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
   cv::Mat difference = cv::abs(in01 - meanImage);
   cv::Mat squaredDif;
   cv::pow(difference, 2, squaredDif); // sqrt(squaredDif
   cv::Mat singleChannelImage;
   cv::transform(squaredDif, singleChannelImage, cv::Matx13f(1,1,1));
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
   result /= 0.5774;
   
   return result;
}
std::tuple<double,double,double> TMONafchi17::calculatePearsonCoeff(cv::Mat &in01,cv::Mat &contrastMap, int width, int height){
   // Calculating corelation between one dimension and contrast map
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
   cout << "Contrast Map at width/2 height/2: " << contrastMap.at<double>(height/2, width/2) << endl;
   
   std::cout << "m0[0]: " << m0 << std::endl;

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

   cout << "dr.size: " << dr.size() << " d1.size: " << d1.reshape(1).size() << endl;
   auto Rho1 = sum(d1.mul(dr))[0] / std::pow(sumd1 * sumdr,0.5);
   auto Rho2 = sum(d1.mul(dg))[0] / std::pow(sumd1 * sumdg,0.5);
   auto Rho3 = sum(d1.mul(db))[0] / std::pow(sumd1 * sumdb,0.5); 
   cout << "Rho1 " << Rho1 << endl;
   cout << "Rho2 " << Rho2 << endl;
   cout << "Rho3 " << Rho3 << endl;
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

   auto minRho = min(rho1, min(rho2, rho3));
   auto maxRho = max(rho1, max(rho2, rho3));

   cv::Mat Gamma = ((matFromTuple - minRho) / (maxRho - minRho)) - 0.5;
   cout<< "tady" << Gamma << endl;
   cv::Mat beta = cv::abs(Gamma);

   beta = beta / sum(beta)[0]; // Normalization

   cv::Mat lambda = beta + min(beta,Gamma);
   lambda = abs(lambda);
   lambda = lambda / sum(lambda); // Normalization

   return std::make_tuple(lambda.at<double>(0,0), lambda.at<double>(0,1), lambda.at<double>(0,2));

}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMONafchi17::Transform()
{
   double r = dParameter; 
   pSrc->Convert(TMO_RGB);
   pDst->Convert(TMO_RGB); 
   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();
 
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   
   // Convert source data to OpenCV Mat
   cv::Mat in01(height, width, CV_64FC3, pSourceData);
   std::vector<cv::Mat> channels;
   cv::split(in01, channels);
   cv::Mat R = channels[0];
   cv::Mat G = channels[1];
   cv::Mat B = channels[2];

   std::cout << "d1 size: " << in01.size() << std::endl;
   // cv::imshow("Source Image", in01);
   cv::Mat Mu = calculateMeanImage(in01, width, height);
   cv::Mat Sigma = calculateStdDevImage(in01, Mu, width, height);
   cv::Mat Q = Mu.mul(Sigma);

   double Rho1a, Rho2a, Rho3a;
   std::tie(Rho1a,Rho2a,Rho3a) = calculatePearsonCoeff(in01, Q, width, height);
   cv::Mat Q_minus1= 1 - Q;
   auto koefs2 = calculatePearsonCoeff(in01, Q_minus1, width, height);
   cout << "Koefs1: " << Rho1a << " " << Rho2a << " " << Rho3a << endl;
   // cout << "Koefs2: " << std::get<0>(koefs2) << " " << std::get<1>(koefs2) << " " << std::get<2>(koefs2) << endl;
   double Lambda1, Lambda2, Lambda3;
   std::tie(Lambda1, Lambda2, Lambda3) = calculateLambda(in01, koefs2, width, height);
   cout << "Lambda1: " << Lambda1 << " " << Lambda2 << " " << Lambda3 << endl;
   cv::Mat result = (R.mul(Lambda1) + G.mul(Lambda2) + B.mul(Lambda3));
   cv::merge(
      std::vector<cv::Mat>{
         result,
         result,
         result,
      },
      result
   );
   cout << "Original Size:" << in01.size() << endl;
   cout << "Result: " << result.size() << endl;
   cout << "Result: " << result.at<double>(20,20) << endl;

   memcpy(pDestinationData, result.data, width*height*3*sizeof(double));
   
	return 0;
}
