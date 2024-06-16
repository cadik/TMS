/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2023-24                                           *
*                                                                              *
*                       Implementation of the TMOWanat14 class                 *
*                       Last part - color retargeting not implemented          *
*                                                                              *
*                       Author: Lukas Macejka (xmacej03)                       *
*                       Mail: xmacej03@vutbr.cz                                *
*                                                                              *
*******************************************************************************/
#include "TMOWanat14.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <stdlib.h>
#include <math.h> 
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include "opencv2/opencv.hpp"
#include <cppad/cppad.hpp>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOWanat14::TMOWanat14()
{
	SetName(L"Wanat14");					  
	SetDescription(L"Simulating and compensating changes in appearance between day and night vision"); 

	dParameter.SetName(L"ParameterName");	
	dParameter.SetDescription(L"ParameterDescription"); 
	dParameter.SetDefault(1);	
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); 
	this->Register(dParameter);
}

TMOWanat14::~TMOWanat14()
{
}



/*
 * ----------------------------CSF----------------------------------
 *
 * Converts TMOImage(input) to cv::Mat
 * Function from file TMOParis11 
 */
cv::Mat TMOImage2Mat(int width, int height, std::vector<double> input)
{
	int rowsCnt, colsCnt;

	rowsCnt = height;
	colsCnt = width;

	cv::Mat srcConvMat(rowsCnt, colsCnt, CV_64FC1);
	
   int i = 0;

	for (int y = 0; y < rowsCnt ; y++)
	{
		for (int x = 0; x <  colsCnt ; x++)
		{
         srcConvMat.at<double>(y,x) = input[i];
		   i++;
      }
	}
	return srcConvMat;
}

/**
 * @brief Converts cv::Mat back to TMOImage format.
 * @param width Width of the image.
 * @param height Height of the image.
 * @param input Input image as cv::Mat.
 * @return Converted image as cv::Mat.
 */
cv::Mat Mat2TMOImage(int width, int height,  cv::Mat input){

   int rowsCnt = height;
	int colsCnt = width;

   cv::Mat outh(width*height, 1, CV_64FC1);

   unsigned i=0;
	for (int y = 0; y < rowsCnt ; y++)
	{
		for (int x = 0; x <  colsCnt ; x++)
		{
         outh.at<double>(i,0) = input.at<double>(y,x);
		   i++;
      }
	}

   return outh;
}

/** 
 *  Linear Interpolation and Contrast Sensitivity Function (CSF)
 */
double singleLinearInterpolate( double x, int pos) {
	if(x<0.002) x=0.002;
	if(x>150)   x=150;
	if(x<0.02){
		double t = (x - 0.002) / (0.02 - 0.002);
		switch(pos){
			case 0:  return 0.0160737 + t * (0.383873 - 0.0160737);
			case 1:  return 0.991265 + t * (0.800889 - 0.991265);
			case 2:  return 3.74038 + t * (3.54104 - 3.74038);
			default: return 0.50722 + t * (0.682505 - 0.50722);
		}
            	
	}
	else
	if(x<0.2){
		double t = (x - 0.02) / (0.2 - 0.02);
		switch(pos){
			case 0:  return 0.383873 + t * (0.929301 - 0.383873);
			case 1:  return 0.800889 + t * (0.476505 - 0.800889);
			case 2:  return 3.54104 + t * (4.37453 - 3.54104);
			default: return 0.682505 + t * (0.750315 - 0.682505);
		}
	}
	else
	if(x<2){
		double t = (x - 0.2) / (2 - 0.2);
		switch(pos){
			case 0:  return 0.929301 + t * (1.29776 - 0.929301);
			case 1:  return 0.476505 + t * (0.405782 - 0.476505);
			case 2:  return 4.37453 + t * (4.40602 - 4.37453);
			default: return 0.750315 + t * (0.935314 - 0.750315);
		}
	}
	else
	if(x<20){
		double t = (x - 2) / (20 - 2);
		switch(pos){
			case 0:  return 1.29776 + t * (1.49222 - 1.29776);
			case 1:  return 0.405782 + t * (0.334278 - 0.405782);
			case 2:  return 4.40602 + t * (3.79542 - 4.40602);
			default: return 0.935314 + t * (1.07327 - 0.935314);
		}
	}
	else{
		double t = (x - 20) / (150 - 20);
		switch(pos){
			case 0:  return 1.49222 + t * (1.46213 - 1.49222);
			case 1:  return 0.334278 + t * (0.394533 - 0.334278);
			case 2:  return 3.79542 + t * (2.7755 - 3.79542);
			default: return 1.07327 + t * (1.16577 - 1.07327);
		}
	}
}

double hdrvdp_joint_rod_cone_sens(double la) {
   if ( la == 0) return 0;
   return 30.162 * pow((4.0627 / la), 1.6596 + 1) * 0.2712;
}

/** 
 *  Custom-fit MTF eye
 */
double hdrvdp_mtf(double rho,  std::vector<double>& mtf_params_a,  std::vector<double>& mtf_params_b) {
   double MTF = 0;
   for (int kk = 0; kk < 4; kk++) {
      MTF += mtf_params_a[kk] * exp(-mtf_params_b[kk] * rho);
   }
   return MTF;
}

/**
 * @brief Computes the Contrast Sensitivity Function (CSF) value for a given spatial frequency and luminance.
 * @param rho Spatial frequency.
 * @param lum Luminance value.
 * @return CSF value.
 */
double csf_hdrvdp(double rho, double lum) {

   double log_lum = log10(lum);
   if (lum <= 0) log_lum = 0;
   std::vector<double> par = {0.061466549455263, 0.99727370023777070};
   std::vector<double> mtf_params_a = {par[1]*0.426, par[1]*0.574, (1-par[1])*par[0], (1-par[1])*(1-par[0])};
   std::vector<double> mtf_params_b = {0.028, 0.37, 37, 360};
   std::vector<double> csf_sa = {30.162, 4.0627, 1.6596, 0.2712};

   std::vector<double> interp_pars(4);
   interp_pars[0] = singleLinearInterpolate( log_lum, 0);
   interp_pars[1] = singleLinearInterpolate( log_lum, 1);
   interp_pars[2] = singleLinearInterpolate( log_lum, 2);
   interp_pars[3] = singleLinearInterpolate( log_lum, 3);

   double S = (interp_pars[3] * 1.0 / (pow(1 + pow(interp_pars[0] * rho, interp_pars[1]), interp_pars[2]) * pow(1.0 - exp(-pow(rho / 7.0, 2.0)), interp_pars[2]))) * hdrvdp_mtf(rho, mtf_params_a, mtf_params_b) * hdrvdp_joint_rod_cone_sens(lum);

   return S;
}

double deriv_hdrvp_joint_rod_cone_sens(double x){
   if ( x == 0) return 0;
	return 905.253 * pow(1/x , 2.6396);
}

double deriv_linear_interpolate(int x,int pos){
	if(x<0.002){
		return 0;
	}
	else
	if(x<0.02){
		if(pos==0)return 0.3677993;
		else if(pos==1)return -1.90376;
		else if(pos==2) return -0.19934;
		else return 0.175285;
	}
	else
	if(x<0.2){
		if(pos == 0)return 0.545428;
		else if(pos==1)return -0.324384;
		else if(pos==2) return 0.83349;
		else return 0.06781;
	}
	else
	if(x<2){
		if(pos == 0)return 0.368459;
		else if(pos==1)return -0.070723;
		else if(pos==2) return 0.03149;
		else return 0.184999;
	}
	else
	if(x<20){
		if(pos==0)return 0.19446;
		else if(pos==1)return -0.071504;
		else if(pos==2) return -0.6106;
		else return 0.137956;
	}
	else
	if(x<150){
		if (pos==0)return -0.03009;
		else if(pos==1)return -0.066255;
		else if(pos==2) return -1.01992;
		else return 0.0925;
	}
	else{
		return 0;
	}
}

double deriv_csf_hdrvdp(double x){
	double f0 = singleLinearInterpolate(x,0);
	double f1 = singleLinearInterpolate(x,1);
	double f2 = singleLinearInterpolate(x,2);
	double f3 = singleLinearInterpolate(x,3);
	double f4 = hdrvdp_joint_rod_cone_sens(x);
   
	double f0_d = deriv_linear_interpolate(x,0);
	double f1_d = deriv_linear_interpolate(x,1);
	double f2_d = deriv_linear_interpolate(x,2);
	double f3_d = deriv_linear_interpolate(x,3);
	double f4_d = deriv_hdrvp_joint_rod_cone_sens(x);

	double ender = pow(pow(2,f1)*pow(f0,f1)+1,-f2);
	double bgn = 0.674818 * pow(0.0783896,f2) * f3 * f4 *
        (-( f2 * (pow(2,f1) * pow(f0,f1) *
                (f1 * f0_d / f0 + log(f0) * f1_d) + f1 * pow(f0,f1) * f1_d)
            / ender)
        - f2_d * log(ender)) * ender;
	double ln1 = 1.71813* pow(0.0783896,f2)*f3*f4*f2_d*ender;
	double ln2 = 0.674818*pow(0.0783896,f2)*f4*f3_d*ender;
	double ln3 = 0.674818*pow(0.0783896,f2)*f3*f4_d*ender;
	return bgn - ln1+ln2+ln3;
}

double deriv_c_thr(double l){
	double c = csf_hdrvdp(2, std::pow(10, l));
	double deriv_c = deriv_csf_hdrvdp(l);
	return -0.125*deriv_c/pow(c,2);
}

double single_diff(double x,double xplus1, double delta){
	return 0.1*((xplus1 - x)/((delta==1)?0:delta) -1);
}

double deriv_diff(double delta){
	return -0.1 / (delta==0)?1:delta;
}

double c_thr(double l) {
    return 1.0 / (8 * csf_hdrvdp(2, std::pow(10, l)));
}

std::vector<double> deriv_objektivna_funkcia( std::vector<double> &l_in, std::vector<double> &l_out){
	double delta = l_in[1] - l_in[0];
	std::vector<double> err(l_in.size() - 1);
	for (int i = 0; i < l_in.size()-1; ++i) {
		double diff = single_diff(l_out[0],l_out[1],delta);
      err[i] = deriv_c_thr(l_in[0])- 2*c_thr(diff)*deriv_diff(delta)*deriv_c_thr(diff);
	}
	return err;
}


/**
 *  linear interpolate ekvivalent from MATLAB dates
 */
double linearInterpolate( std::vector<double>& x,  std::vector<double>& y, double x0) {
   for (size_t i = 0; i < x.size() - 1; i++) {
      if (x0 >= x[i] && x0 <= x[i + 1]) {
         double t = (x0 - x[i]) / (x[i + 1] - x[i]);
         return y[i] + t * (y[i + 1] - y[i]);
      }
   }
   return y.back(); 
}

/** 
 *Get tone curve
 *----------------------------GTC----------------------
*/
using namespace Eigen;

double csf_hdrvdp(double rho, double lum); 

/**
 *  help function for calculate diferencial vector
  */
std::vector<double> diff( std::vector<double>& vec, double delta) {
   std::vector<double> d(vec.size() - 1);
   for (int i = 0; i < vec.size() - 1; ++i) {
      d[i] = 0.1 * ((vec[i + 1] - vec[i])/delta - 1);
   }
   return d;
}

/**
 * Error calculating funcion
  */
double objektivneFunkcia(std::vector<double> &l_in, std::vector<double> &l_out) {

   double delta = l_in[1] - l_in[0];
   double c = 0.1;  
   std::vector<double> err = diff(l_out, delta);
   double f = 0.0;  
   for (int i = 0; i < err.size(); ++i) {
      err[i] += c_thr(l_in[i]) - c_thr(l_out[i]);
      f += (err[i] * err[i]);
   }        
   return f;
}

std::vector <double> refactor(VectorXd eigenVec){ 
   std::vector<double> stdVec;
   stdVec.resize(eigenVec.size()); 

   for(int i = 0; i < eigenVec.size(); ++i) {
      stdVec[i] = eigenVec(i); 
   }
   return stdVec;
}
/**
 * Function to find gradient
 */
std::vector<double> GetGradiant(std::vector<double> &l_in, std::vector<double> &l_out)
{
   double h = 0.00001;
   std::vector<double> pred;
   for ( int i = 0; i < 3; i++)
   {
      std::vector<double> h_vector = l_out;
      h_vector[i] += h;
      
      pred.push_back((objektivneFunkcia(l_in, h_vector) - objektivneFunkcia(l_in, l_out)) / h);
   }
   return pred;
}

VectorXd optimizeToneCurve(VectorXd& l_out,  VectorXd& l_in) {

   std::vector<double> l_in_h = refactor(l_in);
   std::vector<double> l_out_h = refactor(l_out);
   
   double minf = 0.0;

   /** 
    * Compute the gradient for the initial values of x 
    * EDITABLE NUMBERS
    */
   int max_steps = 10; 
   double tolerance = 0.00001;
   double learn_rate = 1;

   for ( int i = 0; i < max_steps; i++)
   {
      std::vector<double> gradient = deriv_objektivna_funkcia(l_in_h, l_out_h);
      for (size_t j = 0; j < l_out_h.size(); ++j) 
      {
         l_out_h[j] -= learn_rate * gradient[j];
      }

      double sum_of_gradients = 0.0;

      for ( double grad : gradient) 
      {
         sum_of_gradients += std::abs(grad);
      }
      if (sum_of_gradients < tolerance) 
      {
         break;
      }
   }
   l_out = Eigen::Map<VectorXd>(l_out_h.data(), l_out_h.size());
   return l_out;
}

 double get_c_localized_broadband_contrast(int x, int y,double sigma,cv::Mat input, cv::Mat gaussian_kernel, cv::Mat g_times_l){
     
   double snd = input.at<double>(x,y) - g_times_l.at<double>(x,y);
   snd *= snd;
   cv::Mat copy = gaussian_kernel.clone();

   for (int i = 0; i < copy.rows; ++i) {
      for (int j = 0; j < copy.cols; ++j) {
         copy.at<double>(i,j) *= snd;
      }
   }
   if ( (copy.at<double>(x,y)) < 0) return 0;

   return std::sqrt(copy.at<double>(x,y));  
}

void get_log10_mat(cv::Mat input, cv::Mat retMat){
   for (int i = 0; i < retMat.rows; ++i) {
        for (int j = 0; j < retMat.cols; ++j) {
            retMat.at<double>(i,j) = std::log10(retMat.at<double>(i,j));
        }
  }
}

int Get_index(int position, int size){
   return size - position;
}

double get_G(double M){
  return std::log10((M+1)/(1-M))/2;
}

double get_rppd(){
  return 360 / (2 * M_PI * 20) / 1920;
}

double get_spatial_frequency(int k){
  return std::pow(2,-1*(k+1))*get_rppd();
}

double get_standard_deviation(int k){
  return 0.5*get_rppd()/get_spatial_frequency(k);
}

cv::Mat generateGaussianKernel(double sigma) {
   int size = sigma /2;
   /** 
    * Create an empty kernel matrix
   */
   
   int kernelSize = 5; 
   sigma = 1.0;

   /** 
    * Create the 1D Gaussian kernel
   */
   cv::Mat gaussianKernel1D = cv::getGaussianKernel(kernelSize, sigma, CV_64F);

   /** 
    * Create the 2D Gaussian kernel by taking the outer product of the 1D kernel with itself
   */

   cv::Mat gaussianKernel2D = gaussianKernel1D * gaussianKernel1D.t(); 
   return gaussianKernel2D;
}
   /** 
    * Main function to calculate Laplacian pyramid
   */
std::tuple<std::vector<cv::Mat>, std::vector<cv::Mat>> Get_laplacian(cv::Mat input)
{
   int levels = 3;
   std::vector<cv::Mat> gaussianPyramid;
   cv::Mat currentLevel = input;
   for (int i = 0; i < levels; ++i) {
      gaussianPyramid.push_back(currentLevel);
      cv::pyrDown(currentLevel, currentLevel);
   }

   /** 
    * Create Laplacian pyramid
   */
   std::vector<cv::Mat> laplacianPyramid;
   for (int i = 0; i < levels - 1; ++i) {
      cv::Mat expanded;
      cv::pyrUp(gaussianPyramid[i + 1], expanded, gaussianPyramid[i].size());
      cv::Mat laplacian;
      subtract(gaussianPyramid[i], expanded, laplacian);
      laplacianPyramid.push_back(laplacian);
   }
   laplacianPyramid.push_back(gaussianPyramid.back());
   return std::make_tuple(gaussianPyramid, laplacianPyramid); 
}

/** 
 * FUnction for calculating M
*/
double get_m(int k, int x, int y,cv::Mat input, cv::Mat kernel,cv::Mat source, cv::Mat log_mat, cv::Mat g_times_l, cv::Mat retar_gaussian, cv::Mat sorce_gaussian){
   double sigma = get_standard_deviation(k);
   double rho = get_spatial_frequency(k);
   double ret;

   double c_localized_broadband_contrast = get_c_localized_broadband_contrast(x,y,sigma,log_mat,kernel, g_times_l);
  
   double M = 1.0 / (8.6 * csf_hdrvdp(rho, sorce_gaussian.at<double>(x,y))); // 
   double Mroof =  1.0 / (8.6 * csf_hdrvdp(rho, retar_gaussian.at<double>(x,y)));

   if ( c_localized_broadband_contrast < 0.000001 )
      c_localized_broadband_contrast = 0.000001;

   ret = ( c_localized_broadband_contrast  - get_G(M) + get_G(Mroof) ) / c_localized_broadband_contrast;

   if ( std::isnan(ret)) 
      ret = 1;

   return ret;
}
//Function is to calculate for every floor in pyramid
cv::Mat edit_floor(int k, cv::Mat retar_floor, cv::Mat source_floor, cv::Mat log_mat, cv::Mat retar_gaussian, cv::Mat sorce_gaussian){

   cv::Mat retMat(source_floor.rows, source_floor.cols, CV_64F, cv::Scalar(0));
   double sigma = get_standard_deviation(k);

   cv::Mat gaussian_kernel =generateGaussianKernel(sigma);
   cv::Mat g_times_l;
   cv::filter2D(log_mat, g_times_l, -1, gaussian_kernel , cv::Point(-1, -1), 0, cv::BORDER_DEFAULT);
   for (int i = 0; i < retar_floor.rows; ++i) {
      for (int j = 0; j < retar_floor.cols; ++j) {
         retMat.at<double>(i,j) = source_floor.at<double>(i,j) * get_m(k, i, j, retar_floor,gaussian_kernel,source_floor, log_mat, g_times_l, retar_gaussian, sorce_gaussian);
      }
      
   }
   return retMat;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOWanat14::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		
	double *pDestinationData = pDst->GetData(); 
												
	double pY, px, py;

   int numPixels = pSrc->GetWidth() * pSrc->GetHeight();
   
    for (int i = 0; i < numPixels; i++) {
        for (int channel = 0; channel < 3; channel++) { 
            double& val = pSourceData[i * 3 + channel]; 
            if (val <= 0.04045) {
                val = val / 12.92;
            } else {
                val = pow((val + 0.055) / 1.055, 2.4);
            }
        }
    }
    
     pSrc->Convert(TMO_Yxy);

//////////////////////////////////

   std::vector<double> lum;   
   Eigen::VectorXd l_in(numPixels);
   Eigen::VectorXd l_out(numPixels);

   double sumLuminance = 0.0;

   for (int i = 0; i < numPixels; i++)
   {
      double Y = pSrc->GetData()[i * 3];
      sumLuminance += Y;
   }

   double avgLuminance = sumLuminance / numPixels;

   for (int i = 0; i < numPixels; i++) {
      double Y = pSrc->GetData()[i * 3]; 
      lum.push_back(std::log10(Y));
      l_in[i] = Y;
      l_out[i] = Y; 
   }

   l_out = optimizeToneCurve(l_out, l_in);

   std::vector<double> source_vec(l_in.data(), l_in.data() + l_in.size());
   std::vector<double> retarget_vec(l_out.data(), l_out.data() + l_out.size());

   cv::Mat source_mat = TMOImage2Mat(pSrc->GetWidth(), pSrc->GetHeight(), source_vec);
   cv::Mat retarget_mat = TMOImage2Mat(pSrc->GetWidth(), pSrc->GetHeight(), retarget_vec);

   
   auto sorce_tuple_pyramid = Get_laplacian(source_mat);
   auto retarget_tuple_pyramid = Get_laplacian(retarget_mat);

   std::vector<cv::Mat> sorce_gaussian = std::get<0>(sorce_tuple_pyramid);
   std::vector<cv::Mat> sorce_laplacian = std::get<1>(sorce_tuple_pyramid);
   std::vector<cv::Mat> retard_gaussian = std::get<0>(retarget_tuple_pyramid);
   std::vector<cv::Mat> retard_laplacian = std::get<1>(retarget_tuple_pyramid);
   int levels = 3;
   for ( int i = 0; i < levels ; i++)
   {
      int index = Get_index(i, levels);
      if (index == 1) continue;
      cv::Mat log_mat(sorce_laplacian[i].rows, sorce_laplacian[i].cols, CV_64F, cv::Scalar(0));
       for (int k = 0; k < log_mat.rows; ++k) {
        for (int j = 0; j < log_mat.cols; ++j) {
            if (sorce_laplacian[i].at<double>(k,j) == 0 ) log_mat.at<double>(k,j) = 0;
            else {
               if (sorce_laplacian[i].at<double>(k,j) < 0 ) log_mat.at<double>(k,j) = -1 * std::log(-1 * sorce_laplacian[i].at<double>(k,j));
               else log_mat.at<double>(k,j) = std::log(sorce_laplacian[i].at<double>(k,j));
            } 
        }
        }
      cv::Mat edited_mat = edit_floor(index, retard_laplacian[i], sorce_laplacian[i], log_mat, retard_gaussian[i], sorce_gaussian[i]);
      retard_laplacian[i] = edited_mat;
   }
    
   cv::Mat reconstructedImage = retard_laplacian.back(); 

   for (int i = retard_laplacian.size() - 2; i >= 0; --i) {
      cv::Mat expanded;
      pyrUp(reconstructedImage, expanded, retard_laplacian[i].size());
      add(expanded, retard_laplacian[i], reconstructedImage);

   }

   cv::Mat xy = Mat2TMOImage(reconstructedImage.cols, reconstructedImage.rows, reconstructedImage);

	pDst->Convert(TMO_Yxy);

   int j,k = 0;

	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); 
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			
			*pDestinationData++ = xy.at<double>(k);
            k++;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}