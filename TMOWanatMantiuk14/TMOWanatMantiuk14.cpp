/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2023-24                                           *
*                                                                              *
*                       Implementation of the TMOWanatMantiuk14 class          *
*                                                                              *
*                       Author: Lukas Macejka (xmacej03)                       *
*                       Mail: xmacej03@vutbr.cz                                *
*                                                                              *
*******************************************************************************/

#include "TMOWanatMantiuk14.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <eigen3/Eigen/Dense>// Pre prácu s vektormi a maticami
#include <stdlib.h>
#include <math.h> 
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
#include <nlopt.hpp>
#include "opencv2/opencv.hpp"
//#include "../../../../usr/local/MATLAB/R2023b/extern/include/MatlabEngine.hpp"
//#include "../../../../usr/local/MATLAB/R2023b/extern/include/MatlabDataArray.hpp"
//#include "Ipopt/src/Interfaces/IpTNLP.hpp"
//#define OPTIM_USE_RCPP_ARMADILLO
//#define OPTIM_ENABLE_ARMA_WRAPPERS
//#include "./optim/header_only_version/optim.hpp"
//#include <armadillo>
#include <eigen-3.4.0/unsupported/Eigen/NumericalDiff>
#include <eigen-3.4.0/unsupported/Eigen/NonLinearOptimization>


/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOWanatMantiuk14::TMOWanatMantiuk14()
{
	SetName(L"WanatMantiuk14");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOWanatMantiuk14::~TMOWanatMantiuk14()
{
}

//----------------------------CSF----------------------------------//

// linearna interpolacia ekvivalentna z MATLAB dat
double linearInterpolate(const std::vector<double>& x, const std::vector<double>& y, double x0) {
    for (size_t i = 0; i < x.size() - 1; i++) {
        if (x0 >= x[i] && x0 <= x[i + 1]) {
            double t = (x0 - x[i]) / (x[i + 1] - x[i]);
            return y[i] + t * (y[i + 1] - y[i]);
        }
    }
    return y.back(); 
}

// Custom-fit MTF oka
double hdrvdp_mtf(double rho, const std::vector<double>& mtf_params_a, const std::vector<double>& mtf_params_b) {
    double MTF = 0;
    for (int kk = 0; kk < 4; kk++) {
        MTF += mtf_params_a[kk] * exp(-mtf_params_b[kk] * rho);
    }
    return MTF;
}

// Join rod-cone sensitivity
double hdrvdp_joint_rod_cone_sens(double la, const std::vector<double>& csf_sa) {
    double cvi_sens_drop = csf_sa[1];
    double cvi_trans_slope = csf_sa[2];
    double cvi_low_slope = csf_sa[3];

    return csf_sa[0] * pow((cvi_sens_drop / la), cvi_trans_slope + 1) * -cvi_low_slope;
}

// Main function for csf
double csf_hdrvdp(double rho, double lum) {
    // Parametre 
    std::vector<std::vector<double>> csf_pars = {
        {0.0160737, 0.991265, 3.74038, 0.50722, 4.46044},
        {0.383873, 0.800889, 3.54104, 0.682505, 4.94958},
        {0.929301, 0.476505, 4.37453, 0.750315, 5.28678},
        {1.29776, 0.405782, 4.40602, 0.935314, 5.61425},
        {1.49222, 0.334278, 3.79542, 1.07327, 6.4635},
        {1.46213, 0.394533, 2.7755, 1.16577, 7.45665}
    };
    std::vector<double> lum_lut = {0.002, 0.02, 0.2, 2, 20, 150};
    double log_lum = log10(lum);
    std::vector<double> par = {0.061466549455263, 0.99727370023777070};
    std::vector<double> mtf_params_a = {par[1]*0.426, par[1]*0.574, (1-par[1])*par[0], (1-par[1])*(1-par[0])};
    std::vector<double> mtf_params_b = {0.028, 0.37, 37, 360};
    std::vector<double> csf_sa = {30.162, 4.0627, 1.6596, 0.2712};

    // Interpolacia
    std::vector<double> interp_pars(4);
    for (int k = 0; k < 4; k++) {
        interp_pars[k] = linearInterpolate(lum_lut, csf_pars[k], std::clamp(log_lum, lum_lut.front(), lum_lut.back()));
    }

    // Vypocet hodnoty S pouzivanim interpolacnych parametrov a funcii
    double S = (interp_pars[3] * 1.0 / (pow(1 + pow(interp_pars[0] * rho, interp_pars[1]), interp_pars[2]) * pow(1 - exp(-pow(rho / 7, 2)), interp_pars[2]))) * hdrvdp_mtf(rho, mtf_params_a, mtf_params_b) * hdrvdp_joint_rod_cone_sens(lum, csf_sa);

    return S;
}

//Get tone curve
//----------------------------GTC----------------------//


using namespace Eigen;

double csf_hdrvdp(double rho, double lum); 

double c_thr(double l) {
    return 1.0 / (8 * csf_hdrvdp(2, std::pow(10, l)));
}

// Pomocná funkcia pre výpočet diferencií vektoru
std::vector<double> diff(const std::vector<double>& vec, double delta) {
    std::vector<double> d(vec.size() - 1);
    for (int i = 0; i < vec.size() - 1; ++i) {
        d[i] = 0.1 * ((vec[i + 1] - vec[i])/delta - 1);
    }
    return d;
}
 
double objektivneFunkcia(const std::vector<double> &l_out, std::vector<double> &grad_out, void *my_func_data) {
//double ArmadilloOperator()(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data) const {
    std::cerr << "zdena"<< std::endl;
    std::vector<double> *l_in_h = static_cast<std::vector<double>*> (my_func_data);
    std::vector<double> &l_in = *l_in_h;
    std::cerr << "zdena"<< std::endl;
    double delta = l_in[1] - l_in[0];
    double c = 0.1;      std::cerr << "zdena"<< std::endl;
    std::vector<double> err = diff(l_out, delta);
    double f = 0.0;      std::cerr << "zdena"<< std::endl;
    for (int i = 0; i < err.size(); ++i) {
        err[i] += c_thr(l_in[i]) - c_thr(l_out[i]);

        f += (err[i] * err[i]);
    }        

    return f;
}

std::vector <double> refactor(VectorXd eigenVec){ 
    std::vector<double> stdVec;
    stdVec.resize(eigenVec.size()); // Zabezpečíme, že má správnu veľkosť

    // Kopírovanie prvkov
    for(int i = 0; i < eigenVec.size(); ++i) {
        stdVec[i] = eigenVec(i); // Priradenie každého prvku
    }
    return stdVec;
}

 VectorXd optimizeToneCurve(VectorXd& l_out, const VectorXd& l_in) {
    double delta = l_in[1] - l_in[0];
    
    // NLopt optimalizace
    nlopt::opt opt(nlopt::LN_COBYLA, l_out.size());

   // Nastavenie dolných a horných hraníc
    std::vector<double> lb(l_in.size(), -10);
    opt.set_lower_bounds(lb);
    std::vector<double> ub(l_in.size(), 10); 
    opt.set_upper_bounds(ub);

    std::vector<double> l_in_h = refactor(l_in);
    std::cerr << l_in_h[299] << std::endl;
    opt.set_min_objective(objektivneFunkcia, static_cast<void*> (& l_in_h));

    std::vector<double> x(l_out.size(), 1.0); // Prevod Eigen::VectorXd na std::vector
    double minf = 0.0;
    std::cerr << "Nalezené fsfrgsrg: " << minf << " lout size "<< l_out.size() << std::endl;
    std::cerr << "x " << x.size() << std::endl;
    opt.set_maxeval(100);

    // segfault - opt.optimize je chybne  

    nlopt::result vysledok = opt.optimize(x, minf);
std::cerr << "Nalezené fsfrgsrgeeee: " << minf << std::endl;
    try {
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "Nalezené minimum: " << minf << std::endl;
        // Aktualizujte l_out podle optimalizovaných výsledků
        l_out = Eigen::Map<VectorXd>(x.data(), x.size());
        return l_out;
    } catch (std::exception &e) {
        std::cerr << "NLopt exception: " << e.what() << std::endl;
        return l_out;
    }
}

/*struct ObjectiveFunction {
   typedef double Scalar;
    typedef Eigen::VectorXd InputType;
    typedef Eigen::VectorXd ValueType;
    typedef Eigen::MatrixXd JacobianType;

    enum {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
    };
   
    int operator()(const Eigen::VectorXd& l_in, const Eigen::VectorXd& l_out, ValueType& E) const {
      // Example objective function (Rosenbrock)
      
      double delta = l_in(1) - l_in(0);
     
      double c = 0.1;
      VectorXd err = c * (diff(l_out) / delta - VectorXd::Ones(l_out.size() - 1));
      double f = 0.0;
      for (int i = 0; i < err.size(); ++i) {
         err(i) += c_thr(l_in(i)) - c_thr(l_out(i));
         
         f += (err(i) * err(i));
      }   
      E = f; 

      return 0;
   }
};*/


double obmedzenie(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    // Telo obmedzenia
    return 0;//x[0] + x[1] + x[2] + x[3] - 10; // Príklad: súčet parametrov musí byť <= 10
}


/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOWanatMantiuk14::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py;

   //pSourceData->Convert(TMO_Yxy); 
   // double *pSourceDataLuminance = pSrc->GetData();

   //sRGB to RGB - gamma correction undone
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

//////////////////////////////////

   //CSF(32,33);

   double* lum = new double[numPixels];

   pSrc->Convert(TMO_Yxy); //RGB -> luminance
   
   Eigen::VectorXd l_in(numPixels);
   Eigen::VectorXd l_out(numPixels);

   std::cerr << "l_in vector" << l_in[55] << "velkost" << l_in.size() << std::endl;

   double sumLuminance = 0.0;
   // Average luminance for global contrast
   for (int i = 0; i < numPixels; i++)
   {
      double Y = pSrc->GetData()[i * 3];
      sumLuminance += Y;
   }
   

   double avgLuminance = sumLuminance / numPixels;
std::cerr << "aga" << numPixels << std::endl;
    for (int i = 0; i < numPixels; i++) {
        double Y = pSrc->GetData()[i * 3]; // Předpokládáme, že Y je první složka v Yxy formátu
        lum[i] =  std::log10(Y);
        l_in[i] = Y;
        l_out[i] = Y; // Inicializace výstupní luminance stejnou hodnotou jako vstupní
    }
std::cerr << "aga" << numPixels << std::endl;
    // Vypočítáme optimální tone curve
   l_out = optimizeToneCurve(l_out, l_in);
   

   double result = csf_hdrvdp(2, * lum);//32, std::pow(10, 33));
   std::cerr <<"fafaf" << result << std::endl;

    //optimizeToneCurve(l_out, l_in);

 std::cerr <<"@@@@@@@2332@@@@" << std::endl;

    // Aplikujeme výslednou tone curve na obrázek
    for (int i = 0; i < numPixels; i++) {
        pDst->GetData()[i * 3] = l_out[i]; // Aktualizace luminance v cílovém obrázku
        // px a py hodnoty zůstanou nezměněné
       // std::cerr <<"fafaf" << l_out[i] << std::endl;
    }

    // Převedeme cílový obrázek zpět do RGB, pokud je to potřeba
    pDst->Convert(TMO_RGB);


   int j = 0;

	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			pY *= dParameter; // Parameters can be used like
							  // simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}