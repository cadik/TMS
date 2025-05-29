/********************************************************************************
*                                                                               *
*                         Brno University of Technology                         *
*                       Faculty of Information Technology                       *
*                                                                               *
*                         Color-to-Grayscale Conversions                        *
*                                                                               *
*             Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]          *
*                                    Brno 2025                                  *
*                                                                               *
*                     Implementation of the TMOZhang08 class                    *
*           A Kernel Based Algorithm for Fast Color-To-Gray Processing          *
*                      https://doi.org/10.1109/CISP.2008.411                    *
*                                                                               *
********************************************************************************/

#include "TMOZhang08.h"
#include <cmath>
#include <fstream>


using namespace std;
using namespace Eigen;


/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOZhang08::TMOZhang08()
{
	SetName(L"Zhang08");
	SetDescription(L"A Kernel Based Algorithm for Fast Color-To-Gray Processing"); 

   HDRParameter.SetName(L"HDR");
	HDRParameter.SetDescription(L"is input image HDR");
	HDRParameter.SetDefault(false);
	HDRParameter = false;

   this->Register(HDRParameter);
}

TMOZhang08::~TMOZhang08()
{
}

/**
 * Stretches LAB values to the full range [0,1] for L and [-1,1] for A, B
 */
void TMOZhang08::stretchLabToFullRange(double* data, int width, int height) {
   double minL = data[0], maxL = data[0];
   double minA = data[1], maxA = data[1];
   double minB = data[2], maxB = data[2];

   // Find min/max values in chanels
   for (int i = 0; i < width * height; ++i) {
      double L = data[i * 3 + 0];
      double a = data[i * 3 + 1];
      double b = data[i * 3 + 2];

      minL = min(minL, L); maxL = max(maxL, L);
      minA = min(minA, a); maxA = max(maxA, a);
      minB = min(minB, b); maxB = max(maxB, b);
   }

   // Linear scaling
   for (int i = 0; i < width * height; ++i) {
       data[i * 3] = (data[i * 3 ] - minL) / (maxL - minL); 
       data[i * 3 + 1] = (data[i * 3 + 1] - minA) / (maxA - minA);
       data[i * 3 + 2] = (data[i * 3 + 2] - minB) / (maxB - minB);
   }
}

/**
 * Gamma correction for sRGB → linear RGB
 */
double TMOZhang08::gammaCorrection(double c) 
{
   return (c <= 0.04045) ? (c / 12.92) : pow((c + 0.055) / 1.055, 2.4);
}

/**
 * Helper function for LAB transformation
 */
double TMOZhang08::fLab(double t) 
{
   return (t > 0.008856) ? pow(t, 1.0 / 3.0) : (7.787 * t + 16.0 / 116.0);
}

/**
 * Converts sRGB to XYZ (D65)
 */
void TMOZhang08::RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z) {
   // sRGB -> linear RGB
   R = gammaCorrection(R);
   G = gammaCorrection(G);
   B = gammaCorrection(B);

   // Matrix fot computing sRGB -> XYZ (D65 reference white)
   X = R * 0.4124564 + G * 0.3575761 + B * 0.1804375;
   Y = R * 0.2126729 + G * 0.7151522 + B * 0.0721750;
   Z = R * 0.0193339 + G * 0.1191920 + B * 0.9503041;
}

/**
 * Converts XYZ to LAB into [0, 1] range
 */
void TMOZhang08::XYZtoLAB(double X, double Y, double Z, double &L, double &a, double &b) {
   // Reference white point D65
   const double Xn = 0.95047, Yn = 1.00000, Zn = 1.08883;

   double fx = fLab(X / Xn);
   double fy = fLab(Y / Yn);
   double fz = fLab(Z / Zn);

   L = (116.0 * fy) - 16.0;
   a = 500.0 * (fx - fy);
   b = 200.0 * (fy - fz);

   // Normalize into [0,1]
   L /= 100.0;
   a = (a / 128.0 + 1) / 2;
   b = (b / 128.0 + 1) / 2;
}

/**
 * Converts an RGB image to LAB - output LAB is in range [0,1]
 */
void TMOZhang08::convertRGBtoLAB(double* data, int width, int height) {

   bool range0to1 = isInRange0to1(data, width * height);

   for (int i = 0; i < width * height; ++i) {
       double R = data[i * 3 + 0];
       double G = data[i * 3 + 1];
       double B = data[i * 3 + 2];
       
      // If format is in range 0-255
      if (!range0to1 && !HDRParameter)
      {
         R /= 255;
         G /= 255;
         B /= 255;
      }

       double X, Y, Z, L, a, b;
       RGBtoXYZ(R, G, B, X, Y, Z);
       XYZtoLAB(X, Y, Z, L, a, b);

       data[i * 3 + 0] = L;
       data[i * 3 + 1] = a;
       data[i * 3 + 2] = b;
   }
   stretchLabToFullRange(data, width, height);
}


/**
 * Kernel function K5
 */
double TMOZhang08::kernelK5(const VectorXd& yi, const VectorXd& yj) 
{
   double a = 0.1, b = 1;
   double size = yi.size();
   double dot_ij = yi.dot(yj);  // Scalar product y_i * y_j
   double dot_ii = yi.dot(yi);  // Scalar product y_i * y_i
   double dot_jj = yj.dot(yj);  // Scalar product y_j * y_j

   double numerator = pow(a * dot_ij + b, 2);
   double denominator = (a * dot_ii + b) * (a * dot_jj + b);

   return numerator / denominator;
}

/**
 * Function for better kernel K6
 */
double TMOZhang08::kernelK6(const VectorXd& yi, const VectorXd& yj, int i, int j) 
{
   double K5_val = kernelK5(yi, yj);
   auto logK5 = log10(1+K5_val);

   if (i != 0 && j != 0) {
      return logK5 * logK5;
   } else {
      return logK5;
   }
}

/**
 * Maps projection values to luminance range (0-1)
 */
vector<double> TMOZhang08::normalizeToLuminance(const VectorXd& projections) 
{
   double minP = projections.minCoeff();
   double maxP = projections.maxCoeff();
   vector<double> luminance(projections.size());

   double range = maxP - minP;
   if (range == 0) {
       // If all values are the same, set everything to 0
       std::fill(luminance.begin(), luminance.end(), 0.0);
       return luminance;
   }

   
   for (size_t i = 0; i < projections.size(); i++) {
      luminance[i] = ((projections[i] - minP) / (range));
   }
   return luminance;
}

// Finds if range is 0-1 or in 0-255
bool TMOZhang08::isInRange0to1(double *pSourceData, int numPix)
{
   for (int i = 0; i < numPix * 3; i++)
   {
      if(pSourceData[i] > 1)
         return false;
   }
   return true;
}

/* --------------------------------------------------------------------------- *
 * This is the main funcion of the Zhang08                                     *
 * --------------------------------------------------------------------------- */

int TMOZhang08::Transform()
{
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int numPixels = width * height;

   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();

   convertRGBtoLAB(pSourceData, width, height);

   // Matrix A (3 × numPixels) 
   MatrixXd A(3, numPixels);
   Matrix3d K;

   // Filling matrix A
   int index = 0;
   for (int j = 0; j < height; j++) 
   {
      for (int i = 0; i < width; i++) 
      {
         (A)(0, index) = *pSourceData++;
         (A)(1, index) = *pSourceData++;   
         (A)(2, index) = *pSourceData++;
         index++;
      }
   }

   // Filling natrix K(3x3)
   for(int i = 0; i < 3 ; i++)
   {
      for(int j = 0; j < 3; j++)
      {        
         K(i, j) = kernelK6(A.row(i), A.row(j), i, j);
      }
   }

   // Matrix I_M with all elements 1/m 
   Matrix3d I_M = Matrix3d::Constant(1.0 / numPixels);

   // Compute K'
   Matrix3d K_prime = K - I_M * K - K * I_M + I_M * K * I_M;


   // Compute eigenvalues and eigenvectors of K'
   SelfAdjointEigenSolver<MatrixXd> eigensolver(K_prime);
   if (eigensolver.info() != Success) {
      std::cerr << "Eigen decomposition failed!" << std::endl;
      return -1;
   }

   // Find first eigenvector
   auto eigenvalues = eigensolver.eigenvalues();

   int maxIndex = 0;
   for (int i = 1; i < eigenvalues.size(); i++) {
      if (abs(eigenvalues[i]) > abs(eigenvalues[maxIndex])) {
        maxIndex = i;
      }
   }

   // First eigenvector
   Vector3d beta = eigensolver.eigenvectors().col(maxIndex); 
   beta.normalize();

   // Compute mean vector mu
   Vector3d mu = A.rowwise().mean();
   
   if (!mu.allFinite()) {
      std::cerr << "Error: mu contains NaN or Inf values!" << std::endl;
      return -1;
   }

   VectorXd result(numPixels);

   for (size_t i = 0; i < numPixels; ++i)
   {
      auto value = (A.col(i) - mu).dot(beta);
      (result)(i) = value;
   }

   // Normalization
   vector<double> luminance = normalizeToLuminance(result);   
   
   for (int j = 0; j < height; j++)
   {
       for (int i = 0; i < width; i++)
       {
         *pDestinationData++ = luminance[j*width + i];
         *pDestinationData++ = luminance[j*width + i];     
         *pDestinationData++ = luminance[j*width + i];  
       }
   }
   
   return 0;
}
