/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOZhang08.h"
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
using namespace Eigen;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOZhang08::TMOZhang08()
{
	SetName(L"Zhang08");					  // TODO - Insert operator name
	SetDescription(L"Zhang08"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);


   gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma Correction Factor [-]");
	gamma.SetDefault(2.2);
	gamma = 2.2;
	gamma.SetRange(1.0e-3, 1.0e+2);
	this->Register(gamma);
}

TMOZhang08::~TMOZhang08()
{
}

void stretchLabToFullRange(double* data, int width, int height) {
   double minL = data[0], maxL = data[0];
   double minA = data[1], maxA = data[1];
   double minB = data[2], maxB = data[2];

   // Find min max values in chanels
   for (int i = 0; i < width * height; ++i) {
       double L = data[i * 3 + 0];
       double a = data[i * 3 + 1];
       double b = data[i * 3 + 2];

       if (L < minL) minL = L;
       if (L > maxL) maxL = L;
       if (a < minA) minA = a;
       if (a > maxA) maxA = a;
       if (b < minB) minB = b;
       if (b > maxB) maxB = b;
   }

   // Linear scaling into [-1,1]
   for (int i = 0; i < width * height; ++i) {
       data[i * 3 + 0] = 2.0 * (data[i * 3 + 0] - minL) / (maxL - minL) - 1.0;  // L
       data[i * 3 + 1] = 2.0 * (data[i * 3 + 1] - minA) / (maxA - minA) - 1.0;  // a
       data[i * 3 + 2] = 2.0 * (data[i * 3 + 2] - minB) / (maxB - minB) - 1.0;  // b
   }
}

void stretchLuminanceToFullRange(vector<double> data, int width, int height) {
   
   double minL = 0.0, maxL = -1.0;

   // // Find min max values in chanel
   for (int i = 0; i < width * height; ++i) {
       double L = data[i ];

       if (L < minL) minL = L;
       if (L > maxL) maxL = L;
   }

   // Linear scaling into [0,1]
   for (int i = 0; i < width * height; ++i) {
       data[i] = 2.0 * (data[i] - minL) / (maxL - minL) - 1.0;  // L
   }
}

// Gamma correction for sRGB → linear RGB
double gammaCorrection(double c) {
   return (c <= 0.04045) ? (c / 12.92) : pow((c + 0.055) / 1.055, 2.4);
}
// Inverze gamma correction (linear RGB -> sRGB)
double inverseGammaCorrection(double c) {
   return (c <= 0.0031308) ? (12.92 * c) : (1.055 * pow(c, 1.0 / 2.4) - 0.055);
}

// Helping function for LAB transformation
double fLab(double t) {
   return (t > 0.008856) ? pow(t, 1.0 / 3.0) : (7.787 * t + 16.0 / 116.0);
}
// LAB -> XYZ helping function
double fInvLab(double t) {
   return (t > 0.206893) ? (t * t * t) : ((t - 16.0 / 116.0) / 7.787);
}

// sRGB to XYZ (D65)
void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z) {
   // sRGB -> linear RGB
   R = gammaCorrection(R);
   G = gammaCorrection(G);
   B = gammaCorrection(B);

   // Matrix fot computing sRGB -> XYZ (D65 reference white)
   X = R * 0.4124564 + G * 0.3575761 + B * 0.1804375;
   Y = R * 0.2126729 + G * 0.7151522 + B * 0.0721750;
   Z = R * 0.0193339 + G * 0.1191920 + B * 0.9503041;
}

// LAB back to XYZ
void LABtoXYZ(double L, double a, double b, double &X, double &Y, double &Z) {
   // Reference white point D65
   const double Xn = 0.95047, Yn = 1.00000, Zn = 1.08883;

   double fy = (L + 16.0) / 116.0;
   double fx = fy + (a / 500.0);
   double fz = fy - (b / 200.0);

   X = Xn * fInvLab(fx);
   Y = Yn * fInvLab(fy);
   Z = Zn * fInvLab(fz);
}

// XYZ to RGB
void XYZtoRGB(double X, double Y, double Z, double &R, double &G, double &B) {
   // Invert matrix XYZ -> sRGB
   R = X *  3.2404542 + Y * -1.5371385 + Z * -0.4985314;
   G = X * -0.9692660 + Y *  1.8760108 + Z *  0.0415560;
   B = X *  0.0556434 + Y * -0.2040259 + Z *  1.0572252;

   // Aplication of inverse gamma correction
   R = inverseGammaCorrection(R);
   G = inverseGammaCorrection(G);
   B = inverseGammaCorrection(B);

   // Clip into range [0, 1]
   R = std::max(0.0, std::min(1.0, R));
   G = std::max(0.0, std::min(1.0, G));
   B = std::max(0.0, std::min(1.0, B));
}

// XYZ to LAB - LAB is in range [-1,1] for L, A and B
void XYZtoLAB(double X, double Y, double Z, double &L, double &a, double &b) {
   // Reference white point D65
   const double Xn = 0.95047, Yn = 1.00000, Zn = 1.08883;

   double fx = fLab(X / Xn);
   double fy = fLab(Y / Yn);
   double fz = fLab(Z / Zn);

   L = (116.0 * fy) - 16.0;
   a = 500.0 * (fx - fy);
   b = 200.0 * (fy - fz);

   // Normalize into range L[0,1] A, B [-1, 1]
   L = (L / 100.0);
   a = a / 128.0;
   b = b / 128.0;

   a = (a+1)/2;
   b = (b+1)/2;
   //Gamma correct L
   //L = pow(L, 2.2);

}

// Main RGB to LAB funciton - output LAB is in range [-1,1]
void convertRGBtoLAB(double* data, int width, int height) {
   for (int i = 0; i < width * height; ++i) {
       double R = data[i * 3 + 0];
       double G = data[i * 3 + 1];
       double B = data[i * 3 + 2];

       double X, Y, Z, L, a, b;
       RGBtoXYZ(R, G, B, X, Y, Z);
       XYZtoLAB(X, Y, Z, L, a, b);

       data[i * 3 + 0] = L;
       data[i * 3 + 1] = a;
       data[i * 3 + 2] = b;
   }
   stretchLabToFullRange(data, width, height);
}

// Main funtion for back conversion LAB (classic range) to RGB [0,1]
void convertLABtoRGB(double* data, int width, int height) {
   for (int i = 0; i < width * height; ++i) {
       double L = data[i * 3 + 0];
       double a = data[i * 3 + 1];
       double b = data[i * 3 + 2];

       double X, Y, Z, R, G, B;
       LABtoXYZ(L, a, b, X, Y, Z);
       XYZtoRGB(X, Y, Z, R, G, B);

       data[i * 3 + 0] = R;
       data[i * 3 + 1] = G;
       data[i * 3 + 2] = B;
   }
}



// Function for cos kernel K5
double kernelK5(const VectorXd& yi, const VectorXd& yj, const double numPixels) 
{
   //double a = 1, b = 0.1;
   double a = 0.1, b = 1;
   double size = yi.size();
   double dot_ij = yi.dot(yj);  // Scalar product y_i * y_j
   double dot_ii = yi.dot(yi);  // Scalar product y_i * y_i
   double dot_jj = yj.dot(yj);  // Scalar product y_j * y_j

   double numerator = pow(a * dot_ij + b, 2);
   double denominator = (a * dot_ii + b) * (a * dot_jj + b);

   return numerator / denominator;
}

// Function for better kernel K6
double kernelK6(const VectorXd& yi, const VectorXd& yj, int i, int j, const double numPixels) 
{

   double K5_val = kernelK5(yi, yj, numPixels);
   auto logK5 = log10(1+K5_val);

   if (i != 0 && j != 0) {
      return logK5 * logK5;
   } else {
      return logK5;
   }
}

double kernelSigmoide(const VectorXd& yi, const VectorXd& yj) 
{
   auto dot = yi.dot(yj) / yi.size();
   auto size = yi.size();
   return std::tanh(0.01 * dot + 1.0); 
}

void saveToFile(const double* data, int height, int width, const std::string& filename = "../../out.txt") {
   std::ofstream file(filename);
   if (!file) {
       std::cerr << "Error: Cannot open file " << filename << " for writing!" << std::endl;
       return;
   }
   file << std::fixed << std::setprecision(1);

   for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++)
      {
         file << data[(i*width + j)*3];
         file << '\t';
         if (j == width - 1) file << '\n'; 
      }
   }
   file << "------------------------------------------------------------------------" << '\n';

   file << '\n'; 
   file.close();
}


void saveToFile(const vector<double> data, int height, int width, const std::string& filename = "../../out.txt") {
   std::ofstream file(filename);
   if (!file) {
       std::cerr << "Error: Cannot open file " << filename << " for writing!" << std::endl;
       return;
   }
   file << std::fixed << std::setprecision(1);

   for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++)
      {
         file << data[(i*width + j)];
         file << '\t';
         if (j == width - 1) file << '\n'; 
      }
   }
   file << "------------------------------------------------------------------------" << '\n';

   file << '\n'; 
   file.close();
}


// Maps projection values to the luminance range (0-100) using min-max normalization  
vector<double> normalizeToLuminance(const VectorXd& projections) 
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

   
   for (int i = 0; i < projections.size(); i++) {
       // Apply min-max normalization to range [0, 1]
       luminance[i] = ((projections[i] - minP) / (range));
       // Gamma correction
       //luminance[i] = inverseGammaCorrection(luminance[i]);
   }
   return luminance;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
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

   double maxL = -std::numeric_limits<double>::max();
   double minL = std::numeric_limits<double>::max();
   double maxA = pSourceData[1];
   double minA = pSourceData[1];
   double maxB = pSourceData[2];
   double minB = pSourceData[2];

   // Filling matrix A
   int index = 0;
   for (int j = 0; j < height; j++) 
   {
      for (int i = 0; i < width; i++) 
      {
         (A)(0, index) = *pSourceData++;     // L 
         maxL = std::max((A)(0, index), maxL);
         minL = std::min((A)(0, index), minL);

         (A)(1, index) = *pSourceData++; // a 
         maxA = std::max((A)(1, index), maxA);
         minA = std::min((A)(1, index), minA);
         
         (A)(2, index) = *pSourceData++;  // b
         maxB = std::max((A)(2, index), maxB);
         minB = std::min((A)(2, index), minB);
         index++;
      }
   }

   Matrix3d L = A * A.transpose();

   // Filling natrix K (3x3)
   for(int i = 0; i <3 ; i++)
   {
      for(int j = 0; j < 3; j++)
      {        
         K(i, j) = kernelK6(A.row(i), A.row(j), i, j, numPixels);
      }
   }

   // Matrix I_M with all elements 1/m 
   Matrix3d I_M = Matrix3d::Constant(1.0 / numPixels);

   // Compute K'
   Matrix3d K_prime = K - I_M * K - K * I_M + I_M * K * I_M;


   // Compute eigenvalues and eigenvectors of K'
   SelfAdjointEigenSolver<MatrixXd> eigensolver(K_prime);//K_prime
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

   //beta = beta * eigenvalues[maxIndex];
   //beta.normalize();

   // Compute mean vector mu
   Vector3d mu = A.rowwise().mean();
   
   if (!mu.allFinite()) {
      std::cerr << "Error: mu contains NaN or Inf values!" << std::endl;
      return -1;
   }

   // Compute projection P_i for each vector
   /*auto projections = (A.colwise() - mu);
   VectorXd result = (projections.transpose() * beta).array();*/

   VectorXd result(numPixels);

   for (size_t i = 0; i < numPixels; ++i)
   {
      auto value = (A.col(i) - mu).dot(beta);

      (result)(i) = value;
   }

   // Normalization
   vector<double> luminance = normalizeToLuminance(result);
   
   if (luminance.size() != numPixels)
      return -2;
   
   stretchLuminanceToFullRange(luminance, width, height);
   
   for (int j = 0; j < height; j++)
   {
       for (int i = 0; i < width; i++)
       {
         *pDestinationData++ = luminance[j*width + i];
         *pDestinationData++ = luminance[j*width + i];     
         *pDestinationData++ = luminance[j*width + i];  
       }
   }
   
   pSrc->ProgressBar(height, height);
   //pDst->CorrectGamma(gamma);


   return 0;
}
