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

// Function for cos kernel K5
double kernelK5(const VectorXd& yi, const VectorXd& yj) 
{
   double a = 0.1, b = 1;
   double dot_ij = yi.dot(yj);  // Scalar product y_i * y_j
   double dot_ii = yi.dot(yi);  // Scalar product y_i * y_i
   double dot_jj = yj.dot(yj);  // Scalar product y_j * y_j

   double numerator = pow(a * dot_ij + b, 2);
   double denominator = (a * dot_ii + b) * (a * dot_jj + b);

   return numerator / denominator;
}

// Function for better kernel K6
double kernelK6(const VectorXd& yi, const VectorXd& yj, int i, int j) 
{

   double K5_val = kernelK5(yi, yj);


   ///////
   //return K5_val;
   ///////

   if (i == 0 || j == 0) {
      return std::log(K5_val);
   } else {
      double logK5 = std::log(K5_val);
      return logK5 * logK5;
   }
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
       // Apply min-max normalization to range [0, 100]
       luminance[i] = 8*((projections[i] - minP) / (maxP - minP));
   }
   
   return luminance;
}


// Maps projection values to the luminance range (0-100) using min-max normalization
void normalizeToLuminance(double* projections, int size) 
{
    if (size <= 0) return;

    double minP = projections[0];
    double maxP = projections[0];

    // Find min and max values
    for (int i = 1; i < size; i++) {
        if (projections[i*3] < minP) minP = projections[i*3];
        if (projections[i*3] > maxP) maxP = projections[i*3];
    }

    double range = maxP - minP;
    if (range == 0) {
        // If all values are the same, set everything to 0
        std::fill(projections, projections + size, 0.0);
        return;
    }

    // Apply min-max normalization to range [0, 100]
    for (int i = 0; i < size; i = i+3) {
      projections[i*3] = ((projections[i*3] - minP) / range) * 100.0;
    }
}


void removeGammaCorrection(double* src, int width, int height, double gamma = 2.2) {
   if (!src) return;

   double invGamma = gamma; // Applying inverse gamma correction

   int numPixels = width * height * 3; // 3 channels per pixel
   for (int i = 0; i < numPixels; i++) {
       src[i] = std::pow(src[i], 2.2);
   }
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */

int TMOZhang08::Transform()
{
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int numPixels = width * height;

   //removeGammaCorrection(pSrc->GetData(), width, height);

   pSrc->Convert(TMO_LAB); // Convert to L*a*b* format
   pDst->Convert(TMO_LAB); // Destination in L*a*b*




   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();

   //saveToFile(pSourceData, height, width);

   //normalizeToLuminance(pSourceData, numPixels);

   //saveToFile(pSourceData, height, width);

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
         (A)(0, index) = *pSourceData++;     // L component centralised arround 0
         maxL = std::max((A)(0, index), maxL);
         minL = std::min((A)(0, index), minL);

         (A)(1, index) = *pSourceData++; // a component centralised arround 0
         maxA = std::max((A)(1, index), maxA);
         minA = std::min((A)(1, index), minA);
         
         (A)(2, index) = *pSourceData++;  // b component centralised arround 0
         maxB = std::max((A)(2, index), maxB);
         minB = std::min((A)(2, index), minB);
         index++;
      }
   }

   // Filling natrix K (3x3)
   for(int i = 0; i <3 ; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         (K)(i, j) = kernelK6(A.row(i), A.row(j), i, j);
      }
   }

   // Matrix I_M with all elements 1/m 
   Matrix3d I_M = Matrix3d::Ones() / numPixels; 

   for (int i = 0; i < 3; i++)
   {
      I_M(i,i) = 1;
   }

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
      if (std::abs(eigenvalues[i]) > std::abs(eigenvalues[maxIndex])) {
        maxIndex = i;
      }
   }

   // First eigenvector
   Vector3d beta = eigensolver.eigenvectors().col(maxIndex); 
   if (beta.size() != 3) {
      std::cerr << "Error: beta size is incorrect! Expected 3, got " << beta.size() << std::endl;
      return -1;
   }


   // Compute mean vector mu
   Vector3d mu = A.rowwise().mean();
   
  if (!mu.allFinite()) {
   std::cerr << "Error: mu contains NaN or Inf values!" << std::endl;
   return -1;
   }


   // Compute projection P_i for each vector
   auto projections = (A.colwise() - mu);
   VectorXd result = (projections.transpose() * beta).array();

   auto r = result.rows();
   auto c = result.cols();
   // Normalization
   vector<double> luminance = normalizeToLuminance(result);

   if (luminance.size() != numPixels)
      return -2;

   
   for (int j = 0; j < height; j++)
   {
       for (int i = 0; i < width; i++)
       {
           *pDestinationData++ = luminance[j*width + i];
           *pDestinationData++ = 0;
           *pDestinationData++ = 0;
       }
   }
   
   pSrc->ProgressBar(height, height);

   //saveToFile(luminance, height, width);
   pDst->Convert(TMO_RGB);
   //pDst->CorrectGamma(gamma);

   return 0;
}
