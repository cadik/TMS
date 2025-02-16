/* --------------------------------------------------------------------------- *
 * TMOYourOperatorName.cpp: implementation of the TMOYourOperatorName class.   *
 * --------------------------------------------------------------------------- */

#include "TMOZhang08.h"
#include <Eigen/Dense>
#include <vector>
#include <cmath>

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
}

TMOZhang08::~TMOZhang08()
{
}

// Function for cos kernel K5
double kernelK5(const VectorXd& yi, const VectorXd& yj) {
   double a = 0.1, b = 1;
   double dot_ij = yi.dot(yj);  // Scalar product y_i * y_j
   double dot_ii = yi.dot(yi);  // Scalar product y_i * y_i
   double dot_jj = yj.dot(yj);  // Scalar product y_j * y_j

   double numerator = pow(a * dot_ij + b, 2);
   double denominator = (a * dot_ii + b) * (a * dot_jj + b);

   return numerator / denominator;
}

// Function for better kernel K6
double kernelK6(const VectorXd& yi, const VectorXd& yj, int i, int j) {
   double K5_val = kernelK5(yi, yj);
   if (i == 0 || j == 0) {
       return std::log(K5_val);
   } else {
       return std::log2(K5_val);
   }
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */

int TMOZhang08::Transform()
{
   pSrc->Convert(TMO_LAB); // Convert to L*a*b* format
   pDst->Convert(TMO_LAB); // Destination in L*a*b*
   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();
   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int numPixels = width * height;

   // Matrix A (3 × numPixels) 
   Eigen::MatrixXd A(3, numPixels);
   Eigen::Matrix3d K;

   // Filling matrix A
   int index = 0;
   for (int j = 0; j < height; j++) 
   {
      for (int i = 0; i < width; i++) 
      {
         (A)(0, index) = pSourceData[index * 3] - 0.5;     // L component centralised arround 0
         (A)(1, index) = pSourceData[index * 3 + 1] - 0.5; // a component centralised arround 0
         (A)(2, index) = pSourceData[index * 3 + 2] - 0.5; // b component centralised arround 0
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
   Eigen::Matrix3d I_M = Eigen::Matrix3d::Ones() / numPixels; 

   // Compute K'
   Eigen::Matrix3d K_prime = K - I_M * K - K * I_M + I_M * K * I_M;


   SelfAdjointEigenSolver<MatrixXd> solver(K_prime);
   VectorXd eigenvalues = solver.eigenvalues();
   MatrixXd eigenvectors = solver.eigenvectors();


   for (int j = 0; j < height; j++)
   {
       for (int i = 0; i < width; i++)
       {
           *pDestinationData++ = 0;
           *pDestinationData++ = 0;
           *pDestinationData++ = 0;
       }
   }
   
   pSrc->ProgressBar(height, height);
   pDst->Convert(TMO_RGB);
   return 0;
}
