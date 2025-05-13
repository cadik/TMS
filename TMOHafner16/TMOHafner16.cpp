/*******************************************************************************
*                                                                              *
*                        Brno University of Technology                         *
*                      Faculty of Information Technology                       *
*                                                                              *
*                        Color-to-Grayscale Conversions                        *
*                                                                              *
*            Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]          *
*                                   Brno 2025                                  *
*                                                                              *
*                     Implementation of the TMOHafner16 class                  *
*             Variational Image Fusion with Optimal Local Contrast             *
*                 https://doi.org/10.1007/978-3-319-18461-6_34                 *
*                                                                              *
*******************************************************************************/

#include "TMOHafner16.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <omp.h>

#include <fstream>

using namespace std;


/*
* Constructor - Initializes the operator's name and description
*/
TMOHafner16::TMOHafner16()
{
	SetName(L"Hafner16");					  
	SetDescription(L"Variational Image Fusion with Optimal Local Contrast"); 

   HDRParameter.SetName(L"HDR");
	HDRParameter.SetDescription(L"is input image HDR");
	HDRParameter.SetDefault(false);
	HDRParameter = false;

   this->Register(HDRParameter);
}

TMOHafner16::~TMOHafner16()
{
}

/*
* Gaussian filter using recursive approximation
*/
double TMOHafner16::gaussianWeight(double x, double y) {

   double dist2 = x * x + y * y;
   return exp(-dist2 / (2.0 * sigma * sigma)) * (inv_2pi / (sigma * sigma));
}

/*
* Approximated psiDerivative using polynomial expansion
*/
double TMOHafner16::psiDerivative(double z) {
   return z / sqrt(z * z + lambda * lambda);
}

/*
* Computing laplacian
*/
void TMOHafner16::computeLaplacian(const vector<double> weights, vector<double>& laplacian, int width, int height) {
   #pragma omp parallel for
   for (int i = width; i < width * (height - 1); ++i) {
      if (i % width == 0 || (i + 1) % width == 0) continue; // Skip borders
      laplacian[i] = weights[i - width] + weights[i + width] +
                     weights[i - 1] + weights[i + 1] -
                     4 * weights[i];
   }
}

/*
* Perform one iteration of gradient descent update
*/
void TMOHafner16::gradientDescentStep(vector<double>& wr, vector<double>& wg, vector<double>& wb, const vector<double>& r, const vector<double>& g, const vector<double>& b, int width, int height)
{
   int N = width * height;

   // Compute u(x) = weighted sum of RGB components
   vector<double> u(N, 0.0);

   #pragma omp parallel for
   for (int i = 0; i < N; i++) {
       u[i] = wr[i] * r[i] + wg[i] * g[i] + wb[i] * b[i];
   }

   // Compute Laplacians for smoothing
   vector<double> laplacian_wr(N, 0.0);
   vector<double> laplacian_wg(N, 0.0);
   vector<double> laplacian_wb(N, 0.0);

   computeLaplacian(wr, laplacian_wr, width, height);
   computeLaplacian(wg, laplacian_wg, width, height);
   computeLaplacian(wb, laplacian_wb, width, height);

   // Gradient Descent Update
   vector<double> new_wr(N), new_wg(N), new_wb(N);

   #pragma omp parallel for collapse(2)
   for (int y = 0; y < height; y++) 
   {
      for (int x = 0; x < width; x++) 
      {
         int i = y * width + x; // Pixel index
         double f_bar = (r[i] + g[i] + b[i]) / 3.0; // Mean intensity
  
         // Compute data term gradients
         double grad = u[i] - f_bar + delta * (u[i] - mu);
         
         // Compute regularization term using 3×3 window
         double reg = 0.0;
         int count = 0;
         for (int dy = -1; dy <= 1; dy++)
         {
            for (int dx = -1; dx <= 1; dx++) 
            {
              int nx = x + dx, ny = y + dy;
              if (nx >= 0 && nx < width && ny >= 0 && ny < height) 
              {
                  int j = ny * width + nx; // neighbour pixel
                  double diff = u[i] - u[j];
                  reg += psiDerivative(diff)*gaussianWeight(dx, dy);
              }
            }
         }
         reg *= gama;
  
         // Update weights
         wr[i] -= tau * (r[i] * (grad - reg) - alpha * laplacian_wr[i]);
         wg[i] -= tau * (g[i] * (grad - reg) - alpha * laplacian_wg[i]);
         wb[i] -= tau * (b[i] * (grad - reg) - alpha * laplacian_wb[i]);
      }
   }
}


/*
* Projects the weights onto a simplex constraint - algorithm 2 
*/
void TMOHafner16::projectOntoSimplex(vector<double>& wr, vector<double>& wg, vector<double>& wb) {
   int n = wr.size();
   std::vector<double> s;
   s.reserve(3*n);

   // Combine all weights into a single vector
   s.insert(s.end(), wr.begin(), wr.end());
   s.insert(s.end(), wg.begin(), wg.end());
   s.insert(s.end(), wb.begin(), wb.end());

   // Sorting values in descending order
   std::sort(s.begin(), s.end(), std::greater<double>());

   // Finding m, sum_s and theta
   double sum_s = 0.0, theta = 0.0;
   int m = 1;
   
   for (int j = 0; j < n; ++j) {
       sum_s += s[j];
       if (s[j] - (sum_s / (j + 1)) > 0) {
           m = j;
       } else {
           break;
       }
   }
   
   theta = (sum_s - 1.0) / m;

   // Projection onto simplex
   for (int i = 0; i < n; ++i) {
       wr[i] = std::max(wr[i] - theta, 0.0);
       wg[i] = std::max(wg[i] - theta, 0.0);
       wb[i] = std::max(wb[i] - theta, 0.0);
   }
}

/*
* Finds if range is 0-1 or in 0-255
*/
bool TMOHafner16::isInRange0to1(double *pSourceData, int numPix)
{
   for (int i = 0; i < numPix * 3; i++)
   {
      if(pSourceData[i] > 1)
         return false;
   }
   return true;
}


/* --------------------------------------------------------------------------- *
 * Main transformation funcion - accorting to algorithm 1                      *
 * --------------------------------------------------------------------------- */
int TMOHafner16::Transform()
{
   auto start = std::chrono::high_resolution_clock::now();
   
	double *pSourceData = pSrc->GetData();	
	double *pDestinationData = pDst->GetData(); 

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int numPix = width * height;

   vector<double> r(numPix), g(numPix), b(numPix);

   bool range0to1 = isInRange0to1(pSourceData, numPix);

   // Compute mean intensity (mu)
   mu = 0;
   for (int i = 0; i < numPix; i ++)
   {
      r[i] = *pSourceData++;
      g[i] = *pSourceData++;
      b[i] = *pSourceData++;


      // If format is in range 0-255
      if (!range0to1 && !HDRParameter)
      {
         r[i] /= 255;
         g[i] /= 255;
         b[i] /= 255;
      }

      if (HDRParameter)
      {
         r[i] /= r[i] + 1;
         g[i] /= g[i] + 1;
         b[i] /= b[i] + 1;
      }
      

      mu += r[i] + g[i] + b[i]; 
   }
   mu /= (numPix * 3);

   // Set sigma based on image size
   sigma = 0.1 * sqrt(height+width);

   // Initialize weights uniformly
   vector<double> wr(numPix, 1.0/3.0), wg(numPix, 1.0/3.0), wb(numPix, 1.0/3.0);
      
   // Iterative optimization process
   for (int i = 0; i < numIter; i++)
   {
      gradientDescentStep(wr, wg, wb, r, g, b, width, height);
      projectOntoSimplex(wr, wg, wb);
   }

   // Normalization of output
   double minOut = std::numeric_limits<double>::max();
   double maxOut = std::numeric_limits<double>::lowest();

   for (int i = 0; i < numPix; i++) {
      double pOut = wr[i] * r[i] + wg[i] * g[i] + wb[i] * b[i];
      minOut = std::min(minOut, pOut);
      maxOut = std::max(maxOut, pOut);
   }
   
   double range = maxOut - minOut;
   if (range < 1e-6) range = 1.0; // Avoid division by zero
   
   for (int i = 0; i < numPix; i++) {
      double pOut = wr[i] * r[i] + wg[i] * g[i] + wb[i] * b[i];

      *pDestinationData++ = pOut;
      *pDestinationData++ = pOut;
      *pDestinationData++ = pOut;
   }

	return 0;
}