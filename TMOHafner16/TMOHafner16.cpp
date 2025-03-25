/* --------------------------------------------------------------------------- *
 * TMOHafner16.cpp: implementation of the TMOHafner16 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOHafner16.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

const int n = 3; // Number of canals
const double lambda = 0.1; 
const double alpha = 1; 
const double gama = 0.25;
const double delta = 1;
const double tau = 0.15;
const double sigma = 0.1; // For gausian weight
double mu = 0.5;
const int numIter = 100; // Number of optimization iterations


/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOHafner16::TMOHafner16()
{
	SetName(L"Hafner16");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOHafner16::~TMOHafner16()
{
}

void RGB2YCbCr(double* data, int numPix)
{
   double r, g, b;
   for(int i = 0; i < numPix*3; i += 3)
   {
      r = data[i];
      g = data[i+1];
      b = data[i+2];

      data[i] = 0.299 * r + 0.587 * g + 0.114 * b; // Y
      data[i+1] = 128 + (-0.168736 * r - 0.33126 * g + 0.5 * b); // Cb
      data[i+2] = 128 + (0.5 * r - 0.418688 * g - 0.081312 * b); // Cr
   }
}

double gaussianWeight(double x, double y, double sigma) {
   double dist2 = x * x + y * y;
   return exp(-dist2 / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
}


double psiDerivative(double z, double lambda) {
   return z / sqrt(z * z + lambda * lambda);
}

//////////////////////////////////
// Correct laplacian? 5 point
///////////////////////////////
void computeLaplacian(const vector<double> weights, vector<double> &laplacian, int width, int height) {
   for (int y = 1; y < height - 1; ++y) { // Eliminate edges
       for (int x = 1; x < width - 1; ++x) {
           int i = y * width + x;
           laplacian[i] = weights[i - width]  // up
                        + weights[i + width]  // down
                        + weights[i - 1]      // left
                        + weights[i + 1]      // right
                        - 4 * weights[i];     // center
       }
   }
}

void gradientDescentStep(vector<double>& wr, vector<double>& wg, vector<double>& wb, const vector<double>& r, const vector<double>& g, const vector<double>& b, int width, int height)
{
   int N = width * height;

   // 1. Compute u(x)
   vector<double> u(N, 0.0);

   for (int i = 0; i < N; i++) {
       u[i] = wr[i] * r[i] + wg[i] * g[i] + wb[i] * b[i];
   }

   // 2. Compute laplacian
   vector<double> laplacian_wr(N, 0.0);
   vector<double> laplacian_wg(N, 0.0);
   vector<double> laplacian_wb(N, 0.0);

   computeLaplacian(wr, laplacian_wr, width, height);
   computeLaplacian(wg, laplacian_wg, width, height);
   computeLaplacian(wb, laplacian_wb, width, height);

   // 3. Gradient Descent Update
   vector<double> new_wr(N), new_wg(N), new_wb(N);
   for (int y = 0; y < height; y++) 
   {
      for (int x = 0; x < width; x++) 
      {
         int i = y * width + x; // pixel index
  
         // Mean final intensity
         double f_bar = (r[i] + g[i] + b[i]) / 3.0;
  
         // Data term
         double grad_r = u[i] - f_bar + delta * (u[i] - mu);
         double grad_g = u[i] - f_bar + delta * (u[i] - mu);
         double grad_b = u[i] - f_bar + delta * (u[i] - mu);
         
         // 3×3 window for regulation
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
                  reg += psiDerivative(diff, lambda)*gaussianWeight(dx, dy, sigma);
              }
            }
         }
         reg *= gama;
  
         // Laplacian smoothing
         double laplacian_r = laplacian_wr[i];
         double laplacian_g = laplacian_wg[i];
         double laplacian_b = laplacian_wb[i];
  
         // Final weights
         new_wr[i] = wr[i] - tau * (r[i] * (grad_r - reg) - alpha * laplacian_r);
         new_wg[i] = wg[i] - tau * (g[i] * (grad_g - reg) - alpha * laplacian_g);
         new_wb[i] = wb[i] - tau * (b[i] * (grad_b - reg) - alpha * laplacian_b);
      }
   }
   // 3. Aktualizace vah
   wr = new_wr;
   wg = new_wg;
   wb = new_wb;
}



void projectOntoSimplex(vector<double>& wr, vector<double>& wg, vector<double>& wb) {
   int n = wr.size();

   std::vector<double> s;
   s.reserve(wr.size() + wg.size() + wb.size());
   s.insert(s.end(), wr.begin(), wr.end());
   s.insert(s.end(), wg.begin(), wg.end());
   s.insert(s.end(), wb.begin(), wb.end());

   // Sorting values
   std::sort(s.begin(), s.end(), std::greater<double>());

   // Finding m, sum s and theta
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

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOHafner16::Transform()
{
	double *pSourceData = pSrc->GetData();	
	double *pDestinationData = pDst->GetData(); 

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();
   int numPix = width * height;

   //RGB2YCbCr(pSourceData, numPix);
	double pR, pG, pB, pOut;
   vector<double> r(numPix), g(numPix), b(numPix);

   for (int i = 0; i < numPix; i ++)
   {
      r[i] = *pSourceData++;
      g[i] = *pSourceData++;
      b[i] = *pSourceData++;
      mu += r[i] + g[i] + b[i]; 
   }

   mu /= (numPix * 3);

   vector<double> wr(numPix, 1.0/3.0), wg(numPix, 1.0/3.0), wb(numPix, 1.0/3.0);

   int numIter = 100; // Number of optimization iterations
   for (int i = 0; i < numIter; i++)
   {
      gradientDescentStep(wr, wg, wb, r, g, b, width, height);
      projectOntoSimplex(wr, wg, wb);
   }


   double minOut = std::numeric_limits<double>::max();
   double maxOut = std::numeric_limits<double>::lowest();
   
   for (int i = 0; i < numPix; i++) {
      double pOut = wr[i] * r[i] + wg[i] * g[i] + wb[i] * b[i];
      minOut = std::min(minOut, pOut);
      maxOut = std::max(maxOut, pOut);
   }
   
   // Zabránit dělení nulou
   double range = maxOut - minOut;
   if (range < 1e-6) range = 1.0;  // Když jsou všechny hodnoty skoro stejné
   
   // Normalizace a uložení výsledku
   for (int i = 0; i < numPix; i++) {
      double pOut = wr[i] * r[i] + wg[i] * g[i] + wb[i] * b[i];
      pOut = (pOut - minOut) / range; // Normalizace do [0,1]
   
      *pDestinationData++ = pOut;
      *pDestinationData++ = pOut;
      *pDestinationData++ = pOut;
   }   
	return 0;
}