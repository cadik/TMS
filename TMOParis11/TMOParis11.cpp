/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2021                                              *
*                                                                              *
*                       Implementation of the TMOParis11 class                 *
*                                                                              *
*******************************************************************************/

#include <opencv2/opencv.hpp>
#include <omp.h>

#include "TMOParis11.h"

#define CHANNELSCNT 3
#define RGB true
#define LUM false
#define LIN true
#define LOG false

using namespace std;
using namespace cv;

TMOParis11::TMOParis11()
{
	SetName(L"Paris11");
	SetDescription(L"Local Laplacian Filters: Edge-aware "
	                "Image Processing with a Laplacian Pyramid");

	/* Inverse tone mapping */
	invToneMp.SetName(L"invToneMp");
	invToneMp.SetDescription(L"Set Inverse Tone Mapping -> if Detail Manipulation is OFF, "
	                          "setting this parameter ON change processing to Inverse Tone Mapping "
	                          "for LDR images. Otherwise is processing Tone Mapping for HDR images.");
	invToneMp.SetDefault(false);
	invToneMp=false;
	this->Register(invToneMp);

	/* Detail manipulation */
	detailMnpl.SetName(L"detailMnpl");
	detailMnpl.SetDescription(L"Set Detail Manipulation -> color remapping: rgb, domain: lin. "
	                           "Otherwise (inverse) tone mapping -> color remapping: lum, domain: log.");
	detailMnpl.SetDefault(false);
	detailMnpl=false;
	this->Register(detailMnpl);

	/* Gamma */
	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma correction value: <1.0e-3,1.0e+2>");
	gamma.SetDefault(2.2);
	gamma=2.2;
	gamma.SetRange(1.0e-3,1.0e+2);
	this->Register(gamma);

	/* Sigma_r parameter */
	sigma_r.SetName(L"sigma_r");
	sigma_r.SetDescription(L"Sigma_r parameter: <0.0, 5.0>.");
	sigma_r.SetDefault(2.5);
	sigma_r=2.5;
	sigma_r.SetRange(0.0, 5.0);
	this->Register(sigma_r);

	/* Beta parameter */
	beta.SetName(L"beta");
	beta.SetDescription(L"Beta parameter: <0.0, 5.0>. Tone mapping: 0.0 <= beta < 1.0. "
	                     "Detail enhancement: beta = 1.0. Inverse tone mapping: beta > 1.0");
	beta.SetDefault(0.0);
	beta=0.0;
	beta.SetRange(0.0, 5.0);
	this->Register(beta);

	/* Alpha parameter */
	alpha.SetName(L"alpha");
	alpha.SetDescription(L"Alpha parameter: <0.0, 5.0>. Detail smoothing: alpha > 1.0. "
	                      "Detail enhancement: 0.0 < alpha < 1.0");
	alpha.SetDefault(0.25);
	alpha=0.25;
	alpha.SetRange(0.0, 5.0);
	this->Register(alpha);
}

TMOParis11::~TMOParis11()
{
}

/*
 * Convert TMOImage to cv::Mat
 */
Mat TMOImage2Mat(TMOImage* pSrc)
{
	double* pSourceData;
	Vec3d* ptrMat;
	int rowsCnt, colsCnt;

	pSourceData = pSrc->GetData();
	rowsCnt = pSrc->GetHeight();
	colsCnt = pSrc->GetWidth();

	Mat srcConvMat(rowsCnt, colsCnt, CV_64FC3);
	
	/*
	 * If data are continuous, it is possible
	 * to work in one for loop
	 */
	if (srcConvMat.isContinuous())
	{
		colsCnt *= rowsCnt;
		rowsCnt = 1;
	}
	
	#pragma omp parallel for collapse(2)
	for (int y = 0; y < rowsCnt; y++)
	{
		ptrMat = srcConvMat.ptr<Vec3d>(y);

		for (int x = 0; x < colsCnt; x++)
		{
			ptrMat[x] = Vec3d(pSourceData[2],
			                  pSourceData[1],
			                  pSourceData[0]);
			
			/* Add count of channels (RGB) to pointer */
			pSourceData += CHANNELSCNT;
		}
	}

	return srcConvMat;
}

/*
 * Create Gaussian pyramid
 */
vector<Mat> gaussianPyramid(Mat srcMat, int nlev)
{
	vector<Mat> gaussPyr(nlev);			
	
	gaussPyr[0] = srcMat;
	
	for (int i = 1; i < nlev; i++)	
	{
		pyrDown(gaussPyr[i-1], gaussPyr[i]);		
	}

	return gaussPyr;
}

/*
 * Create Laplacian pyramid
 */
vector<Mat> laplacianPyramid(Mat I, int nlev)
{
	vector<Mat> lapPyr(nlev);
	Mat J, down, up;	

	J = I;

	for (int i = 0; i < nlev-1; i++)
	{
		pyrDown(J, down);				
		pyrUp(down, up, J.size());		
		lapPyr[i] = J - up;
		J = down;
	}
		
	lapPyr[nlev-1] = J;
	
	return lapPyr;
}

/*
 * Reconstruct laplacian pyramid (collapse)
 */
Mat reconstructLaplacianPyramid(vector<Mat> lapPyr, int nlev)
{
	Mat rec, tmp;
		
	rec = lapPyr[nlev-1];

	for(int i = nlev-2; i >= 0; i--) {		
		pyrUp(rec, tmp, lapPyr[i].size());
		rec = tmp + lapPyr[i];
	}

	return rec;		
}

/**********************************/
/* Helper functions for remapping */
/**********************************/

/*
 * Source: official MATLAB implementation
 * https://people.csail.mit.edu/sparis/publi/2011/siggraph/
 */

double smoothStep(double xmin, double xmax, double x)
{
	double y;
	
	y = (x - xmin)/(xmax-xmin);
	y = max(0.0, min(1.0,y));

	return (pow(y, 2)*pow(y-2,2));
}

double fd(double d, double alpha, double sigma_r)
{
	double tau, out;
	double noiseLevel = 0.01;

	out = pow(d, alpha);

	if (alpha < 1.0)
	{
		tau = smoothStep(noiseLevel, 2*noiseLevel, d*sigma_r);
		out = tau*out + (1-tau)*d; 
	}
	
	return out;
}

double fe(double a, double beta)
{
	return beta*a;
}

/**********************************/

double sign(double d)
{
	if (d > 0.0) return 1.0;
	if (d < 0.0) return -1.0;
	return 0.0;
}

/*
 * Remapping function for color vectors (RGB)
 * used for detail manipulation
 */
Mat remapColor(Mat& subMat, Vec3d g0, double alpha, double beta, double sigma_r)
{
	int height, width;
	double normV;
	Vec3d  v;
	Vec3d* ptrSubMat;
	Vec3d* ptrRemapped;

	height = subMat.rows;
	width = subMat.cols;

	Mat remapped(height, width, CV_64FC3);	

	if (remapped.isContinuous() &&
	    subMat.isContinuous())
	{
		width *= height;
		height = 1;
	}

	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrSubMat = subMat.ptr<Vec3d>(y);
		ptrRemapped = remapped.ptr<Vec3d>(y);

		for (int x = 0; x < width; x++)
		{
			v = ptrSubMat[x] - g0;
			normV = norm(v);

			if (normV != 0.0)
			{
				v /= normV;
			}

			if (normV > sigma_r)
			{
				/* r_e */
				ptrRemapped[x] = g0 + v * (fe(normV - sigma_r, beta) + sigma_r);
			}
			else
			{
				/* r_d */
				ptrRemapped[x] = g0 + v * sigma_r * fd(normV / sigma_r, alpha, sigma_r);
			}			
		}
	}	

	return remapped;
}

/*
 * Remapping function for luminance
 * used for tone mapping and inverse 
 * tone mapping
 */
Mat remapGray(Mat& subMat, double g0, double alpha, double beta, double sigma_r)
{
	int height, width;
	double dnrm, dsgn;
	double rd, re;
	double* ptrSubMat;
	double* ptrRemapped;

	height = subMat.rows;
	width = subMat.cols;
	
	Mat remapped(height, width, CV_64F);

	if (remapped.isContinuous() &&
	    subMat.isContinuous())
	{
		width *= height;
		height = 1;
	}

	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrSubMat = subMat.ptr<double>(y);
		ptrRemapped = remapped.ptr<double>(y);

		for (int x = 0; x < width; x++)
		{
			dnrm = abs(ptrSubMat[x]-g0);			
			dsgn = sign(ptrSubMat[x]-g0);			

			if (dnrm > sigma_r)
			{
				/* r_e */
				ptrRemapped[x] = g0 + dsgn*(fe(dnrm-sigma_r, beta)+sigma_r);
			}
			else
			{
				/* r_d */
				ptrRemapped[x] = g0 + dsgn*sigma_r*fd(dnrm/sigma_r, alpha, sigma_r);				
			}			
		}
	}
	
	return remapped;
}

/*
 * Calculate number of levels for input image
 */
int numLevels(int height, int width)
{
	double min_d = (double)min(height, width);
	int nlev = 1;

	while (min_d > 1.0)
	{
		nlev++;
		min_d = floor((min_d+1.0)/2.0);
	}

	return nlev;	
}

/*
 * Source:        TMS/TMOAubry14/TMOAubry14.cpp
 * Original code: https://github.com/daikiyamanaka/L0-gradient-smoothing
 */
void cvMat2Vec(const Mat &mat, vector<double> &vec){
	int rows = mat.rows;
	int cols = mat.cols;
	vec.resize(rows*cols);
	
	#pragma omp parallel for collapse(2)
	for(int i=0; i<rows; i++){
		double *ptr = reinterpret_cast<double*>(mat.data+mat.step*i);
		for(int j=0; j<cols; j++){
			vec[i*cols+j] = *ptr;
			++ptr;
		}
	}
}

/*
 * Source: TMS/TMOAubry14/TMOAubry14.cpp
 */
double prctileNearestRank(const vector<double> &vector, double percentile) {
	if(percentile < 0.0 || percentile > 100.0) {
		cerr << "Percentile not in range [0,100], setting it to 50\n";
		percentile = 50.0;
	}
	double ordinalRank = percentile/100.0 * vector.size();
	return vector[ceil(ordinalRank)-1];
}

/*
 * Convert color image (RGB) to intensity and color ratio
 */
void color2Intensity(Mat* srcMat, Mat* colorRatio, const double eps)
{
	int height, width;
	double* ptrIntChannel;
	Vec3d* ptrSrcMat;
	Vec3d* ptrColRatio;

	height = srcMat->rows;
	width = srcMat->cols;

	Mat intensityChannel(height, width, CV_64F);

	if (intensityChannel.isContinuous() &&
	    srcMat->isContinuous())
	{
		width *= height;
		height = 1;
	}

	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrIntChannel = intensityChannel.ptr<double>(y);
		ptrSrcMat = srcMat->ptr<Vec3d>(y);

		for (int x = 0; x < width; x++)
		{
			/* 
			 * Data are stored in BGR order
			 *
			 * Formula for intensity: 			 
			 * I_i = 1/61 * (20I_r + 40I_g + I_b)	
			 * 		 
			 */

			ptrIntChannel[x] = (ptrSrcMat[x][2]*20.0
			                 +  ptrSrcMat[x][1]*40.0
			                 +  ptrSrcMat[x][0])/61.0;			
		}
	}

	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrIntChannel = intensityChannel.ptr<double>(y);
		ptrSrcMat = srcMat->ptr<Vec3d>(y);
		ptrColRatio = colorRatio->ptr<Vec3d>(y);

		for (int x = 0; x < width; x++)
		{
			/*			 
			 * Formula for color ratio:
			 * (p_r, p_g, p_b) = 1/I_i * (I_r, I_g, I_b)
			 * 			  
			 */ 
			double divider = ptrIntChannel[x]+eps;
			ptrColRatio[x][2] = ptrSrcMat[x][2] / divider;
			ptrColRatio[x][1] = ptrSrcMat[x][1] / divider;
			ptrColRatio[x][0] = ptrSrcMat[x][0] / divider;
		}
	}

	srcMat->release();

	*srcMat = intensityChannel;	
}

/*
 * Convert intensity and color ratio back to color (BGR order)
 */
void intensity2Color(Mat* result, Mat* colorRatio, Mat* rec, const double eps)
{
	int height, width;
	double multiplier;
	double* ptrRec;
	Vec3d* ptrResult;
	Vec3d* ptrColRatio;

	height = rec->rows;
	width = rec->cols;

	if (result->isContinuous() &&
	    colorRatio->isContinuous() &&
	    rec->isContinuous())
	{
		width *= height;
		height = 1;
	}
	
	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrResult = result->ptr<Vec3d>(y);
		ptrColRatio = colorRatio->ptr<Vec3d>(y);
		ptrRec = rec->ptr<double>(y);

		for (int x = 0; x < width; x++)
		{		
			multiplier = ptrRec[x];
			ptrResult[x][2] = multiplier * ptrColRatio[x][2];
			ptrResult[x][1] = multiplier * ptrColRatio[x][1];
			ptrResult[x][0] = multiplier * ptrColRatio[x][0];			
		}
	}			
}

/*
 * Estimation of robust maximum and minimum
 * with the 99.5th and 0.5th percentiles
 * 
 * Source: official MATLAB implementation
 * https://people.csail.mit.edu/sparis/publi/2011/siggraph/
 */
void processPercentilesOutput(Mat *rec)
{
	double* ptrRec;
	double division;
	double DRDesired = 100.0;
	double prcClip = 0.5;	

	vector<double> pixelsVector;
	cvMat2Vec(*rec, pixelsVector);
	sort(pixelsVector.begin(), pixelsVector.end());	

	double maxClip = prctileNearestRank(pixelsVector, 100.0-prcClip);
	double minClip = prctileNearestRank(pixelsVector, prcClip);
	double DRClip = maxClip / minClip;
	double exponent = log(DRDesired) / log(DRClip);

	int height = rec->rows;
	int width = rec->cols;

	if (rec->isContinuous())
	{
		width *= height;
		height = 1;
	}

	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrRec = rec->ptr<double>(y);
		for (int x = 0; x < width; x++)
		{
			division = ptrRec[x] / maxClip;
			ptrRec[x] = (division > 0.0) ? pow(division, exponent) : 0.0;			
		}
	}
	
	pixelsVector.clear();		
}

/*
 * Map to the displayable range [0, 1]
 */
void normalizeMat(Mat* srcMat)
{
	int height, width;
	Vec3d* ptrSrcMat;

	height = srcMat->rows;
	width = srcMat->cols;

	if (srcMat->isContinuous())
	{
		width *= height;
		height = 1;
	}

	#pragma omp parallel for collapse(2)
	for (int y = 0; y < height; y++)
	{
		ptrSrcMat = srcMat->ptr<Vec3d>(y);

		for (int x = 0; x < width; x++)
		{
			if (ptrSrcMat[x][2] < 0.0)
			{
				ptrSrcMat[x][2] = 0.0;
			}
			else if (ptrSrcMat[x][2] > 1.0)
			{
				ptrSrcMat[x][2] = 1.0;
			}
			
			if (ptrSrcMat[x][1] < 0.0)
			{
				ptrSrcMat[x][1] = 0.0;
			}
			else if (ptrSrcMat[x][1] > 1.0)
			{
				ptrSrcMat[x][1] = 1.0;
			}

			if (ptrSrcMat[x][0] < 0.0)
			{
				ptrSrcMat[x][0] = 0.0;
			}
			else if (ptrSrcMat[x][0] > 1.0)
			{
				ptrSrcMat[x][0] = 1.0;
			}	
		}	
	}
}

/*
 * Core function for processing Laplacian filter for color vectors (BGR)
 * Working with Vec3d values
 */
Mat lapFilterColor(Mat& srcMat,
                   int nlev,				   
                   double alpha,
                   double beta,
                   double sigma_r)
{
	int gaussPyrHeight, gaussPyrWidth, origGaussWidth;	
	bool isCont;
	Vec3d* ptrGaussPyr;
	Vec3d* ptrLapPyr;	
	
	int height = srcMat.rows;	
	int width  = srcMat.cols;

	isCont = false;

	/* GAUSSIAN PYRAMID */
	vector<Mat> gaussPyr = gaussianPyramid(srcMat, nlev);

	/* LAPLACIAN PYRAMID */	
	int dim = srcMat.channels();
	
	// Allocate space for result
	vector<Mat> lapPyr = laplacianPyramid(Mat::zeros(srcMat.size(), 
	                                      CV_MAKETYPE(CV_64FC3, dim)),
	                                      nlev);
	
	#pragma omp parallel for	
	for (int lev = 0; lev < nlev-1; lev++)
	{
		gaussPyrHeight = gaussPyr[lev].rows;
		gaussPyrWidth = gaussPyr[lev].cols;

		if (gaussPyr[lev].isContinuous() &&
		    lapPyr[lev].isContinuous())
		{
			origGaussWidth = gaussPyrWidth;
			gaussPyrWidth *= gaussPyrHeight;			
			gaussPyrHeight = 1;			
			isCont = true;
		}
				
		int hw = 3*pow(2,lev+1)-2;	

		#pragma omp parallel for collapse(2)
		for (int y = 0; y < gaussPyrHeight; y++)
		{
			ptrGaussPyr = gaussPyr[lev].ptr<Vec3d>(y);
			ptrLapPyr = lapPyr[lev].ptr<Vec3d>(y);

			for (int x = 0; x < gaussPyrWidth; x++)
			{
				int xf = (isCont) ? (x%origGaussWidth)*pow(2, lev) : x*pow(2, lev);
				int yf = (isCont) ? (x/origGaussWidth)*pow(2, lev) : y*pow(2, lev);										

				int yrng1 = max(0, yf-hw); 
				int yrng2 = min(height-1, yf+hw);
				int xrng1 = max(0, xf-hw);
				int xrng2 = min(width-1, xf+hw);							
								
				Mat subMat = srcMat(Rect(xrng1, yrng1, xrng2-xrng1+1, yrng2-yrng1+1));				
								
				Vec3d g0 = ptrGaussPyr[x];
												
				Mat remapMat = remapColor(subMat, g0, alpha, beta, sigma_r);
				vector<Mat> remapLap = laplacianPyramid(remapMat, lev+2);		

				Mat tmp;

				for(int sublev=0; sublev<lev; sublev++) {
					pyrUp(remapLap[lev], tmp);					
					remapLap[lev] = tmp;
				}
												
				ptrLapPyr[x] = remapLap[lev].at<Vec3d>(yf-yrng1, xf-xrng1);
			}
		}
	}
	
	lapPyr.back() = gaussPyr.back();

	Mat rec = reconstructLaplacianPyramid(lapPyr, nlev);
	
	return rec;
}

/*
 * Core function for processing Laplacian filter for luminance
 * Working with double values
 */
Mat lapFilterToneMapping(Mat& srcMat, int nlev, 
                         double alpha, 
                         double beta, 
                         double sigma_r)
{
	int gaussPyrHeight, gaussPyrWidth, origGaussWidth;	
	double* ptrGaussPyr;
	double* ptrLapPyr;
	bool isCont;

	int height = srcMat.rows;	
	int width  = srcMat.cols;

	isCont = false;
	
	/* GAUSSIAN PYRAMID */
	vector<Mat> gaussPyr = gaussianPyramid(srcMat, nlev);

	/* LAPLACIAN PYRAMID */	
	int dim = srcMat.channels();

	// Allocate space for result
	vector<Mat> lapPyr = laplacianPyramid(Mat::zeros(srcMat.size(),
	                                      CV_MAKETYPE(CV_64F, dim)),
	                                      nlev);	

	#pragma omp parallel for	
	for (int lev = 0; lev < nlev-1; lev++)
	{
		gaussPyrHeight = gaussPyr[lev].rows;
		gaussPyrWidth = gaussPyr[lev].cols;

		if (gaussPyr[lev].isContinuous() &&
		    lapPyr[lev].isContinuous())
		{
			origGaussWidth = gaussPyrWidth;
			gaussPyrWidth *= gaussPyrHeight;			
			gaussPyrHeight = 1;			
			isCont = true;
		}
				
		int hw = 3*pow(2,lev+1)-2;

		#pragma omp parallel for collapse(2)
		for (int y = 0; y < gaussPyrHeight; y++)
		{			
			ptrGaussPyr = gaussPyr[lev].ptr<double>(y);
			ptrLapPyr = lapPyr[lev].ptr<double>(y);

			for (int x = 0; x < gaussPyrWidth; x++)
			{
				int xf = (isCont) ? (x%origGaussWidth)*pow(2, lev) : x*pow(2, lev);
				int yf = (isCont) ? (x/origGaussWidth)*pow(2, lev) : y*pow(2, lev);										

				int yrng1 = max(0, yf-hw); 
				int yrng2 = min(height-1, yf+hw);
				int xrng1 = max(0, xf-hw);
				int xrng2 = min(width-1, xf+hw);							
								
				Mat subMat = srcMat(Rect(xrng1, yrng1, xrng2-xrng1+1, yrng2-yrng1+1));
						
				double g0 = ptrGaussPyr[x];
								
				Mat remapMat = remapGray(subMat, g0, alpha, beta, sigma_r);								
				vector<Mat> remapLap = laplacianPyramid(remapMat, lev+2);		

				Mat tmp;

				for(int sublev=0; sublev<lev; sublev++) {
					pyrUp(remapLap[lev], tmp);					
					remapLap[lev] = tmp;
				}
					
				ptrLapPyr[x] = remapLap[lev].at<double>(yf-yrng1, xf-xrng1);
			}
		}
	}
	
	lapPyr.back() = gaussPyr.back();

	Mat rec = reconstructLaplacianPyramid(lapPyr, nlev);

	return rec;
}

/*
 * Core function for remapping 
 * and calling main operator algorithm
 */
Mat laplacianFilter(TMOImage* pSrc,
                    double alpha,
                    double beta,
                    double sigma_r,
                    double gamma,
                    bool detailMnpl,
                    bool invToneMp)
{
	/* Convert TMOImage to cv::Mat */
	Mat srcMat = TMOImage2Mat(pSrc);		

	const double eps = pow(2,-52);	

	int height = srcMat.rows;
	int width  = srcMat.cols;
	/* Get number of levels for input image */
	int nlev = numLevels(height, width);			

	bool colorRemapping, domain;

	/*
	 * Set colorRemapping and domain for Detail Manipulation,
	 * Tone Mapping or Inverse Tone Mapping
	 */
	if (detailMnpl)
	{
		colorRemapping = RGB;
		domain = LIN;

		/* Input is LDR image */
		srcMat.convertTo(srcMat, CV_64FC3, 1.0/255.0);
	}
	else
	{
		colorRemapping = LUM;
		domain = LOG;

		if (invToneMp)
		{
			/* Input is LDR image */
			srcMat.convertTo(srcMat, CV_64FC3, 1.0/255.0);
			pow(srcMat, gamma, srcMat);
		}		
	}

	/************************************/
	/* Color Remapping for Tone Mapping */
	/************************************/
	
	Mat colorRatio(height, width, CV_64FC3);	

	if (colorRemapping == LUM)
	{
		color2Intensity(&srcMat, &colorRatio, eps);
	}
	else
	{
		colorRatio.release();
	}

	if (domain == LOG)
	{
		sigma_r = log(sigma_r);		
		log(srcMat + eps, srcMat);			
	}				


	/************************************/	
	/*          Lapfilter Core          */
	/************************************/	

	Mat rec;

	if (colorRemapping == LUM && domain == LOG)
	{
		rec = lapFilterToneMapping(srcMat, nlev, alpha, beta, sigma_r);
	}
	else
	{		
		rec = lapFilterColor(srcMat, nlev, alpha, beta, sigma_r);	
	}


	/***************************************/
	/* Color mapping back for Tone Mapping */
	/***************************************/

	if (domain == LOG)
	{
		exp(rec, rec);		
		rec -= eps;	
	}	

	if (domain == LOG && beta <= 1.0)
	{
		processPercentilesOutput(&rec);
	}

	Mat result(height, width, CV_64FC3);

	if (colorRemapping == LUM)
	{
		intensity2Color(&result, &colorRatio, &rec, eps);
	}
	else
	{
		result = rec;
	}

	if(domain == LOG && beta <= 1.0)
	{
		/* Map to range [0, 1] */
		normalizeMat(&result);
		/* Gamma correction */		
		pow(result, 1/gamma, result);
	}		

	return result;		
}

/* --------------------------------------------------------------------------------------- */

/*
 * Main method for Local Laplacian Filters operator
 */
int TMOParis11::Transform()
{
	
	int height, width;
	bool pDetailMnpl, pInvToneMp;
	double pAlpha, pBeta, pSigmaR, pGamma;
	double* pDestData;
	Vec3d* ptrResult;
	Mat result;	

	/***********************/
	/* Get user parameters */
	/***********************/
	
	pAlpha  = alpha.GetDouble();
	pBeta   = beta.GetDouble();
	pSigmaR = sigma_r.GetDouble(); 
	pGamma  = gamma.GetDouble();
	pDetailMnpl = detailMnpl.GetBool();
	pInvToneMp  = invToneMp.GetBool();


	/***********************/
	/* Call main algorithm */
	/***********************/
	
	pDestData = pDst->GetData();
	result = laplacianFilter(pSrc, pAlpha, pBeta, pSigmaR, pGamma, pDetailMnpl, pInvToneMp);		

	/***********************/	

	width = result.cols;
	height = result.rows;	

	if (result.isContinuous())
	{
		width *= height;
		height = 1;
	}	

	int y = 0;

	/*
	 * Save result to the destination image
	 */
	#pragma omp parallel for collapse(2)
	for (; y < height; y++)	
	{
		pSrc->ProgressBar(y, height);
		
		ptrResult = result.ptr<Vec3d>(y);

		for (int x = 0; x < width; x++)		
		{
			*pDestData++ = ptrResult[x][2];	
			*pDestData++ = ptrResult[x][1];
			*pDestData++ = ptrResult[x][0];			
		}
	}

	pSrc->ProgressBar(y, height);

	return 0;
}
