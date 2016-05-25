/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                 diploma thesis                               *
*             Author: Petr Pospisil [xpospi68 AT stud.fit.vutbr.cz]            *
*                                    Brno 2016                                 *
*                                                                              *
*******************************************************************************/

#include "TMO.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <limits>

// this is multiplier for hitoogram levels to avoid posterize effect, has to be integer
#define HISTOGRAM_BOOSTER 100
#define HISTOGRAM_BOOSTER_F 100.0

// number of histogram levels
#define HISTOGRAM_LEVELS (256 * HISTOGRAM_BOOSTER)
#define HISTOGRAM_LEVELS_F (256.0 * HISTOGRAM_BOOSTER_F)

// maximal acceptable value in color range
#define MAX_VALUE_IN_RANGE 255.0

// comulative histogram will be normalised to this value
#define HISTOGRAM_NORMALISATION 100000000

// window size for bilateral filter for LDR
#define WINDOW_SIZE_LDR 5

// window size for bilateral filter for HDR
#define WINDOW_SIZE_HDR 3

// defines how many times is bigger sigmaS for textureness bilateral filter than normal bilateral filter
#define SIGMA_S_TEX_MULTIPLIER 8.0

// maximal detail multiplier
#define MAX_RHO 5.0

// ignoring rho ratio over this
// #define RHO_IGNORED_MAX 1000.0

// value to scale from [0, 100] range to [0, 255] range
#define SCALE_LAB_TO_RGB 2.55

// sigma_r for LDR, for simplicity is used constant
#define SIGMA_R_LDR 15.0

// sigma_r for HDR, for simplicity is used constant
#define SIGMA_R_HDR 1.0

// library for array 2d
#include "../tmolib/pfstmo.h"

class TMOBae06 : public TMO  
{
public:
	TMOBae06();
	virtual ~TMOBae06();
	virtual int Transform();

protected:
	
private: 	
	TMOImage * model;
	double sigmaR;
	double sigmaS;
	TMOBool verbose;
	TMOBool hdr;
	TMOString modelFileNameParam;
	
	virtual double RgbToGray(double, double, double);
	virtual void InitialiseHistogram(int *);
	virtual void PrintHistogram(int *, std::string);
	virtual void FillHistogram(pfstmo::Array2D *, int *);
	//virtual std::string GetModelFilename(std::string);
	virtual void ComputeComulativeHistogram(int *, int *);
	virtual int FindClosestShade(int, int *);
	virtual void NormaliseHistogram(int *, int);
	virtual void BilateralFilter(pfstmo::Array2D *, pfstmo::Array2D *);
	virtual double BilateralFilterWeight(double, double, double, double, pfstmo::Array2D *);
	//virtual void HighPassFilter(pfstmo::Array2D, pfstmo::Array2D, double);
	virtual void CreateGrayscale(pfstmo::Array2D, TMOImage *);
	virtual void GetDetailFromBase(pfstmo::Array2D, pfstmo::Array2D, pfstmo::Array2D);
	virtual double ComputeSigmaS(int, int);
	virtual void HistogramMatching(int *, int *, pfstmo::Array2D *);
	virtual void FillRho(pfstmo::Array2D, pfstmo::Array2D, pfstmo::Array2D, pfstmo::Array2D);
	virtual void FillTextureness(pfstmo::Array2D *, pfstmo::Array2D *, bool);
	virtual void HighPassFilterV2(pfstmo::Array2D *);
	virtual void CrossBilateralFilter(pfstmo::Array2D *, pfstmo::Array2D *, pfstmo::Array2D *);
	virtual double CrossBilateralFilterWeight(double, double, double, double, pfstmo::Array2D *, pfstmo::Array2D *, bool);
};
