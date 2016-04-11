#include "TMO.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <limits>

// number of histogram levels
#define HISTOGRAM_LEVELS 256

// maximal acceptable value in color range
#define MAX_VALUE_IN_RANGE 255.0

// every input filename has to contain this
#define INPUT_FILENAME_SUBSTR "_input."

// every model filename has to contain this
#define MODEL_FILENAME_SUBSTR "_model."

// comulative histogram will be normalised to this value
#define HISTOGRAM_NORMALISATION 100000

// window size for bilateral filter
#define WINDOW_SIZE 5

// defines how many times is bigger sigmaS for textureness bilateral filter than normal bilateral filter
#define SIGMA_S_TEX_MULTIPLIER 8.0

// maximal detail multiplier
#define MAX_RHO 5.0

// ignoring rho ratio over this
#define RHO_IGNORED_MAX 1000.0

// library for array 2d
//#include "../pftools/pfstmo.h"
#include "../tmolib/pfstmo.h"

class TMOBae06 : public TMO  
{
public:
	TMOBae06();
	virtual ~TMOBae06();
	virtual int Transform();

protected:
	const char * modelFileName;

private: 	
	TMOImage * model;
	double sigmaR;
	double sigmaS;
	TMOBool verbose;
	virtual double RgbToGray(double, double, double);
	virtual void InitialiseHistogram(int *);
	virtual void PrintHistogram(int *, std::string);
	virtual void FillHistogram(pfstmo::Array2D *, int *);
	virtual std::string GetModelFilename(std::string);
	virtual void ComputeComulativeHistogram(int *, int *);
	virtual int FindClosestShade(int, int *);
	virtual void NormaliseHistogram(int *, int);
	virtual void BilateralFilter(pfstmo::Array2D *, pfstmo::Array2D *);
	virtual double BilateralFilterWeight(double, double, double, double, pfstmo::Array2D *);
	virtual void HighPassFilter(pfstmo::Array2D, pfstmo::Array2D, double);
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
