#include "TMO.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

// number of histogram levels
#define HISTOGRAM_LEVELS 256

// every input filename has to contain this
#define INPUT_FILENAME_SUBSTR "_input."

// every model filename has to contain this
#define MODEL_FILENAME_SUBSTR "_model."

// comulative histogram will be normalised to this value
#define HISTOGRAM_NORMALISATION 100000

// window size for bilateral filter
#define WINDOW_SIZE 5

// library for array 2d
#include "../pftools/pfstmo.h"

class TMO2scaleToneManagement : public TMO  
{
public:
	TMO2scaleToneManagement();
	virtual ~TMO2scaleToneManagement();
	virtual int Transform();

protected:
	const char * modelFileName;

private: 	
	TMOImage model;
	virtual void InitialiseHistogram(int *);
	virtual void PrintHistogram(int *, std::string);
	virtual void FillHistogram(TMOImage *, int *);
	virtual std::string GetModelFilename(std::string);
	virtual void ComputeComulativeHistogram(int *, int *);
	virtual int FindClosestShade(int, int *);
	virtual void NormaliseHistogram(int *, int);
	virtual void BilateralFilter(pfstmo::Array2D, pfstmo::Array2D, double, double);
	virtual double BilateralFilterWeight(double, double, double, double, double, double, pfstmo::Array2D);
	virtual void CreateGrayscale(pfstmo::Array2D, TMOImage *);
	virtual void GetDetailFromBase(pfstmo::Array2D, pfstmo::Array2D, pfstmo::Array2D);
};
