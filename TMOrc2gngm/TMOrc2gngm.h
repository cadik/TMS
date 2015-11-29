#ifndef TMORC2GNGM_H
#define TMORC2GNGM_H

#include "TMO.h"
#include <iostream>
#include <math.h>

// degree of freedom
#define N 4

// TODO presunout do private?

class TMOrc2gngm : public TMO  
{
private:
	TMODouble alpha;
protected:	
	double Gradient(double*, double*);	
	double FunctionF(double theta);
	double Hk_effect_predictor(double *);
	void PixelLabToLuv(double*, double*);
	void PixelLchToLab(double*, double*);
public:
	TMOrc2gngm();
	virtual ~TMOrc2gngm();
	virtual int Transform();
};

#endif
