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

#ifndef TMORC2GNGM_H
#define TMORC2GNGM_H

#include "TMO.h"
#include <iostream>
#include <math.h>

// degree of freedom
#define N 4

// h in lch color space of white color
#define LCH_H_WHITE 296.812926236627

// lab to luv conversion rate
#define LAB_TO_LUV 2.55

class TMOKim09 : public TMO  
{
private:
	TMODouble alpha;
	TMOBool verbose;
protected:	
	double Gradient(double*, double*);	
	double FunctionF(double theta);
	double HkEffectPredictor(double *);
	void PixelLabToLuv(double*, double*);
	void PixelLchToLab(double*, double*);
public:
	TMOKim09();
	virtual ~TMOKim09();
	virtual int Transform();
};

#endif
