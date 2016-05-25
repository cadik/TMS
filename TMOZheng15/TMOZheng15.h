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
#include <limits>
#include <math.h>

//#include "../matrix.h"							// matrix library

// constant c is used to ensure that the largest chromatic difference will not 
// be completely scaled down
#define CONSTANT_C 2.0

class TMOZheng15 : public TMO  
{	
public:
	TMOZheng15();
	virtual ~TMOZheng15();
	virtual int Transform();
	
	double AttenuationFunction(double);
	double ModulatedColorDifference(double, double, double);
	double Sign(double, double, double, double, double);
	double ChromaticGradientComponent(double, double, double, double);
	double GradientFiledComponent(double, double, double, double, double, double);	
	
private: 
	double x_max;
	TMODouble alpha;
	TMODouble beta;
	TMODouble gamma;
	TMODouble theta;	
	TMOBool verbose;
};

