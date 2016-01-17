#include "TMO.h"
#include <limits>
#include <math.h>

//#include "../matrix.h"								// matrix library

class TMOEfficientC2GConversionForDIInGD : public TMO  
{	
public:
	TMOEfficientC2GConversionForDIInGD();
	virtual ~TMOEfficientC2GConversionForDIInGD();
	virtual int Transform();
	// TODO presunout do private
	double AttenuationFunction(double);
	double ModulatedColorDifference(double, double, double);
	double Sign(double, double, double, double, double);
	double Chromatic_gradient_component(double, double, double, double);
	double Gradient_filed_component(double, double, double, double, double, double);
	//void Get_grayscale_field_simple(Matrix, Matrix);
	
private: 
	double x_max;
	TMODouble alpha;
	TMODouble beta;
	TMODouble gamma;
	TMODouble theta;	
};

