#include "TMO.h"
#include <limits>

class TMOEfficientC2GConversionForDIInGD : public TMO  
{
public:
	TMOEfficientC2GConversionForDIInGD();
	virtual ~TMOEfficientC2GConversionForDIInGD();
	virtual int Transform();
	double AttenuationFunction(double);
	double ModulatedColorDifference(double, double, double);
	double Sign(double, double, double, double, double);
	double Chromatic_gradient_component(double, double, double, double);
	double Gradient_filed_component(double, double, double, double, double, double);
	
	
private: 
	double x_max;
	TMODouble alpha;
	TMODouble beta;
	TMODouble gamma;
	TMODouble theta;
};
