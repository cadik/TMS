#include "TMO.h"

class TMOEisemann04 : public TMO  
{
public:
	TMOEisemann04();
	virtual ~TMOEisemann04();
	virtual int Transform();
protected:
	double* ComputeIntensity(TMOImage *inputImage, TMOImage *weightImage);
	TMOImage ComputeColor(TMOImage *inputImage, double *intensity);
	TMOString flashImagePathParameter;
	TMOImage flashImage;
};
