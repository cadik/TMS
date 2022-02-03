#include "TMO.h"

class TMODrago03 : public TMO
{
public:
	TMODrago03();
	virtual ~TMODrago03();
	virtual int Transform();

protected:
	TMODouble bias;
	TMODouble exposure;
	TMODouble gamma;
	TMODouble kernel;
	TMOBool gammaForm;
	TMOBool center;
	TMOInt centerX;
	TMOInt centerY;
};
