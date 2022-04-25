#include "TMO.h"

class TMOGastal11 : public TMO
{
public:
	TMOGastal11();
	virtual ~TMOGastal11();
	virtual int Transform();

protected:
    TMOInt numIter;
    TMOInt filterType;
    TMODouble sigma_s;
    TMODouble sigma_r;
};
