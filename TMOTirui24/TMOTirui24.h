#include "TMO.h"
#include <numeric>

class TMOTirui24 : public TMO
{
public:
	TMOTirui24();
	virtual ~TMOTirui24();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
