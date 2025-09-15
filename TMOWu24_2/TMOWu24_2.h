#include "TMO.h"
#include <numeric>

class TMOWu24_2 : public TMO
{
public:
	TMOWu24_2();
	virtual ~TMOWu24_2();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
