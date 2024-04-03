#include "TMO.h"

class TMONafchi17 : public TMO
{
public:
	TMONafchi17();
	virtual ~TMONafchi17();
	virtual int Transform();

protected:
	TMODouble dParameter;
	TMOBool bParameter;
};
