#include "TMO.h"

class TMOWanatMantiuk14 : public TMO
{
public:
	TMOWanatMantiuk14();
	virtual ~TMOWanatMantiuk14();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
