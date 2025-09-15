#include "TMO.h"

class TMOWu24 : public TMO
{
public:
	TMOWu24();
	virtual ~TMOWu24();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
