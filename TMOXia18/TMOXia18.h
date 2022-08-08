#include "TMO.h"

class TMOXia18 : public TMO
{
public:
	TMOXia18();
	virtual ~TMOXia18();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
