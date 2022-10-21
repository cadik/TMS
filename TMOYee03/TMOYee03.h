#include "TMO.h"

class TMOYee03 : public TMO
{
public:
	TMOYee03();
	virtual ~TMOYee03();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
