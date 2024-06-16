#include "TMO.h"

class TMOKirk11 : public TMO
{
public:
	TMOKirk11();
	virtual ~TMOKirk11();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
