#include "TMO.h"

class TMONguy15 : public TMO  
{
public:
	TMONguy15();
	virtual ~TMONguy15();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
