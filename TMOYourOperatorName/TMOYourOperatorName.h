#include "TMO.h"

class TMOYourOperatorName : public TMO  
{
public:
	TMOYourOperatorName();
	virtual ~TMOYourOperatorName();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
