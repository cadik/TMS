#include "TMO.h"

class TMOMikamo14 : public TMO
{
public:
	TMOMikamo14();
	virtual ~TMOMikamo14();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
