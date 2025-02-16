#include "TMO.h"

class TMOThompson02 : public TMO
{
public:
	TMOThompson02();
	virtual ~TMOThompson02();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
