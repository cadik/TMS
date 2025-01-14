#include "TMO.h"

class TMOCyriac15 : public TMO
{
public:
	TMOCyriac15();
	virtual ~TMOCyriac15();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
