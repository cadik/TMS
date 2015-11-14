#include "TMO.h"

class TMOrc2gngm : public TMO  
{
public:
	TMOrc2gngm();
	virtual ~TMOrc2gngm();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
