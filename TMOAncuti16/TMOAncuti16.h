#include "TMO.h"
#include <stdio.h>
#include <stdlib.h>

class TMOAncuti16: public TMO  
{
public:
	TMOAncuti16();
	virtual ~TMOAncuti16();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
