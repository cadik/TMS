//(c)Martin Cadik
//see [tumblin99two], pp. 73
//11/2004

#include "TMO.h"
#include <stdio.h>
#include <math.h>
#include <string>


class TMOWard94 : public TMO  
{
public:
	TMOWard94();
	virtual ~TMOWard94();
	virtual int Transform();

protected:
	TMODouble Ld_max;
	TMODouble gamma; //gamma correction factor
};
