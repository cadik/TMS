#include "TMO.h"

class TMOAlsam20 : public TMO
{
public:
	TMOAlsam20();
	virtual ~TMOAlsam20();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
