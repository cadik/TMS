#include "TMO.h"

class TMOAlsam09 : public TMO  
{
public:
	TMOAlsam09();
	virtual ~TMOAlsam09();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
