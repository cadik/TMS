#include "TMO.h"

class TMOQueiroz06 : public TMO
{
public:
	TMOQueiroz06();
	virtual ~TMOQueiroz06();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
