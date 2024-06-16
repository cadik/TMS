#include "TMO.h"

class TMOYang18 : public TMO
{
public:
	TMOYang18();
	virtual ~TMOYang18();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
