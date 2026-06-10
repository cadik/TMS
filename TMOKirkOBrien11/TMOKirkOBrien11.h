#include "TMO.h"

class TMOKirkOBrien11 : public TMO
{
public:
	TMOKirkOBrien11();
	virtual ~TMOKirkOBrien11();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
