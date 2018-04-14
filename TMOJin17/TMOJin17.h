#include "TMO.h"

class TMOJin17: public TMO  
{
public:
	TMOJin17();
	virtual ~TMOJin17();
	virtual int Transform();

protected:
	TMODouble sigma;
};
