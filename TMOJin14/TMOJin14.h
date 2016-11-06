#include "TMO.h"

class TMOJin14 : public TMO  
{ 
public:
	TMOJin14();
	virtual ~TMOJin14();
	virtual int Transform();
protected: 
	TMODouble tau;
	TMODouble mu;
	TMOBool rescale;
};
