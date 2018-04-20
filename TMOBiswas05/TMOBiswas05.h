#include "TMO.h"

class TMOBiswas05 : public TMO  
{
public:
	TMOBiswas05();
	virtual ~TMOBiswas05();
	virtual int Transform();
	
protected:
	TMODouble dParameter;
	TMOInt iParameter;
};
