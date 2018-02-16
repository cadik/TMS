#include "TMO.h"

class TMOBiswas05 : public TMO  
{
public:
	TMOBiswas05();
	virtual ~TMOBiswas05();
	virtual int Transform();
	
private:
	//should be always odd
	const int MEDIAN_DIMENSION = 9;

protected:
	TMODouble dParameter;
};
