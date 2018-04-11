#include "TMO.h"

class TMOXiong17 : public TMO  
{
public:
	TMOXiong17();
	virtual ~TMOXiong17();
	virtual int Transform();

protected:
	TMOInt maximum_of_iterations;
	TMOInt order;
};
