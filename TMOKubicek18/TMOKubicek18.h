#include "TMO.h"

class TMOKubicek18 : public TMO  
{
public:
	TMOKubicek18();
	virtual ~TMOKubicek18();
	virtual int Transform();

protected:
	TMOInt order;
};