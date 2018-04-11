#include "TMO.h"

class TMOZhongping15 : public TMO  
{
public:
	TMOZhongping15();
	virtual ~TMOZhongping15();
	virtual int Transform();

protected:
	TMODouble sigma;
};
