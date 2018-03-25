#include "TMO.h"

class TMOAlsam06 : public TMO  
{
public:
	TMOAlsam06();
	virtual ~TMOAlsam06();
	virtual int Transform();
	
private:
	//should be always odd
	const int MEDIAN_DIMENSION = 9;

protected:
	TMODouble dParameter;
};
