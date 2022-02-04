#include "TMO.h"

class TMOParis11 : public TMO  
{
public:
	TMOParis11();
	virtual ~TMOParis11();
	virtual int Transform();

protected:
	TMODouble alpha;
	TMODouble beta;
	TMODouble sigma_r;
	TMODouble gamma;
	TMOBool detailMnpl;
	TMOBool invToneMp;
};
