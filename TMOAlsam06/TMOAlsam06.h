#include "TMO.h"

class TMOAlsam06 : public TMO  
{
public:
	TMOAlsam06();
	virtual ~TMOAlsam06();
	virtual int Transform();
	
private:
	TMODouble alpha;
	TMODouble beta;
	TMODouble gamma;

	virtual void RecalculateWeights();
};
