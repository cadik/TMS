#include "TMO.h"

#include <Eigen/Dense>

class TMOAlsam06 : public TMO  
{
public:
	TMOAlsam06();
	virtual ~TMOAlsam06();
	virtual int Transform();
private:
	TMOInt kernelSizeParameter;
	TMOInt standardDeviationParameter;
};
