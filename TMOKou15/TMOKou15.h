#include "TMO.h"
#include <algorithm>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "L0minimization.h"

#include <Eigen/Sparse>
#undef EPS

class TMOKou15 : public TMO  
{
public:
	TMOKou15();
	virtual ~TMOKou15();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
