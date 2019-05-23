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
// @eta :   Parameter controlling how many times of detail adding to the
//          input image
// @lambda: Smoothing parameter controlling the degree of enhance. 
//          Typically it is within the range [1e-3, 1e-1], 2e-2 by default.
// @kappa : Parameter that controls the rate. 
//          Small kappa results in more iteratioins and with sharper edges.   
//          We select kappa in (1, 2].    
//          kappa = 2 is suggested for natural images.  
	TMODouble etaParameter;
	TMODouble lambdaParameter;
	TMODouble kappaParameter;
};
