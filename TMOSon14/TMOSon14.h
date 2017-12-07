#include "TMO.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
// #include "L0minimization.h"

#include <Eigen/Sparse>
#undef EPS


class TMOSon14: public TMO  
{
public:
	TMOSon14();
	virtual ~TMOSon14();
	virtual int Transform();
	

protected:
	TMODouble kappa;
};
