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

#undef EPS


class TMOAncuti16: public TMO  
{
public:
	TMOAncuti16();
	virtual ~TMOAncuti16();
	virtual int Transform();
	

protected:
	TMODouble dParameter;
};
