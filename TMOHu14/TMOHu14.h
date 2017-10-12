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


class TMOHu14: public TMO  
{
public:
	TMOHu14();
	virtual ~TMOHu14();
	virtual int Transform();
	virtual cv::Mat getEdgeMat(cv::Mat channel);
	

protected:
	TMODouble dParameter;
};
