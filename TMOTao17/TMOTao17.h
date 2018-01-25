#include "TMO.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>

#include <map>
#include <iostream>



#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#undef EPS

#include <dlib/optimization.h>

#define LAMBDA 0.9
#define TAU 15
#define DELTA 10.0

typedef dlib::matrix<double,0,1> column_vector;

class TMOTao17: public TMO  
{
public:
	TMOTao17();
	virtual ~TMOTao17();
	virtual int Transform();
	static double rosen (const column_vector& m);
	void rgb2xyz(double *data);
	void xyz2lab(double *data);
	std::vector<double> getPixelDifferences(std::vector<cv::Vec3d> labVector, int pixelCount);
	
	

protected:
	TMODouble dParameter;
};
