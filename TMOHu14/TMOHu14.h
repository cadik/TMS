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

struct lessVec3b
{
    bool operator()(const cv::Vec3f& lhs, const cv::Vec3f& rhs) const {
        return (lhs[0] != rhs[0]) ? (lhs[0] < rhs[0]) : ((lhs[1] != rhs[1]) ? (lhs[1] < rhs[1]) : (lhs[2] < rhs[2]));
    }
};

class TMOHu14: public TMO  
{
public:
	TMOHu14();
	virtual ~TMOHu14();
	virtual int Transform();
	virtual cv::Mat getEdgeMat(cv::Mat channel);
	 std::map<cv::Vec3f, int, lessVec3b> getPalette(const cv::Mat& src);
	

protected:
	TMODouble dParameter;
};
