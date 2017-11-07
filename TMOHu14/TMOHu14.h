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

#define LAMBDA 0.9
#define TAU 15
#define DELTA 10

struct lessVec3b
{
    bool operator()(const cv::Vec3b& lhs, const cv::Vec3b& rhs) const {
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
	void kmeansColorQuantization(const cv::Mat3b& src, cv::Mat3b& dst);
	 void getPalette(std::map<cv::Vec3d, float, lessVec3b>& paletteRGB, cv::Mat& src);
	 std::map<cv::Vec3d, float, lessVec3b> getDominantColorDescriptor(std::map<cv::Vec3d, float, lessVec3b> palette);
	 cv::Vec3d rgb2Luv(cv::Vec3i bgrVector);
	 cv::Vec3d Luv2rgb(cv::Vec3d luvVector);
	 cv::Vec3d xyz2bgr(cv::Vec3d xyzVector);
	 cv::Vec3d getBestWeightsCandidate(std::map<cv::Vec3d, float, lessVec3b> luvPalette,std::map<cv::Vec3d, float, lessVec3b> bgrPalette,cv::Mat redMat,cv::Mat greenMat ,cv::Mat blueMat);
	 std::map<int, float> getGrayscalePalette (float weight_r, float weight_g, float weight_b, std::map<cv::Vec3d, float, lessVec3b> bgrPalette);
	 double getXiMetric(std::map<cv::Vec3d, float, lessVec3b> luvPalette, std::map<int, float> grayscalePalette);
	 double getHMetric(std::map<cv::Vec3d, float, lessVec3b> luvPalette, float percentage1, float percentage2);
	

protected:
	TMODouble dParameter;
};
