#include "TMO.h"
#include "TMOv.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>

#include <map>
#include <iostream>



/*#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"


#undef EPS*/

#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

#define LAMBDA 0.9
#define TAU 15
#define DELTA 10.0
#define SUBSAMPLE 32

//typedef dlib::matrix<double,0,1> column_vector;


class TMOTao17: public TMOv  
{
  private:
	
	TMODouble beta;
public:
  //static 
	TMOTao17();
	virtual ~TMOTao17();
	virtual int Transform();
	virtual int TransformVideo();
	
	void rgb2xyz(double *data);
	void xyz2lab(double *data);
	std::vector<float> getPixelDifferences(std::vector<cv::Vec3f> labVector, int pixelCount);
	
	
	virtual int getWeights(float &wr,float &wg, float &wb);
	virtual int getWeightsTemp(float &wr,float &wg, float &wb);
	
	virtual int getProximityValues(cv::Mat diffMat,std::vector<cv::Point2f> &proximityVals);
	virtual int getDifferentialRefinementFrame(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat &DFrame);
	virtual int getCoherenceRefinementFrame(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat previousGray, cv::Mat &CFrame);
	virtual int getOpticalFlow(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat &flow);
	
	virtual int getLowProximityFrame(cv::Mat currentFrame, cv::Mat previousFrame,cv::Mat previousGray, cv::Mat &currentGray);
	virtual int classifier(std::vector<cv::Point2f> &proximityVals);
	
	

protected:
	TMODouble dParameter;
};
