#include "TMO.h"

class TMOSlomp12 : public TMO
{
public:
	TMOSlomp12();
	virtual ~TMOSlomp12();
	cv::Mat TMOImageToLogLuminanceMat();
	void logLuminanceImage(cv::Mat srcMat);
	cv::Mat mipmap(cv::Mat srcMat, int levels);
	void scaleLuminance(cv::Mat *srcMat, double keyValue, double alpha);
	void scaledLuminanceImage(cv::Mat srcMat);
	double redResponseValue(double illuminance);
	double arithLuminanceAverage();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
