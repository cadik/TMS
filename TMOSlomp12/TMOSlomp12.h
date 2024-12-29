#include "TMO.h"

class TMOSlomp12 : public TMO
{
public:
	TMOSlomp12();
	virtual ~TMOSlomp12();
	cv::Mat TMOImageToLogLuminanceMat();
	void logLuminanceImage(cv::Mat luminanceMat);
	cv::Mat mipmap(cv::Mat mat, int levels);
	void scaleLuminance(cv::Mat *luminanceMat, double keyValue);
	void scaledLuminanceImage(cv::Mat luminanceMat);
	double boxFilter(cv::Mat *SAT, int x, int y, int s);
	double getNormalizedDifference(double conv0, double conv1, int s);
	int getMaxScale(cv::Mat *SAT, int x, int y);
	double redResponseValue(double illuminance);
	double arithLuminanceAverage();
	virtual int Transform();

protected:
	TMOBool local;
	TMOBool mesopic;
	TMOBool varying;

	double alpha = 0.18;
	double phi = 8.;
	double epsilon = 0.025;
};
