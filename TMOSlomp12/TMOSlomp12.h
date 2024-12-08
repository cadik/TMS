#include "TMO.h"

class TMOSlomp12 : public TMO
{
public:
	TMOSlomp12();
	virtual ~TMOSlomp12();
	cv::Mat TMOImageToLogLuminanceMat(TMOImage *pSrc);
	void logLuminanceImage(cv::Mat srcMat);
	cv::Mat mipmap(cv::Mat srcMat, int levels);
	virtual int Transform();

protected:
	TMODouble dParameter;
};
