#include "TMO.h"

class TMOSlomp12 : public TMO
{
public:
	TMOSlomp12();
	virtual ~TMOSlomp12();
	cv::Mat TMOImageToLogLuminanceMat(TMOImage *pSrc);
	virtual int Transform();

protected:
	TMODouble dParameter;
};
