#include "TMO.h"
#include "TMOv.h"

class TMOAncuti11 : public TMOv
{
public:
	TMOAncuti11();
	virtual ~TMOAncuti11();
	virtual int Transform();
	virtual int TransformVideo();

	cv::Mat decolorization(cv::Mat &input, double eta);

protected:
	TMODouble dParameter;
};
