#include "TMO.h"

class TMOLu14 : public TMO
{
public:
	TMOLu14();
	virtual ~TMOLu14();
	virtual int Transform();
	cv::Mat calculatePolyGrad(const cv::Mat& img, int order);
protected:
	cv::Mat powElement(const cv::Mat& base, int exp);
	TMODouble sigmaParameter;
};
