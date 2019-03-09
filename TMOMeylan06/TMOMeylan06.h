#include "TMO.h"

class TMOMeylan06 : public TMO
{
public:
	TMOMeylan06();
	virtual ~TMOMeylan06();
	virtual int Transform();

protected:
	TMODouble dParameter;
	cv::Mat RGBToPCA(cv::Mat rgb);
	cv::Mat PCAToRGB(cv::Mat pca);
};
