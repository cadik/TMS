#include "TMO.h"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

class TMOMeylan06 : public TMO
{
public:
	TMOMeylan06();
	virtual ~TMOMeylan06();
	virtual int Transform();

protected:
	TMODouble dParameter;
	int numberOfPixels;
	cv::PCA pca;
	cv::Mat RGBToPCA(double *rgbSourceData);
	cv::Mat PCAToRGB(cv::Mat &PCAProjection);
};
