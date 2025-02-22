#include "TMO.h"

class TMOThompson02 : public TMO
{
public:
	typedef struct matAndDouble
	{
		cv::Mat mat;
		double d;
	} matAndDouble;

	TMOThompson02();
	virtual ~TMOThompson02();
	matAndDouble getLuminanceMat();
	double getDisplayableLuminance(double L, double maxLuminance);
	void mapLuminance(cv::Mat &luminanceMat, double maxLuminance);
	cv::Mat getScotopicLuminanceMat();
	double getMesopicFactor(double L);
	virtual int Transform();

protected:
	TMODouble dParameter;
	std::vector<double> bluishGreyRGB = {1.05, 0.97, 1.27};
};
