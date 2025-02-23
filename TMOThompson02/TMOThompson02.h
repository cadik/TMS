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
	cv::Mat TMO2mat();
	void mat2TMO(cv::Mat &input);
	cv::Mat applyGaussianBlur(cv::Mat &input, double sigma);
	cv::Mat applyNightFilter(cv::Mat &input);
	cv::Mat addGaussianNoise(cv::Mat &input);
	virtual int Transform();

protected:
	TMODouble sigmaBlur = 1.6;
	TMODouble gammaEdge = 1.25;
	TMODouble sigmaNoise = 0.0125;
	std::vector<double> bluishGreyRGB = {1.05, 0.97, 1.27};
};
