#include "TMO.h"
#include "wavelib.h"
#include <opencv2/opencv.hpp>

struct YCbCr
{
	double Y;
	double Cb;
	double Cr;
};

class TMOQueiroz06 : public TMO
{
public:
	TMOQueiroz06();
	virtual ~TMOQueiroz06();
	virtual int Transform();

	YCbCr RGBToYCbCr(double fr, double fg, double fb);

protected:
	TMODouble dParameter;
};
