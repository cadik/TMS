#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <algorithm>
#include <random>
#include <cmath>

using namespace cv;

class TMOLiu15 : public TMO
{
public:
	TMOLiu15();
	virtual ~TMOLiu15();
	virtual int Transform();
	Mat wei();

protected:
	TMODouble Lpp;
};
