#include "TMO.h"

#include <vector>
#include <cmath>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#undef EPS

class TMOAubry14 : public TMO
{
public:
	TMOAubry14();
	virtual ~TMOAubry14();
	virtual int Transform();

protected:
	TMODouble sigmaParameter;
	TMOInt factParameter;
	TMOInt NParameter;
	TMOBool HDRParameter;
private:
	std::vector<double> linspace(double min, double max, int n);
};
