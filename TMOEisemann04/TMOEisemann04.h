#include "TMO.h"

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif

#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <cmath>
#include <vector>
using namespace std;
using namespace cv;

class TMOEisemann04 : public TMO  
{
public:
	TMOEisemann04();
	virtual ~TMOEisemann04();
	virtual int Transform();
protected:
	TMOString flashImagePathParameter;
	TMOBool shadowCorrectionParameter;
	TMOImage flashImage;
};
