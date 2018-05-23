#include "TMO.h"
#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#undef EPS

class TMOKuhn08 : public TMO  
{
public:
	TMOKuhn08();
	virtual ~TMOKuhn08();
	virtual int Transform();

private:
	cv::Mat convertToMat(TMOImage* image);
	void copy(TMOImage* in, TMOImage* out);
protected:
	TMOInt iK;
	TMODouble dTime;
	TMOInt iAttempts;
	TMOBool bColor;
	TMOBool bOriginal;
	TMOBool bInterpolate;
	TMOBool bChrominance;
};


