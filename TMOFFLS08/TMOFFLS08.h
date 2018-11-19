#include "TMO.h"

#include <vector>
#include <cmath>

#include <iostream>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "opencv2/ximgproc/edge_filter.hpp"

#undef EPS

class TMOFFLS08 : public TMO
{
public:
	TMOFFLS08();
	virtual ~TMOFFLS08();
	virtual int Transform();

protected:
	TMODouble sigmaParameter;
	TMODouble lambdaParameter;
	TMODouble multiplyParameter;
};
