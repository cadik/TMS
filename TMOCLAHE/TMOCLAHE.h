#include "TMO.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "clahe.h"
#undef EPS

class TMOCLAHE : public TMO
{
public:
	TMOCLAHE();
	virtual ~TMOCLAHE();
	virtual int Transform();

protected:
	TMOInt gridRegions;
	TMODouble cl;
};
