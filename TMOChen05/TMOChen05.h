#include "TMO.h"
#include "stdio.h"
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
using namespace std;

#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgproc/types_c.h"

class TMOChen05 : public TMO
{
public:
	TMOChen05();
	virtual ~TMOChen05();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
