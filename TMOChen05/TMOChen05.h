#include "TMO.h"
#include "stdio.h"
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <map>
#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/imgproc/types_c.h"
#include "opencv2/imgcodecs.hpp"
#include <eigen3/Eigen/Dense>
using namespace std;




class TMOChen05 : public TMO
{
public:
	TMOChen05();
	virtual ~TMOChen05();
	virtual int Transform();

protected:
	TMODouble Theta;
	TMODouble Delta;
};
