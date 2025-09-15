#include "TMO.h"
#include <cmath>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/ximgproc.hpp>

class TMOSong14 : public TMO
{
public:
	TMOSong14();
	virtual ~TMOSong14();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
