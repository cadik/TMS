#include "TMO.h"
#include "stdio.h"
#include <math.h>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

class TMOYee03 : public TMO
{
public:
	TMOYee03();
	virtual ~TMOYee03();
	virtual int Transform();

protected:
	TMODouble bin_size1;
	TMODouble bin_size2;
	TMODouble small_threshold;
	TMODouble big_threshold;
	TMODouble max_layers;
};
