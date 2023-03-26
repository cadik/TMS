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
class TMOMantiuk08 : public TMO
{
public:
	TMOMantiuk08();
	virtual ~TMOMantiuk08();
	virtual int Transform();

protected:
	TMODouble dParameter;
};
