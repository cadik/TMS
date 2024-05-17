#include "TMO.h"
#include "stdio.h"
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include "cqp/gsl_cqp.h"
using namespace std;
class TMOMantiuk08 : public TMO
{
public:
	TMOMantiuk08();
	virtual ~TMOMantiuk08();
	virtual int Transform();

protected:
};
