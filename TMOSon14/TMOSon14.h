/* -------------------------------------------------------------------------- *
 * TMOSon14.h: C++ implementation of "Art-Photographic Detail Enhancement"    *
 *             Minjung Son, Yunjin Lee, Henry Kang, Seungyong Lee             *
 *                        Computer Graphics Forum 2014                        *
 * Authors of code: Pavel Sedlar (2018)                                       *
 *                  - initial code                                            *
 *                  Tomas Hudziec (2019)                                      *
 *                  - added debug code                                        *
 *                  - replaced QP library qpOASES with OSQP                   *
 *                  - code cleanup                                            *
 * -------------------------------------------------------------------------- */
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
#undef EPS

#include <Eigen/Sparse>

#include "L0minimization.h"
#include "baseAndDetailDecomposition.h"
#include "osqpOptimization.h"
// #include "qpOasesOptimization.h"
// #include "pallasOptimization.h"
// #include "cgalOptimization.h"

class TMOSon14 : public TMO
{
public:
	TMOSon14();
	virtual ~TMOSon14();
	virtual int Transform();

protected:
	TMOInt optim1Iteration;
	TMODouble mu;
	TMOBool debugFlag;
};
