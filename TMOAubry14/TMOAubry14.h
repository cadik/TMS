/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*             Fast Local Laplacian Filters: Theory and Applications (2014)     *
* by Mathieu Aubry, Sylvain Paris, Samuel W. Hasinoff, Jan Kautz, Fredo Durand *
*                         ACM Transactions on Graphics                         *
*                                                                              *
*             Author: Tomas Hudziec [xhudzi01 AT stud.fit.vutbr.cz]            *
*         Term project for Computational Photography course - 2018             *
*       Part of master thesis (HDR support, code reorganization) - 2019        *
*                                                                              *
*******************************************************************************/

#include "TMO.h"

#include <vector>
#include <cmath>

#ifdef EPS
#undef EPS
#define EPS EPS2
#endif
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#undef EPS

#include "FastLocalLaplFilt.h"

class TMOAubry14 : public TMO
{
public:
	TMOAubry14();
	virtual ~TMOAubry14();
	virtual int Transform();

protected:
	TMODouble sigmaParameter;
	TMOInt factParameter;
	TMOInt NParameter;
	TMOBool HDRParameter;
};
