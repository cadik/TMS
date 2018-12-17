#include "TMO.h"

#include <vector>
#include <cmath>

#include <iostream>

// tmolib contains definition for EPS, so does openCV,
// so we undefine it and then put it back after includes
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
	TMOBool bFineP, bMediumP, bCoarseP;

	// TODO
	// TMODouble dFineVal0P, dFineVal1P, dFineVal2P;
	// TMODouble dMediumVal0P, dMediumVal1P, dMediumVal2P;
	// TMODouble dCoarseVal0P, dCoarseVal1P, dCoarseVal2P;

	TMODouble dFineExposureP, dFineGammaP, dFineSaturationP;
	TMODouble dMediumExposureP, dMediumGammaP, dMediumSaturationP;
	TMODouble dCoarseExposureP, dCoarseGammaP, dCoarseSaturationP;

private:
	cv::Mat sigmoid(cv::Mat x, double a);
	cv::Mat tonemapLAB(cv::Mat Lab, cv::Mat L0, cv::Mat L1,
											double val0, double val1, double val2,
											double exposure, double gamma, double saturation);
};
