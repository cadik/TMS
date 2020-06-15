/**
 * A Two-Stage Contrast Enhancement Algorithm for Digital Images
 * https://doi.org/10.1109/CISP.2008.400
 *
 * Author: Michal Vlnas
 */

#include "TMO.h"

enum Predictor : uint8_t
{
    BETA_UNDEFINED = 0,
    BETA_DEFINED   = 1,
};

class TMOTai08 : public TMO  
{
public:
	TMOTai08();
	~TMOTai08() override = default;
	int Transform() override;

private:
    double _ApproxKneeCurve(double value);
    double _ApproxKneeCurveDifferential(double value);

    double GetAlpha(int x, int y, double value, double avg);
    void BuildBetaMap();

    double GetMeanAroundPixel(int x, int y, double value);
    double GetBetaAroundPixel(int x, int y, double value);
    double NL(double value);

    double t;
    double k;
    double m;

    double a;
    double b;
    double c;
    double d;

    double gmax;

    cv::Mat lum_avg;
    cv::Mat grad;

    cv::Mat beta;
    cv::Mat predictor;

protected:
	TMODouble tc;
	TMODouble kc;
	TMODouble gamma;
	TMODouble th;
};
