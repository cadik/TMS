#include "TMO.h"
#include <Eigen/Dense>
#include <vector>

class TMOZhang08 : public TMO
{
public:
    TMOZhang08();
    virtual ~TMOZhang08();
    virtual int Transform();

private:
    TMODouble dParameter;
    TMODouble gamma;

    void stretchLabToFullRange(double* data, int width, int height);
    double gammaCorrection(double c);
    double fLab(double t);
    void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z);
    void XYZtoLAB(double X, double Y, double Z, double &L, double &a, double &b);
    void convertRGBtoLAB(double* data, int width, int height);
    double kernelK5(const Eigen::VectorXd& yi, const Eigen::VectorXd& yj);
    double kernelK6(const Eigen::VectorXd& yi, const Eigen::VectorXd& yj, int i, int j);
    std::vector<double> normalizeToLuminance(const Eigen::VectorXd& projections);
};
