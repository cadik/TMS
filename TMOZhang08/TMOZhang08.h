/********************************************************************************
*                                                                               *
*                         Brno University of Technology                         *
*                       Faculty of Information Technology                       *
*                                                                               *
*                         Color-to-Grayscale Conversions                        *
*                                                                               *
*             Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]          *
*                                    Brno 2025                                  *
*                                                                               *
*                     Implementation of the TMOZhang08 class                    *
*           A Kernel Based Algorithm for Fast Color-To-Gray Processing          *
*                      https://doi.org/10.1109/CISP.2008.411                    *
*                                                                               *
********************************************************************************/

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

    /**
    * Stretches LAB values to the full range [0,1] for L, A and B
    */
    void stretchLabToFullRange(double* data, int width, int height);

    /**
    * Gamma correction for sRGB → linear RGB
    */
    double gammaCorrection(double c);

    /**
    * Helper function for LAB transformation
    */
    double fLab(double t);

    /**
    * Converts sRGB to XYZ (D65)
    */
    void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z);

    /**
    * Converts XYZ to LAB into [0, 1] range
    */
    void XYZtoLAB(double X, double Y, double Z, double &L, double &a, double &b);

    /**
    * Converts an RGB image to LAB - output LAB is in range [0,1]
    */
    void convertRGBtoLAB(double* data, int width, int height);

    /**
    * Kernel function K5
    */
    double kernelK5(const Eigen::VectorXd& yi, const Eigen::VectorXd& yj);

    /**
    * Function for better kernel K6
    */
    double kernelK6(const Eigen::VectorXd& yi, const Eigen::VectorXd& yj, int i, int j);

    /**
    * Maps projection values to luminance range (0-1)
    */
    std::vector<double> normalizeToLuminance(const Eigen::VectorXd& projections);

    /*
	* Finds if range is 0-1 or in 0-255
	*/
	bool isInRange0to1(double *pSourceData, int numPix);

protected:
	TMOBool HDRParameter;
};
