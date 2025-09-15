/*******************************************************************************
*                                                                              *
*                        Brno University of Technology                         *
*                      Faculty of Information Technology                       *
*                                                                              *
*                        Color-to-Grayscale Conversions                        *
*                                                                              *
*            Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]          *
*                                   Brno 2025                                  *
*                                                                              *
*                     Implementation of the TMOHafner16 class                  *
*             Variational Image Fusion with Optimal Local Contrast             *
*                 https://doi.org/10.1007/978-3-319-18461-6_34                 *
*                                                                              *
*******************************************************************************/

#include "TMO.h"
#include <vector>

class TMOHafner16 : public TMO
{
public:
TMOHafner16();
	virtual ~TMOHafner16();
	virtual int Transform();

protected:
	TMODouble dParameter;

private:
    // Perform one iteration of gradient descent update
    void gradientDescentStep(std::vector<double>& wr, std::vector<double>& wg, std::vector<double>& wb,
                             const std::vector<double>& r, const std::vector<double>& g, const std::vector<double>& b,
                             int width, int height);
    
    //Projects the weights onto a simplex constraint - algorithm 2 
    void projectOntoSimplex(std::vector<double>& wr, std::vector<double>& wg, std::vector<double>& wb);

    // Gaussian filter using recursive approximation
    double gaussianWeight(double x, double y);

    // Approximated psiDerivative using polynomial expansion
    double psiDerivative(double z);

    // Computing laplacian
    void computeLaplacian(const std::vector<double> weights, std::vector<double>& laplacian, int width, int height);

    // Finds if range is 0-1 or in 0-255
	bool isInRange0to1(double *pSourceData, int numPix);
         

    // Constants for optimization according to the article
    static constexpr double lambda = 0.1;
    static constexpr double alpha = 0.05;
    static constexpr double gama = 0.25;
    static constexpr double delta = 1.0;
    static constexpr double tau = 0.2;
    static constexpr int numIter = 100;
    static constexpr double inv_2pi = 1.0 / (2.0 * 3.141592653589793);

    double sigma = 0.1;
    double mu = 0.5;
   
protected:
	TMOBool HDRParameter;
};
