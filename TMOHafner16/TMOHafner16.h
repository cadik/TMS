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
    void computeLaplacian(const std::vector<double> weights, std::vector<double>& laplacian, int width, int height);
    void gradientDescentStep(std::vector<double>& wr, std::vector<double>& wg, std::vector<double>& wb,
                             const std::vector<double>& r, const std::vector<double>& g, const std::vector<double>& b,
                             int width, int height);
    void projectOntoSimplex(std::vector<double>& wr, std::vector<double>& wg, std::vector<double>& wb);
    double gaussianWeight(double x, double y);
    double psiDerivative(double z);

    // Constants for optimization
    static constexpr double lambda = 0.1;
    static constexpr double alpha = 0.05;
    static constexpr double gama = 0.25;
    static constexpr double delta = 1.0;
    static constexpr double tau = 0.15;
    static constexpr int numIter = 100;
    static constexpr double inv_2pi = 1.0 / (2.0 * 3.141592653589793);

    double sigma = 0.1;
    double mu = 0.5;
};
