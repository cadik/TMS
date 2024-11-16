#include "TMO.h"
#include <vector>
#include <memory>

class TMOYu21 : public TMO
{
public:
	TMOYu21();
	virtual ~TMOYu21();
	virtual int Transform();
protected:
	struct SImageStats
	{
		double covRG, covGB, covBR;
		double stddevR, stddevG, stddevB;
		double meanR, meanG, meanB;
		double Krg, Kgb, Kbr;
	};

	std::array<double, 3> computeK();

	SImageStats computeCorrelationCoefficient();
	std::unique_ptr<std::vector<double>> createContrastImage(const SImageStats &correlationCoefficients);
protected:
	TMODouble dParameter;
};
