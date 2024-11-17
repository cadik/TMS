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

	struct CImagePlusStats
	{
		std::unique_ptr<std::vector<double>> contrastPicture;
		double stddevC;
		double meanC;
	};

	std::array<double, 3> computeK(const SImageStats &imageStats);

	CImagePlusStats createContrastImage(const SImageStats &imageStatistics);
	SImageStats computeCorrelationCoefficient();
	std::array<double, 3> computeCovContrastRGB(const SImageStats &imageStatistics, const CImagePlusStats &contrastImageStat);
	std::array<double, 3> computeSSIM(const SImageStats &imageStatistics, const CImagePlusStats &contrastImageStat);


protected:
	TMODouble dParameter;
};
