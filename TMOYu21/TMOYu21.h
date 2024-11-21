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
	inline double getPixel(const double* data, int width, int x, int y, int channel);
	inline void setPixel(double* data, int width, int x, int y, int channel, double value) ;
	std::unique_ptr<double[]> resizeImage(const double* input, int srcWidth, int srcHeight, int destWidth, int destHeight);
	std::vector<double> computeContrastDifferences(const double* image64, const double* image32, int channel);



protected:
	TMODouble dParameter;
};
