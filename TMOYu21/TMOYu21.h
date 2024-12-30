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
	std::shared_ptr<std::vector<double>> computeContrastDifferences(const std::vector<std::pair<int, int>> &pairs, const double* image64, const double* image32, int channel);
	void logDebug(const std::string& message);
	std::array<double, 3> computeWeights(const std::vector<double> &allIr, const std::vector<double> &allIg,
				const std::vector<double> &allIb, const std::array<double, 3> &kr_kg_kb);

	double 	computeColorEnergy(const std::array<double, 3> &w, double k, const std::array<std::vector<double>, 3> &I, size_t colorIndex);
	double 	computeColorEnergy2(const std::array<double, 3> &w, const std::array<double, 3> &k, const std::array<std::vector<double>, 3> &I);

	std::vector<std::pair<int, int>> findRandomPairs(int size) const;

	static void normalizeGrayscaleImage(TMOImage &image);

	static std::unique_ptr<TMOImage> createImage(const double *data, int width, int height);
	static std::unique_ptr<TMOImage> createImageFromIntenzities(const double *data, int width, int height);

protected:
	TMODouble dParameter;
};
