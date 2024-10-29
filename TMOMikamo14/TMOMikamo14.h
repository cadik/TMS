#include "TMO.h"

class TMOMikamo14 : public TMO
{
public:
	struct OpponentColor
	{
    	double achromatic;
    	double redGreen;
    	double yellowBlue;
	};
	TMOMikamo14();
	virtual ~TMOMikamo14();
	virtual int Transform();
	double computeAdaptedLuminance(int x, int y);
	std::vector<double> computeSigmoidParams(double I);
	double getIntensityFromWavelength(std::vector<double> spectrum, double lambda);
	double getWavelengthFromIntensity(std::vector<double> spectrum, double intensity);
	OpponentColor applyTwoStageModel(std::vector<double> spectrum, double I);
	std::vector<double> RGBtoSpectrum(int x, int y);


	// 10 bin spectral representation from 380nm to 720nm. The bins are all equal size. source: https://www.researchgate.net/profile/Shahab-Askarian/post/What_is_a_good_way_to_convert_a_RGB_pixel_to_a_wavelength/attachment/59d62ba6c49f478072e9d9a7/AS%3A273535619010561%401442227382042/download/An+RGB+to+Spectrum+Conversion+for+Reflectances.pdf
	const static int bins = 10;
	const static int colors = 7;
	enum color { White, Cyan, Magenta, Yellow, Red, Green, Blue };

	double color_data[bins][colors] = {
        {1.0000, 0.9710, 1.0000, 0.0001, 0.1012, 0.0000, 1.0000},  // Bin 1
        {1.0000, 0.9426, 1.0000, 0.0000, 0.0515, 0.0000, 1.0000},  // Bin 2
        {0.9999, 1.0007, 0.9685, 0.1088, 0.0000, 0.0273, 0.8916},  // Bin 3
        {0.9993, 1.0007, 0.2229, 0.6651, 0.0000, 0.7937, 0.3323},  // Bin 4
        {0.9992, 1.0007, 0.0000, 1.0000, 0.0000, 1.0000, 0.0000},  // Bin 5
        {0.9998, 1.0007, 0.0458, 1.0000, 0.0000, 0.9418, 0.0000},  // Bin 6
        {1.0000, 0.1564, 0.8369, 0.9996, 0.8325, 0.1719, 0.0003},  // Bin 7
        {1.0000, 0.0000, 1.0000, 0.9586, 1.0149, 0.0000, 0.0369},  // Bin 8
        {1.0000, 0.0000, 1.0000, 0.9685, 1.0149, 0.0000, 0.0483},  // Bin 9
        {1.0000, 0.0000, 0.9959, 0.9840, 1.0149, 0.0025, 0.0496}   // Bin 10
    };

protected:
	TMODouble ri;
};
