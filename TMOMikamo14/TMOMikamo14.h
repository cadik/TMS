#include "TMO.h"

class TMOMikamo14 : public TMO
{
public:
	struct LMS
	{
		std::vector<double> L;
		std::vector<double> M;
		std::vector<double> S;
	};
	TMOMikamo14();
	virtual ~TMOMikamo14();
	virtual int Transform();
	double computeAdaptedLuminance();
	std::vector<double> computeSigmoidParams(double I);
	std::vector<double> lambdaAdjust(int cone, double lambdaDiff);
	cv::Mat applyTwoStageModel(std::vector<double> spd, double I);
	std::vector<double> RGBtoSpectrum(double red, double green, double blue);

	// 10 bin spectral representation from 380nm to 720nm. The bins are all equal size. source: https://www.researchgate.net/profile/Shahab-Askarian/post/What_is_a_good_way_to_convert_a_RGB_pixel_to_a_wavelength/attachment/59d62ba6c49f478072e9d9a7/AS%3A273535619010561%401442227382042/download/An+RGB+to+Spectrum+Conversion+for+Reflectances.pdf
	const static int bins = 10;
	const static int colors = 7;
	const static int binWidth = 34;
	enum color
	{
		White,
		Cyan,
		Magenta,
		Yellow,
		Red,
		Green,
		Blue
	};

	double color_data[bins][colors] = {
		{1.0000, 0.9710, 1.0000, 0.0001, 0.1012, 0.0000, 1.0000}, // Bin 1
		{1.0000, 0.9426, 1.0000, 0.0000, 0.0515, 0.0000, 1.0000}, // Bin 2
		{0.9999, 1.0007, 0.9685, 0.1088, 0.0000, 0.0273, 0.8916}, // Bin 3
		{0.9993, 1.0007, 0.2229, 0.6651, 0.0000, 0.7937, 0.3323}, // Bin 4
		{0.9992, 1.0007, 0.0000, 1.0000, 0.0000, 1.0000, 0.0000}, // Bin 5
		{0.9998, 1.0007, 0.0458, 1.0000, 0.0000, 0.9418, 0.0000}, // Bin 6
		{1.0000, 0.1564, 0.8369, 0.9996, 0.8325, 0.1719, 0.0003}, // Bin 7
		{1.0000, 0.0000, 1.0000, 0.9586, 1.0149, 0.0000, 0.0369}, // Bin 8
		{1.0000, 0.0000, 1.0000, 0.9685, 1.0149, 0.0000, 0.0483}, // Bin 9
		{1.0000, 0.0000, 0.9959, 0.9840, 1.0149, 0.0025, 0.0496}  // Bin 10
	};

	double LMSsensitivities[bins][3] = {
		{1.48756E-03, 1.37367E-03, 3.39832E-02}, // Bin 1
		{2.93303E-02, 4.17546E-02, 8.24818E-01}, // Bin2
		{8.06894E-02, 1.44541E-01, 7.38268E-01}, // Bin3
		{2.76587E-01, 4.12260E-01, 1.29937E-01}, // Bin4
		{8.04084E-01, 9.56733E-01, 9.68268E-03}, // Bin5
		{9.97895E-01, 8.51173E-01, 3.75940E-04}, // Bin6
		{8.22859E-01, 3.19843E-01, 1.68551E-05}, // Bin7
		{3.27864E-01, 4.44879E-02, 0.0}, // Bin8
		{5.31896E-02, 3.69803E-03, 0.0}, // Bin9
		{4.74086E-03, 2.93607E-04, 0.0} // Bin10
	};

protected:
	TMODouble ri; // adapted retinal illuminance
	TMODouble al; // average luminance
};
