#include "TMO.h"

class LuminanceLUT {
private:
    std::vector<std::pair<double, int>> lut;

public:
    LuminanceLUT(int bins, int truncation, const std::vector<double>& luminance,
                 int lum_len, double min_luminance, double max_luminance);

    [[nodiscard]] std::tuple<double, int, double, int> getValue(double key) const;

    [[maybe_unused]] void printLUT() const;
};

class TMOKhan20 : public TMO
{
public:
	TMOKhan20();
	virtual ~TMOKhan20();
	virtual int Transform();

protected:
    TMOInt binParameter;
    TMOInt truncationParameter;

    static void ToneMap(std::vector<double> &luminance, int lum_len, LuminanceLUT &lut);
    std::tuple<double, double> ApplyPerceptualQuantizer(std::vector<double> &luminance);
    void ApplyTransformationToDstImage(std::vector<double> luminance);
};
