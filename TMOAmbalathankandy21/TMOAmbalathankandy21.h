/* --------------------------------------------------------------------------- *
 * TMOAmbalathankandy.h														   *
 *																			   *
 * Author: Ludmila Krejčová													   *
 * --------------------------------------------------------------------------- */


#include <TMO.h>

class TMOAmbalathankandy21 : public TMO
{
public:
	TMOAmbalathankandy21();
	
	
	virtual ~TMOAmbalathankandy21();
	virtual int Transform() override;
	void normalizeGrayscaleImage(TMOImage &image);
	std::pair<double, double> getImageMinMax(TMOImage &image);

protected:
	TMODouble dParameter;
};
