/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*             Author: Ludmila Krejcova [xkrejc85 AT stud.fit.vutbr.cz]         *
*                                    Brno 2025                                 *
*                                                                              *
*                 Implementation of the TMOAmbalathankandy21 class             *
*                                                                              *
*               A color temperature-based high-speed decolorization:           *
*               an empirical approach for tone mapping applications            * 	
*                    https://doi.org/10.48550/arXiv.2108.13656                 *
*                                                                              *
*******************************************************************************/


#include <TMO.h>

class TMOAmbalathankandy21 : public TMO
{
public:
	TMOAmbalathankandy21();
	
	
	virtual ~TMOAmbalathankandy21();
	virtual int Transform() override;

	// Normalizes the grayscale image (not part of the original article, added for compariso
	void normalizeGrayscaleImage(TMOImage &image);

	// Computes the minimum and maximum values in the image for normalization
	std::pair<double, double> getImageMinMax(TMOImage &image);

	// Finds if range is 0-1 or in 0-255
	bool isInRange0to1(double *pSourceData, int numPix);

protected:
	TMOBool HDRParameter;
};
