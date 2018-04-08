/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                                 Gamut mapping                       		   *
*                                                                              *
*                                   VYF project                                *
*            Author: Martin Kubicek [xkubic34 AT stud.fit.vutbr.cz]            *
*                                    Brno 2018                                 *
*                                                                              *
*******************************************************************************/

/* --------------------------------------------------------------------------- *
 * TMOKubicek18.cpp: implementation of the Kubicek18 class.                    *
 *                                 Gamut mapping                               *
 * --------------------------------------------------------------------------- */



#include "TMOKubicek18.h"
#include <fftw3.h>

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */

TMOKubicek18::TMOKubicek18()
{
	SetName(L"Kubicek18");
	SetDescription(L"Gamut mapping");

	method.SetName(L"Method");
	method.SetDescription(L"Method for mapping; 1:Perceptual, 2:Saturation, 3:AbsoluteColor, 4: RelativeColor"); 
	method.SetDefault(1);
	method = 2;		
	method.SetRange(1,4);	
	this->Register(method);

	/*
	Nastaveni toho jaky bude prevod 
	*/
}

TMOKubicek18::~TMOKubicek18()
{
}


/*
Spatne navrzeno predelat!!!!!
Rozdelit upravovani RGB, kazda hodnota bude mit pevne nadefinovane body, bude stacit jen posilat double hodnotu.
u HSV, CMYK a dalsich nezalezi na profilu viz poznamky

void PerceptualRGBToRGB(ColorSpace *RGB, ColorSpace *RGB){
//Sekce na výpočet koeficientu, pokud je druhý model větší první zůstavají koeficienty stejné.
Pokud je druhý menší vypočítá se jeho poměrná část 

	double koeficientRed = 1.0;
	double koeficientGreen = 1.0;
	double koeficientBlue = 1.0;

	if(RGB1.RedMax > RGB2.RedMax){
		koeficientRed = (RGB2.RedMax/RGB1.RedMax);
	}

	if(RGB1.GreenMax > RGB2.GreenMax){
		koeficientGreen = (RGB2.GreenMax/RGB1.GreenMax);
	}

	if(RGB1.BlueMax > RGB2.BlueMax){
		koeficientBlue = (RGB2.BlueMax/RGB1.BlueMax);
	}

	RGB2.Red = RGB1.Red * koeficientRed;
	RGB2.Green = RGB1.Green * koeficientGreen;
	RGB2.Blue = RGB1.Blue * koeficientBlue;
}

void PerceptualRGBToCMYK(ColorSpace *RGB, ColorSpace *CMYK){
// Prevod RGB na hodnotu CMYK a nasledne aplikovani prevodu 

	double tempCyan = (1 - RGB.RedMax / 255.0 );
	double tempMagenta = (1 - RGB.GreenMax / 255.0 );
	double tempYellow = (1 - RGB.BlueMax / 255.0 );

	double koeficientCyan = 1.0;
	double koeficientMagenta = 1.0;
	double koeficientYellow = 1.0;


	if(tempCyan > CMYK.CyanMax){
		koeficientCyan = (CMYK.CyanMax/tempCyan);
	}

	if(tempMagenta > CMYK.Magenta){
		koeficientMagenta = (CMYK.Magenta/tempMagenta);
	}

	if(tempYellow > CMYK.Yellow){
		koeficientYellow = (CMYK.Yellow/tempYellow);
	}

	CMYK.Cyan = (1 - RGB.Red / 255.0 ) * koeficientCyan;
	CMYK.Magenta = (1 - RGB.Green / 255.0 ) * koeficientMagenta;
	CMYK.Yellow = (1 - RGB.Blue / 255.0 ) * koeficientYellow;
}

void PerceptualRGBToXYZ(ColorSpace *RGB, ColorSpace *XYZ){
// Prevod RGB na hodnotu XYZ a nasledne aplikovani prevodu 

	double tempX = (RGB.RedMax/255.0)*0.4124564 + (RGB.GreenMax/255.0)*0.3575761 + (RGB.BlueMax/255.0)*0.1804375;
	double tempY = (RGB.RedMax/255.0)*0.2126729 + (RGB.GreenMax/255.0)*0.7151522 + (RGB.BlueMax/255.0)*0.0721750;
	double tempZ = (RGB.RedMax/255.0)*0.0193339 + (RGB.GreenMax/255.0)*0.1191920 + (RGB.BlueMax/255.0)*0.9503041;

	double koeficientX = 1.0;
	double koeficientY = 1.0;
	double koeficientZ = 1.0;


	if(tempX > XYZ.XMax){
		koeficientX = (XYZ.XMax/tempX);
	}

	if(tempY > XYZ.YMax){
		koeficientY = (XYZ.YMax/tempY);
	}

	if(tempZ > XYZ.ZMax){
		koeficientZ = (XYZ.ZMax/tempZ);
	}

	XYZ.X = ((RGB.Red/255.0)*0.4124564 + (RGB.Green/255.0)*0.3575761 + (RGB.Blue/255.0)*0.1804375) * koeficientX;
	XYZ.Y = ((RGB.Red/255.0)*0.2126729 + (RGB.Green/255.0)*0.7151522 + (RGB.Blue/255.0)*0.0721750) * koeficientY;
	XYZ.Z = ((RGB.Red/255.0)*0.0193339 + (RGB.Green/255.0)*0.1191920 + (RGB.Blue/255.0)*0.9503041) * koeficientZ;
}

void PerceptualRGBToHSV(ColorSpace *RGB, ColorSpace *HSV){
// TODO 
}

void PerceptualRGBToHSL(ColorSpace *RGB, ColorSpace *HSL){
// TODO 
}

void PerceptualRGBToLAB(ColorSpace *RGB, ColorSpace *LAB){
// TODO - převod na XYZ 
}


void SaturationRGBToRGB(ColorSpace *RGB1, ColorSpace *RGB2){
// Saturační model se mapuje do druhého. Je tedy nutne vypocitat vzdy koeficient 

	double koeficientRed = (RGB2.RedMax/RGB1.RedMax);
	double koeficientGreen = (RGB2.GreenMax/RGB1.GreenMax);
	double koeficientBlue = (RGB2.BlueMax/RGB1.BlueMax);

	RGB2.Red = RGB1.Red * koeficientRed;
	RGB2.Green = RGB1.Green * koeficientGreen;
	RGB2.Blue = RGB1.Blue * koeficientBlue;
}

void SaturationRGBToCMYK(ColorSpace *RGB, ColorSpace *CMYK){
// Saturační model RGB -> CMYK

	double koeficientCyan = (CMYK.CyanMax/(1 - RGB.RedMax / 255.0 ));
	double koeficientMagenta = (CMYK.Magenta/(1 - RGB.GreenMax / 255.0 ));
	double koeficientYellow = (CMYK.Yellow/(1 - RGB.BlueMax / 255.0 ));

	CMYK.Cyan = (1 - RGB.Red / 255.0 ) * koeficientCyan;
	CMYK.Magenta = (1 - RGB.Green / 255.0 ) * koeficientMagenta;
	CMYK.Yellow = (1 - RGB.Blue / 255.0 ) * koeficientYellow;
}

void SaturationRGBToXYZ(ColorSpace *RGB, ColorSpace *XYZ){
// Prevod RGB na hodnotu XYZ a nasledne aplikovani prevodu 

	double koeficientX = (XYZ.XMax/((RGB.RedMax/255.0)*0.4124564 + (RGB.GreenMax/255.0)*0.3575761 + (RGB.BlueMax/255.0)*0.1804375));
	double koeficientY = (XYZ.YMax/((RGB.RedMax/255.0)*0.2126729 + (RGB.GreenMax/255.0)*0.7151522 + (RGB.BlueMax/255.0)*0.0721750));
	double koeficientZ = (XYZ.ZMax/((RGB.RedMax/255.0)*0.0193339 + (RGB.GreenMax/255.0)*0.1191920 + (RGB.BlueMax/255.0)*0.9503041));

	XYZ.X = ((RGB.Red/255.0)*0.4124564 + (RGB.Green/255.0)*0.3575761 + (RGB.Blue/255.0)*0.1804375) * koeficientX;
	XYZ.Y = ((RGB.Red/255.0)*0.2126729 + (RGB.Green/255.0)*0.7151522 + (RGB.Blue/255.0)*0.0721750) * koeficientY;
	XYZ.Z = ((RGB.Red/255.0)*0.0193339 + (RGB.Green/255.0)*0.1191920 + (RGB.Blue/255.0)*0.9503041) * koeficientZ;
}

void SaturationRGBToHSV(ColorSpace *RGB, ColorSpace *HSV){
// TODO 
}

void SaturationRGBToHSL(ColorSpace *RGB, ColorSpace *HSL){
// TODO 
}

void SaturationRGBToLAB(ColorSpace *RGB, ColorSpace *LAB){
// TODO 
}

*/

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOKubicek18::Transform()
{
	double* data;
	int size = pSrc->getHeight()*pSrc->getWidth();
	double *rgb[] = {fftw_alloc_real(size),fftw_alloc_real(size),fftw_alloc_real(size)};
	// Pridat zjisteni profilu v podobne 

	data=pSrc->getData();
	for (int i = 0; i < size; ++i)
	{
		rgb[0][i] = *data++;
		rgb[1][i] = *data++;
		rgb[2][i] = *data++;
	}


	if(method == 1){
	/* Percepcni metoda */
			
		switch(gamutSrc){
			case 'RGB':
				switch(gamutSrc){

					case 'RGB':
						PerceptualRGBToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualRGBToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						PerceptualRGBToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						PerceptualRGBToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						PerceptualRGBToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						PerceptualRGBToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'CMYK':
				switch(gamutSrc){

					case 'RGB':
						PerceptualCMYKoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualCMYKToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						PerceptualCMYKToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						PerceptualCMYKToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						PerceptualCMYKToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						PerceptualCMYKToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'XYZ':
				switch(gamutSrc){

					case 'RGB':
						PerceptualXYZoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualXYZToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						PerceptualXYZToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						PerceptualXYZToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						PerceptualXYZToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						PerceptualXYZToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSV':
				switch(gamutSrc){

					case 'RGB':
						PerceptualHSVToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualHSVToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						PerceptualHSVToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						PerceptualHSVToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						PerceptualHSVToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						PerceptualHSVToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSL':
				switch(gamutSrc){

					case 'RGB':
						PerceptualHSLToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualHSLToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						PerceptualHSLToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						PerceptualHSLToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						PerceptualHSLToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						PerceptualHSLToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'LAB':
				switch(gamutSrc){

					case 'RGB':
						PerceptualLABToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualLABToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						PerceptualLABToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						PerceptualLABToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						PerceptualLABToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						PerceptualLABToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			default:
				//TODO error
		}

	}else if(method == 2){
	/* Saturacni metoda */

		switch(gamutSrc){
			case 'RGB':
				switch(gamuDst){

					case 'RGB':
						SaturationRGBToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						SaturationRGBToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						SaturationRGBToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						SaturationRGBToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						SaturationRGBToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						SaturationRGBToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'CMYK':
				switch(gamutDst){

					case 'RGB':
						SaturationCMYKoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						PerceptualCMYKToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						SaturationCMYKToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						SaturationCMYKToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						SaturationCMYKToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						SaturationCMYKToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'XYZ':
				switch(gamutDst){

					case 'RGB':
						SaturationXYZoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						SaturationXYZToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						SaturationXYZToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						SaturationXYZToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						SaturationXYZToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						SaturationXYZToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSV':
				switch(gamutDst){

					case 'RGB':
						SaturationHSVToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						SaturationHSVToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						SaturationHSVToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						SaturationHSVToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						SaturationHSVToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						SaturationHSVToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSL':
				switch(gamutDst){

					case 'RGB':
						SaturationHSLToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						SaturationHSLToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						SaturationHSLToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						SaturationHSLToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						SaturationHSLToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						SaturationHSLToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'LAB':
				switch(gamutDst){

					case 'RGB':
						SaturationLABToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						SaturationLABToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						SaturationLABToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						SaturationLABToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						SaturationLABToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						SaturationLABToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			default:
				//TODO error
		}

	}else if(method == 3){
	/* Absolute colorimetric metoda */

		switch(gamutSrc){
			case 'RGB':
				switch(gamuDst){

					case 'RGB':
						AbsoluteColorimetricRGBToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						AbsoluteColorimetricRGBToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						AbsoluteColorimetricRGBToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						AbsoluteColorimetricRGBToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						AbsoluteColorimetricRGBToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						AbsoluteColorimetricRGBToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'CMYK':
				switch(gamutDst){

					case 'RGB':
						AbsoluteColorimetricCMYKoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						AbsoluteColorimetricCMYKToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						AbsoluteColorimetricCMYKToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						AbsoluteColorimetricCMYKToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						AbsoluteColorimetricCMYKToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						AbsoluteColorimetricCMYKToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'XYZ':
				switch(gamutDst){

					case 'RGB':
						AbsoluteColorimetricXYZoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						AbsoluteColorimetricXYZToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						AbsoluteColorimetricXYZToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						AbsoluteColorimetricXYZToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						AbsoluteColorimetricXYZToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						AbsoluteColorimetricXYZToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSV':
				switch(gamutDst){

					case 'RGB':
						AbsoluteColorimetricHSVToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						AbsoluteColorimetricHSVToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						AbsoluteColorimetricHSVToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						AbsoluteColorimetricHSVToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						AbsoluteColorimetricHSVToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						AbsoluteColorimetricHSVToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSL':
				switch(gamutDst){

					case 'RGB':
						AbsoluteColorimetricHSLToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						AbsoluteColorimetricHSLToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						AbsoluteColorimetricHSLToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						AbsoluteColorimetricHSLToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						AbsoluteColorimetricHSLToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						AbsoluteColorimetricHSLToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'LAB':
				switch(gamutDst){

					case 'RGB':
						AbsoluteColorimetricLABToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						AbsoluteColorimetricLABToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						AbsoluteColorimetricLABToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						AbsoluteColorimetricLABToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						AbsoluteColorimetricLABToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						AbsoluteColorimetricLABToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			default:
				//TODO error
		}

	}else{
	/* Relative colorimetric metoda */

		switch(gamutSrc){
			case 'RGB':
				switch(gamuDst){

					case 'RGB':
						RelativeColorimetricRGBToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						RelativeColorimetricRGBToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						RelativeColorimetricRGBToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						RelativeColorimetricRGBToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						RelativeColorimetricRGBToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						RelativeColorimetricRGBToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'CMYK':
				switch(gamutDst){

					case 'RGB':
						RelativeColorimetricCMYKoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						RelativeColorimetricCMYKToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						RelativeColorimetricCMYKToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						RelativeColorimetricCMYKToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						RelativeColorimetricCMYKToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						RelativeColorimetricCMYKToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'XYZ':
				switch(gamutDst){

					case 'RGB':
						RelativeColorimetricXYZoRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						RelativeColorimetricXYZToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						RelativeColorimetricXYZToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						RelativeColorimetricXYZToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						RelativeColorimetricXYZToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						RelativeColorimetricXYZToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSV':
				switch(gamutDst){

					case 'RGB':
						RelativeColorimetricHSVToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						RelativeColorimetricHSVToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						RelativeColorimetricHSVToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						RelativeColorimetricHSVToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						RelativeColorimetricHSVToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						RelativeColorimetricHSVToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'HSL':
				switch(gamutDst){

					case 'RGB':
						RelativeColorimetricHSLToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						RelativeColorimetricHSLToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						RelativeColorimetricHSLToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						RelativeColorimetricHSLToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						RelativeColorimetricHSLToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						RelativeColorimetricHSLToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			case 'LAB':
				switch(gamutDst){

					case 'RGB':
						RelativeColorimetricLABToRGB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'CMYK':
						RelativeColorimetricLABToCMYK(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'XYZ':
						RelativeColorimetricLABToXYZ(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSV':
						RelativeColorimetricLABToHSV(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'HSL':
						RelativeColorimetricLABToHSL(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					case 'LAB':
						RelativeColorimetricLABToLAB(rgb[0][i],rgb[1][i],rgb[2][i]);
						break;
					default:
					//TODO error
				}
				break;
			default:
				//TODO error
		}

	}


	data=pDst->GetData();
	for (int i = 0; i < size; ++i)
	{
		*data++ = rgb[0][i];
		*data++ = rgb[1][i];
		*data++ = rgb[2][i];
	}

	pDst->Convert(TMO_RGB);

return 0;	
}