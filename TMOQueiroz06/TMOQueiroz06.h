/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*       Color to Gray and Back: Color Embedding Into Textured Gray Images      *
*                                                                              *
*             Author: Peter Zdravecký [xzdrav00 AT stud.fit.vutbr.cz]          *
*                                    Brno 2024                                 *
*                                                                              *
*******************************************************************************/
#include "TMO.h"
#include "wavelib.h"
#include <opencv2/opencv.hpp>

/**
 * @brief Structure for YCbCr color space
 * 
 */
struct YCbCr
{
	double Y;
	double Cb;
	double Cr;
};

class TMOQueiroz06 : public TMO
{
public:
	TMOQueiroz06();
	virtual ~TMOQueiroz06();
	virtual int Transform();

	YCbCr RGBToYCbCr(double fr, double fg, double fb);

protected:
	TMODouble dParameter;
};
