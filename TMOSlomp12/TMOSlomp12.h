/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   Fast Local Tone Mapping, Summed-Area Tables                *
 *                          and Mesopic Vision Simulation                       *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMO.h"

class TMOSlomp12 : public TMO
{
public:
	TMOSlomp12();
	virtual ~TMOSlomp12();
	cv::Mat TMOImageToLogLuminanceMat();
	void logLuminanceImage(cv::Mat *luminanceMat);
	double fullMipmap(cv::Mat *mat);
	void scaleLuminance(cv::Mat *luminanceMat, double keyValue);
	void scaledLuminanceImage(cv::Mat *luminanceMat);
	double boxFilter(cv::Mat *SAT, int x, int y, int s);
	double getNormalizedDifference(double conv0, double conv1, int s);
	int getMaxScale(cv::Mat *SAT, int x, int y);
	double redResponseValue(double illuminance);
	double arithLuminanceAverage();
	virtual int Transform();

protected:
	TMOBool local;
	TMOBool mesopic;
	TMOBool varying;

	double alpha = 0.72;
	double phi = 8.;
	double epsilon = 0.025;
};
