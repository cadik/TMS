/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   A Tone Mapping Operator Based on Neural and                *
 *                    Psychophysical Models of Visual Perception                *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMO.h"

class TMOCyriac15 : public TMO
{
public:
	TMOCyriac15();
	virtual ~TMOCyriac15();
	cv::Mat normalizeAndLuminance();
	void addToCumulativeHistogram(std::vector<double> *cumulativeHistogram, double value);
	virtual int Transform();

protected:
	int bins = std::pow(2, 16);
	TMODouble dParameter;
};
