/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                       Author: Jakub Krystufek                                     *
*                       Brno 2024                                                   *
*                                                                                   *
*                       Implementation of the SPDecolor method                      *
*                       Credits to Matlab implementation                            *
*                       https://github.com/yqx7150/SPDecolor                        *
*                                                                                   *
************************************************************************************/

#include "TMO.h"

class TMOLiu16 : public TMO
{
public:
	TMOLiu16();
	virtual ~TMOLiu16();
	virtual int Transform();
	cv::Mat calculatePolyGrad(const cv::Mat& img, int order);
protected:
	cv::Mat powElement(const cv::Mat& base, int exp);
	TMODouble sigmaParameter;
};
