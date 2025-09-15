/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio	                                       *
*                                                                                   *
*                       Project to subject Computational photography (VYF)          *
*                       Author: Boris Strbak [xstrba05 AT stud.fit.vutbr.cz]        *
*                       Brno 2024                                                   *
*                                                                                   *
*                      Converting color images to grayscale using correlation       *
*                      method by Nafchi, H. Z., Shahkolaei, A.,                     *
*                      Hedjam, R., and Cheriet, M.                                  *
*                                                                                   *
************************************************************************************/
/**
 * @file TMONafchi17_2.h
 * @brief Converting color images to grayscale using correlation method by Nafchi, H. Z., Shahkolaei, A., Hedjam, R., and Cheriet, M.
 * @author Boris Strbak
 * @class TMONafchi17_2.h
 */

#include "TMO.h"

class TMONafchi17_2 : public TMO
{
public:
	TMONafchi17_2();
	virtual ~TMONafchi17_2();
	virtual int Transform();

protected:
	TMODouble invCorrParam;
	TMOBool comParam;
	TMOInt downSampleParam;
};
