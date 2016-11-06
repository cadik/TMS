/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                bachelor thesis                               *
*             Author: Martin Molek [xmolek00 AT stud.fit.vutbr.cz]             *
*                                    Brno 2017                                 *
*                                                                              *
*******************************************************************************/

#ifndef TMORC2GNGM_H
#define TMORC2GNGM_H

#include "TMO.h"
#include <iostream>
#include <math.h>

class TMOMM16 : public TMO  
{
private:
	TMODouble red;
	TMODouble green;
	TMODouble blue;
	TMOBool negative;
protected:	

public:
	TMOMM16();
	virtual ~TMOMM16();
	virtual int Transform();
};

#endif
