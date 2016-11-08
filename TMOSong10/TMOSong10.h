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

class TMOSong10 : public TMO  
{
private:
	TMODouble k;
	TMODouble alpha;
protected:	

public:
	TMOSong10();
	virtual ~TMOSong10();
	virtual int Transform();
};

#endif
