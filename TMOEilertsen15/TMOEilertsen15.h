/*******************************************************************************
*                                                                              *
*                        Brno University of Technology                         *
*                      Faculty of Information Technology                       *
*                                                                              *
*                              Tone Mapping Methods                            *
*                                                                              *
*            Author: Hugo Bohácsek [xbohach00 AT stud.fit.vutbr.cz]            *
*                                   Brno 2025                                  *
*                                                                              *
*                    Implementation of the TMOEilertsen15 header               *
*                         Real-time noise-aware tone mapping                   *
*                     https://computergraphics.on.liu.se/rntm/                 *
*                                                                              *
*******************************************************************************/

#ifndef TMO_EILERTSEN15_H
#define TMO_EILERTSEN15_H

#include "TMO.h"
#include <vector>
#include <cmath>

class TMOEilertsen15 : public TMO {
    public:
        TMOEilertsen15();
        virtual ~TMOEilertsen15();
        virtual int Transform();

    protected:
        // Parameters
        TMODouble dPeakLuminance;      // Peak display luminance [cd/m2]
        TMODouble dBlackLevel;         // Display black level [cd/m2]
        TMODouble dAmbientLight;       // Ambient light [lux]
        TMODouble dReflectivity;       // Display reflectivity [0-1%]
        TMODouble dDetailScaling;      // Detail enhancement factor
        TMODouble dNoiseControl;       // Noise visibility control
        TMODouble dTonePriority;       // Tone priority <-1; 1>
        TMOBool bLocalToneCurves;      // Use local tone curves
        TMOInt iFilterIterations;      // Number of diffusion iterations
};

#endif // TMO_EILERTSEN15_H