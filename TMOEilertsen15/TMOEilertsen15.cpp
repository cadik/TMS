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
*                     Implementation of the TMOEilertsen15 class               *
*                         Real-time noise-aware tone mapping                   *
*                     https://computergraphics.on.liu.se/rntm/                 *
*                                                                              *
*******************************************************************************/

#include "TMOEilertsen15.h"
#include <algorithm>
#include <cmath>
#include <cstring>

TMOEilertsen15::TMOEilertsen15(){
   
   SetName(L"Eilertsen15");
   SetDescription(L"Real-time noise-aware tone mapping");

   // Display parameters
   dPeakLuminance.SetName(L"Peak Luminance");
   dPeakLuminance.SetDescription(L"Peak display luminance in cd/m2");
   dPeakLuminance.SetDefault(200.0);
   dPeakLuminance.SetRange(50.0, 10000.0);
   this->Register(dPeakLuminance);

   dBlackLevel.SetName(L"Black Level");
   dBlackLevel.SetDescription(L"Display black level in cd/m2");
   dBlackLevel.SetDefault(0.5);
   dBlackLevel.SetRange(0.01, 10.0);
   this->Register(dBlackLevel);

   dAmbientLight.SetName(L"Ambient Light");
   dAmbientLight.SetDescription(L"Ambient illuminance in lux");
   dAmbientLight.SetDefault(100.0);
   dAmbientLight.SetRange(0.0, 5000.0);
   this->Register(dAmbientLight);

   dReflectivity.SetName(L"Reflectivity");
   dReflectivity.SetDescription(L"Display reflectivity in %");
   dReflectivity.SetDefault(0.8);
   dReflectivity.SetRange(0.1, 5.0);
   this->Register(dReflectivity);

   // Processing parameters
   dDetailScaling.SetName(L"Detail Scaling");
   dDetailScaling.SetDescription(L"Detail enhancement factor");
   dDetailScaling.SetDefault(1.0);
   dDetailScaling.SetRange(0.0, 4.0);
   this->Register(dDetailScaling);

   dNoiseControl.SetName(L"Noise Control");
   dNoiseControl.SetDescription(L"Noise visibility control");
   dNoiseControl.SetDefault(1.0);
   dNoiseControl.SetRange(0.0, 2.0);
   this->Register(dNoiseControl);

   dTonePriority.SetName(L"Tone Priority");
   dTonePriority.SetDescription(L"Priority: -1=bright, 0=neutral, 1=dark");
   dTonePriority.SetDefault(0.0);
   dTonePriority.SetRange(-1.0, 1.0);
   this->Register(dTonePriority);

   bLocalToneCurves.SetName(L"Local Tone Curves");
   bLocalToneCurves.SetDescription(L"Use local tone curves");
   bLocalToneCurves.SetDefault(true);
   this->Register(bLocalToneCurves);

   iFilterIterations.SetName(L"Filter Iterations");
   iFilterIterations.SetDescription(L"Number of diffusion iterations");
   iFilterIterations.SetDefault(12);
   iFilterIterations.SetRange(1, 20);
   this->Register(iFilterIterations);
}

TMOEilertsen15::~TMOEilertsen15()
{
}

