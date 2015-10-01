// TMOGUIAdjustValues.cpp: implementation of the TMOGUIAdjustValues class.
//
//////////////////////////////////////////////////////////////////////

#include "TMOGUIAdjustValues.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMOGUIAdjustValues::TMOGUIAdjustValues(
	double RMinimum,
	double GMinimum,
	double BMinimum,
	double RMaximum,
	double GMaximum,
	double BMaximum,
	double RGamma,
	double GGamma,
	double BGamma):
	dRMinimum(RMinimum),
	dGMinimum(GMinimum),
	dBMinimum(BMinimum),
	dRMaximum(RMaximum),
	dGMaximum(GMaximum),
	dBMaximum(BMaximum),
	dRGamma(RGamma),
	dGGamma(GGamma),
	dBGamma(BGamma),
	dMinimum(0),
	dMaximum(1)
{
}

TMOGUIAdjustValues::TMOGUIAdjustValues(
	double Minimum,
	double Maximum,
	double Gamma):
	dRMinimum(Minimum),
	dGMinimum(Minimum),
	dBMinimum(Minimum),
	dRMaximum(Maximum),
	dGMaximum(Maximum),
	dBMaximum(Maximum),
	dRGamma(Gamma),
	dGGamma(Gamma),
	dBGamma(Gamma),
	dMinimum(Minimum),
	dMaximum(Maximum)
{
}

TMOGUIAdjustValues::TMOGUIAdjustValues():
	dRMinimum(0),
	dGMinimum(0),
	dBMinimum(0),
	dRMaximum(1),
	dGMaximum(1),
	dBMaximum(1),
	dRGamma(1),
	dGGamma(1),
	dBGamma(1),
	dMinimum(0),
	dMaximum(1)
{
}

TMOGUIAdjustValues& TMOGUIAdjustValues::operator=(TMOGUIAdjustValues& v)
{
	dRMinimum = v.dRMinimum;
	dGMinimum = v.dGMinimum;
	dBMinimum = v.dBMinimum;
	dRMaximum = v.dRMaximum;
	dGMaximum = v.dGMaximum;
	dBMaximum = v.dBMaximum;
	dRGamma = v.dRGamma;
	dGGamma = v.dGGamma;
	dBGamma = v.dBGamma;
	dMaximum = v.dMaximum;
	dMinimum = v.dMinimum;
	return *this;
}

void TMOGUIAdjustValues::reset()
{
	dRMinimum = 0;
	dGMinimum = 0;
	dBMinimum = 0;
	dRMaximum = 1;
	dGMaximum = 1;
	dBMaximum = 1;
	dRGamma = 1;
	dGGamma = 1;
	dBGamma = 1;
}

void TMOGUIAdjustValues::UpdateValues()
{
	emit valueschanged();
}

void TMOGUIAdjustValues::RenderValues()
{
	emit render();
}
