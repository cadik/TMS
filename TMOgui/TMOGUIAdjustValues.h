// TMOGUIAdjustValues.h: interface for the TMOGUIAdjustValues class.
//
//////////////////////////////////////////////////////////////////////
#ifndef TMOGUIADJUSTVALUES_H
#define TMOGUIADJUSTVALUES_H

#include <qobject.h>
//#include <iostream>

class TMOGUIAdjustValues : public QObject
{
	Q_OBJECT
public:
	void UpdateValues();
	void RenderValues();
	TMOGUIAdjustValues(
		double RMinimum,
		double GMinimum,
		double BMinimum,
		double RMaximum,
		double GMaximum,
		double BMaximum,
		double RGamma,
		double GGamma,
		double BGamma);
	TMOGUIAdjustValues(
		double Minimum,
		double Maximum,
		double Gamma);
	TMOGUIAdjustValues();
	TMOGUIAdjustValues& operator=(TMOGUIAdjustValues& v);
	void reset();
	double dRMinimum;
	double dGMinimum;
	double dBMinimum;
	double dRMaximum;
	double dGMaximum;
	double dBMaximum;
	double dRAverage;
	double dGAverage;
	double dBAverage;
	double dRGamma;
	double dGGamma;
	double dBGamma;
	double dMaximum;
	double dMinimum;
	double dAverage;
	double dExtreme;
/*
void print(int i=0)
	{
std::cout << i << ": "
<< dRMinimum
<< " "
<< dGMinimum
<< " "
<< dBMinimum
<< " "
<< dMinimum
<< " - "
<< dRMaximum
<< " "
<< dGMaximum
<< " "
<< dBMaximum
<< " "
<< dMaximum
<< " "
<< std::endl;

	}*/
signals:
	void valueschanged();
	void render();
};

#endif //TMOGUIADJUSTVALUES_H
