// TMOParameter.cpp: implementation of the TMOParameter class.
//
//////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <wchar.h>
#include "TMOParameter.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
unsigned int TMOParameter::iCounter = 0;
double TMODouble::dInitial = (double)0xc0c0c0c0;
int TMOInt::iInitial = (int)0xc0c0c0c0;

TMOParameter::TMOParameter()
{
	iType = TMO_VOID;
	iCounter++;
	sName = new wchar_t[10];
#ifndef LINUX
        sName = _itow(iCounter, 0, 10); //TODO: najit nahradu _itow
#endif
	sDescription = 0;
	
}

TMOParameter::~TMOParameter()
{
	if (sName) 
	{
		delete[] sName;
		sName = 0;
	}
	if (sDescription) 
	{
		delete[] sDescription;
		sDescription = 0;
	}
	iCounter--;
}

int TMOParameter::SetName(const wchar_t* name)
{
	int iLength = wcslen(name);
	if (sName) 
	{
		delete[] sName;
		sName = 0;
	}
	sName = new wchar_t[iLength + 1];
	wcscpy(sName, name);
	return iLength;
}

int TMOParameter::SetDescription(const wchar_t* description)
{
	int iLength = wcslen(description);
	if (sDescription) 
	{
		delete[] sDescription;
		sDescription = 0;
	}
	sDescription = new wchar_t[iLength + 1];
	wcscpy(sDescription, description);
	return iLength;
}

//////////////////////////////////////////////////////////////////////
// TMODouble definitions
//////////////////////////////////////////////////////////////////////

TMODouble::TMODouble():TMOParameter()
{
	dMinimum = dMaximum = dValue = 0.;
	dDefault = dInitial;
	iType = TMO_DOUBLE;
}

TMODouble::TMODouble(double value):TMOParameter()
{
	dMinimum = dMaximum = 0.;
	dValue = value;
	dDefault = dInitial;
	iType = TMO_DOUBLE;
}

TMODouble::~TMODouble()
{

}

TMOParameter& TMODouble::operator=(double value)
{
	if (dMaximum==dMinimum) 
	{
		dValue = value;
		if (dDefault == dInitial) dDefault = dValue;
		return *this;
	}
	if (value < dMinimum) 
	{
		dValue = dMinimum;
		if (dDefault == dInitial) dDefault = dValue;
		return *this;
	}
	dValue = value;
	if (value > dMaximum) dValue = dMaximum;
	if (dDefault == dInitial) dDefault = dValue;
	return *this;
}

TMOParameter& TMODouble::operator=(int value)
{
	return *this = static_cast<double>(value);
}

TMOParameter& TMODouble::operator=(bool value)
{
	return *this = static_cast<bool>(value);
}

double TMODouble::GetDouble()
{
	return dValue;
}

TMODouble::operator double()
{
	return dValue;
}

int TMODouble::GetInt()
{
	return static_cast<int>(dValue);
}

bool TMODouble::GetBool()
{
	return dValue > 0;
}

TMODouble& TMODouble::operator+=(double v)
{
	double value = dValue + v;
	if (dMaximum==dMinimum)
	{
		dValue = value;
		return *this;
	}
	if (value < dMinimum) 
	{
		dValue = dMinimum;
		return *this;
	}
	dValue = value;
	if (value > dMaximum) dValue = dMaximum;
	return *this;
}

TMODouble& TMODouble::operator*=(double v)
{
	double value = dValue * v;
	if (dMaximum==dMinimum)
	{
		dValue = value;
		return *this;
	}
	if (value < dMinimum) 
	{
		dValue = dMinimum;
		return *this;
	}
	dValue = value;
	if (value > dMaximum) dValue = dMaximum;
	return *this;
}

TMODouble& TMODouble::operator/=(double v)
{
	double value = dValue / v;
	if (dMaximum==dMinimum)
	{
		dValue = value;
		return *this;
	}
	if (value < dMinimum) 
	{
		dValue = dMinimum;
		return *this;
	}
	else dValue = value;
	if (value > dMaximum) dValue = dMaximum;
	return *this;
}


int TMODouble::SetRange(double minimum, double maximum)
{
	dMinimum = minimum;
	dMaximum = maximum;
	return 0;
}

int TMODouble::GetRange(double & minimum, double & maximum)
{
	minimum = dMinimum;
	maximum = dMaximum;
	return 0;
}

int TMODouble::SetDefault(double value)
{
	if (dMaximum==dMinimum) 
	{
		dDefault = value;
		return 0;
	}
	if (value < dMinimum) 
	{
		dDefault = dMinimum;
		return 0;
	}
	dDefault = value;
	if (value > dMaximum) dDefault = dMaximum;
	return 0;
}

int TMODouble::Reset()
{
	if (dDefault == dInitial) return -1;
	dValue = dDefault;
	return 0;
}

//////////////////////////////////////////////////////////////////////
// TMOInt definitions
//////////////////////////////////////////////////////////////////////

TMOInt::TMOInt():TMOParameter()
{
	iMinimum = iMaximum = iValue = 0;
	iDefault = iInitial;
	iType = TMO_INT;
}

TMOInt::TMOInt(int value):TMOParameter()
{
	iMinimum = iMaximum = 0;
	iValue = value;
	iDefault = iInitial;
	iType = TMO_INT;
}

TMOInt::~TMOInt()
{

}

TMOParameter& TMOInt::operator=(int value)
{
	if (iMaximum==iMinimum) 
	{
		iValue = value;
		if (iDefault == iInitial) iDefault = iValue;
		return *this;
	}
	if (value < iMinimum) 
	{
		iValue = iMinimum;
		if (iDefault == iInitial) iDefault = iValue;
		return *this;
	}
	iValue = value;
	if (value > iMaximum) iValue = iMaximum;
	if (iDefault == iInitial) iDefault = iValue;
	return *this;
}

TMOParameter& TMOInt::operator=(double value)
{
	return *this = static_cast<int>(value);
}

TMOParameter& TMOInt::operator=(bool value)
{
	return *this = static_cast<int>(value);
}

int TMOInt::GetInt()
{
	return iValue;
}

TMOInt::operator int()
{
	return iValue;
}

double TMOInt::GetDouble()
{
	return static_cast<double>(iValue);
}

bool TMOInt::GetBool()
{
	return iValue > 0;
}

TMOInt& TMOInt::operator+=(int v)
{
	int value = iValue + v;
	if (iMaximum==iMinimum)
	{
		iValue = value;
		return *this;
	}
	if (value < iMinimum) 
	{
		iValue = iMinimum;
		return *this;
	}
	iValue = value;
	if (value > iMaximum) iValue = iMaximum;
	return *this;
}

TMOInt& TMOInt::operator*=(int v)
{
	int value = iValue * v;
	if (iMaximum==iMinimum)
	{
		iValue = value;
		return *this;
	}
	if (value < iMinimum) 
	{
		iValue = iMinimum;
		return *this;
	}
	iValue = value;
	if (value > iMaximum) iValue = iMaximum;
	return *this;
}

TMOInt& TMOInt::operator/=(int v)
{
	int value = iValue / v;
	if (iMaximum==iMinimum)
	{
		iValue = value;
		return *this;
	}
	if (value < iMinimum) 
	{
		iValue = iMinimum;
		return *this;
	}
	else iValue = value;
	if (value > iMaximum) iValue = iMaximum;
	return *this;
}


int TMOInt::SetRange(int minimum, int maximum)
{
	iMinimum = minimum;
	iMaximum = maximum;
	return 0;
}

int TMOInt::GetRange(int & minimum, int & maximum)
{
	minimum = iMinimum;
	maximum = iMaximum;
	return 0;
}

int TMOInt::SetDefault(int value)
{
	if (iMaximum==iMinimum) 
	{
		iDefault = value;
		return *this;
	}
	if (value < iMinimum) 
	{
		iDefault = iMinimum;
		return *this;
	}
	iDefault = value;
	if (value > iMaximum) iDefault = iMaximum;
	return *this;
}

int TMOInt::Reset()
{
	if (iDefault == iInitial) return -1;
	iValue = iDefault;
	return 0;
}

//////////////////////////////////////////////////////////////////////
// TMOBool definitions
//////////////////////////////////////////////////////////////////////

TMOBool::TMOBool():TMOParameter()
{
	bValue = false;
	bDefault = false;
	bInitial = true;
	iType = TMO_BOOL;
}

TMOBool::TMOBool(bool value):TMOParameter()
{
	bValue = value;
	bDefault = false;
	bInitial = true;
	iType = TMO_BOOL;
}

TMOBool::~TMOBool()
{

}

TMOParameter& TMOBool::operator=(bool value)
{
	bValue = value;
	if (bInitial) 
	{
		bDefault = bValue;
		bInitial = false;
	}
	return *this;
}

TMOParameter& TMOBool::operator=(int value)
{
	return *this = static_cast<int>(value);
}

TMOParameter& TMOBool::operator=(double value)
{
	return *this = static_cast<double>(value);
}

bool TMOBool::GetBool()
{
	return bValue;
}

TMOBool::operator bool()
{
	return bValue;
}

int TMOBool::GetInt()
{
	return static_cast<int>(bValue);
}

double TMOBool::GetDouble()
{
	return static_cast<double>(bValue);
}

int TMOBool::SetDefault(bool value)
{
	bDefault = value;
	bInitial = false;
	return *this;
}

int TMOBool::Reset()
{
	if (bInitial) return -1;
	bValue = bDefault;
	return 0;
}

//////////////////////////////////////////////////////////////////////
// TMOString definitions
//////////////////////////////////////////////////////////////////////

TMOString::TMOString():TMOParameter()
{
	iType = TMO_STRING;
	sValue = "";
}

TMOString::TMOString(std::string value):TMOParameter()
{
	iType = TMO_STRING;
	sValue = value;
}

TMOString::~TMOString()
{

}

TMOParameter& TMOString::operator=(std::string value)
{	
	sValue = value;		
	return *this;
}

int TMOString::SetDefault(std::string value)
{
	sValue = "";
	return 0;
}

int TMOString::Reset()
{	
	sValue = "";
	return 0;
}

TMOString::operator std::string()
{
	return sValue;
}

std::string TMOString::GetString()
{
	return sValue;
}