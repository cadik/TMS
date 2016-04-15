// TMOParameter.h: interface for the TMOParameter class.
//
//////////////////////////////////////////////////////////////////////
#include "string.h"
#include "TMOExports.h"
#include <iostream>

//typedef unsigned short wchar_t;

enum TMOLIB_API TMO_PARAMETER
{
	TMO_VOID = 0,
	TMO_DOUBLE = 1,
	TMO_INT = 2,
	TMO_BOOL = 4,
	TMO_STRING = 5,
};

class TMOLIB_API TMOParameter  
{
public:
	TMOParameter();
	virtual ~TMOParameter();

	virtual TMOParameter& operator=(double value) {return *this;}
	virtual TMOParameter& operator=(int value) {return *this;}
	virtual TMOParameter& operator=(bool value) {return *this;}
	virtual TMOParameter& operator=(std::string value) {return *this;}

	virtual double GetDouble () {return 0.;}
	virtual int GetInt () {return 0;}
	virtual bool GetBool () {return false;}

	virtual const wchar_t* GetName() { return sName; }
	virtual const wchar_t* GetDescription() { return sDescription; }
	virtual int SetName(const wchar_t* name);
	virtual int SetDescription(const wchar_t* description);
	//virtual bool Is(int type) {return (type&iType)>0;}
	virtual bool Is(int type) {return (type == iType);}

	virtual int Reset() {return 0;}
protected:
	wchar_t* sName;
	wchar_t* sDescription;
	int iType;
private:
	static unsigned int iCounter;
};

class TMOLIB_API TMODouble : public TMOParameter  
{
public:
	virtual int Reset();
	virtual int SetDefault(double value);
	virtual double GetDefault(){return dDefault;};
	virtual int SetRange(double minimum, double maximum);
	virtual int GetRange(double & minimum, double & maximum);

	TMODouble();
	TMODouble(double value);
	virtual ~TMODouble();

	virtual TMOParameter& operator=(double value);
	virtual TMOParameter& operator=(int value);
	virtual TMOParameter& operator=(bool value);

	virtual TMODouble& operator+=(double value);
	virtual TMODouble& operator*=(double value);
	virtual TMODouble& operator/=(double value);

	virtual bool operator> (double value) {return dValue > value;}
	virtual bool operator>=(double value) {return dValue >=value;}
	virtual bool operator< (double value) {return dValue < value;}
	virtual bool operator<=(double value) {return dValue <=value;}
	virtual bool operator==(double value) {return dValue == value;}
	virtual bool operator!=(double value) {return dValue != value;}


	virtual double GetDouble ();
	virtual operator double();
	virtual int GetInt ();
	virtual bool GetBool ();

protected:
	double dDefault;
	double dMinimum;
	double dMaximum;
	double dValue;
	static double dInitial;
};

class TMOLIB_API TMOInt : public TMOParameter  
{
public:
	virtual int Reset();
	virtual int SetDefault(int value);
	virtual int GetDefault(){return iDefault;};
	virtual int SetRange(int minimum, int maximum);
	virtual int GetRange(int & minimum, int & maximum);

	TMOInt();
	TMOInt(int value);
	virtual ~TMOInt();

	virtual TMOParameter& operator=(double value);
	virtual TMOParameter& operator=(int value);
	virtual TMOParameter& operator=(bool value);

	virtual TMOInt& operator+=(int value);
	virtual TMOInt& operator*=(int value);
	virtual TMOInt& operator/=(int value);

	virtual bool operator> (int value) {return iValue > value;}
	virtual bool operator>=(int value) {return iValue >=value;}
	virtual bool operator< (int value) {return iValue < value;}
	virtual bool operator<=(int value) {return iValue <=value;}
	virtual bool operator==(int value) {return iValue == value;}
	virtual bool operator!=(int value) {return iValue != value;}


	virtual double GetDouble ();
	virtual operator int();
	virtual int GetInt ();
	virtual bool GetBool ();

protected:
	int iDefault;
	int iMinimum;
	int iMaximum;
	int iValue;
	static int iInitial;
};

class TMOLIB_API TMOBool : public TMOParameter  
{
public:
	virtual int Reset();
	virtual int SetDefault(bool value);
	virtual bool GetDefault(){return bDefault;};

	TMOBool();
	TMOBool(bool value);
	virtual ~TMOBool();

	virtual TMOParameter& operator=(double value);
	virtual TMOParameter& operator=(int value);
	virtual TMOParameter& operator=(bool value);

	virtual double GetDouble ();
	virtual operator bool();
	virtual int GetInt ();
	virtual bool GetBool ();

protected:
	bool bDefault;
	bool bValue;
	bool bInitial;
};

class TMOLIB_API TMOString : public TMOParameter
{
public:
	virtual int Reset();
	virtual int SetDefault(std::string value);
	
	TMOString();
	TMOString(std::string value);	
	virtual ~TMOString();
	
	virtual TMOParameter& operator=(std::string value);

	virtual operator std::string();	
	virtual std::string GetString ();
protected:
	std::string sValue;
};
