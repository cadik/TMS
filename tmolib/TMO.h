// TMO.h: interface for the CTMO class.
//
//////////////////////////////////////////////////////////////////////
#ifndef TMO_LIB
#define TMO_LIB

#include <string.h>
#include "TMOExports.h"
#include "TMOImage.h"
#include "TMOParameter.h"

enum TMOLIB_API TMO_EXCEPTION
{
	TMO_EOK = 0,
	TMO_EPARAM
};


struct vect2D {
	double x;
	double y;
};


class TMOLIB_API TMO  
{
private:
	class TMOParameterList
	{
	public:
		TMOParameterList(TMOParameter *param, TMOParameterList *next):pParam(param),pNext(next) {}
		TMOParameter* pParam;
		TMOParameterList* pNext;
	};

public:
	virtual int EnumParameters(TMOParameter *buffer[]);
	virtual int GetParameterCount();
	virtual TMOParameter& GetParameter(const wchar_t* s);
	virtual int GetImages(TMOImage **src, TMOImage **dest);
	virtual int SetImages(TMOImage& src, TMOImage& dest);
	virtual int GetImage(TMOImage **image);
	virtual int SetImage(TMOImage& image);
	virtual TMOImage& GetSource(TMOImage **p);
	virtual TMOImage& GetDestination(TMOImage **p);
	virtual int SetSource(TMOImage& p);
	virtual int SetDestination(TMOImage& p);
	virtual const wchar_t* GetName() { return sName; }
	virtual const wchar_t* GetDescription() { return sDescription; }
	virtual int SetName(const wchar_t* name);
	virtual int SetDescription(const wchar_t* description);

	TMO (TMOImage& source, TMOImage& dest);
	TMO();
	virtual ~TMO();

	virtual int Transform() = 0;

protected:
	int iParamCount;
	virtual int Register(TMOParameter& p);

	TMOParameterList* pParameters;	// registered parameters
	TMOImage* pSrc;					// source picture
	TMOImage* pDst;					// transformed picture
	int iWidth, iHeight;			// image dimensions
	wchar_t* sName;					// operator name
	wchar_t* sDescription;			// description (author, year, technique)
	int (*WriteLine)(const wchar_t* text);
};

#endif

