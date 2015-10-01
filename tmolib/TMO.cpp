// TMO.cpp: implementation of the TMO class.
//
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "TMO.h"

#ifdef LINUX
#include <wchar.h>
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMO::TMO()
{
	pSrc = 0;
	pDst = 0;
	pParameters = 0;
	sName = 0;
	sDescription = 0;
	iParamCount = 0;
	pParameters = 0;
}

TMO::~TMO()
{
	TMOParameterList *temp;
	
	while (pParameters) 
	{
		temp = pParameters;
		pParameters = pParameters->pNext;
		delete temp;
	}
}

TMO::TMO(TMOImage& source, TMOImage& dest)
{
	pSrc = &source;
	pDst = &dest;
	pParameters = 0;
	iParamCount = 0;
	sName = 0;
	sDescription = 0;
}

int TMO::SetSource(TMOImage& p)
{
	pSrc = &p;
	return 0;
}

int TMO::SetDestination(TMOImage& p)
{
	pDst = &p;
	return 0;
}

TMOImage& TMO::GetSource(TMOImage **p)
{
	if (p) *p = pSrc;
	return *pSrc;
}

TMOImage& TMO::GetDestination(TMOImage **p)
{
	if (p) *p = pDst;
	return *pDst;
}

int TMO::SetImages(TMOImage& src, TMOImage& dest)
{
	pSrc = &src;
	pDst = &dest;
	iWidth = pSrc->GetWidth();
	iHeight = pSrc->GetHeight();
	return 0;
}

int TMO::GetImages(TMOImage **src, TMOImage **dest)
{
	*src = pSrc;
	*dest = pDst;
	return 0;
}

int TMO::SetImage(TMOImage& src)
{
	pSrc = &src;
	pDst = &src;
	iWidth = pSrc->GetWidth();
	iHeight = pSrc->GetHeight();
	return 0;
}

int TMO::GetImage(TMOImage **src)
{
	*src = pSrc;
	return 0;
}



int TMO::Register(TMOParameter &p)
{
	TMOParameterList *temp;

	for (temp = pParameters; temp; temp = temp->pNext)
		if (temp->pParam->GetName() == p.GetName()) return -1;

	pParameters = new TMOParameterList(&p, pParameters);
	iParamCount++;
	return 0;
}

TMOParameter& TMO::GetParameter(const wchar_t* s)
{
	TMOParameterList *temp;

	for (temp = pParameters; temp; temp = temp->pNext)
		if (wcscmp(temp->pParam->GetName(), s) == 0) return *(temp->pParam);
	throw TMO_EPARAM;
}

int TMO::GetParameterCount()
{
	return iParamCount;
}

int TMO::EnumParameters(TMOParameter *buffer[])
{
	TMOParameterList *temp;
	int index = 0;

	for (temp = pParameters; temp; temp = temp->pNext)
	{
		buffer[index++] = temp->pParam;
	}
	return iParamCount;
}

int TMO::SetName(const wchar_t* name)
{
	int iLength = wcslen(name);
	if (sName) delete[] sName;
	sName = new wchar_t[iLength + 1];
	wcscpy(sName, name);
	return iLength;
}

int TMO::SetDescription(const wchar_t* description)
{
	int iLength = wcslen(description);
	if (sDescription) delete[] sDescription;
	sDescription = new wchar_t[iLength + 1];
	wcscpy(sDescription, description);
	return iLength;
}

