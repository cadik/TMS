#include <stdio.h>
#include <math.h>
//#include "TMO.h"
#include "TMOv.h"






TMOv::TMOv()
{
	vSrc = 0;
	vDst = 0;
	pParameters = 0;
	sName = 0;
	sDescription = 0;
	iParamCount = 0;
	pParameters = 0;
	
	
}

TMOv::~TMOv()
{
	TMOParameterList *temp;
	
	while (pParameters) 
	{
		temp = pParameters;
		pParameters = pParameters->pNext;
		delete temp;
	}
}

int TMOv::SetVideos(TMOVideo& src,TMOVideo& dst)
{
	vSrc = &src;
	vDst = &dst;
	return 0;
}
int TMOv::SetOutVideo(TMOVideo& dst)
{
	
	vDst = &dst;
	return 0;
}
int TMOv::Register(TMOParameter &p)
{
	TMOParameterList *temp;

	for (temp = pParameters; temp; temp = temp->pNext)
		if (temp->pParam->GetName() == p.GetName()) return -1;

	pParameters = new TMOParameterList(&p, pParameters);
	iParamCount++;
	return 0;
}

