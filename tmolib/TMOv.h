
#ifndef TMOv_LIB
#define TMOv_LIB



#include "TMO.h"








class TMOLIB_API TMOv : public TMO
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
	
	
	virtual int SetVideos(TMOVideo& src, TMOVideo& dst);
	virtual int SetOutVideo(TMOVideo& dst);

	TMOv();
	virtual ~TMOv();

	
	virtual int TransformVideo()=0;
	

protected:
	int iParamCount;
	virtual int Register(TMOParameter& p);

	TMOParameterList* pParameters;	// registered parameters
	
	TMOVideo* vSrc;
	TMOVideo* vDst;
};

#endif