
//#ifndef TMO_LIB
//#define TMO_LIB



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
	
	
	//virtual int SetVideos(TMOVideo& src, TMOVideo& dst);
	//virtual int SetOutVideo(TMOVideo& dst);

	//TMO (TMOImage& source, TMOImage& dest);
	TMOv();
	virtual ~TMOv();

	
	virtual int TransformVideo()=0;
	

protected:
	int iParamCount;
	virtual int Register(TMOParameter& p);

	TMOParameterList* pParameters;	// registered parameters
					// transformed picture
	/*int iWidth, iHeight;			// image dimensions
	wchar_t* sName;					// operator name
	wchar_t* sDescription;			// description (author, year, technique)
	int (*WriteLine)(const wchar_t* text);
	TMOVideo* vSrc;
	TMOVideo* vDst;*/
};

//#endif