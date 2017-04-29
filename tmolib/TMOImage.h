/*
 * TMOImage.h
 *
 * Techniky mapování tónù
 * Autor: Odrej Hajdok
 *
 * Popisuje tridu TMOImage, ktera slouzi pro praci s vysoce kontrastnimi obrazky.
 * Trida poskytuje zakladni sadu funkci pro nacitani a ukladani obrazku
 * ulozenych ve forme TIFF, HDR, OpenEXR.
 *
 * Definice funkci : TMOImage.cpp
 * Zavislosti : knihovna libtiff (http://www.libtiff.org)
 */
#include "tiffio.h"
#include "TMOExports.h"
#include <assert.h>

//#define ERIKLIB
//#undef ERIKLIB

/*#ifdef ERIKLIB
extern "C" {
	#include "defines.h"
}
#endif*/

#define TAKE_LOG(a)  (log((double)(0.00001+(a))))
#define TAKE_LOG2(a)  (log((double)(0.00001+(a)))/log(2.))
#define TAKE_LOG10(a)  (log10((double)(0.00001+(a))))

#define EPS (216.0/24389.0)
#define KAPPA (24389.0/27.0)

// XYZ reference white
#define XYZ_WHITE_X 95.047
#define XYZ_WHITE_Y 100.0
#define XYZ_WHITE_Z 108.8830

// Luv reference white
#define LUV_WHITE_U ((4 * XYZ_WHITE_X) / (XYZ_WHITE_X + 15 * XYZ_WHITE_Y + 3 * XYZ_WHITE_Z))
#define LUV_WHITE_V ((9 * XYZ_WHITE_Y) / (XYZ_WHITE_X + 15 * XYZ_WHITE_Y + 3 * XYZ_WHITE_Z))

// adapting luminance
#define L_A 20.0

enum TMOLIB_API TMO_IMAGEFORMAT
{
	TMO_RGB = 0,
	TMO_XYZ = 1,
	TMO_Yxy = 2,
	TMO_LCH = 3,
	TMO_LAB = 4,
	TMO_Y = 256,
	TMO_NOTSPEC,
};

//MC - output file format
enum TMOLIB_API TMO_FILEFORMAT
{
	TMO_TIFF_8	= 0,
	TMO_TIFF_32 = 1,
	//TMO_RAW_16	= 2, // not implemented
	TMO_RAW_32	= 3,
	TMO_EXR_16	= 4,
	//TMO_EXR_32	= 5, // not implemented
	TMO_PFM_32	= 6,
	TMO_HDR_32	= 7,
	TMO_JPEG_8	= 8,
	//TMO_JPEG_32	= 8, // not implemented
};


/*#ifdef ERIKLIB
int TMOLIB_API SaveImageset(imageset1_t *imageset, const char* filename, bool pyramid=true);
#endif*/


class TMOLIB_API TMOImage
{
protected:
// promenne
	double dStonits;		// input conversion to nits
	double* pData;			// pixel data
	float fXres;			// image resolution X
	float fYres;			// image resolution Y
	int iWidth;				// image width
	int iHeight;			// image height
	char* pName;			// image name
	int iFormat;			// pixel data format (enumerated in TMO_IMAGEFORMAT)
	unsigned short iBitsPerSample;	// bits per sample -> byte => 8, float => 32...
	unsigned short iPhotometric;	// photometric type
	static double XYZ2RGB[][3];		// conversion matrix
	static double RGB2XYZ[][3];		// conversion matrix
//funkce	
	virtual int ReadHeader(TIFF*);
	virtual int WriteHeader(TIFF*,bool);
	virtual int ReadData(TIFF*);
	virtual int WriteData(TIFF*,bool);
	virtual int WriteDataRAW(FILE*);
	virtual int OpenTIFF_8_32();
	virtual int OpenPFM_32();
	virtual int OpenEXR_16();
	virtual int OpenJPEG_32();
	virtual int OpenRAW_32(); //filename_XRESxYRES.hdrraw !
	virtual int OpenRAW_16(); //filename_XRESxYRES.raw ! - canon 3xSHORT
	virtual int SaveTIFF_32();
	virtual int SaveTIFF_8();
	virtual int SaveRAW_32();
	virtual int SaveEXR_16_pData_is_RGB();
	virtual int SaveEXR_16();
	virtual int SavePFM_32();
	virtual int SaveHDR_32();
	virtual int SaveJPEG_8(int quality=90);
	virtual int OpenHDR_32();
	int (*pProgressBar)(TMOImage*, int part, int all);
	int (*pWriteLine)(TMOImage*, const wchar_t* text);
	static int DefaultProgressBar(TMOImage*, int part, int all);
	static int DefaultWriteLine(TMOImage*, const wchar_t* text);
	virtual int Clear();
	double InverseSrgbCompanding(double);
	double SrgbCompanding(double);
	
public:
	inline double* GetOffset(int offset)
	{
	 return (iFormat < 0xFF) ? &pData[3 * offset] : &pData[offset];
	}
	inline double* GetPixel(int x, int y)
	{
	 return (iFormat < 0xFF) ? &pData[3 * (y * iWidth + x)] : &pData[y * iWidth + x];
	}
	virtual double GetStonits();
	virtual int Convert(int format, bool fast=false);
	virtual double* GetData();
	virtual int SaveWithSuffix(const char *suffix, int fileFormat=TMO_TIFF_8);
	virtual int Close();
	virtual int New(int width, int height, int format = TMO_RGB);
	virtual int New(const TMOImage &ref, int format = TMO_NOTSPEC, bool bCopy = false);
	virtual int SaveAs(const char* filename, int fileFormat=TMO_TIFF_8);
	virtual int Save(int fileFormat=TMO_TIFF_8);
	virtual int Open(const char *filename);
	virtual int GetWidth() {return iWidth;}
	virtual int GetHeight() {return iHeight;}
	virtual operator double*() {return(GetData());}
	virtual double** ExportLuminance();
	virtual int ImportLuminance(double **luminance, int width, int height);
	virtual TMOImage& operator=(const TMOImage& image);
	virtual int SetProgress(int(*)(TMOImage*, int, int));
	virtual int SetWriteLine(int(*)(TMOImage*, const wchar_t*));
	virtual int ProgressBar(int part, int all);
	virtual int WriteLine(const wchar_t* text);
	virtual double GetLuminance(int x, int y, int r = 0);
	virtual void GetStatisticsCircleRange(int x, int y, int r = 0);
	virtual void GetStatisticsSquareRange(int x, int y, int size = 0);
	virtual double GetLuminanceYxy(int x, int y, int r = 0);
	virtual double SetLuminance(int x, int y, double L);
	virtual int GetMinMaxAvg(double *minimum, double *maximum, double *average);
	virtual int GetMinMaxAvgLog10(double *minimum, double *maximum, double *average);
	virtual int GetDimensions(int *x, int *y);
	virtual int SetDimensions(int x, int y) {iWidth = x; iHeight = y; return 0;}

	virtual const char* GetFilename() {return(pName);};
	virtual int SetFilename(const char* filename);
	virtual int CorrectGamma(double gamma); 
	virtual int CorrectGammaYxy(double gamma); //just for testing, not reasonable
	
	struct TMOStatistics
	{
		float avgLum;
		float minLum;
		float maxLum;
		float curLum;
		double avgColor[3];		
		double color[3];		
	};
	TMOStatistics statistics;
	TMOImage();
	TMOImage(const char *filename);
	virtual ~TMOImage();
	static double RadiansToDegrees(double);
	static double DegreesToRadians(double);
	static void LabToXyz(double, double, double, double *, double *, double *);
	static void XyzToLuv(double, double, double, double *, double *, double *);	
};


