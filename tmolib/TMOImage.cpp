/*
 * TMOImage.cpp
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "TMOImage.h"
#include "TMORadiance.h"
//-OpenEXR-
#ifndef LINUX
 #define OPENEXR_DLL
#else
 #define _stricmp strcasecmp
#endif

#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>

using namespace std;
using namespace Imf;
using namespace Imath;

//~OpenEXR~

//-JPEG-HDR-
#include "freejpeghdr.h"

extern "C" {
#define XMD_H // fixing error C2371: 'INT32' : redefinition; different basic types
#include "jpeglib.h"
} // extern "C"
//~JPEG-HDR~

// tone mapping
#include "TMO_Tonemap.h"

// original version
/*double TMOImage::XYZ2RGB[3][3] = 
{
	{ 2.5651,   -1.1665,   -0.3986},
	{-1.0217,    1.9777,    0.0439},
	{ 0.0753,   -0.2543,    1.1892}
};

double TMOImage::RGB2XYZ[3][3] = 
{ 
	{0.5141364, 0.3238786, 0.16036376},
	{0.265068, 0.67023428, 0.06409157},
	{0.0241188, 0.1228178, 0.84442666}
};*/

// sRGB version
double TMOImage::XYZ2RGB[3][3] = 
{
	{ 3.2404542,   -1.5371385,   -0.4985314},
	{-0.9692660,    1.8760108,    0.0415560},
	{ 0.0556434,   -0.2040259,    1.0572252}
};

double TMOImage::RGB2XYZ[3][3] = 
{ 
	{0.4124564, 0.3575761, 0.1804375},
	{0.2126729, 0.7151522, 0.0721750},
	{0.0193339, 0.1191920, 0.9503041}
};

double TMOImage::RadiansToDegrees(double radians){
	return (radians > 0) ? (radians / M_PI) * 180.0 : 360 - (abs(radians) / M_PI) * 180.0;	
}

double TMOImage::DegreesToRadians(double degrees){
	return degrees * M_PI / 180.0;
}

// konstruktor
TMOImage::TMOImage()
{
	SetProgress(DefaultProgressBar);
	SetWriteLine(DefaultWriteLine);
	Clear();
}

TMOImage::TMOImage(const char *filename)
{
	SetProgress(DefaultProgressBar);
	SetWriteLine(DefaultWriteLine);
	Clear();
	Open(filename);
}

// destruktor
TMOImage::~TMOImage( )
{
	//printf("TMOImage destruction\n");
	if (pData) delete[] pData;
	pData = 0;
	if (pName) delete[] pName;
	pName = 0;
}

int TMOImage::Open(const char *filename)
{
	int length = strlen(filename), i;

	if (pName) delete[] pName;
	pName = new char[length + 6];
	strcpy(pName, filename);

	if (!length) return -1;
	i = length - 1;
	while (i)
	{
		if (pName[i] == '.') break;
		i--;
	}
	i++;
	
	if (_stricmp(&pName[i],"hdr") == 0) return OpenHDR_32();
	if (_stricmp(&pName[i],"pic") == 0) return OpenHDR_32();
	if (_stricmp(&pName[i],"pfm") == 0) return OpenPFM_32();
	if (_stricmp(&pName[i],"exr") == 0) return OpenEXR_16();
	if (_stricmp(&pName[i],"raw") == 0) return OpenRAW_16();
	if (_stricmp(&pName[i],"hdrraw") == 0) return OpenRAW_32();
	if (_stricmp(&pName[i],"tif") == 0) return OpenTIFF_8_32();
	if (_stricmp(&pName[i],"jpg") == 0) return OpenJPEG_32();
	if (_stricmp(&pName[i],"jpeg") == 0) return OpenJPEG_32();
	
	return OpenTIFF_8_32();
}

int TMOImage::OpenTIFF_8_32()
{
	TIFF *pFile;			// TIFF file pointer
	TIFFSetErrorHandler(0);

	// otevreni vstupniho souboru TIFF
	if (( pFile = TIFFOpen(pName, "r")) == NULL) 
	{
		strcat(pName, ".tif");
		if (( pFile = TIFFOpen(pName, "r")) == NULL) 
		{
			strcat(pName, "f");
			if (( pFile = TIFFOpen(pName, "r")) == NULL) 
			{
				delete[] pName;
				pName = 0;
				throw -1;
			}
		}
	}
	ReadHeader(pFile);
	ReadData(pFile);
	TIFFClose(pFile);
	return 0;
}

int TMOImage::OpenHDR_32()
{
	FILE *f;
	COLR *scanin;
	COLOR col;
	int x, y;

	if ((f = fopen(pName, "rb")) == 0) throw -1;

	if (TMORadiance::checkheader(f, (char *) COLRFMT, NULL) < 0 || TMORadiance::fgetresolu(&iWidth, &iHeight, f) < 0)
		throw -2;

	dStonits = 179.0 / TMORadiance::dExposure;

	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];
	iFormat = TMO_RGB;

	if ((scanin = (COLR*) new unsigned char[iWidth*4]) == 0) throw -3;

	for (y = 0; y < iHeight; y++) 
	{
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		if (TMORadiance::freadcolrs(scanin, iWidth, f) < 0) throw -4;
		for (x = 0; x < iWidth; x++) 
		{
			TMORadiance::colr_color(col, scanin[x]);
			GetPixel(x,y)[0] = col[0];
			GetPixel(x,y)[1] = col[1];
			GetPixel(x,y)[2] = col[2];
		}
		if (ferror(f)) throw -5;
			
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
					
	delete[] scanin;
	fclose(f);
	return 0;
}

// ***************************************************************************
// Open PFM image. Just color (no greyscale) images are opened.  
// Scale of the image isn't taken into consideration.
// Implemented version with header elements separated by newline character (!)
// ***************************************************************************
int TMOImage::OpenPFM_32() {
	FILE* pfmFile;
    if (( pfmFile = fopen(pName, "rb")) == NULL) {
        throw -10;
    }
	// *** read file header ***	
    char buf[1000];
    fgets(buf, 100, pfmFile);
    if(buf[1]=='F') {
        iFormat = TMO_RGB; // color image
    } else {
	// don't open greyscale image
	// code for working with greyscale images can be implemented later
        iFormat = TMO_Y;
		throw -100;
	}
	fgets(buf, 100, pfmFile);
	// skip commentary (if present)
	// commentary = line beginning with '#'
	while(buf[0]=='#') {
		fgets(buf, 100, pfmFile);
	}
	// read width and height of image
	sscanf(buf,"%d %d",&iWidth, &iHeight);

	fgets(buf, 100, pfmFile);
	atof(buf); // read scale - to be (maybe) implemented, isn't necessary
	// **** end of header***********************

	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];
	
	// *** read data *****
	
	int row_nr = 0;
	float fr, fg, fb;
	char data_r[4], data_g[4], data_b[4];
	int y;
	for	(y = 0; y < iHeight ; y++) {	
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		for(int x = 0; x < iWidth; x++) {
			// read B, R, G, convert it to float*
			fread(data_b,4,1,pfmFile);
			fread(data_r,4,1,pfmFile);
			fread(data_g,4,1,pfmFile);
			fr = *(float*)data_r; 
			fg = *(float*)data_g; 
			fb = *(float*)data_b; 
			// write image data to our buffer
			// write from bottom left
			double *output_pixel = GetPixel(x,iHeight -1 - y);
			output_pixel[0]=fb;
			output_pixel[1]=fr;
			output_pixel[2]=fg;
		}		
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
	fclose(pfmFile);
	pfmFile = 0;
	
	return 0;
}

int TMOImage::OpenEXR_16()
{
	Array2D<Rgba> pixels;
	
	RgbaInputFile *file = NULL;
	try {
		file = new RgbaInputFile(pName);
	} catch(...) {
		if (file) delete file;
		throw -2; // TODO: Co mam nahodit pri spatnem formatu souboru????
	}
	
    Box2i dw = file->dataWindow();
	
	iWidth = dw.max.x - dw.min.x + 1;
    iHeight = dw.max.y - dw.min.y + 1;
    pixels.resizeErase(iHeight, iWidth);
	
    file->setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * iWidth, 1, iWidth);
    file->readPixels(dw.min.y, dw.max.y);
	delete file;
	
	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];
	iFormat = TMO_RGB; // Je nutne nastavit pred volanim GetPixel!!!
	
	Rgba *act_pixel_input;
	double *act_pixel_output;
	ProgressBar(0, iHeight);
	for (int y = 0; y < iHeight; y++) {
		for (int x = 0; x < iWidth; x++) {
			act_pixel_input = &pixels[y][x];
			act_pixel_output = GetPixel(x, y);
			act_pixel_output[0] = act_pixel_input->r;
			act_pixel_output[1] = act_pixel_input->g;
			act_pixel_output[2] = act_pixel_input->b;
		}
		ProgressBar(y, iHeight);
	}
	
	dStonits = 179.0;
	fXres = 72.0f;
	fYres = 72.0f;
	iPhotometric = 0;
	
	return 0;
}

int TMOImage::OpenJPEG_32()
{
	FILE *file;
	struct jpeg_decompress_struct cinfo;
	
	FILE *subband_file;
	
	struct jpeg_decompress_struct subband_cinfo;
	
	JSAMPARRAY buffer;
	JSAMPARRAY subband_buffer;

	int row_stride; // physical row width in buffer
	int subband_row_stride; // physical row width in subband_buffer
	
	// subband vars
	float ln0, ln1, alpha, beta, samples2nits;
	float sb_dec[256];
	
	int i; // for-loop var
	
	float R, G, B, Y, Cb, Cr;
	
	double *out_data;
	
	int prograss_bar_part;
	
	int hdr;
	
	file = fopen(pName, "rb");
	if ( ! (hdr=JPEG_HDR_prepare_reading(file, &cinfo))) {
		fclose(file);
		throw -2;
	}
	
	iWidth = cinfo.output_width;
	iHeight = cinfo.output_height;
	iFormat = TMO_RGB;

	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];
	out_data = pData;

	if(hdr==1){//hdr jpeg
		try {
			subband_file = tmpfile();
			//subband_file = fopen("subband.jpg", "w+b");
			if ( ! JPEG_HDR_write_subband_image(subband_file, &cinfo, &ln0, &ln1, &alpha, &beta, &samples2nits)) {
				fprintf(stderr, "Error while reading JPEG-HDR: Cannot write tmp subband image or read JPEG-HDR subband header.\n");
				jpeg_destroy_decompress(&cinfo);
				
				jpeg_destroy_decompress(&subband_cinfo);
				fclose(subband_file);
				throw 1;
			}
			
			dStonits = samples2nits;
		
			JPEG_HDR_generate_sb_dec((float *)sb_dec, ln0, ln1);
		
			fseek(subband_file, 0, SEEK_SET);
			
			if ( ! JPEG_HDR_open_JPEG_file_and_prepare_to_reading(subband_file, &subband_cinfo, 0)) {
				fprintf(stderr, "Error while reading JPEG-HDR: Can't uncompress subband image.\n");
				jpeg_destroy_decompress(&cinfo);
				
				jpeg_destroy_decompress(&subband_cinfo);
		//			fclose(subband_file);
				throw 1;
			}
			
			if (cinfo.output_width != subband_cinfo.output_width || cinfo.output_height != subband_cinfo.output_height) {
				fprintf(stderr, "Error while reading JPEG-HDR: Resolution of imaga and subband image is not equal.");
				printf("Error while reading JPEG-HDR: Resolution of imaga and subband image is not equal.");
				jpeg_destroy_decompress(&cinfo);
				
				jpeg_destroy_decompress(&subband_cinfo);
				fclose(subband_file);
				throw 1;
			}
			
			row_stride = cinfo.output_width * cinfo.output_components;
			subband_row_stride = subband_cinfo.output_width * subband_cinfo.output_components;
			
			buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
			subband_buffer = (*subband_cinfo.mem->alloc_sarray) ((j_common_ptr) &subband_cinfo, JPOOL_IMAGE, subband_row_stride, 2);
			
			/* Step 6: while (scan lines remain to be read) */
			/*           jpeg_read_scanlines(...); */
		
			/* Here we use the library's state variable cinfo.output_scanline as the
			* loop counter, so that we don't have to keep track ourselves.
			*/
			ProgressBar(0, cinfo.output_height);
		
			while (cinfo.output_scanline < cinfo.output_height) {
				JSAMPROW buf_iter;
				JSAMPROW subband_buf_iter;
				
				/* jpeg_read_scanlines expects an array of pointers to scanlines.
				* Here the array is only one element long, but you could ask for
				* more than one scanline at a time if that's more convenient.
				*/
				(void) jpeg_read_scanlines(&cinfo, buffer, 1);
				(void) jpeg_read_scanlines(&subband_cinfo, subband_buffer, 1);
				/* Assume put_scanline_someplace wants a pointer and sample count. */
				
				buf_iter = buffer[0];
				subband_buf_iter = subband_buffer[0];
				for (i = 0; i < cinfo.output_width; i++) {
					Y = buf_iter[0] / 256.;
					Cb = buf_iter[1] / 256. - 0.5;
					Cr = buf_iter[2] / 256. - 0.5;
					
					R = sb_dec[*subband_buf_iter] * JPEG_HDR_inverse_gamma(Y                + 1.40200 * Cr);
					G = sb_dec[*subband_buf_iter] * JPEG_HDR_inverse_gamma(Y - 0.34414 * Cb - 0.71414 * Cr);
					B = sb_dec[*subband_buf_iter] * JPEG_HDR_inverse_gamma(Y + 1.77200 * Cb);
		
					JPEG_HDR_gamutDecompanding(R, G, B, alpha, beta, &out_data[0], &out_data[1], &out_data[2]);
		
					out_data += 3;
					buf_iter += 3;
					subband_buf_iter++;
				}
				
				ProgressBar(cinfo.output_scanline, cinfo.output_height);
			}
			ProgressBar(cinfo.output_height, cinfo.output_height);
		
			/* Step 7: Finish decompression */
		
			(void) jpeg_finish_decompress(&subband_cinfo);
			/* We can ignore the return value since suspension is not possible
			* with the stdio data source.
			*/
		
			/* Step 8: Release JPEG decompression object */
		
			/* This is an important step since it will release a good deal of memory. */
			jpeg_destroy_decompress(&subband_cinfo);
			/* After finish_decompress, we can close the input file.
			* Here we postpone it until after no more JPEG errors are possible,
			* so as to simplify the setjmp error logic above.  (Actually, I don't
			* think that jpeg_destroy can do an error exit, but why assume anything...)
			*/
		
			/* At this point you may want to check to see whether any corrupt-data
			* warnings occurred (test whether jerr.pub.num_warnings is nonzero).
			*/
			
			fclose(subband_file);
		} catch (int i_val) {
			throw -2;
		}
	}else{ //normal jpeg
		row_stride = cinfo.output_width * cinfo.output_components;
		buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
		
		/* Step 6: while (scan lines remain to be read) */
		/*           jpeg_read_scanlines(...); */
	
		/* Here we use the library's state variable cinfo.output_scanline as the
		* loop counter, so that we don't have to keep track ourselves.
		*/
		ProgressBar(0, cinfo.output_height);
	
		while (cinfo.output_scanline < cinfo.output_height) {
			JSAMPROW buf_iter;
			
			/* jpeg_read_scanlines expects an array of pointers to scanlines.
			* Here the array is only one element long, but you could ask for
			* more than one scanline at a time if that's more convenient.
			*/
			(void) jpeg_read_scanlines(&cinfo, buffer, 1);
			/* Assume put_scanline_someplace wants a pointer and sample count. */
			
			buf_iter = buffer[0];
			for (i = 0; i < cinfo.output_width; i++) {
				R = buf_iter[0] / 256.;
				G = buf_iter[1] / 256.;
				B = buf_iter[2] / 256.;
				
				*out_data++ = R;
				*out_data++ = G;
				*out_data++ = B;
				buf_iter += 3;
			}
			
			ProgressBar(cinfo.output_scanline, cinfo.output_height);
		}
		ProgressBar(cinfo.output_height, cinfo.output_height);
	}
	fXres = 72.0f;
	fYres = 72.0f;
	iPhotometric = 0;
	
	/* Step 7: Finish decompression */

	(void) jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);
	/* After finish_decompress, we can close the input file.
	* Here we postpone it until after no more JPEG errors are possible,
	* so as to simplify the setjmp error logic above.  (Actually, I don't
	* think that jpeg_destroy can do an error exit, but why assume anything...)
	*/

	/* At this point you may want to check to see whether any corrupt-data
	* warnings occurred (test whether jerr.pub.num_warnings is nonzero).
	*/

	fclose(file);

	return 0;
}

int TMOImage::OpenRAW_32()
{
	//determine the resolution from the filename first
	//filename format is: filename_XRESxYRES.hdrraw !
	int length=strlen(pName);
	char *resolution = new char[128];
	int underscore_pos=0, i=0, j=0; 

	for (underscore_pos = length - 1; underscore_pos > 0; underscore_pos--)
		if (pName[underscore_pos] == '_') break;
	if (!underscore_pos) throw -1;
	for (i = underscore_pos+1, j=0; i < length; i++, j++) {
		if(pName[i]!='.') resolution[j] = pName[i];//resolution[j] = tolower(pName[i]);
		else break;
	}
	resolution[j] = 0;

	sscanf(resolution, "%dx%d", &iWidth, &iHeight);
    //printf("Opening RAW file %s, Width=%d, Height=%d\n", pName, iWidth, iHeight);
	delete [] resolution;
	iFormat = TMO_RGB;

	FILE *f;
	unsigned char *scanline;
	int scanline_size=iWidth*3*4; //imageWidth * R,G,B * sizeof(float) = iWidth*3*4
	int x, y;

	if ((f = fopen(pName, "rb")) == 0) throw -1;

	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];

	
	if ((scanline = new unsigned char[scanline_size]) == 0) throw -3;

	for (y = 0; y < iHeight; y++) 
	{
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		//read scanline
		if ( fread(scanline, scanline_size, 1, f) != 1 ) throw -11;
		for (x = 0; x < iWidth; x++) 
		{
			GetPixel(x,y)[0] = (double) ((float*)scanline)[3*x];
			GetPixel(x,y)[1] = (double) ((float*)scanline)[3*x+1];
			GetPixel(x,y)[2] = (double) ((float*)scanline)[3*x+2];
		}
		if (ferror(f)) throw -5;
			
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
					
	delete[] scanline;
	fclose(f);
	return 0;
}


//Canon raw format (unsigned)
int TMOImage::OpenRAW_16()
{
	//determine the resolution from the filename first
	//filename format is: filename_XRESxYRES.raw !
	int length=strlen(pName);
	char *resolution = new char[128];
	int underscore_pos=0, i=0, j=0; 

	for (underscore_pos = length - 1; underscore_pos > 0; underscore_pos--)
		if (pName[underscore_pos] == '_') break;
	if (!underscore_pos) throw -1;
	for (i = underscore_pos+1, j=0; i < length; i++, j++) {
		if(pName[i]!='.') resolution[j] = pName[i];//resolution[j] = tolower(pName[i]);
		else break;
	}
	resolution[j] = 0;

	sscanf(resolution, "%dx%d", &iWidth, &iHeight);
    //printf("Opening RAW file %s, Width=%d, Height=%d\n", pName, iWidth, iHeight);
	delete [] resolution;
	iFormat = TMO_RGB;

	FILE *f;
	unsigned char *scanline;
	int scanline_size=iWidth*3*2; //imageWidth * R,G,B * sizeof(unsigned short) = iWidth*3*2
	int x, y;

	if ((f = fopen(pName, "rb")) == 0) throw -1;

	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];

	
	if ((scanline = new unsigned char[scanline_size]) == 0) throw -3;

	for (y = 0; y < iHeight; y++) 
	{
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		//read scanline
		if ( fread(scanline, scanline_size, 1, f) != 1 ) throw -11;
		for (x = 0; x < iWidth; x++) 
		{
			GetPixel(x,y)[0] = (double) ((unsigned short*)scanline)[3*x];
			GetPixel(x,y)[1] = (double) ((unsigned short*)scanline)[3*x+1];
			GetPixel(x,y)[2] = (double) ((unsigned short*)scanline)[3*x+2];
		}
		if (ferror(f)) throw -5;
			
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
					
	delete[] scanline;
	fclose(f);
	return 0;
}//raw16


// Reads information about the tiff file
int TMOImage::ReadHeader(TIFF* pFile)
{
	unsigned short iCompression = 0;		// compression type
	unsigned short iOrientation = 1;		// visual orientation (TIFF spec.)

	if (!TIFFGetField(pFile, TIFFTAG_COMPRESSION, &iCompression))
	{
		iCompression = COMPRESSION_NONE;
	}

	TIFFGetField(pFile, TIFFTAG_XRESOLUTION, &fXres);
	TIFFGetField(pFile, TIFFTAG_YRESOLUTION, &fYres);

	TIFFGetFieldDefaulted(pFile, TIFFTAG_ORIENTATION, &iOrientation);

	if (!TIFFGetFieldDefaulted(pFile, TIFFTAG_PHOTOMETRIC, &iPhotometric)) throw -2;
  
	switch (iPhotometric)
	{
    case PHOTOMETRIC_LOGLUV:
		TIFFSetField( pFile, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT );
		TIFFSetField( pFile, TIFFTAG_SGILOGENCODE, SGILOGENCODE_RANDITHER);
		break;
    case PHOTOMETRIC_LOGL:
		TIFFSetField( pFile, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		TIFFSetField( pFile, TIFFTAG_SGILOGENCODE, SGILOGENCODE_NODITHER);
		break;
	case PHOTOMETRIC_RGB:
		break;
    default:
		throw -3;
	}  

	if (!TIFFGetField(pFile, TIFFTAG_BITSPERSAMPLE, &iBitsPerSample)) {
		iBitsPerSample = 8;
	}
	
	if (!TIFFGetField(pFile, TIFFTAG_IMAGEWIDTH, &iWidth) ||
		!TIFFGetField(pFile, TIFFTAG_IMAGELENGTH, &iHeight)) throw -5;

	if ( !TIFFGetField(pFile, TIFFTAG_STONITS, &dStonits) )
		dStonits = 1.;

	if (pData) delete[] pData;
	pData = new double[3 * iWidth * iHeight];

	return 0;
}


// convert TIFF image to color parts
int TMOImage::ReadData(TIFF *pFile)
{ 
	int x, y, i, offset;
	float *scanline = new float[3 * iWidth];
	unsigned char *scanlineb = new unsigned char[3 * iWidth];
	// initialize conversion
	// covert whole picture - all scanlines
	
	switch (iPhotometric)
	{
	case PHOTOMETRIC_LOGLUV:
	case PHOTOMETRIC_LOGL:
		iFormat = TMO_XYZ;
		//scanline = new float[3 * iWidth];		
		break;
	case PHOTOMETRIC_RGB:
		iFormat = TMO_RGB;
		//scanlineb = new unsigned char[3 * iWidth];		
		break;
	}

	switch (iBitsPerSample)
	{
	case 32:		
		for (y = 0; y < iHeight; y++)
		{
			if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
			offset = y * iWidth;
			if ( TIFFReadScanline(pFile, scanline, y, 0) < 0 ) throw -8;
			for ( x = iWidth; x--; )
			{
				for ( i = 0; i < 3; i++ )
				{
					GetOffset(offset+x)[i] = scanline[3*x+i];
				}
			}
		}
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		break;
	case 8:				
		for (y = 0; y < iHeight; y++)
		{
			if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
			offset = y * iWidth;
			if ( TIFFReadScanline(pFile, scanlineb, y, 0) < 0 ) throw -8;
			for ( x = iWidth; x--; )
			{
				GetOffset(offset+x)[0] = (double)scanlineb[3*x]/255;
				GetOffset(offset+x)[1] = (double)scanlineb[3*x + 1]/255;
				GetOffset(offset+x)[2] = (double)scanlineb[3*x + 2]/255;
			}
		}
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		break;
	}
	if (scanline) delete[] scanline;
	if (scanlineb) delete[] scanlineb;
	return 0;
}

// Writes information about the TIFF file
int TMOImage::WriteHeader(TIFF* pFile, bool HDR)
{	
  TIFFSetField(pFile, TIFFTAG_PLANARCONFIG, 1 );
  TIFFSetField(pFile, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
  if(HDR) {
	TIFFSetField(pFile, TIFFTAG_BITSPERSAMPLE, 32, 32, 32);
  } else {
	TIFFSetField(pFile, TIFFTAG_BITSPERSAMPLE, 8, 8, 8);
  }
  TIFFSetField(pFile, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(pFile, TIFFTAG_XRESOLUTION, fXres);
  TIFFSetField(pFile, TIFFTAG_YRESOLUTION, fYres);
  TIFFSetField(pFile, TIFFTAG_ORIENTATION, 1);
  TIFFSetField(pFile, TIFFTAG_PHOTOMETRIC, 2);	//PHOTOMETRIC_RGB
  TIFFSetField(pFile, TIFFTAG_IMAGEWIDTH, iWidth);
  TIFFSetField(pFile, TIFFTAG_IMAGELENGTH, iHeight);
  return 0;
}

int TMOImage::WriteData(TIFF* pFile, bool HDR)
{ 
	Convert(TMO_RGB,false);

	int x, y, offset;
	double pixel[3];
	
	float *scanline = NULL;			// 32b HDR
	unsigned char *scanlineb = NULL;	// 8b LDR
	
	if(HDR) { // write 32b floats
		scanline = new float[3 * iWidth];
		for (y = 0; y < iHeight; y++)
		{
			if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
			offset = y * iWidth;
			for ( x = iWidth; x--; )
			{
				if (iFormat < 256)
				{
					pixel[0] = GetOffset(offset+x)[0];
					pixel[1] = GetOffset(offset+x)[1];
					pixel[2] = GetOffset(offset+x)[2];
					scanline[3 * x    ] = pixel[0];
					scanline[3 * x + 1] = pixel[1];
					scanline[3 * x + 2] = pixel[2];
				}
				else 
				{
					pixel[0] = GetOffset(offset+x)[0];
					if (pixel[0] < 0.0) pixel[0] = 0.0;
					if (pixel[0] > 1.0) pixel[0] = 1.0;
					scanline[3 * x    ] = 
						scanline[3 * x + 1] = 
						scanline[3 * x + 2] = (unsigned char)(255*pixel[0]);
				}
			}
			if ( TIFFWriteScanline(pFile, scanline, y, 0) < 0 ) throw -11;
		}
	} else {
		scanlineb = new unsigned char[3 * iWidth];
		for (y = 0; y < iHeight; y++)
		{
			if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
			offset = y * iWidth;
			for ( x = iWidth; x--; )
			{
				if (iFormat < 256)
				{
				pixel[0] = GetOffset(offset+x)[0];
				pixel[1] = GetOffset(offset+x)[1];
				pixel[2] = GetOffset(offset+x)[2];
				if (pixel[0] < 0.0) pixel[0] = 0.0;
				if (pixel[1] < 0.0) pixel[1] = 0.0;
				if (pixel[2] < 0.0) pixel[2] = 0.0;
				if (pixel[0] > 1.0) pixel[0] = 1.0;
				if (pixel[1] > 1.0) pixel[1] = 1.0;
				if (pixel[2] > 1.0) pixel[2] = 1.0;
				scanlineb[3 * x    ] = (unsigned char)(255*pixel[0]);
				scanlineb[3 * x + 1] = (unsigned char)(255*pixel[1]);
				scanlineb[3 * x + 2] = (unsigned char)(255*pixel[2]);
				}
				else 
				{
				pixel[0] = GetOffset(offset+x)[0];
				if (pixel[0] < 0.0) pixel[0] = 0.0;
				if (pixel[0] > 1.0) pixel[0] = 1.0;
				scanlineb[3 * x    ] = 
				scanlineb[3 * x + 1] = 
				scanlineb[3 * x + 2] = (unsigned char)(255*pixel[0]);
				}
			}
			if ( TIFFWriteScanline(pFile, scanlineb, y, 0) < 0 ) throw -11;
		}
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
        if (pFile) TIFFClose(pFile);
	pFile = 0;  
	if(scanline) delete[] scanline;
	if(scanlineb) delete[] scanlineb;
	return 0;
}


int TMOImage::WriteDataRAW(FILE* pFile)
{ 
	int x, y, offset;
	float *scanline;
	double pixel[3];
	
	scanline = new float[3 * iWidth];
	
	// convert whole picture - all scanlines
	for (y = 0; y < iHeight; y++)
	{
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		offset = y * iWidth;
		for ( x = iWidth; x--; )
		{
			if (iFormat < 256)//"RGB" formaty
			{
				pixel[0] = GetOffset(offset+x)[0];
				pixel[1] = GetOffset(offset+x)[1];
				pixel[2] = GetOffset(offset+x)[2];
				scanline[3 * x    ] = (float)pixel[0];
				scanline[3 * x + 1] = (float)pixel[1];
				scanline[3 * x + 2] = (float)pixel[2];
			}
			else //Y (samotny) format
			{
				pixel[0] = GetOffset(offset+x)[0];
				scanline[3 * x    ] = 
				scanline[3 * x + 1] = 
				scanline[3 * x + 2] = (float)pixel[0];
			}
		}
		if ( fwrite((char *)scanline, 4/*sizeof(float)*/, 3*iWidth, pFile) != (3*iWidth)) throw -11;
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
    if (pFile) fclose(pFile);
	pFile = 0;  
	delete[] scanline;
	return 0;
}//writeDataRAW

int TMOImage::SaveWithSuffix(const char *suffix, int fileFormat)
{
	int length = strlen(pName);
	int suffixlength = strlen(suffix);
	int point;
	int i;
	char *filename = new char[length + suffixlength + 1 + 128];
	char *oldName = new char[length+1];
	char *resolution = new char[128];

	strcpy(oldName, pName);

	for (point = length - 1; point > 0; point--)
		if (pName[point] == '.') break;
	if (!point) point = length; 
	for (i = 0; i < point; i++) filename[i] = pName[i];
	for (i = 0; i < suffixlength; i++) filename[point+i] = suffix[i];
	filename[point+i] = 0;
	//strcat(filename, &pName[point]);
	switch(fileFormat){
		case TMO_RAW_32: 
			sprintf(resolution, "_%dx%d.hdrraw", iWidth, iHeight);
			strcat(filename, resolution);
			break;
		case TMO_EXR_16:
			strcat(filename, ".exr");
			break;
		case TMO_PFM_32:
			strcat(filename, ".pfm");
			break;
		case TMO_HDR_32:
			strcat(filename, ".hdr");
			break;
		case TMO_JPEG_8:
			strcat(filename, ".jpeg");
			break;
		case TMO_TIFF_32: 
		case TMO_TIFF_8:
		default:
			strcat(filename, ".tif");
			break;
	}

	SaveAs(filename, fileFormat);
		
	if (pName) delete[] pName;
	pName = new char[strlen(oldName)+1];
	strcpy(pName, oldName);

	delete[] filename;
	delete[] oldName;
	delete[] resolution;
	return 0;
}


int TMOImage::SaveAs(const char *filename, int fileFormat)
{
	int length = strlen(filename);

	if (pName) delete[] pName;
	pName = new char[length+1];
	strcpy(pName, filename);

	Save(fileFormat);
	return 0;
}


int TMOImage::Save(int fileFormat)
{
	switch(fileFormat){
		case TMO_RAW_32:
			SaveRAW_32();
			break;
		case TMO_EXR_16:
			SaveEXR_16();
			break;
		case TMO_PFM_32:
			SavePFM_32();
			break;
		case TMO_HDR_32:
			SaveHDR_32();
			break;
		case TMO_JPEG_8:
			SaveJPEG_8();
			break;
		case TMO_TIFF_8: 
			SaveTIFF_8();
			break;
		case TMO_TIFF_32: 
			SaveTIFF_32();
			break;
		default:
			SaveTIFF_32();
			break;
	}
	return 0;
}


int TMOImage::SaveRAW_32()
{
	FILE* rawFile;
    if (( rawFile = fopen(pName, "wb")) == NULL) 
	{
		throw -10;
	}

	WriteDataRAW(rawFile);
	return 0;
}

// ***************************************************************************
// Save HDR TIFF image.
// Image is saved in RGB color space using LZW lossless compression.
// ***************************************************************************
int TMOImage::SaveTIFF_32()
{
	TIFF* pFile;
    if (( pFile = TIFFOpen(pName, "w")) == NULL) 
	{
		throw -10;
	}

	WriteHeader(pFile,true);
	WriteData(pFile,true);
	return 0;
}

// ***************************************************************************
// Save LDR TIFF image.
// Image is saved in RGB color space using LZW lossless compression.
// ***************************************************************************
int TMOImage::SaveTIFF_8()
{
	TIFF* pFile;
    if (( pFile = TIFFOpen(pName, "w")) == NULL) 
	{
		throw -10;
	}

	WriteHeader(pFile,false);
	WriteData(pFile,false);
	return 0;
}


// ***************************************************************************
// Save PFM image. Just color (no greyscale) images are considered.
// If image data isn't in RGB, backup is made, image is converted and saved,
// then original format is restored.
// ***************************************************************************
int TMOImage::SavePFM_32() {
	FILE* pfmFile;
    if (( pfmFile = fopen(pName, "wb")) == NULL) 
	{
		throw -10;
	}
	// ********** write header ***********
	// write 'PF'
	char wChar[100] = "PF\n";
	fwrite(wChar,1,strlen(wChar),pfmFile);
	// write width, length and scale/endianess
	sprintf(wChar,"%i",iWidth);
	strcat(wChar," ");
	fwrite(wChar,1,strlen(wChar),pfmFile);
	sprintf(wChar,"%i",iHeight);
	strcat(wChar,"\r\n");
	strcat(wChar, "-1.000000\n"); // scale/endianess
    fwrite(wChar, 1, strlen(wChar), pfmFile);
	// ********** end of header **********
	
	// ************************************************************
	// check current image format and (if necessary) convert to RGB
	double* pDataXYZ = NULL;
	int oldFormat = iFormat;
	// if format = RGB, do nothing, else backup data in original format and convert pData to RGB	
	if(iFormat!=TMO_RGB) {
		pDataXYZ = pData;
		pData = new double[3 * iWidth * iHeight];
		memcpy(pData, pDataXYZ, (3 * iWidth * iHeight) * sizeof(double));
		try {
			Convert(TMO_RGB, false);
		} catch (...) {
			delete[] pData;
			pData = pDataXYZ;
			iFormat = oldFormat;
			throw -10;
		}
	}
	// *************** end of format check ************************

	// ******************* write image data ***********************
	float dr, dg, db;
	int y;
	for	(y = 0; y < iHeight ; y++) {
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		for (int x = 0; x < iWidth; x++) {
			double *input_pixel = GetPixel(x,iHeight -1 - y);
			
			dr = (float)(input_pixel[0]);
			dg = (float)(input_pixel[1]);
			db = (float)(input_pixel[2]);
			
			fwrite(&dr, 4, 1, pfmFile);
			fwrite(&dg, 4, 1, pfmFile);
            fwrite(&db, 4, 1, pfmFile);
		}	
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;

	// restore data to original format
	if(oldFormat!=TMO_RGB) {
		delete[] pData;
		pData = pDataXYZ;
		iFormat = oldFormat;
	}

	fclose(pfmFile);
	pfmFile = 0;
	return 0;
}

// ***************************************************************************
// Save HDR image.
// Header: Just format and dimensions are written.
// Data: Written in run-length encoding.
// If image data isn't in RGB, backup is made, image is converted and saved,
// then original format is restored.
// ***************************************************************************
int TMOImage::SaveHDR_32() {
	FILE* hdrFile;
    if (( hdrFile = fopen(pName, "wb")) == NULL) {
		throw -10;
	}
	// ********** write header ***********
	char wChar[100] = "FORMAT=32-bit_rle_rgbe\n\n";
	fwrite(wChar,1,strlen(wChar),hdrFile);
	char wHeight[100] = "-Y ";
	sprintf(wChar,"%i ",iHeight);
	strcat(wHeight,wChar);
	fwrite(wHeight,1,strlen(wHeight),hdrFile);
	char wWidth[100] = "+X ";
	sprintf(wChar,"%i\n",iWidth);
	strcat(wWidth,wChar);
	fwrite(wWidth,1,strlen(wWidth),hdrFile);
	// ********** end of header **********

	// ************************************************************
	// check format and (if necessary) convert to RGB
	double* pDataXYZ = NULL;
	int oldFormat = iFormat;
	// if format = RGB, do nothing, else backup data in original format and convert pData to RGB	
	if(iFormat!=TMO_RGB) {
		pDataXYZ = pData;
		pData = new double[3 * iWidth * iHeight];
		memcpy(pData, pDataXYZ, (3 * iWidth * iHeight) * sizeof(double));
		try {
			Convert(TMO_RGB, false);
		} catch (...) {
			delete[] pData;
			pData = pDataXYZ;
			iFormat = oldFormat;
			throw -10;
		}
	}
	// *************** end of format check ************************

	// ******************* write image data ***********************
	COLR *scanin;
	COLOR col;
	int x, y;
	
	if ((scanin = (COLR*) new unsigned char[iWidth*4]) == 0) throw -3;

	for (y = 0; y < iHeight; y++) 
	{	
		if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
		for (x = 0; x < iWidth; x++) 
		{
			col[0] = GetPixel(x,y)[0];
			col[1] = GetPixel(x,y)[1];
			col[2] = GetPixel(x,y)[2];
			TMORadiance::setcolr(scanin[x], col[0], col[1], col[2]); // convert to hdr color format (r,g,b,e)
		}
		if (TMORadiance::fwritecolrs(scanin, iWidth, hdrFile) < 0) throw -4;
		if (ferror(hdrFile)) throw -5;		
	}
	if (y%10 == 0) if (ProgressBar(y, iHeight)==1) throw -19;
	
	// restore data to original format
	if(oldFormat!=TMO_RGB) {
		delete[] pData;
		pData = pDataXYZ;
		iFormat = oldFormat;
	}

	fclose(hdrFile);
	hdrFile = 0;
	delete[] scanin;
	return 0;
}

// ***********************************************************************************************
// Save JPEG LDR image.
// Image is tonemapped by drago03 algorithm (Drago, F., Myszkowski, K., Annen, T., Chiba, N.:
// Adaptive Logarithmic Mapping For Displaying High Contrast Scenes). 
// Concrete implementation is taken from PFSTmo package (http://www.mpi-inf.mpg.de/resources/tmo/).
// Implementation of saving jpeg itself according to http://www.smalleranimals.com/jpegfile.htm
// Default quality of saved image is set to '90'
// ***********************************************************************************************
int TMOImage::SaveJPEG_8(int quality) {
		
	// ************************************************************
	// need to tonemap in XYZ
	// backup orignal picture data into pDataOrig - !necessary! - data will be modified by tone mapping
	// the original picture data are restored after jpeg is saved
	double* pDataOrig = NULL;
	int oldFormat = iFormat;	
	pDataOrig = pData;
	pData = new double[3 * iWidth * iHeight];
	memcpy(pData, pDataOrig, (3 * iWidth * iHeight) * sizeof(double));
	// convert picture to XYZ
	try {
		Convert(TMO_XYZ, false);
	} catch (...) {
		delete[] pData;
		pData = pDataOrig;
		iFormat = oldFormat;
		throw -10;
	}
	// *************** end of backup and conversion ***************
		
	//***************** tonemap************************************
	// if tone mapping should not be used, comment out these lines
	// (from here to 'end of tonemap')
	// collect luminance information for tone mapping into Y
	float * Y = new float[iWidth * iHeight];
	for (int y = 0; y < iHeight; y++) {
		for(int x = 0; x < iWidth; x++) {
			Y[x+y*iWidth] = GetPixel(x,y)[1];
		}
	}	
	// output tonemapped luminance L
	float * L = new float[iWidth * iHeight];
	float maxLum,avLum;
	float biasValue = 0.85f;

	// calculate max and avg luminance for tone mapping
	calculateLuminance(Y, iHeight, iWidth, avLum, maxLum );
	// tonemap
	tmo_drago03(Y, L, iHeight, iWidth, maxLum, avLum, biasValue);
	
	// apply result of tonemapping on our image
	for (int y = 0; y < iHeight; y++) {
		for(int x = 0; x < iWidth; x++) {
			float scale = L[x+iWidth*y] / Y[x+iWidth*y];
			GetPixel(x,y)[0] *= scale;
			GetPixel(x,y)[1] = L[x+y*iWidth];
			GetPixel(x,y)[2] *= scale;
		}
	}
	
	// free memory
	L = NULL;
	Y = NULL;
	//***************** end of tonemap ****************************

	// save jpeg in RGB
	Convert(TMO_RGB, false);
	
	// convert floats to unsigned chars
	unsigned char *dataBuf = new unsigned char[3 * iWidth * iHeight];
	float pixel[3];
	for (int y = 0; y < iHeight; y++) {
		for(int x = 0; x < iWidth; x++) {
				pixel[0] = GetPixel(x,y)[0];
				pixel[1] = GetPixel(x,y)[1];
				pixel[2] = GetPixel(x,y)[2];
				if (pixel[0] < 0.0) pixel[0] = 0.0;
				if (pixel[1] < 0.0) pixel[1] = 0.0;
				if (pixel[2] < 0.0) pixel[2] = 0.0;
				if (pixel[0] > 1.0) pixel[0] = 1.0;
				if (pixel[1] > 1.0) pixel[1] = 1.0;
				if (pixel[2] > 1.0) pixel[2] = 1.0;
				dataBuf[3 * (y * iWidth + x)]	= (unsigned char)(255*pixel[0]);
				dataBuf[3 * (y * iWidth + x)+1] = (unsigned char)(255*pixel[1]);
				dataBuf[3 * (y * iWidth + x)+2] = (unsigned char)(255*pixel[2]);
		}
	}
	

	jpeg_error_mgr err;	
	struct jpeg_compress_struct cinfo;
	/* More stuff */
	FILE * outfile = NULL;			/* target file */
	int row_stride;			/* physical row widthPix in image buffer */
	

	/* Step 1: allocate and initialize JPEG compression object */
	cinfo.err = jpeg_std_error(&err);

	/* Now we can initialize the JPEG compression object. */
	jpeg_create_compress(&cinfo);

	/* Step 2: specify data destination (eg, a file) */
	/* Note: steps 2 and 3 can be done in either order. */

	if ((outfile = fopen(pName, "wb")) == NULL) {				
		throw -10;
	}

	jpeg_stdio_dest(&cinfo, outfile);

	/* Step 3: set parameters for compression */

	/* First we supply a description of the input image.
	* Four fields of the cinfo struct must be filled in:
	*/
	cinfo.image_width = iWidth; 	/* image widthPix and height, in pixels */
	cinfo.image_height = iHeight;
	cinfo.input_components = 3;		/* # of color components per pixel */
	cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */


	/* Now use the library's routine to set default compression parameters.
	* (You must set at least cinfo.in_color_space before calling this,
	* since the defaults depend on the source color space.)
	*/

	jpeg_set_defaults(&cinfo);
	/* Now you can set any non-default parameters you wish to.
	* Here we just illustrate the use of quality (quantization table) scaling:
	*/
	jpeg_set_quality(&cinfo, quality, true /* limit to baseline-JPEG values */);

	/* Step 4: Start compressor */

	/* TRUE ensures that we will write a complete interchange-JPEG file.
	* Pass TRUE unless you are very sure of what you're doing.
	*/
	jpeg_start_compress(&cinfo, true);

	/* Step 5: while (scan lines remain to be written) */
	/*           jpeg_write_scanlines(...); */

	/* Here we use the library's state variable cinfo.next_scanline as the
	* loop counter, so that we don't have to keep track ourselves.
	* To keep things simple, we pass one scanline per call; you can pass
	* more if you wish, though.
	*/
	row_stride = iWidth * 3;	/* JSAMPLEs per row in image_buffer */

	while (cinfo.next_scanline < cinfo.image_height) {
		/* jpeg_write_scanlines expects an array of pointers to scanlines.
		* Here the array is only one element long, but you could pass
		* more than one scanline at a time if that's more convenient.
		*/
		//LPBYTE outRow;
		unsigned char * outRow;

		outRow = dataBuf + (cinfo.next_scanline * iWidth * 3);


		(void) jpeg_write_scanlines(&cinfo, &outRow, 1);
	}

	/* Step 6: Finish compression */

	jpeg_finish_compress(&cinfo);

	/* After finish_compress, we can close the output file. */
	fclose(outfile);

	/* Step 7: release JPEG compression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_compress(&cinfo);

	
	
	// restore data to original format
	delete[] pData;
	pData = pDataOrig;
	iFormat = oldFormat;
	// free memory
	dataBuf = NULL;

	return 0;
}

int TMOImage::SaveEXR_16_pData_is_RGB() {
	Rgba *pixels = NULL;
		
	double *input_pixel;
	Rgba *output_pixel;
	
	int progress = 0;
	int progress_max = 2 * iHeight;
	ProgressBar(0, progress_max);
	try {
		pixels = new Rgba[iWidth * iHeight];
		
		output_pixel = &pixels[0];
		for	(int y = 0; y < iHeight; y++) {
			for (int x = 0; x < iWidth; x++) {
				double *input_pixel = GetPixel(x, y);
				
				output_pixel->r = input_pixel[0];
				output_pixel->g = input_pixel[1];
				output_pixel->b = input_pixel[2];
				output_pixel->a = 0;
				
				output_pixel++;
			}
			progress++;
			ProgressBar(progress, progress_max);
		}
		
		RgbaOutputFile file(pName, iWidth, iHeight, WRITE_RGBA);
		
		file.setFrameBuffer(pixels, 1, iWidth);
		
		for (int y = 10; y <= iHeight; y += 10) {
			file.writePixels(10);
			progress += 10;
			ProgressBar(progress, progress_max);
		}
		
		if (iHeight % 10 != 0) {
			file.writePixels(iHeight % 10);
		}
	} catch (...) {
		throw -10;
	}
	
	ProgressBar(progress_max, progress_max);

	if (pixels != NULL) {
		delete[] pixels;
	}
	
	return 0;
}

int TMOImage::SaveEXR_16()
{
	double* pDataXYZ = NULL;

	switch (iFormat) {
	case TMO_RGB:
		SaveEXR_16_pData_is_RGB();
		break;
	case TMO_XYZ:
		pDataXYZ = pData;
		
		pData = new double[3 * iWidth * iHeight];
		memcpy(pData, pDataXYZ, (3 * iWidth * iHeight) * sizeof(double));
		try {
			Convert(TMO_RGB, false);
			SaveEXR_16_pData_is_RGB();
		} catch (...) {
			delete[] pData;
			pData = pDataXYZ;
			iFormat = TMO_XYZ;
			
			throw -10;
		}
		
		delete[] pData;
		pData = pDataXYZ;
		iFormat = TMO_XYZ;
		
		break;
	default:
		throw -10;
	}

	return 0;
}

/*#ifdef ERIKLIB
int SaveImageset(imageset1_t *imageset, const char* filename, bool pyramid)
{
	int i, y, x, tmp_y, _xmax, _ymax;
	TMOImage temp;
	_xmax=imageset->width;
	_ymax=imageset->height;
	temp.New(_xmax, _ymax, TMO_Y);
	temp.SetFilename(filename);
	double *data=temp.GetData();

	char str[10];
	for(i=0; i<imageset->levels; i++){
		if(pyramid){
			_xmax=imageset->width>>i;
			_ymax=imageset->height>>i;
			temp.SetDimensions(_xmax, _ymax); 
		}
		for (y = 0; y < _ymax; y++){
			tmp_y=y*_ymax;
            for (x = 0; x < _xmax; x++)
				data[tmp_y + x]=imageset->images[i][y][x];
		}
		sprintf(str, "_level_%.2d", i);
		temp.SaveWithSuffix(str);
		temp.SaveWithSuffix(str, TMO_RAW);
	}
	temp.SetDimensions(imageset->width, imageset->height);
	return 0;
}//SaveImageset
#endif*/


int TMOImage::New(int width, int height, int format)
{
	Close();
	iFormat = format;
	iWidth = width;
	iHeight = height;

	if (pData) delete[] pData;
	pData=0;

	try {
		switch(format){
			case TMO_Y:
				pData = new double[iWidth * iHeight];
				break;
			default:
				pData = new double[3 * iWidth * iHeight];
				break;
		}
	} 
	catch(...) {
		printf("Memory allocation error!\n");
		exit(-1);
	}
	if(!pData){
		printf("Memory allocation error!\n");
		exit(-1);
	}
	return 0;
}

int TMOImage::New(const TMOImage &ref, int format, bool bCopy)
{
	int length; 
	int i;

	Close();
	iWidth = ref.iWidth;
	iHeight = ref.iHeight;
	fXres = ref.fXres;
	fYres = ref.fYres;
	dStonits = ref.dStonits;
	pProgressBar = ref.pProgressBar;
	pWriteLine = ref.pWriteLine;
	if (pName) 
	{
		delete[] pName;
		length = strlen(ref.pName);
		pName = new char[length+1];
		strcpy(pName, ref.pName);
	}

	switch(format)
	{
	case TMO_NOTSPEC:
		iFormat = ref.iFormat;
		if (pData) delete[] pData;
		if (iFormat < 256) length = 3 * iWidth * iHeight;
		else length = iWidth * iHeight;
		pData = new double[length];
		if (bCopy) for (i = 0; i < length; i++) pData[i] = ref.pData[i];
		break;
	case TMO_Y:
		iFormat = TMO_Y;
		if (pData) delete[] pData;
		length = iWidth * iHeight;
		pData = new double[length];
		for (i = 0; i < length; i++)
		{
			switch (ref.iFormat)
			{
			case TMO_Yxy:
				pData[i] = ref.pData[3*i];
				break;
			case TMO_XYZ:
				pData[i] = ref.pData[3*i+1];
				break;
			case TMO_RGB:
				pData[i] = 
					RGB2XYZ[1][0] * ref.pData[3 * i + 0] +
					RGB2XYZ[1][1] * ref.pData[3 * i + 1] +
					RGB2XYZ[1][2] * ref.pData[3 * i + 2];
				break;
			case TMO_Y:
				pData[i] = ref.pData[i];
				break;
			}
		}
		break;
	default:
		iFormat = ref.iFormat;
		if (pData) delete[] pData;
		if (iFormat < 256) length = 3 * iWidth * iHeight;
		else length = iWidth * iHeight;
		pData = new double[length];
		for (i = 0; i < length; i++) pData[i] = ref.pData[i];
		Convert(format);
		break;
	}
	return 0;
}

int TMOImage::Close()
{
	if (pName) delete[] pName;
	if (pData) delete[] pData;
	Clear();
	return 0;
}

double* TMOImage::GetData()
{
	return pData;
}

/**
 * converts one pixel in Lab to XYZ
 * 
 * @param L - L component in LAB (input)
 * @param a - a component in LAB (input)
 * @param b - b component in LAB (input)
 * @param x - x component in XYZ (output)
 * @param y - y component in XYZ (output)
 * @param z - z component in XYZ (output)
 */
void TMOImage::LabToXyz(double L, double a, double b, double * x, double * y, double * z){
	double fy = (L + 16.0) / 116.0;
	double fx = a / 500.0 + fy;
	double fz = fy - b / 200.0;
	
	double xr = (pow(fx, 3) > EPS) ? pow(fx, 3) : (116.0 * fx - 16.0) / KAPPA;
	double yr = (L > KAPPA * EPS) ? pow((L + 16.0) / 116.0, 3) : L / KAPPA;
	double zr = (pow(fz, 3) > EPS) ? pow(fz, 3) : (116.0 * fz - 16.0) / KAPPA;
	
	*x = xr * XYZ_WHITE_X;
	*y = yr * XYZ_WHITE_Y;
	*z = zr * XYZ_WHITE_Z;
}

/**
 * converts one pixel in XYZ to LUV
 * 
 * @param x - x component in XYZ (input)
 * @param y - y component in XYZ (input)
 * @param z - z component in XYZ (input)
 * @param L - L component in Luv (output)
 * @param u - u component in Luv (output)
 * @param v - v component in Luv (output)
 */
void TMOImage::XyzToLuv(double x, double y, double z, double * L, double * u, double * v){
	double u_line = (4 * x) / (x + 15 * y + 3 * z);
	double v_line = (9 * y) / (x + 15 * y + 3 * z);
	
	double yr = y / XYZ_WHITE_Y;
	
	double l = (yr > EPS) ? (116 * pow(yr, 1.0/3.0) - 16.0) : KAPPA * yr;
	*u = 13 * l * (u_line - LUV_WHITE_U);
	*v = 13 * l * (v_line - LUV_WHITE_V);	
	*L = l;
}

	/*TMO_RGB = 0,
	TMO_XYZ = 1,
	TMO_Yxy = 2,
	TMO_LCH = 3,
	TMO_LAB = 4,
	TMO_Y = 256,
	TMO_NOTSPEC,*/	

/**
 * Inverse sRGB Companding
 * 
 * @param color - one component of color
 * @return converted color
 */
double TMOImage::InverseSrgbCompanding(double color){
	return (color <= 0.04045) ? (color / 12.92) : pow((color + 0.055) / 1.055, 2.4);		
}

/**
 * sRGB Companding
 * 
 * @param color - one component of color
 * @return converted color
 */
double TMOImage::SrgbCompanding(double color){
	return (color <= 0.0031308) ? (12.92 * color) : (1.055 * pow(color, 1.0/2.4) - 0.055);		
}
	
int TMOImage::Convert(int format, bool fast)
{
	//std::cerr << "In TMOImage::Convert, format (dest): " << format << ", iFormat (src): " << iFormat << ", fast is: " << fast << std::endl; 
  
	int i, j, k, tmp_y;
	double pixel[3], W, X, Y, Z, EPSILON = 1e-06;

	if (fast) 
	{
		iFormat = format;
		return 0;
	}		

	if (iFormat == TMO_RGB)
	{
		if (format == TMO_RGB) return 1;
		if (format == TMO_XYZ) 
		{			
			double r, g, b;
			for (i = 0; i < iHeight; i++)
			{
				if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
				for (j = 0; j < iWidth; j++)
				{
					pixel[0] = pixel[1] = pixel[2] = 0.0;
					for (k = 0; k < 3; k++)
					{
						r = InverseSrgbCompanding(GetPixel(j,i)[0]);
						g = InverseSrgbCompanding(GetPixel(j,i)[1]);
						b = InverseSrgbCompanding(GetPixel(j,i)[2]);
						
						pixel[k] += RGB2XYZ[k][0] * r; 
						pixel[k] += RGB2XYZ[k][1] * g; 
						pixel[k] += RGB2XYZ[k][2] * b;
					}
					
					/*GetPixel(j,i)[0] = pixel[0];
					GetPixel(j,i)[1] = pixel[1];
					GetPixel(j,i)[2] = pixel[2];*/
					
					GetPixel(j,i)[0] = pixel[0] * 100;
					GetPixel(j,i)[1] = pixel[1] * 100;
					GetPixel(j,i)[2] = pixel[2] * 100;
					
					//std::cerr << "x: " << GetPixel(j,i)[0] << ", y: " << GetPixel(j,i)[1] << ", z: " << GetPixel(j,i)[2] << std::endl;
				}
			}
			if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
			iFormat = TMO_XYZ;						
			
			return 0;
		}
		if (format == TMO_Yxy)
		{
			Convert(TMO_XYZ);
			Convert(TMO_Yxy);
			return 0;
		}
		if (format == TMO_LCH)						// RGB to LCH
		{			
			Convert(TMO_LAB);					// RGB to LAB ... LAB to LCH
			double h, a, b, c;
			
			for (i = 0; i < iHeight; i++){				
				for (j = 0; j < iWidth; j++){			
					a = GetPixel(j,i)[1];
					b = GetPixel(j,i)[2];
					
					h = atan2(b, a);
					h = RadiansToDegrees(h);
					
					if (h < 0){
						h += 360.0;
					} else if (h > 360){
						h -= 360.0;
					}
					
					c = sqrt(pow(a, 2) + pow(b, 2));
				
					//GetPixel(j,i)[0] = GetPixel(j,i)[0];					
					GetPixel(j,i)[1] = (c > 100.0) ? 100.0 : c;					
					GetPixel(j,i)[2] = h;
				}
			}
			
			iFormat = TMO_LCH;
			return 0;			
		}
		if (format == TMO_LAB)							// RGB to LAB
		{			
			Convert(TMO_XYZ);						// RGB to XYZ ... XYZ to LAB
			double L, a, b, fx, fy, fz, x, y, z, xr, yr, zr;						
			
			for (i = 0; i < iHeight; i++){				
				for (j = 0; j < iWidth; j++){
					// get xyz values
					x = GetPixel(j,i)[0];
					y = GetPixel(j,i)[1];
					z = GetPixel(j,i)[2];
										
					xr = x / XYZ_WHITE_X;
					yr = y / XYZ_WHITE_Y;
					zr = z / XYZ_WHITE_Z;
					
					fx = (xr > EPS) ? pow(xr, 1.0 / 3.0) : ((KAPPA * xr + 16.0) / 116.0);
					fy = (yr > EPS) ? pow(yr, 1.0 / 3.0) : ((KAPPA * yr + 16.0) / 116.0);
					fz = (zr > EPS) ? pow(zr, 1.0 / 3.0) : ((KAPPA * zr + 16.0) / 116.0);
										
					L = (116.0 * fy - 16.0) > 0.0 ? 116.0 * fy - 16.0 : 0.0;
			        	a = 500 * (fx - fy);
					b = 200 * (fy - fz);
					
					// set LAB values
					GetPixel(j,i)[0] = L;
					GetPixel(j,i)[1] = a;
					GetPixel(j,i)[2] = b;
					
					// debug
					//std::cerr << "XYZ2LAB x: " << x << ", y: " << y << ", z: " << z << "L: " << L << ", a: " << a << ", b: " << b << std::endl;
				}
			}			
			
			iFormat = TMO_LAB;
			return 0;	
			
		}
	}
	if (iFormat == TMO_XYZ)
	{
		if (format == TMO_XYZ) return 1;
		if (format == TMO_RGB) 
		{			
			for (i = 0; i < iHeight; i++)
			{
				if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
				for (j = 0; j < iWidth; j++)
				{
					pixel[0] = pixel[1] = pixel[2] = 0.0;
					for (k = 0; k < 3; k++)
					{
						pixel[k] += XYZ2RGB[k][0] * (GetPixel(j,i)[0] / 100.0); 
						pixel[k] += XYZ2RGB[k][1] * (GetPixel(j,i)[1] / 100.0); 
						pixel[k] += XYZ2RGB[k][2] * (GetPixel(j,i)[2] / 100.0); 
					}
					
					//std::cerr << "r:" << pixel[0] << ", g:" << pixel[1] << ", b:" << pixel[2] << std::endl;
					
					GetPixel(j,i)[0] = SrgbCompanding(pixel[0]);
					GetPixel(j,i)[1] = SrgbCompanding(pixel[1]);
					GetPixel(j,i)[2] = SrgbCompanding(pixel[2]);										
				}
			}
			if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
			iFormat = TMO_RGB;
			return 0;
		}
		if (format == TMO_Yxy)
		{
			for (i = 0; i < iHeight; i++)
			{
				if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
				for (j = 0; j < iWidth; j++)
				{
					tmp_y = i * iWidth;
					pixel[0] = GetOffset(tmp_y+j)[0];
					pixel[1] = GetOffset(tmp_y+j)[1];
					pixel[2] = GetOffset(tmp_y+j)[2];
					if ((W = pixel[0] + pixel[1] + pixel[2]) > 0.) 
					{ 
						GetOffset(tmp_y+j)[0] = pixel[1];     // Y 
						GetOffset(tmp_y+j)[1] = pixel[0] / W;	// x 
						GetOffset(tmp_y+j)[2] = pixel[1] / W;	// y 	
					} 
					else
						GetOffset(tmp_y+j)[0] = 
						GetOffset(tmp_y+j)[1] = 
						GetOffset(tmp_y+j)[2] = 0.;
				}
			}
			if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
			iFormat = TMO_Yxy;
			return 0;
		}
	}
	if (iFormat == TMO_Yxy)
	{
		if (format == TMO_Yxy) return 1;
		if (format == TMO_XYZ)
		{
			for (i = 0; i < iHeight; i++)
			{
				if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
				for (j = 0; j < iWidth; j++)
				{
					tmp_y = i * iWidth;
					Y = GetOffset(tmp_y+j)[0];				//Y
					pixel[1] = GetOffset(tmp_y+j)[1];		//x
					pixel[2] = GetOffset(tmp_y+j)[2];		//y
					if ((Y > EPSILON) && (pixel[1] > EPSILON) && (pixel[2] > EPSILON)) 
					{ 
				        X = (pixel[1] * Y) / pixel[2];
						Z = (X / pixel[1]) - X - Y;
					} 
					else
						X = Z = EPSILON;
					GetOffset(tmp_y+j)[0] = X;
					GetOffset(tmp_y+j)[1] = Y;
					GetOffset(tmp_y+j)[2] = Z;
				}
			}
			if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
			iFormat = TMO_XYZ;
			return 0;
		}
		if (format == TMO_RGB)
		{
			Convert(TMO_XYZ);
			Convert(TMO_RGB);
			return 0;
		}
	}
	if (iFormat == TMO_LCH)						// LCH to ...
	{
		if (format == TMO_LCH) return 1;
		if (format == TMO_RGB){
			Convert(TMO_LAB);
			Convert(TMO_XYZ);
			Convert(TMO_RGB);
			return 0;
		}
		if (format == TMO_LAB){					// LCH to LAB			
			double h, c;
			for (i = 0; i < iHeight; i++){
				if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
				for (j = 0; j < iWidth; j++){	
					h = GetPixel(j,i)[2] * M_PI / 180.0;			// h in radians					
					c = GetPixel(j,i)[1];
					
					GetPixel(j,i)[0] = GetPixel(j,i)[0];			// L
					GetPixel(j,i)[1] = cos(h) * c;			// a
					GetPixel(j,i)[2] = sin(h) * c;			// b										
				}
			}
			if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
			iFormat = TMO_LAB;
			return 0;					
		}
		
	}	
	if (iFormat == TMO_LAB)							// LAB to ...
	{
		if (format == TMO_LAB) return 1;
		if (format == TMO_XYZ){						// LAB to XYZ
			double fx, fy, fz, L, a, b, xr, yr, zr, x, y, z;
			
			for (i = 0; i < iHeight; i++){				
				for (j = 0; j < iWidth; j++){												
					L = GetPixel(j,i)[0];
					a = GetPixel(j,i)[1];
					b = GetPixel(j,i)[2];
					
					LabToXyz(L, a, b, &x, &y, &z);
					
					GetPixel(j,i)[0] = x;
					GetPixel(j,i)[1] = y;
					GetPixel(j,i)[2] = z;
				}
			}	
			
			iFormat = TMO_XYZ;
			return 0;			
		}
		if (format == TMO_LCH){						// LAB to LCH
			double h, a, b, c;
			
			for (i = 0; i < iHeight; i++){				
				for (j = 0; j < iWidth; j++){			
					a = GetPixel(j,i)[1];
					b = GetPixel(j,i)[2];
					
					h = atan2(b, a);
					h = RadiansToDegrees(h);
					
					if (h < 0){
						h += 360.0;
					} else if (h > 360){
						h -= 360.0;
					}
					
					c = sqrt(pow(a, 2) + pow(b, 2));
				
					//GetPixel(j,i)[0] = GetPixel(j,i)[0];					
					GetPixel(j,i)[1] = c;
					GetPixel(j,i)[2] = h;
				}
			}			
			
			iFormat = TMO_LCH;
			return 0;			
		}
		if (format == TMO_RGB){						// LAB to RGB
			Convert(TMO_XYZ);
			Convert(TMO_RGB);
			return 0;
		}
	}		
		
	std::cerr << "Error, color space conversion implemented yet." << std::endl;
	return -1;
}

double TMOImage::GetStonits()
{
	return dStonits;
}

TMOImage& TMOImage::operator=(const TMOImage& image)
{
	int i, length;
	if (image.iFormat == iFormat)
	{
		if (iFormat < 256) length = 3 * iWidth * iHeight;
		else length = iWidth * iHeight;
		for (i = 0; i < length; i++) pData[i] = image.pData[i];
		return *this;
	}
	if (image.iFormat == TMO_Y)
	{
		length = iWidth * iHeight;
		for (i = 0; i < length; i++) 
		{
			switch(iFormat)
			{
			case TMO_Yxy:
				pData[3*i] = image.pData[i];
				break;
			case TMO_XYZ:
				pData[3*i+1] = image.pData[i];
				break;
			}
		}
	}
	New(image);
	return *this;
}


double** TMOImage::ExportLuminance()
{
	int y, x, tmp_y, format=iFormat;
	double **data;

	data = (double **) malloc (iHeight * sizeof (double*));
	for (y = 0; y < iHeight; y++)
		data[y] = (double *) malloc (iWidth * sizeof (double));

	if (iFormat < 256){
		Convert(TMO_Yxy);
		for (y = 0; y < iHeight; y++){
			tmp_y=y*iWidth;
			for (x = 0; x < iWidth; x++)
				data[y][x] = pData[3*(tmp_y + x)];
		}
		Convert(format);//convert back
	}
	else{
		for (y = 0; y < iHeight; y++){
			tmp_y=y*iWidth;
			for (int x = 0; x < iWidth; x++)
				data[y][x] = pData[tmp_y + x];
		}
	}
	return data;
}//ExportLuminance


int TMOImage::ImportLuminance(double **data, int width, int height)
{
	int y, x, tmp_y;

	New(width, height, TMO_Y);

	for (y = 0; y < iHeight; y++){
		tmp_y=y*iWidth;
		for (x = 0; x < iWidth; x++)
			pData[tmp_y + x]=data[y][x];
	}

	return 0;
}

/*
double* TMOImage::GetPixel(int x, int y)
{
	return (iFormat < 0xFF) ? &pData[3 * (y * iWidth + x)] : &pData[y * iWidth + x];
}

double* TMOImage::GetOffset(int offset)
{
	return (iFormat < 0xFF) ? &pData[3 * offset] : &pData[offset];
}
*/
int TMOImage::SetProgress(int(*pFunction)(TMOImage*, int, int))
{
	pProgressBar = pFunction;
	return 0;
}

int TMOImage::SetWriteLine(int(*pFunction)(TMOImage*, const wchar_t*))
{
	pWriteLine = pFunction;
	return 0;
}

int TMOImage::DefaultProgressBar(TMOImage*, int part, int all)
{
	if (part==all) wprintf (L"100%%\n");
	else wprintf (L"%2i%%\b\b\b", 100*part/all);
	return 0;
}

int TMOImage::DefaultWriteLine(TMOImage*, const wchar_t* text)
{
	wprintf(L"%s\n", text);
	return 0;
}

int TMOImage::Clear()
{
	dStonits = 179.0;
	fXres = 72.0f;
	fYres = 72.0f;
	iWidth = 0;
	iHeight = 0;
	iFormat = TMO_NOTSPEC;
	pData = 0;
	pName = new char[11];
	iPhotometric = 0;		
	strcpy(pName, "output.tif");
	return 0;
}

int TMOImage::ProgressBar(int part, int all)
{
	if (pProgressBar) return pProgressBar(this, part, all);
	return 0;
}
int TMOImage::WriteLine(const wchar_t* text)
{
	if (pWriteLine) return pWriteLine(this, text);
	return 0;
}


/////////////////////////////////////////////////////////////////////
//MC

int TMOImage::CorrectGamma(double gamma=2.2){
	if(gamma==1.0) return(0);
	assert(iFormat == TMO_RGB);

	printf("Correcting gamma... ");
	double invGamma=1/gamma;
    double R=0., G=0., B=0.;
	int tmp_y=0;
	for (int i = 0; i < iHeight; i++)
	{
		if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
		for (int j = 0; j < iWidth; j++)
		{
			tmp_y = i * iWidth;
			R = GetOffset(tmp_y+j)[0];
			G = GetOffset(tmp_y+j)[1];
			B = GetOffset(tmp_y+j)[2];
			GetOffset(tmp_y+j)[0] = pow(R, invGamma);
			GetOffset(tmp_y+j)[1] = pow(G, invGamma);
			GetOffset(tmp_y+j)[2] = pow(B, invGamma);
		}
	}
	ProgressBar(1, 1);
	return 0;
}

//it does not make sense to gamma correct just luminance - just for testing!!
int TMOImage::CorrectGammaYxy(double gamma=2.2){
	if(gamma==1.0) return(0);
	assert(iFormat == TMO_Yxy);

	printf("Correcting gamma... ");
	double invGamma=1/gamma;
    double Y=0.;
	int tmp_y=0;
	for (int i = 0; i < iHeight; i++)
	{
		if (i%10 == 0) if (ProgressBar(i, iHeight)==1) throw -19;
		for (int j = 0; j < iWidth; j++)
		{
			tmp_y = i * iWidth;
			Y = GetOffset(tmp_y+j)[0];
			GetOffset(tmp_y+j)[0] = pow(Y, invGamma);
		}
	}
	ProgressBar(1, 1);
	return 0;
}


int TMOImage::SetFilename(const char *filename)
{
	int length = strlen(filename);

	if (pName) delete[] pName;
	pName = new char[length+1];
	strcpy(pName, filename);

	return 0;
}

int TMOImage::GetMinMaxAvg(double *minimum, double *maximum, double *average)
{
	double suma = 0., tmp = 0.;
	int tmp_y, i, j;

	Convert(TMO_RGB);

	*minimum = 1.7E+308;
	*maximum = 0.265068 * GetOffset(0)[0] + 0.67023428 * GetOffset(0)[1] + 0.06409157 * GetOffset(0)[2];
	suma = 0;
			
	for ( i = 0; i < iHeight; i++ )
	{
		tmp_y = i*iWidth;
		for ( j = 0; j < iWidth; j++ )
		{
			tmp  = 0.265068 * GetOffset(j + tmp_y)[0];
			tmp += 0.67023428 * GetOffset(j + tmp_y)[1];
			tmp += 0.06409157 * GetOffset(j + tmp_y)[2];
			suma += tmp;
			if ( tmp > *maximum ) *maximum = tmp;
			if ( tmp < *minimum ) *minimum = tmp;
		}
	}
	*average = suma / (iWidth*iHeight);
	return 0;
}


int TMOImage::GetMinMaxAvgLog10(double *minimum, double *maximum, double *average)
{
	double suma = 0., tmp = 0.;
	int tmp_y, i, j;

	Convert(TMO_RGB);

	*minimum = 1.7E+308;
	*maximum = TAKE_LOG10(0.265068 * GetOffset(0)[0] + 0.67023428 * GetOffset(0)[1] + 0.06409157 * GetOffset(0)[2]);
	suma = 0;
			
	for ( i = 0; i < iHeight; i++ )
	{
		tmp_y = i*iWidth;
		for ( j = 0; j < iWidth; j++ )
		{
			tmp  = 0.265068 * GetOffset(j + tmp_y)[0];
			tmp += 0.67023428 * GetOffset(j + tmp_y)[1];
			tmp += 0.06409157 * GetOffset(j + tmp_y)[2];
			tmp = TAKE_LOG10(tmp);
			suma += tmp;
			if ( tmp > *maximum ) *maximum = tmp;
			if ( tmp < *minimum ) *minimum = tmp;
		}
	}
	*average = suma / (iWidth*iHeight);
	return 0;
}
double TMOImage::GetLuminance(int x, int y, int r)
{
	int i, j, xp, yp, count = 0, offset;
	double retval = .0;

	if (!r) 
	{
		offset = y * iWidth + x;
		if ((x < 0) || (x >= iWidth)) return .0;
		if ((y < 0) || (y >= iHeight)) return .0;
		retval += 0.299 * GetOffset(offset)[0];
		retval += 0.587 * GetOffset(offset)[1];
		retval += 0.114 * GetOffset(offset)[2];
		return retval;
	}
	
	for (i = -r; i <= r; i++)
		for (j = -r; j <= r; j++)
		{
			xp = x + j;
			yp = y + i;
			offset = yp * iWidth + xp;
			if ((xp < 0) || (xp >= iWidth)) continue;
			if ((yp < 0) || (yp >= iHeight)) continue;
			retval += 0.299 * GetOffset(offset)[0];
			retval += 0.587 * GetOffset(offset)[1];
			retval += 0.114 * GetOffset(offset)[2];
			count++;
		}
	return retval / count;
}

void TMOImage::GetStatisticsCircleRange(int x, int y, int r)
{
	int i, j, xp, yp, count = 0, offset, tmp;
	double retval = .0, curLum = .0, red = .0, green = .0, blue = .0, min = 1000000.0, max = .00000000001;

	offset = y * iWidth + x;
	if (!(((x < 0) || (x >= iWidth)) || ((y < 0) || (y >= iHeight))))
	{
		curLum = 0.299 * GetOffset(offset)[0];
		curLum += 0.587 * GetOffset(offset)[1];
		curLum += 0.114 * GetOffset(offset)[2];
		statistics.curLum = curLum;
		statistics.color[0] = GetPixel(x, y)[0];
		statistics.color[1] = GetPixel(x, y)[1];
		statistics.color[2] = GetPixel(x, y)[2];
	}
	double r2 = r * r;
	for (i = 0; i <= r; i++)
	{
		int r1D = sqrt( r2 - i * i);
		for (j = -r; j <= r; j++)
		{
			xp = x + j;
			yp = y + i;
			tmp = abs(j);
			if ((xp < 0) || (xp >= iWidth) || ( tmp > r1D)) continue;
			if ((yp < 0) || (yp >= iHeight)) continue;
			offset = yp * iWidth + xp;
			red += GetOffset(offset)[0];
			curLum = 0.299 * GetOffset(offset)[0];
			green += GetOffset(offset)[1];
			curLum += 0.587 * GetOffset(offset)[1];
			blue += GetOffset(offset)[2];
			curLum += 0.114 * GetOffset(offset)[2];
			retval += curLum;
			if(curLum < min)
				min = curLum;
			if(curLum > max)
				max = curLum;
			count++;			
			if(i != 0)
			{
				yp = y - i;
				if ((yp < 0) || (yp >= iHeight)) continue;
				offset = yp * iWidth + xp;
				red += GetOffset(offset)[0];
				curLum = 0.299 * GetOffset(offset)[0];
				green += GetOffset(offset)[1];
				curLum += 0.587 * GetOffset(offset)[1];
				blue += GetOffset(offset)[2];
				curLum += 0.114 * GetOffset(offset)[2];
				retval += curLum;
				if(curLum < min)
					min = curLum;
				if(curLum > max)
					max = curLum;
				count++;			
			}
		}
	}
	statistics.avgLum = retval / count;
	statistics.avgColor[0] = red / count;
	statistics.avgColor[1] = green / count;
	statistics.avgColor[2] = blue / count;
	statistics.minLum = min;
	statistics.maxLum = max;
}

void TMOImage::GetStatisticsSquareRange(int x, int y, int size)
{
	int i, j, xp, yp, count = 0, offset;
	double retval = .0, curLum = .0, red = .0, green = .0, blue = .0, min = 1000000.0, max = .00000000001;

	offset = y * iWidth + x;
	if (!(((x < 0) || (x >= iWidth)) || ((y < 0) || (y >= iHeight))))
	{
		curLum = 0.299 * GetOffset(offset)[0];
		curLum += 0.587 * GetOffset(offset)[1];
		curLum += 0.114 * GetOffset(offset)[2];
		statistics.curLum = curLum;
		statistics.color[0] = GetPixel(x, y)[0];
		statistics.color[1] = GetPixel(x, y)[1];
		statistics.color[2] = GetPixel(x, y)[2];
	}
	for (i = -size; i <= size; i++)
		for (j = -size; j <= size; j++)
		{
			xp = x + j;
			yp = y + i;
			if ((xp < 0) || (xp >= iWidth)) continue;
			if ((yp < 0) || (yp >= iHeight)) continue;
			offset = yp * iWidth + xp;
			red += GetOffset(offset)[0];
			curLum = 0.299 * GetOffset(offset)[0];
			green += GetOffset(offset)[1];
			curLum += 0.587 * GetOffset(offset)[1];
			blue += GetOffset(offset)[2];
			curLum += 0.114 * GetOffset(offset)[2];
			retval += curLum;
			if(curLum < min)
				min = curLum;
			if(curLum > max)
				max = curLum;
			count++;			
	}
	statistics.avgLum = retval / count;
	statistics.avgColor[0] = red / count;
	statistics.avgColor[1] = green / count;
	statistics.avgColor[2] = blue / count;
	statistics.minLum = min;
	statistics.maxLum = max;
}

double TMOImage::GetLuminanceYxy(int x, int y, int r)
{
	int i, j, xp, yp, count = 0, offset;
	double retval = .0;

	if (!r) 
	{
		offset = y * iWidth + x;
		if ((x < 0) || (x >= iWidth)) return .0;
		if ((y < 0) || (y >= iHeight)) return .0;
		retval = GetOffset(offset)[0];
		return retval;
	}
	
	for (i = -r; i <= r; i++)
		for (j = -r; j <= r; j++)
		{
			xp = x + j;
			yp = y + i;
			offset = yp * iWidth + xp;
			if ((xp < 0) || (xp >= iWidth)) continue;
			if ((yp < 0) || (yp >= iHeight)) continue;
			retval += GetOffset(offset)[0];
			count++;
		}
	return retval / count;
}

double TMOImage::SetLuminance(int x, int y, double L)
{
	double retval = GetLuminance(x, y);
	double r;
	int offset;

	if ((x < 0)||(x > iWidth)) return 0;
	if ((y < 0)||(y > iHeight)) return 0;

	r = L / retval;
	offset = y * iWidth + x;

	GetOffset(offset)[0] *= r;
	GetOffset(offset)[1] *= r;
	GetOffset(offset)[2] *= r;

	return retval;	
}

int TMOImage::GetDimensions(int *x, int *y)
{
	if (x) *x = iWidth;
	if (y) *y = iHeight;
	return 0;
}
