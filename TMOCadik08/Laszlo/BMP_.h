
//BMP Header Write

//BGR a szinsorrend egy pixeln belul

//BGR a pixeleken beluli sorrend
//soronkent a sorvegen dummy szamu tetszoleges toltelek karakter


typedef struct tagBITMAPFILEHEADER {    
	     //UINT    bfType;
         DWORD   bfSize;
         UINT    bfReserved1;
         UINT    bfReserved2;
         DWORD   bfOffBits;
       } BITMAPFILEHEADER;
/*
         The BITMAPFILEHEADER structure contains information about the type, size, and
         layout of a device-independent bitmap (DIB) file.

         Member          Description

         bfType          Specifies the type of file. This member must be BM.
         bfSize          Specifies the size of the file, in bytes.
         bfReserved1     Reserved; must be set to zero.
         bfReserved2     Reserved; must be set to zero.
         bfOffBits       Specifies the byte offset from the BITMAPFILEHEADER structure
         to the actual bitmap data in the file.


typedef struct tagBITMAPFILEHEADER1
{
//  UINT1    bfType;
//  char    bfType[2];      // = "BM"
    DWORD   bfSize;         // file size
    UINT1    bfReserved1;    // = 0
    UINT1    bfReserved2;    // = 0
    DWORD   bfOffBits;      //  offset from here to pixel bytes
} BITMAPFILEHEADER1;
*/

typedef struct tagBITMAPINFOHEADER1
{
    DWORD   biSize;	      // = sizeof(BITMAPINFOHEADER)
    LONG    biWidth;
    LONG    biHeight;
    WORD    biPlanes;	      // = 1
    WORD    biBitCount;       // = 24
    DWORD   biCompression;    // = 0
    DWORD   biSizeImage;
    LONG    biXPelsPerMeter;  // X pixels per m
    LONG    biYPelsPerMeter;  // Y pixels per m
    DWORD   biClrUsed;	      // = 0
    DWORD   biClrImportant;   // = 0
} BITMAPINFOHEADER1;

// --------------------------- into SOURCE file -------------------------------
void BMP_WRITE_HEADER_24(long x_res,
						 long y_res,
						 long * dummy,
						 FILE * fout)
{
//24 bpp bmp header
//to write dummy times empty (arbitrary) char after each rows

  char    bfType[2];	    // = "BM"
  BITMAPFILEHEADER bmphead;
  BITMAPINFOHEADER1 bmpinfo;

    *dummy = (4 - (3 * x_res) %4 ) % 4;

	bfType[0] = 'B';
	bfType[1] = 'M';

	//bmphead.bfType = "BM";

	bmphead.bfSize =	54 + (x_res*3+*dummy)*y_res;
						// [y_res * (soronkenti byte szam)]
	bmphead.bfReserved1 = 0;
	bmphead.bfReserved2 = 0;
	bmphead.bfOffBits   = 54; //sizeof(BITMAPFILEHEADER1)+sizeof(BITMAPINFOHEADER1);

	bmpinfo.biSize = 40; //sizeof(BITMAPINFOHEADER1);
	bmpinfo.biWidth  = x_res;
	bmpinfo.biHeight = y_res;
	bmpinfo.biPlanes = 1;
	bmpinfo.biBitCount = 24; //color_bit_num; // ff:1  //8 bit : 8
	bmpinfo.biCompression = 0;
	bmpinfo.biSizeImage = (x_res*3+*dummy)*y_res;
	bmpinfo.biXPelsPerMeter = 12500;
	bmpinfo.biYPelsPerMeter = 12500;
	bmpinfo.biClrUsed = 0;
	bmpinfo.biClrImportant = 0;

	fwrite(&bfType[0],1,1,fout);
	fwrite(&bfType[1],1,1,fout);
	fwrite(&bmphead, 12, 1, fout);
	fwrite(&bmpinfo, 40, 1, fout);

}