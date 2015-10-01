// freejpeghdr.h: interface for JPEG-HDR images.
//
//////////////////////////////////////////////////////////////////////

#ifndef _FREEJPEGHDR_H
#define _FREEJPEGHDR_H

int JPEG_HDR_open_JPEG_file_and_prepare_to_reading(FILE *infile, struct jpeg_decompress_struct *pcinfo, int with_APP11_marker);
int JPEG_HDR_write_subband_image(FILE *outfile, struct jpeg_decompress_struct *pcinfo,
						float *pln0, float *pln1, float *palpha, float *pbeta, float *psamples2nits);
void JPEG_HDR_gamutDecompanding(double r, double g, double b_, double alpha, double beta, double *red, double *green, double *blue);
double JPEG_HDR_inverse_gamma(double val);
void JPEG_HDR_generate_sb_dec(float (*sb_dec),float ln0, float ln1);

int JPEG_HDR_prepare_reading(FILE *infile, struct jpeg_decompress_struct *pcinfo);

#endif /* _FREEJPEGHDR_H */
