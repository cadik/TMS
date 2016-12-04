// freejpeghdr.cpp: implementation of JPEG-HDR images.
//
//////////////////////////////////////////////////////////////////////

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <jpeglib.h>
#include <setjmp.h>
} // extern "C"
#include "freejpeghdr.h"


struct my_error_mgr {
    struct jpeg_error_mgr pub; /* "public" fields */
    jmp_buf setjmp_buffer;	/* for return to caller  */
};

typedef struct my_error_mgr *my_error_ptr;

struct my_error_mgr jerr;

METHODDEF(void) my_error_exit(j_common_ptr cinfo) {
    /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
    my_error_ptr myerr = (my_error_ptr) cinfo->err;

    /* Always display the message.
     * We could postpone this until after returning, if we chose.
	 */
    (*cinfo->err->output_message) (cinfo);

    /* Return control to the setjmp point */
    longjmp(myerr->setjmp_buffer, 1);
}

void JPEG_HDR_generate_sb_dec(float (*sb_dec),float ln0, float ln1) {
	int i;
	for (i  = 0; i < 256; i++) {
		sb_dec[i] = exp(((((float) i) + 0.5) / 256) * (ln1 - ln0) + ln0);
	}
}

/**
 * Reads JPEG-HDR subband header.
 * @param str string with subband header.
 * @param pln0 ln0 value from header.
 * @param pln1 ln1 value from header.
 * @param palpha alp value from header.
 * @param pbeta bet value from header.
 * @param psample2nits s2n value from header.
 * @retval 0 Error.
 * @retval 1 Ok.
 */
int read_JPEG_HDR_header(char *str, float *pln0, float *pln1, float *palpha, float *pbeta, float *psamples2nits) {
	char *pos;
	char *str_arr[] = { "ln0=",  "ln1=", "alp=", "bet=", "s2n=" };

	int len = sizeof(str_arr) / sizeof(char *);
	int i;
	for (i = 0; i < len; i++) {
		pos = strstr(str, str_arr[i]);
		
		if (pos == NULL) return 0;
		
		pos += strlen(str_arr[i]);
		
		switch (i) {
		case 0:
			*pln0 = atof(pos);
			break;
		case 1:
			*pln1 = atof(pos);
			break;
		case 2:
			*palpha = atof(pos);
			break;
		case 3:
			*pbeta = atof(pos);
			break;
		case 4:
			*psamples2nits = atof(pos);
			break;
		default:
			return 0;
		}
	}
	
	return 1;
}

 /**
 * Reads from marker_list in pcinfo and writes JPEG-HDR subband image into
 * opened file outfile.
 * @param outfile opened output file for subband image.
 * @param pcinfo
 * @param pln0 ln0 value from header.
 * @param pln1 ln1 value from header.
 * @param palpha alp value from header.
 * @param pbeta bet value from header.
 * @param psample2nits s2n value from header.
 * @retval 0 Error.
 * @retval 1 Ok.
 * @see read_JPEG_HDR_header()
 */
int JPEG_HDR_write_subband_image(FILE *outfile, struct jpeg_decompress_struct *pcinfo,
						float *pln0, float *pln1, float *palpha, float *pbeta, float *psamples2nits) {
	struct jpeg_marker_struct *marker;
	int text_data_len;

	if (outfile == NULL) {
		return 0;
    }
	
	for (marker = pcinfo->marker_list; marker != NULL; marker = marker->next) {

		if (strncmp((char *) marker->data, "HDR_RI ver=", 11) == 0) {
			if ( ! read_JPEG_HDR_header((char *)marker->data, pln0, pln1, palpha, pbeta, psamples2nits)) {
				return 0;
			}
		}

		text_data_len = 1 + strlen((char *)marker->data);
		
		fwrite(marker->data + text_data_len, sizeof(JOCTET), marker->data_length - text_data_len, outfile);
	}
	
	return 1;
}


/**
 * Prepare input file for decompression of image.
 * @param infile opened input file.
 * @param pcinfo
 * @param with_APP11_marker If should be read APP11 marker. For reading HDR images must be set to true, else
 * HDR image will be read as LDR.
 * @retval 0 Error.
 * @retval 1 Ok, hdr image.
 * @retval 2 Ok, ldr image.
 */
int JPEG_HDR_open_JPEG_file_and_prepare_to_reading(FILE *infile, struct jpeg_decompress_struct *pcinfo, int with_APP11_marker) {
	/* We use our private extension JPEG error handler.
     * Note that this struct must live as long as the main JPEG parameter
     * struct, to avoid dangling-pointer problems.
     */
    
	
	/* In this example we want to open the input file before doing anything else,
     * so that the setjmp() error recovery below can assume the file is open.
     * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
     * requires it in order to read binary files.
     */
    if (infile == NULL) {
		return 0;
    }
	
	
    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    pcinfo->err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
    if (setjmp(jerr.setjmp_buffer)) {
		/* If we get here, the JPEG code has signaled an error.
		 * We need to clean up the JPEG object, close the input file, and return.
		 */
		jpeg_destroy_decompress(pcinfo);
		fclose(infile);
		return 0;
    }
	
    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(pcinfo);
	
	
    /* Step 2: specify data source (eg, a file) */
	
    jpeg_stdio_src(pcinfo, infile);
	
	if (with_APP11_marker) {
		/* Setting up APP11 makrker to be read */
		jpeg_save_markers(pcinfo, JPEG_APP0 + 11, 0xFFFF);
	}
	
    /* Step 3: read file parameters with jpeg_read_header() */
	
    (void) jpeg_read_header(pcinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.doc for more info.
     */
	
	bool hdr = false;
	if (with_APP11_marker) {
		struct jpeg_marker_struct *marker;
		for (marker = pcinfo->marker_list; marker != NULL; marker = marker->next) {
			if (strncmp((char *) marker->data, "HDR_RI ver=", 11) == 0){
				hdr = true;
				pcinfo->out_color_space = JCS_YCbCr;
			}
		}
	}

	
    /* Step 4: set parameters for decompression */

    /* In this example, we don't need to change any of the defaults set by
     * jpeg_read_header(), so we do nothing here.
     */

    /* Step 5: Start decompressor */

    (void) jpeg_start_decompress(pcinfo);
    /* We can ignore the return value since suspension is not possible
     * with the stdio data source.
     */

	return (hdr)?1:2;
}

/**
 * Prepare input file for decompression of HDR image.
 * @param infile opened input file.
 * @param pcinfo
 * @param with_APP11_marker If should be read APP11 marker. For reading HDR images must be set to true, else
 * HDR image will be read as LDR.
 * @retval 0 Error.
 * @retval 1 Ok, hdr image.
 * @retval 2 Ok, ldr image.
 * @see JPEG_HDR_open_JPEG_file_and_prepare_to_reading()
 */
int JPEG_HDR_prepare_reading(FILE *infile, struct jpeg_decompress_struct *pcinfo) {
	return JPEG_HDR_open_JPEG_file_and_prepare_to_reading(infile, pcinfo, 1);
}

double JPEG_HDR_inverse_gamma(double val) {
	double abs_val = abs(val);
	
	if (abs_val <= 0.04045) {
        return val / 12.92;
	}
	
	return pow(((abs_val + 0.055) / 1.055), 2.4) * (val < 0 ? -1 : 1);
}

void JPEG_HDR_gamutDecompanding(double r, double g, double b_, double alpha, double beta, double *red, double *green, double *blue) {
	double a, b, c, y;
	double *a_out, *b_out, *c_out;
	
	
	if (r == g && g ==b_) {
		*red = r;
		*green = g;
		*blue = b_;
		return;
	}
	

	y = (0.2126 * r + 0.7152 * g + 0.0722 * b_);
	
	/* determine minimum */
	if (r <= g && r <= b_) {
		a = r;
		a_out = red;

		b = g;
		b_out = green;

		c = b_;
		c_out = blue;
	} else if (g <= r && g <= b_) {
		a = g;
		a_out = green;

		b = r;
		b_out = red;

		c = b_;
		c_out = blue;
	} else {
		a = b_;
		a_out = blue;

		b = r;
		b_out = red;
		
		c = g;
		c_out = green;
	}

	*a_out = (y - y * pow((y - a) / (alpha * y), 1 / beta));
	*b_out = (y - ((y - b) / alpha) * pow((1 - *a_out / y), 1 - beta));
	*c_out = (y - ((y - c) / alpha) * pow((1 - *a_out / y), 1 - beta));
}


