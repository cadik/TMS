/*
Editor : Sung-Jun Yoon
E-mail : sungjunyoon@kaist.ac.kr
Title : header.h
Version : 1701241333
*/

#ifdef __cplusplus
extern "C" {
#endif

#define EPSILON 2.204e-16
#define NumOfBit 8
#define MAX_INT pow(2, NumOfBit) - 1 // 2^8-1

typedef struct {
	size_t row;
	size_t col;
	byte_buffer buf_B;
	byte_buffer buf_G;
	byte_buffer buf_R;
} picture_BMP;

typedef struct {
	size_t row;
	size_t col;
} option_size;

typedef struct {
	size_t row;
	size_t col;
} img_sz;

typedef struct {
	double row;
	double col;
} scaling_ratio;

typedef struct {
	double row;
	double col;
} interp_kernel_taps;

typedef enum {
	nearest, 
	bilinear, 
	bicubic,
	lanczos2, 
	lanczos3
} interp_kernel;

typedef enum {
	circular,
	replicate,
	symmetric
} how_to_pad;

double interp_kernel_func(double x, interp_kernel interp_kernel_name)
{
	double res = 0.0;
	double absx = fabs(x);
	double absx2 = absx * absx;
	double absx3 = absx * absx2;

	switch (interp_kernel_name) {
	case nearest:
		res = ((-0.5 <= x) && (x < 0.5));
		break;
	case bilinear:
		res = (x + 1)*((-1 <= x) && (x < 0)) + (-x + 1)*((0 <= x) && (x <= 1));
		break;
	case bicubic:
		/* See Keys, "Cubic Convolution Interpolation for Digital Image
		   Processing, " IEEE Transactions on Acoustics, Speech, and Signal
		   Processing, Vol.ASSP - 29, No. 6, December 1981, p. 1155. */
		res = (1.5*absx3 - 2.5*absx2 + 1.0) * (absx <= 1) + (-0.5*absx3 + 2.5*absx2 - 4.0 * absx + 2.0) * ((1 < absx) && (absx <= 2));
		break;
	case lanczos2:
		/* See Graphics Gems, Andrew S.Glasser(ed), Morgan Kaufman, 1990,
		   pp. 156 - 157. */
		res = (sin(PI*x)*sin(PI*x / 2.0) + EPSILON) / (((pow(PI, 2.0) * pow(x, 2.0)) / 2.0) + EPSILON);
		res *= (absx < 2);
		break;
	case lanczos3:
		/* See Graphics Gems, Andrew S.Glasser(ed), Morgan Kaufman, 1990,
		   pp. 157 - 158. */
		res = (sin(PI*x)*sin(PI*x / 3.0) + EPSILON) / (((pow(PI, 2.0) * pow(x, 2.0)) / 3.0) + EPSILON);
		res *= (absx < 3);
		break;
	default:
		error("The function interp_kernel_func may be error.");
	}

	return res;
}

double iminterp1(double *ptr_window, double idx_fraction, double scaling_ratio, interp_kernel interp_kernel_name, double interp_kernel_taps, bool anti_aliasing_filter)
{
	size_t i = 0;
	size_t kernel_taps = (size_t)ceil(interp_kernel_taps);
	double *interp_kernel_window = (double *)calloc(kernel_taps, sizeof(double));
	double cum_sum_interp_kernel_window = 0.0;
	double interesting_pst = 0.0;
	double initial_pst = 0.0;
	double res = 0.0;

	initial_pst = (idx_fraction - 1.0);
	for (i = 0; i < kernel_taps; i++) {
		interesting_pst = initial_pst - (double)i;
		
		if (anti_aliasing_filter) {
			interesting_pst *= scaling_ratio;
		}
		interp_kernel_window[i] = interp_kernel_func(interesting_pst, interp_kernel_name);
		if (anti_aliasing_filter) {
			interp_kernel_window[i] *= scaling_ratio;
		}
		cum_sum_interp_kernel_window += interp_kernel_window[i];
	}
	
	for (i = 0; i < kernel_taps; i++) {
		res += ptr_window[i] * (interp_kernel_window[i] / cum_sum_interp_kernel_window);
	}

	free(interp_kernel_window);

	return res;
}

double_buffer DB_padarray(double_buffer *Ptr_ori_buf, option_size *Ptr_pad_pair_number, how_to_pad option)
{
	size_t i, j, k = 0;
	double_buffer res = DB_mem_alloc_2_with_size(Ptr_ori_buf->row + 2 * Ptr_pad_pair_number->row, Ptr_ori_buf->col + 2 * Ptr_pad_pair_number->col);

	switch (option) {
	case replicate: // replicate, mirroring
		for (i = 0; i < Ptr_ori_buf->row; i++) {
			memcpy(res.buf[i + Ptr_pad_pair_number->row] + Ptr_pad_pair_number->col, Ptr_ori_buf->buf[i], Ptr_ori_buf->col*sizeof(double));
		}

		for (i = 0; i < Ptr_pad_pair_number->row; i++) {
			memcpy(res.buf[i], res.buf[Ptr_pad_pair_number->row], res.col*sizeof(double));
			memcpy(res.buf[i + Ptr_ori_buf->row + Ptr_pad_pair_number->row], res.buf[Ptr_ori_buf->row + Ptr_pad_pair_number->row - 1], res.col*sizeof(double));
		}

		for (j = 0; j < Ptr_pad_pair_number->col; j++) {
			for (i = 0; i < res.row; i++) {
				res.buf[i][j] = res.buf[i][Ptr_pad_pair_number->col];
				res.buf[i][j + Ptr_ori_buf->col + Ptr_pad_pair_number->col] = res.buf[i][Ptr_ori_buf->col + Ptr_pad_pair_number->col - 1];
			}
		}
		break;

	case symmetric:
		for (i = 0; i < Ptr_ori_buf->row; i++) {
			memcpy(res.buf[i + Ptr_pad_pair_number->row] + Ptr_pad_pair_number->col, Ptr_ori_buf->buf[i], Ptr_ori_buf->col*sizeof(double));
		}

		for (i = Ptr_pad_pair_number->row, j = 0; i >= 1; i--, j += 2) {
			memcpy(res.buf[i - 1], res.buf[i + j], res.col*sizeof(double));
			memcpy(res.buf[res.row - i], res.buf[i + Ptr_ori_buf->row - 1], res.col*sizeof(double));
		}

		for (j = Ptr_pad_pair_number->col, k = 0; j >= 1; j--, k += 2) {
			for (i = 0; i < res.row; i++) {
				res.buf[i][j - 1] = res.buf[i][j + k];
				res.buf[i][res.col - j] = res.buf[i][j + Ptr_ori_buf->col - 1];
			}
		}
		break;

	case circular: // @To do
		error("The option of padarray function has not been implemented yet");
		break;

	default:
		printf("The option of padarray function is not correct.\n");
		assert(option != 's');
	}

	return res;
}

void free_mem_alloc_2(void **dest_buf, size_t row) {
	size_t i = 0;

	for (i = 0; i < row; i++) {
		free(dest_buf[i]);
	}
	free(dest_buf);
}

double_buffer imresize(double_buffer *ptr_ori, img_sz *ptr_out_img_sz, interp_kernel interp_kernel_name, bool anti_aliasing_filter)
{
	interp_kernel_taps kernel_taps;
	option_size pad_sz;
	scaling_ratio scaling_ratio_val;
	size_t kernel_tap = 0;
	size_t i, j, n = 0;
	size_t ori_interp_kernel_taps = 0;
	int idx_window = 0;
	int left = 0;
	int indices = 0;
	double idx_src = 0.0;
	double *window;
	double_buffer temp_pad;	
	double_buffer res_temp = DB_mem_alloc_2_with_size(ptr_ori->row, ptr_out_img_sz->col);
	double_buffer res = DB_mem_alloc_2_with_size(ptr_out_img_sz->row, ptr_out_img_sz->col);
	
	scaling_ratio_val.row = (double)ptr_out_img_sz->row / (double)ptr_ori->row;
	scaling_ratio_val.col = (double)ptr_out_img_sz->col / (double)ptr_ori->col;

	switch (interp_kernel_name) {
	case nearest:
		ori_interp_kernel_taps = 1;
		break;
	case bilinear:
		ori_interp_kernel_taps = 2;
		break;
	case bicubic:
		ori_interp_kernel_taps = 4;
		break;
	case lanczos2:
		ori_interp_kernel_taps = 4;
		break;
	case lanczos3:
		ori_interp_kernel_taps = 6;
		break;
	default:
		break;
	}

	kernel_taps.row = (double)ori_interp_kernel_taps;
	kernel_taps.col = (double)ori_interp_kernel_taps;

	if (scaling_ratio_val.col >= 1.0 || interp_kernel_name == nearest) {
		anti_aliasing_filter = false;
	}
	if (anti_aliasing_filter && (scaling_ratio_val.col < 1.0)) {
		kernel_taps.col /= scaling_ratio_val.col;
	}
	
	pad_sz.row = 0;
	pad_sz.col = (size_t) ceil((double)kernel_taps.col / 2.0);
	temp_pad = DB_padarray(ptr_ori, &pad_sz, symmetric);

	kernel_tap = (size_t)ceil(kernel_taps.col);
	window = (double *)calloc(kernel_tap, sizeof(double));

	for (i = 0; i<ptr_ori->row; i++) {
		for (j = 0; j<ptr_out_img_sz->col; j++) {
			idx_src = ((double)(j + 1) / scaling_ratio_val.col) + 0.5 * (1.0 - 1.0 / scaling_ratio_val.col); // u
			left = (int)floor(idx_src - ((double)kernel_taps.col / 2.0));
			indices = left - 1;
			
			for (n = 0; n < kernel_tap; n++) {
				idx_window = (left + (int)n) + (int)pad_sz.col;
				window[n] = temp_pad.buf[i][idx_window];
			}
			res_temp.buf[i][j] = iminterp1(window, (idx_src - indices - 1), scaling_ratio_val.col, interp_kernel_name, kernel_taps.col, anti_aliasing_filter); //idx_src - floor(idx_src)
		}
	}

	if (scaling_ratio_val.row >= 1.0 || interp_kernel_name == nearest) {
		anti_aliasing_filter = false;
	}
	if (anti_aliasing_filter && (scaling_ratio_val.row < 1.0)) {
		kernel_taps.row /= scaling_ratio_val.row;
	}

	pad_sz.row = (size_t)ceil((double)kernel_taps.row / 2.0);
	pad_sz.col = 0;
	temp_pad = DB_padarray(&res_temp, &pad_sz, symmetric);
	
	kernel_tap = (size_t)ceil(kernel_taps.row);
	window = (double *)calloc(kernel_tap, sizeof(double));	

	for (i = 0; i < ptr_out_img_sz->row; i++) {
		idx_src = ((double)(i + 1) / scaling_ratio_val.row) + 0.5 * (1.0 - 1.0 / scaling_ratio_val.row); // u
		left = (int)floor(idx_src - ((double)kernel_taps.row / 2.0));
		indices = left - 1;

		for (j = 0; j<ptr_out_img_sz->col; j++) {
			for (n = 0; n < kernel_tap; n++) {
				idx_window = (left + (int)n) + (int)pad_sz.row;
				window[n] = temp_pad.buf[idx_window][j];
			}
			res.buf[i][j] = iminterp1(window, (idx_src - indices - 1), scaling_ratio_val.row, interp_kernel_name, kernel_taps.row, anti_aliasing_filter);
		}
	}

	free(window);
	free_mem_alloc_2((void**)temp_pad.buf, temp_pad.row);
	free_mem_alloc_2((void**)res_temp.buf, res_temp.row);

	return res;	
}


#ifdef __cplusplus
}
#endif
