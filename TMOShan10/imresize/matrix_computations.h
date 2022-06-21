/*
Editor : Sung-Jun Yoon
E-mail : sungjunyoon@kaist.ac.kr
Title : header.h
Version : 1701241333
*/

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char BYTE;

typedef struct byte_matrix_size {
	BYTE **buf;
	size_t row;
	size_t col;
} byte_buffer;

typedef struct double_matrix_size {
	double **buf;
	size_t row;
	size_t col;
} double_buffer;

typedef struct int_matrix_size {
	int **buf;
	size_t row;
	size_t col;
} int_buffer;

typedef struct double_3d_matrix_size {
	double ***buf;
	size_t row;
	size_t col;
	size_t page;
} double_buffer3;

typedef struct int_3d_matrix_size {
	int ***buf;
	size_t row;
	size_t col;
	size_t page;
} int_buffer3;

double_buffer DB_mem_alloc_2_with_size(size_t row, size_t col)
{
	// Initialize array with 0 using calloc function because memset function doesn't work for double array
	size_t i = 0;
	double_buffer res;
	res.row = row;
	res.col = col;
	
	res.buf = (double **)calloc(res.row, sizeof(double*));
	for (i = 0; i<res.row; i++) {
		res.buf[i] = (double *)calloc(res.col, sizeof(double));
	}

	return res;
}

#ifdef __cplusplus
}
#endif
