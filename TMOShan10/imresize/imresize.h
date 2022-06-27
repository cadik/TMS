/*
Editor : Sung-Jun Yoon
E-mail : sungjunyoon@kaist.ac.kr
Title : header.h
Version : 1701241333
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <memory.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>

//#include "matrix_computations.h"
//#include "image_processing.h"
// #include "fruc.h"

#define PI 3.1415926535
//#define EPS 2.204e-16

#define COMMAND_SIZE 256
#define COMMAND_CONSOLE_LINES 80
#define COMMAND_CONSOLE_COLS 220


struct TIME {
	time_t startTime, endTime;
	double gapTime;
	int hours, minutes, seconds;
};

void error(char *error_statement);

#ifdef __cplusplus
}
#endif
