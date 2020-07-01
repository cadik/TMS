#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void color2grey(  double **fpI,          // Input
                  double **fpO,          // Output
                  int d_c, int d_w,int d_h);

struct complex_vec{
    double real_x;
    double real_y;
    double imag_x;
    double imag_y;
};

double vc(double val);
double vl(double val);
double uc(double val);
double ul(double val);

double calculate_D(double** fpI,int w,int h ,int index, int d_w, int d_h);
void amplitude_mod(double *imaginery_part, int w, int h, int chanel);
double frekv_mod(double ** Lab_origin,int pos, int chanel);