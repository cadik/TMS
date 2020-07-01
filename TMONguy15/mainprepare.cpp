#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "color2grey.h"
#include "io_png.h"


double* mainprepare(double* d_v, int d_w, int d_h) {

    int d_c = 3;
    int d_wh = d_w * d_h;
    int d_whc = d_c * d_w * d_h;

    size_t nx,ny,nc;
    double *nd_v = new double[d_whc];

    const int nd_w = (int) nx;
    const int nd_h = (int) ny;
    const int nd_c = (int) nc-1;

    int index = 0;
    for(int w = 0; w < (d_w*d_h) ; w += d_w){
        for(int h = 0; h < d_w ; h++){
            for(int c = 0; c < 3; c++) {
                nd_v[(d_wh*c)+(w+h)] = d_v[w+h+index+c]; // converting array to suitable representation
            }
            index+=2;
        }
    }


    double **fpI = new double*[d_c]; //input
    double **fpO = new double*[d_c]; //output
    double *greyscale = new double[d_whc];

    for (int ii=0; ii < d_c; ii++) {
        fpI[ii] = &nd_v[ii * d_wh];
        fpO[ii] = &greyscale[ii * d_wh];
    }



    color2grey(fpI, fpO, d_c, d_w, d_h); // main funtion transfers input to greyscale and store it into output

    index = 0;
    for(int w = 0; w < (d_w*d_h) ; w += d_w){ //returning to tmolib array convetion
        for(int h = 0; h < d_w ; h++){
            for(int c = 0; c < 3; c++) {
                nd_v[w+h+index+c] = greyscale[(d_wh*c)+(w+h)];
            }
            index+=2;
        }
    }

    return nd_v;
}