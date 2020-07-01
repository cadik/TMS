#include "color2grey.h"
#include "rgb2lab.h"
#include "optimalizer.h"
#include "transform.h"
#include <cmath>


void amplitude_mod(double *imaginery_part, int w, int h, int chanel){
    hilbert(imaginery_part,2); //aplication of hilbert
}


double unsigned_distance(int pix1, int pix2, double** frekv_matrx){
    double L =  frekv_matrx[0][pix1] - frekv_matrx[0][pix2];
    double a =  frekv_matrx[1][pix1] - frekv_matrx[1][pix2];
    double b =  frekv_matrx[2][pix1] - frekv_matrx[2][pix2];
    return sqrt(pow(L,2) + pow(a,2) + pow(b,2));
}

int sgn(double prom)
{
    return (0 < prom) - (prom < 0);
}


double vc(double **fpI, int pix){
    double X;
    double Y;
    double Z;
    RGBtoXYZ(fpI[0][pix],fpI[1][pix],fpI[2][pix],X,Y,Z);
    return (9*Y)/(X+15*Y+3*Z);
}

double vl(double **fpI, int pix){
    double X = 0.9504;
    double Y = 1.0000;
    double Z = 1.0888;
    return (9*Y)/(X+15*Y+3*Z);
}

double uc(double **fpI, int pix){
    double X;
    double Y;
    double Z;
    RGBtoXYZ(fpI[0][pix],fpI[1][pix],fpI[2][pix],X,Y,Z);
    return (4*X)/(X+15*Y+3*Z);
}

double ul(double **fpI, int pix){
    double X = 0.9504;
    double Y = 1.0000;
    double Z = 1.0888;
    return (4*X)/(X+15*Y+3*Z);
}

double signed_distance(int pix1, int pix2, double** Lab_modulation, double** frekv_matrx, double **fpI){

    double dL = Lab_modulation[0][pix1] - Lab_modulation[0][pix2];
    double da = Lab_modulation[1][pix1] - Lab_modulation[1][pix2];
    double db = Lab_modulation[2][pix1] - Lab_modulation[2][pix2];

    double vc1 = vc( fpI, pix1);
    double vc2 = vc( fpI, pix2);
    double vl1 = vl( fpI, pix1);
    double vl2 = vl( fpI, pix2);
    double uc1 = uc( fpI, pix1);
    double uc2 = uc( fpI, pix2);
    double ul1 = ul( fpI, pix1);
    double ul2 = ul( fpI, pix2);

    double angle1 = atan((vc1 - vl1) / (uc1 - ul1));
    double angle2 = atan((vc2 - vl2) / (uc2 - ul2));

    double q1 = - 0.01585 - 0.03017*cos(angle1) - 0.04556*cos(2*angle1) - 0.02667*cos(3*angle1) - 0.00295*cos(4*angle1) + 0.14592*sin(angle1) + 0.05084*sin(2*angle1) - 0.01900*sin(3*angle1) - 0.00764*sin(4*angle1);
    double q2 = - 0.01585 - 0.03017*cos(angle2) - 0.04556*cos(2*angle2) - 0.02667*cos(3*angle2) - 0.00295*cos(4*angle2) + 0.14592*sin(angle2) + 0.05084*sin(2*angle2) - 0.01900*sin(3*angle2) - 0.00764*sin(4*angle2);

    double k = 0.2717*((6.469 + 6.362 * pow(20,0.4495)) / (6.469 + pow(20,0.4495)));

    double s1 = 13 * sqrt(pow(uc1 - ul1,2) + pow(vc1 - vl1,2));
    double s2 = 13 * sqrt(pow(uc2 - ul2,2) + pow(vc2 - vl2,2));

    double dLvac1 = Lab_modulation[0][pix1] + ((-0.1340 * q1) + (0.0872 * k)) * s1 * Lab_modulation[0][pix1];
    double dLvac2 = Lab_modulation[0][pix2] + ((-0.1340 * q2) + (0.0872 * k)) * s2 * Lab_modulation[0][pix1];

    double result = dLvac1 - dLvac2;

    double signum_val = 0;

    if(result != 0){
        signum_val = sgn(result);
    } else if(result == 0 && dL != 0){
        signum_val = sgn(dL);
    } else{
        signum_val = sgn(pow(da,3) + pow(db,3));
    }
    return signum_val * unsigned_distance(pix1,pix2,frekv_matrx);
}

double calculate_d(int w,int h ,int index, int d_w, int d_h, double** Lab_modulation ,double **frekv_matrx, double **fpI){
    int rotat = 0;
    double result = 0;
    switch (index){
        case 0:
            rotat = w+(h-1);
            break;
        case 1:
            rotat = h+w+d_h;
            break;
        case 2:
            rotat = w+h+1;
            break;
        case 3:
            rotat = h+(w-d_h) ;
            break;
    }
    //fprintf(stderr,"%d = w, %d  = h, %d = index \n",w,h,index);
    const int max_h = (d_h-1);
    const int max_w = ((d_w*d_h)-d_h);
    if (w == 0){
            if (h == 0) {
                if (index == 1 || index == 2) {
                    result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
                }
            } else if(h == max_h) {
                if (index == 0 || index == 1) {
                    result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
                }
            } else {
                if (index >= 0 && index < 3) {
                    result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
                }
            }
    } else if( max_w == w){
        if (h == 0) {
            if (index > 1) {
                result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
            }
        } else if(h == max_h) {
            if ((index == 0 || index == 3)) {
                result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
            }
        } else {
            if(index == 0 || index == 3 || index == 2){
                result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
            }
        }
    } else {
        if (h == 0) {
            if ((index > 0)) {
                result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
            }
        } else if(h == max_h) {
            if (index == 0 || index == 1 || index == 3) {
                result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
            }
        } else {
            result = signed_distance(w+h, rotat ,Lab_modulation, frekv_matrx, fpI);
        }
    }
    return result;
}

double calculate_D(double** fpI,int w,int h ,int index, int d_w, int d_h){
    int c = 0;
    int rotat = 0;
    double result = 0;
    switch (index){
        case 0:
            c = 0;
            rotat = w+(h-1);
            break;
        case 1:
            c = 1;
            rotat = w+(h-1);
            break;
        case 2:
            c = 2;
            rotat = w+(h-1);
            break;
        case 3:
            c = 0;
            rotat = h+w+d_h;
            break;
        case 4:
            c = 1;
            rotat = h+w+d_h;
            break;
        case 5:
            c = 2;
            rotat = h+w+d_h;
            break;
        case 6:
            c = 0;
            rotat = w+h+1;
            break;
        case 7:
            c = 1;
            rotat = w+h+1;
            break;
        case 8:
            c = 2;
            rotat = w+h+1;
            break;
        case 9:
            c = 0;
            rotat = h+(w-d_h) ;
            break;
        case 10:
            c = 1;
            rotat = h+(w-d_h) ;
            break;
        case 11:
            c = 2;
            rotat = h+(w-d_h) ;
            break;
    }

    const int max_h = (d_h-1);
    const int max_w = ((d_w*d_h)-d_h);
    if (w == 0){
        if (h == 0) {
            if (index > 2 && index < 9) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        } else if(h == max_h) {
            if (index >= 0 && index < 6) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        } else {
            if (index >= 0 && index < 9) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        }
    } else if( max_w == w){
        if (h == 0) {
            if (index > 5) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        } else if(h == max_h) {
            if ((index >= 0 && index < 3) || (index > 8)) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        } else {
            if((index >= 0 && index < 3) || (index > 5)){
                result = fpI[c][w+h] - fpI[c][rotat];
            }
        }
    } else {
        if (h == 0) {
            if ((index > 2)) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        } else if(h == max_h) {
            if ((index >= 0 && index < 6) || (index > 8)) {
                result = fpI[c][w + h] - fpI[c][rotat];
            }
        } else {
            result = fpI[c][w+h] - fpI[c][rotat];
        }
    }
    return result;
}


void color2grey(  double **fpI, double **fpO,  int d_c, int d_w,int d_h) {


    double **Lab_origin = new double*[d_c];
    for (int ii=0; ii < d_c; ii++) Lab_origin[ii] = new double[d_w*d_h];



    /*
     * STEP 1 calculation of Lab space
     */
    for(int w = 0; w < d_w*d_h; w += d_h){
        for(int h = 0; h < d_h; h++){
            RGBtoLab(fpI[0][w + h],fpI[1][w + h],fpI[2][w + h],Lab_origin[0][w + h],Lab_origin[1][w + h],Lab_origin[2][w + h]);
            //if(fpI[0][w + h]!=0 && fpI[1][w + h]!=0 && fpI[2][w + h]!=0 && fpO[0][w + h]!=0 && fpO[1][w + h]!=0 && fpO[2][w + h]);
        }
    }
    double **real_part = new double*[d_c];
    double **imaginery_part = new double*[d_c];
    complex_vec **complex_vector = new complex_vec*[d_c];
    for (int ii=0; ii < d_c; ii++){ // calculate hilbert transofrmation and reprezentation in complex numbers for computing AM-FM modulation
        real_part[ii] = new double[d_w*d_h];
        imaginery_part[ii] = new double[d_w*d_h];
        complex_vector[ii] = new complex_vec[d_w*d_h];

        for(int w = 0; w < d_w*d_h; w += d_h){
            for(int h = 0; h < d_h; h++){
                real_part[ii][w+h] = Lab_origin[ii][w+h];
                imaginery_part[ii][w+h] = Lab_origin[ii][w+h];
            }
        }

        amplitude_mod(imaginery_part[ii],d_w,d_h,ii); // get amplitude information

        double i, powi,bot_real,bot_imag;
        for(int w = 0; w < (d_w*d_h)  ; w += d_h){ // calculating partial derivatives from complex pictures
            for(int h = 0; h < d_h ; h++){

                complex_vec pom_vec;
                if (bot_real != 0 && bot_imag != 0 && h != 0 && h != (d_h -1) && w != 0 && w != ((d_w*d_h) - 1)) {
                    pom_vec.real_x = (real_part[ii][(w+h)-d_h] - real_part[ii][(w+h)+d_h])/d_w;
                    pom_vec.imag_x = (imaginery_part[ii][(w+h)-d_h] - imaginery_part[ii][(w+h)+d_h])/d_w;
                    pom_vec.real_y = (real_part[ii][(w+h)-1] - real_part[ii][(w+h)+1])/d_h;
                    pom_vec.imag_y = (imaginery_part[ii][(w+h)-1] - imaginery_part[ii][(w+h)+1])/d_h;
                    //printf("%f, %f",pom_vec.real_y,pom_vec.real_x);
                    bot_real = real_part[ii][w+h];
                    bot_imag = -imaginery_part[ii][w+h];

                    pom_vec.real_x = ((pom_vec.real_x * bot_real) + (-(-imaginery_part[ii][w + h] * pom_vec.imag_x))) /
                                     (bot_real * bot_real + bot_imag * bot_imag);
                    pom_vec.imag_x = ((pom_vec.imag_x * bot_real) + (-imaginery_part[ii][w + h] * pom_vec.real_x)) /
                                     (bot_real * bot_real + bot_imag * bot_imag);

                    pom_vec.real_y = ((pom_vec.real_y * bot_real) + (-(-imaginery_part[ii][w + h] * pom_vec.imag_y))) /
                                     (bot_real * bot_real + bot_imag * bot_imag);
                    pom_vec.imag_y = ((pom_vec.imag_y * bot_real) + (-imaginery_part[ii][w + h] * pom_vec.real_y)) /
                                     (bot_real * bot_real + bot_imag * bot_imag);
                } else {
                    pom_vec.real_x = 1;
                    pom_vec.imag_y = 1;
                    pom_vec.real_y = 1;
                    pom_vec.imag_x = 1;
                }

               // printf("%f, %f",pom_vec.real_y,pom_vec.real_x);
                complex_vector[ii][w+h] = pom_vec;
            }
        }
    }

    /*
     * STEP 2 prepering numbers used for Distances of pixels and color and modulation
     */

    double **Lab_modulation = new double*[d_c];
    for (int ii=0; ii < d_c; ii++) Lab_modulation[ii] = new double[d_w*d_h];
    double total_L = 0; double total_a = 0; double total_b = 0;
    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            total_L += Lab_origin[0][w + h];
            total_a += Lab_origin[1][w + h];
            total_b += Lab_origin[2][w + h];
        }
    }

    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            Lab_modulation[0][w + h] = sqrtf( pow(imaginery_part[0][w+h],2) + pow(real_part[0][w+h],2) )*(cos(sqrtf(complex_vector[0][w+h].real_x*complex_vector[0][w+h].real_x + complex_vector[0][w+h].real_y*complex_vector[0][w+h].real_y))) + total_L/(d_w*d_h);
            Lab_modulation[1][w + h] = sqrtf( pow(imaginery_part[0][w+h],2) + pow(real_part[0][w+h],2) )*(cos(sqrtf(complex_vector[1][w+h].real_x*complex_vector[1][w+h].real_x + complex_vector[1][w+h].real_y*complex_vector[1][w+h].real_y))) + total_a/(d_w*d_h);
            Lab_modulation[2][w + h] = sqrtf( pow(imaginery_part[0][w+h],2) + pow(real_part[0][w+h],2) )*(cos(sqrtf(complex_vector[2][w+h].real_x*complex_vector[2][w+h].real_x + complex_vector[2][w+h].real_y*complex_vector[2][w+h].real_y))) + total_b/(d_w*d_h);
        }
    }

    /*
     * STEP 3 more preperations
     */
    double *amplitde_vec = new double[d_w*d_h];
    double **feature_vec = new double*[d_w*d_h];
    for (int ii=0; ii < d_c; ii++) feature_vec[ii] = new double[d_w*d_h];


    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            amplitde_vec[w + h] = sqrtf(pow(sqrtf( pow(imaginery_part[0][w+h],2) + pow(real_part[0][w+h],2) ),2) + pow(sqrtf(pow(imaginery_part[1][w+h],2) + pow(real_part[1][w+h],2)),2) + pow(sqrtf(pow(imaginery_part[2][w+h],2) + pow(real_part[2][w+h],2)),2));
            if(amplitde_vec[w + h] == 0)
                amplitde_vec[w + h] = 0.0000001;

        }
    }

    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            feature_vec[0][w + h] = Lab_origin[0][w+h] / 100;
            feature_vec[1][w + h] = Lab_origin[1][w+h] / amplitde_vec[w + h];
            feature_vec[2][w + h] = Lab_origin[2][w+h] / amplitde_vec[w + h];
            //printf("%f, \n",feature_vec[1][w + h]);
        }
    }

    /*
     * MAIN CYCLE STEP get pixel and color distance
     */

    double *D_vec = new double[d_w*d_h*d_c*4];
    double **alpha_vec = new double*[d_c];
    for (int ii=0; ii < d_c; ii++) alpha_vec[ii] = new double[d_w*d_h];
    double *dist_vec = new double[d_w*d_h*4];

    int index = 0;
    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            for(int c = 0; c < d_c*4; c++){
                D_vec[index] = calculate_D(fpI,w,h,index % 12,d_w,d_h);
                index++;
            }
        }
    }
    index = 0;
    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            for(int c = 0; c < 4; c++) {
                dist_vec[index] = calculate_d( w, h, c, d_w, d_h, Lab_origin, feature_vec, fpI);
                //printf("%f, ",dist_vec[index]);
                index++;
            }
        }
    }


    /*
     * RETURNING GREYSCALE IMAGE5 using optimalizer to find minimum for every pixel taking in consideration other pixels
     */
    int counter = 0;
    int Dindex = 0;
    int dindex = 0;
    for(int w = 0; w < (d_w*d_h) ; w += d_h){
        for(int h = 0; h < d_h ; h++){
            /*if(counter == 23957){
                printf("%f R, %f, G, %f, B, %f D",D_vec[Dindex+3],D_vec[Dindex+4],D_vec[Dindex+5],dist_vec[dindex+1]);
            }*/
                //if(isnan(dist_vec[dindex]) || isnan(dist_vec[dindex+1]) || isnan(dist_vec[dindex+2]) || isnan(dist_vec[dindex+3])) {
                if(fpI[0][w + h] == 0 || fpI[1][w + h] == 0 || fpI[2][w + h] == 0 ){
                    fpO[0][w + h] = 0;
                    fpO[1][w + h] = 0;
                    fpO[2][w + h] = 0;
                } else {
                    //printf("\n%f, %f, %f, %f, \n", dist_vec[dindex], dist_vec[dindex+1], dist_vec[dindex+2], dist_vec[dindex+3]);

                    optimalize(D_vec[Dindex], D_vec[Dindex + 1], D_vec[Dindex + 2], D_vec[Dindex + 3],
                               D_vec[Dindex + 4], D_vec[Dindex + 5], D_vec[Dindex + 6], D_vec[Dindex + 7],
                               D_vec[Dindex + 8], D_vec[Dindex + 9], D_vec[Dindex + 10], D_vec[Dindex + 11],
                               fpO[0][w + h], fpO[1][w + h], fpO[2][w + h], dist_vec[dindex], dist_vec[dindex + 1],
                               dist_vec[dindex + 2], dist_vec[dindex + 3]);
                    fpO[0][w + h] = fpI[0][w + h]*fpO[0][w + h] + fpI[1][w + h]*fpO[1][w + h] + fpI[2][w + h]*fpO[2][w + h];
                    fpO[1][w + h] = fpO[0][w + h];
                    fpO[2][w + h] = fpO[0][w + h];
                    fprintf(stderr,"pixel x: %d, y: %d -- ",counter,h);

                }

                Dindex +=12;
                dindex +=4;

        }
        counter++;
    }
}

