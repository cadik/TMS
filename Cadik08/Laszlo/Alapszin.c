//COLOR TO GRAY

#include <io.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <dos.h>
#include <string.h>

#include "_common.h"

void main(void){
	    
		double X1,Y1,Z1, X2,Y2,Z2, A_hue,T,V, grad_LUMINANCE, dA,dT,dV;

		
		READ_SPECTRUM_DATAE();
		READ_COLOROID_DATAE();
		READ_COLOR2GRAY_DATAE();
		_7_basic_fi_computation();
//		regression_hue_saturation_weigths(); it was one times necessary to calibration

		
		//the next function defines the grad value, starting from (X1,Y1,Z1) and (X2,Y2,Z2) neighbor pixels
		//grad_LUMINANCE is measured in Coloroid luminance V
		//the integral of consitent image will be obtained in Coloroid V luminance
		//CIE Y = 0.01 * (V * V) is the formula to Y !!!!!
		//pipeline: from Y the rgb and after the gamma-RGB

		//51 45 46   blue
		//12 50 82   yellow  A, T, V
		X1 = 25.87;  Y1 = 21.16;  Z1 = 79.7;
		X2 = 61.81;  Y2 = 67.24;  Z2 = 23.2;

		grad_LUMINANCE = LUMINANCE_GRAD(X1, Y1, Z1,
							            X2, Y2, Z2,
		    				           &dA, &dT, &dV);

		printf("\n\t   dA = %g  dT = %g   dV = %g  grad = %g", dA, dT, dV, grad_LUMINANCE);
		getch();
		//A , T , V  are the separately not important hue, saturation and luminance parts

		//CURRENTLY:    #define WEIGHT_LIGHTNESS_CHROMINANCE 3.0
		//this weight regultes the luminance chrominance ratio
		//for higher chrominance effect it has to be increased the above default value


	exit(0);


}
