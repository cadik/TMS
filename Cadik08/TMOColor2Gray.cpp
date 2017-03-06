//(c)Martin Cadik
//03--05/2007 - Prague

//POZOR - neni dodelane:
//	double ZETA [psi_xmax][psi_ymax];
//	double PSI [psi_xmax][psi_ymax];
// -jsou "natvrdo" maximalni dimenze - pro mcdesk-test.hdr

//poisson solver based on CATAM,
//http://www.maths.cam.ac.uk/undergrad/tripos/catam/


/* --------------------------------------------------------------------------- *
 * TMOColor2Gray.cpp: implementation of the TMOColor2Gray class.                       *
 * --------------------------------------------------------------------------- */

//const bool NEW_FORMULA=true;

#include "./TMOColor2Gray.h"

#undef CCATSLWIN

extern "C" {
	#include "catam.h"	
}

//#include "RGB_Gamma.cpp"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "Laszlo/_common.h"
#include "common/text_loader.h"

double Cdisplay_01_Clinear(double Cdisplay)
{
//both are normalized to [0,1]
//default 709 gamma = 2.2

	double Clinear;

	Clinear =  Cdisplay <= 0.081 ?  Cdisplay / 4.5 : 
	pow((Cdisplay + 0.099) / 1.099, 2.2); 

	return(Clinear);

}

double** alloc_image1(const unsigned w, const unsigned h)
{
	double** img = new double*[w];
	for (unsigned i{}; i < w; ++i)
		img[i] = new double[h];
	return img;
}

void free_image1(double** const img, const unsigned w, const unsigned h)
{
	for (unsigned i{}; i < w; ++i)
		delete[] img[i];
	delete[] img;
}

void multiply_image1(double** const img, const unsigned w, const unsigned h, const double v)
{
	for (unsigned i{}; i < w; ++i)
		for (unsigned j{}; j < h; ++j)
			img[i][j] *= v;
}

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOColor2Gray::TMOColor2Gray() :
	host{},
	gpu{host, CL_DEVICE_TYPE_CPU}, //gpu{host, CL_DEVICE_TYPE_GPU},
	env{{CL_CONTEXT_PLATFORM,
	     reinterpret_cast<cl_context_properties>
	     ((const cl_platform_id&) host), 0}, gpu},
	step{gpu, env},
	exe{env, com::text_loader{"Cadik08/resources/kernels/color2gray.cl"}().c_str(),
	    gpu, "-ICadik08/resources/kernels"},
	lms{gpu.info<cl_ulong>(CL_DEVICE_LOCAL_MEM_SIZE)},
	wgs{gpu.info<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE)}
{
	SetName(L"Color2Gray");
	SetDescription(L"Color2Gray in gradient domain - version 0.1");

	normalize.SetName(L"normalize");
	normalize.SetDescription(L"Normalize output luminance (using 1/Lo_max) to be in <0, 1>");
	normalize.SetDefault(true);
	normalize = true;
	this->Register(normalize);

	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma Correction Factor [-]");
	gamma.SetDefault(2.2);
	gamma = 2.2;
	gamma.SetRange(1.0e-3, 1.0e+2);
	this->Register(gamma);

	s.SetName(L"s");
	s.SetDescription(L"overprojection step");
	s.SetDefault(1.99);
	s = 1.99;
	s.SetRange(0., 2.);
	this->Register(s);
	
	eps.SetName(L"eps");
	eps.SetDescription(L"epsilon for gradient correction");
	eps.SetDefault(.1);
	eps = .1;
	eps.SetRange(0., 100.);
	this->Register(eps);
}

TMOColor2Gray::~TMOColor2Gray()
{
}
/*

	const int psi_xmax=513;//2272;//400;
	const int psi_ymax=513;//2049;//513;
	double ZETA [psi_xmax][psi_ymax];
	double PSI [psi_xmax][psi_ymax];



double myPow(double num, double power){
	if(fabs(num)<EPS_8) return 0.;
	if(num>0) return (pow(num, power));
	else return(-pow(fabs(num), power));
}//myPow


double TMOColor2Gray::formula(double* data, long y1, long x1, long y2, long x2, double* thresholddata){
	const double wa_green=0.3;
	const double wa_red=0.3;
	const double wb_blue=0.2;
	const double wb_yellow=0.4;
	const double threshold=1.;
	const double wa=0.25;
	const double wb=0.33;

	long tmp_ind_1 = y1*xmax+x1;
	long tmp_ind_2 = y2*xmax+x2;

	double dL= data[3*tmp_ind_1]-data[3*tmp_ind_2];
	double da= data[3*tmp_ind_1+1]-data[3*tmp_ind_2+1];
	double db= data[3*tmp_ind_1+2]-data[3*tmp_ind_2+2];

	double T1=(sqrt(data[3*tmp_ind_1+1]*data[3*tmp_ind_1+1]+data[3*tmp_ind_1+2]*data[3*tmp_ind_1+2]))/data[3*tmp_ind_1];
	double T2=(sqrt(data[3*tmp_ind_2+1]*data[3*tmp_ind_2+1]+data[3*tmp_ind_2+2]*data[3*tmp_ind_2+2]))/data[3*tmp_ind_2];

	T1=(sqrt(data[3*tmp_ind_1+1]*data[3*tmp_ind_1+1]+data[3*tmp_ind_1+2]*data[3*tmp_ind_1+2]))/data[3*tmp_ind_1];
	T2=(sqrt(data[3*tmp_ind_2+1]*data[3*tmp_ind_2+1]+data[3*tmp_ind_2+2]*data[3*tmp_ind_2+2]))/data[3*tmp_ind_2];


	//if(	(T1<threshold) || (T2<threshold) ){
	//		thresholddata[3*tmp_ind_1]=data[3*tmp_ind_1]*0.01;
	//		thresholddata[3*tmp_ind_1+1]=0;
	//		thresholddata[3*tmp_ind_1+2]=0;
	//		return (dL);
	//	}
	//else{
		thresholddata[3*tmp_ind_1]=data[3*tmp_ind_1]*0.01;
		thresholddata[3*tmp_ind_1+1]=data[3*tmp_ind_1]*0.01;
		thresholddata[3*tmp_ind_1+2]=data[3*tmp_ind_1]*0.01;

		if(0){//4 weights
			if(da>0) da*=wa_red;
			else da*=wa_green;
			if(db>0) db*=wb_yellow;
			else db*=wb_blue;
			return (myPow( myPow(dL,3.) + myPow(da,3.) + myPow(db,3.), 1./3.));
		}
		else{//2 weights
			da*=wa;
			db*=wb;
			return (myPow( myPow(dL,3.) + myPow(da,3.) + myPow(db,3.), 1./3.));
		}
		
        //return (myPow( myPow(dL,3.) + myPow(wa*da,3.) + myPow(wb*db,3.), 1./3.));
	//}


////	if(	((data[3*tmp_ind_1]-data[3*tmp_ind_2])>0) && (myPow( myPow(dL,3.) + myPow(wa*da,3.) + myPow(wb*db,3.), 1./3.)<=0) ){
//	if(	(data[3*tmp_ind_1]>80) || (data[3*tmp_ind_2]>80) ){
//			thresholddata[3*tmp_ind_1]=data[3*tmp_ind_1]*0.01;
//			thresholddata[3*tmp_ind_1+1]=0;
//			thresholddata[3*tmp_ind_1+2]=0;
//			return (dL);
//		}
//	else{
//		thresholddata[3*tmp_ind_1]=data[3*tmp_ind_1]*0.01;
//		thresholddata[3*tmp_ind_1+1]=data[3*tmp_ind_1]*0.01;
//		thresholddata[3*tmp_ind_1+2]=data[3*tmp_ind_1]*0.01;
//        return (myPow( myPow(dL,3.) + myPow(wa*da,3.) + myPow(wb*db,3.), 1./3.));
//	}


	//return myPow(dL+w*da+w*db,1./3.);
	//return (dL+w*da+w*db);
	
}//formula*/

double TMOColor2Gray::formulaColoroid(const double* const data,
                                      const long y1,
                                      const long x1,
                                      const long y2,
                                      const long x2,
                                      double* const thresholdData)
{
	const double threshold = 1.;
	double X1, Y1, Z1,
	       X2, Y2, Z2,
	       A_hue, T, V,
	       gradLuminance,
	       dA, dT, dV;

	long tmp_ind_1 = y1 * xmax + x1,
	     tmp_ind_2 = y2 * xmax + x2;

	X1 = data[3 * tmp_ind_1];
	X2 = data[3 * tmp_ind_2];
	Y1 = data[3 * tmp_ind_1 + 1];
	Y2 = data[3 * tmp_ind_2 + 1];
	Z1 = data[3 * tmp_ind_1 + 2];
	Z2 = data[3 * tmp_ind_2 + 2];

	thresholdData[3 * tmp_ind_1] = data[3 * tmp_ind_1] * .01;
	thresholdData[3 * tmp_ind_1 + 1] = data[3 * tmp_ind_1] * .01;
	thresholdData[3 * tmp_ind_1 + 2] = data[3 * tmp_ind_1] * .01;

	//the next function defines the grad value, starting from (X1,Y1,Z1) and (X2,Y2,Z2) neighbor pixels
	//grad_LUMINANCE is measured in Coloroid luminance V
	//the integral of consitent image will be obtained in Coloroid V luminance
	//CIE Y = 0.01 * (V * V) is the formula to Y !!!!!
	//pipeline: from Y the rgb and after the gamma-RGB

	//51 45 46   blue
	//12 50 82   yellow  A, T, V
	//X1 = 25.87;  Y1 = 21.16;  Z1 = 79.7;
	//X2 = 61.81;  Y2 = 67.24;  Z2 = 23.2;

	gradLuminance = LUMINANCE_GRAD(X1, Y1, Z1,
	                               X2, Y2, Z2,
	                               &dA, &dT, &dV);

	//printf("\n\t   dA = %g  dT = %g   dV = %g  grad = %g", dA, dT, dV, grad_LUMINANCE);
	//getch();
	//A , T , V  are the separately not important hue, saturation and luminance parts

	//CURRENTLY:    #define WEIGHT_LIGHTNESS_CHROMINANCE 3.0 (color_datae.h)
	//this weight regultes the luminance chrominance ratio
	//for higher chrominance effect it has to be increased the above default value

	return gradLuminance;
}

/* --------------------------------------------------------------------------- *
 * Gradient inconsistency correction -- parallel                                          *
 * --------------------------------------------------------------------------- */
/*
void TMOColor2Gray::inconsistencyCorrection(TMOImage& G_image,
                                            TMOImage& divG_image,
                                            const double eps)
{
	long xmax = G_image.GetWidth(),
	     ymax = G_image.GetHeight();
	long tmp_y, tmp_ind, i, j, iter = 0;
	double* pG_image = G_image.GetData(),
	      * pdivG_image = divG_image.GetData();
	double E = 0, tmp_divG, maxE = 0, avgE = 0;

	do {
		avgE = maxE = 0;
		for (i = 0; i < ymax; ++i) {
			tmp_y = i * xmax;
			for (j = 0; j < xmax; ++j) {
				tmp_ind = j + tmp_y;
				E = pG_image[3 * tmp_ind] - //Gx
				    pG_image[3 * tmp_ind + 1] + //Gy
				    ((j + 1< xmax) ? (pG_image[3 * (tmp_ind + 1) + 1]) : 0.) - //Gy(i+1,j)
				    ((i + 1< ymax) ? (pG_image[3 * (tmp_ind + xmax)]) : 0.); //Gx(i,j+1)

				if (fabs(E) > maxE)
					maxE = fabs(E);
				avgE += fabs(E);

				E *= .25 * s;

				pG_image[3 * tmp_ind] -= E;	//Gx
				pG_image[3 * tmp_ind+1] += E;	//Gy
				if(j + 1 < xmax)
					pG_image[3 * (tmp_ind + 1) + 1] -= E;		//Gy(i+1,j)
				if(i + 1 < ymax)
					pG_image[3 * (tmp_ind + xmax)] += E;	//Gx(i,j+1)

				tmp_divG = (((j - 1) < 0) ? (pG_image[3 * tmp_ind]) :
				           (pG_image[3 * tmp_ind] - pG_image[3 * (tmp_ind - 1)]));
				tmp_divG= tmp_divG + (((i - 1) < 0) ? (pG_image[3 * tmp_ind + 1]) :
				          (pG_image[3 * tmp_ind + 1] - pG_image[3 * (tmp_ind - xmax) + 1]));
				pdivG_image[3 * tmp_ind] = tmp_divG;
				pdivG_image[3 * tmp_ind + 1] = tmp_divG;
				pdivG_image[3 * tmp_ind + 2] = tmp_divG;
			}
		}
		iter++;

		if(!(iter % 100))
			printf("Iteration: %d,\t maxE=%8.2g, \tavgE=%8.2g \n", iter, maxE, avgE);
	} while (maxE > eps);
}
*/

void TMOColor2Gray::correct_grad(TMOImage& g, const double eps)
{
	const unsigned rows = g.GetHeight(),
	               cols = g.GetWidth();
	const cl::buffer grad{env, CL_MEM_READ_WRITE, rows * cols * 3 * sizeof(double)},
	                 err{env, CL_MEM_READ_WRITE, rows * cols * sizeof(double)};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   g.GetData())};

	double e_max;
	do {
		e_max = 0.;

		exe["calc_error"].set_args(grad, err, rows, cols);
		status = step.ndrange_kernel(exe["calc_error"], {0, 0},
		                             {com::math::ceil_multiple(rows,
		                                                       32),
		                              com::math::ceil_multiple(cols,
		                                                       32)},
		                             {32, 32},
		                             {status});

		for (unsigned char i{}; i < 4; ++i) {
			exe["correct_grad"].set_args(grad, err, s.GetDouble(),
			                             i, rows, cols);
			status = step.ndrange_kernel(exe["correct_grad"], {0, 0},
			                             {com::math::ceil_multiple(rows,
			                                                       32),
			                              com::math::ceil_multiple(cols,
			                                                       32)},
			                             {32, 32},
			                             {status});
		}
		status = reduce_max(err, rows * cols, e_max, {status});
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * 3 * sizeof(double),
	                          g.GetData());
}

//______________________________________________________________________________
cl::event TMOColor2Gray::reduce_max(const cl::buffer& in, unsigned n,
                                    double& out,
                                    const std::vector<cl::event> pending)
{
	// number of work-groups needed to reduce the problem
	unsigned m = std::ceil(static_cast<float>(com::math::ceil_multiple(n,
	                       2 * wgs)) / (2.f * wgs));
	// partial maximas
	cl::buffer maximas{env, CL_MEM_READ_WRITE, m * sizeof(double)};

	exe["reduce_max"].set_args(in, cl::local_mem{2 * wgs * sizeof(double)},
	                           maximas, n);
	cl::event status{step.ndrange_kernel(exe["reduce_max"], {0},
	                                     {m * wgs}, {wgs}, pending)};
	while (m > 1) {
		n = m;
		m = std::ceil(static_cast<float>(com::math::ceil_multiple(n,
		                                 2 * wgs)) / (2 * wgs));
		const cl::buffer tmp{env, CL_MEM_READ_WRITE, m * sizeof(double)};
		exe["reduce_max"].set_args(maximas, cl::local_mem{2 * wgs *
		                           sizeof(double)}, tmp, n);
		status = step.ndrange_kernel(exe["reduce_max"], {0},
		                             {m * wgs}, {wgs}, {status});
		maximas = tmp;
	}
	status = step.read_buffer(maximas, 0, sizeof(double), &out, {status});

	return status;
}
/*
void TMOColor2Gray::GFintegration(TMOImage& G_image, TMOImage& Dst_image)
{
	long xmax = Dst_image.GetWidth(),
	     ymax = Dst_image.GetHeight();
	long tmp_y, tmp_ind, i, j;
	double* pG_image = G_image.GetData(),
	      * pDst_image = Dst_image.GetData();

	pDst_image[0] = pDst_image[1] = pDst_image[2] = 0.;
	for (i = 0; i < ymax; ++i) {
		tmp_y = i*xmax;
		if (i > 0)
			pDst_image[3 * tmp_y] = pDst_image[3 * tmp_y + 1] =
			pDst_image[3*tmp_y+2] = (pDst_image[3 * (tmp_y - xmax)] +
			                        pG_image[3 * (tmp_y - xmax) + 1]);
			//neboli: OUTPUT_BW[0][y] = OUTPUT_BW[0][y-1] + Grad_Y[0][y-1];

		for (j = 1; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			pDst_image[3 * tmp_ind] = pDst_image[3 * tmp_ind + 1] =
			pDst_image[3 * tmp_ind + 2] = (pDst_image[3 * (tmp_ind - 1)] +
			                              pG_image[3 * (tmp_ind - 1)]);
			//neboli: OUTPUT_BW[x][y] = OUTPUT_BW[x-1][y] + Grad_X[x-1][y]; 
		}
	}
}*/

//______________________________________________________________________________
void TMOColor2Gray::integrate(TMOImage& G_image, TMOImage& Dst_image)
{
	const unsigned rows = Dst_image.GetHeight(),
	               cols = Dst_image.GetWidth();
	const cl::buffer grad{env, CL_MEM_READ_WRITE, rows * cols * 3 *
	                      sizeof(double)},
	                 out{env, CL_MEM_READ_WRITE, rows * cols * 3 *
	                     sizeof(double)};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   G_image.GetData())};

	exe["integrate2x"].set_args(grad, out, rows, cols);
	status = step.ndrange_kernel(exe["integrate2x"], {0},
	                             {com::math::ceil_multiple(rows, wgs)},
	                             {wgs}, {status});

	status = step.read_buffer(out, 0, 3 * rows * cols * sizeof(double),
	                          Dst_image.GetData(), {status});
}

/* --------------------------------------------------------------------------- *
 * Gradient inconsistency correction                                           *
 * --------------------------------------------------------------------------- */
void TMOColor2Gray::inconsistencyCorrection(TMOImage& G_image,
                                            TMOImage& divG_image,
                                            const double eps)
{
	long xmax = G_image.GetWidth(),
	     ymax = G_image.GetHeight();
	long tmp_y, tmp_ind, i, j, iter = 0;
	double* pG_image = G_image.GetData(),
	      * pdivG_image = divG_image.GetData();
	double E = 0, tmp_divG, maxE = 0, avgE = 0;

	do {
		avgE = maxE = 0;
		for (i = 0; i < ymax; ++i) {
			tmp_y = i * xmax;
			for (j = 0; j < xmax; ++j) {
				tmp_ind = j + tmp_y;
				E = pG_image[3 * tmp_ind] - //Gx
				    pG_image[3 * tmp_ind + 1] + //Gy
				    ((j + 1< xmax) ? (pG_image[3 * (tmp_ind + 1) + 1]) : 0.) - //Gy(i+1,j)
				    ((i + 1< ymax) ? (pG_image[3 * (tmp_ind + xmax)]) : 0.); //Gx(i,j+1)

				if (fabs(E) > maxE)
					maxE = fabs(E);
				avgE += fabs(E);

				E *= .25 * s;

				pG_image[3 * tmp_ind] -= E;	//Gx
				pG_image[3 * tmp_ind+1] += E;	//Gy
				if(j + 1 < xmax)
					pG_image[3 * (tmp_ind + 1) + 1] -= E;		//Gy(i+1,j)
				if(i + 1 < ymax)
					pG_image[3 * (tmp_ind + xmax)] += E;	//Gx(i,j+1)

				tmp_divG = (((j - 1) < 0) ? (pG_image[3 * tmp_ind]) :
				           (pG_image[3 * tmp_ind] - pG_image[3 * (tmp_ind - 1)]));
				tmp_divG= tmp_divG + (((i - 1) < 0) ? (pG_image[3 * tmp_ind + 1]) :
				          (pG_image[3 * tmp_ind + 1] - pG_image[3 * (tmp_ind - xmax) + 1]));
				pdivG_image[3 * tmp_ind] = tmp_divG;
				pdivG_image[3 * tmp_ind + 1] = tmp_divG;
				pdivG_image[3 * tmp_ind + 2] = tmp_divG;
			}
		}
		iter++;
		std::cout << "maxE: " << maxE << std::endl;

		if(!(iter % 100))
			printf("Iteration: %d,\t maxE=%8.2g, \tavgE=%8.2g \n", iter, maxE, avgE);
	} while (maxE > eps);
}

/* --------------------------------------------------------------------------- *
 * Gradient field integration                                                  *
 * --------------------------------------------------------------------------- */
void TMOColor2Gray::GFintegration(TMOImage& G_image, TMOImage& Dst_image)
{
	long xmax = Dst_image.GetWidth(),
	     ymax = Dst_image.GetHeight();
	long tmp_y, tmp_ind, i, j;
	double* pG_image = G_image.GetData(),
	      * pDst_image = Dst_image.GetData();

	pDst_image[0] = pDst_image[1] = pDst_image[2] = 0.;
	for (i = 0; i < ymax; ++i) {
		tmp_y = i*xmax;
		if (i > 0)
			pDst_image[3 * tmp_y] = pDst_image[3 * tmp_y + 1] =
			pDst_image[3*tmp_y+2] = (pDst_image[3 * (tmp_y - xmax)] +
			                        pG_image[3 * (tmp_y - xmax) + 1]);
			//neboli: OUTPUT_BW[0][y] = OUTPUT_BW[0][y-1] + Grad_Y[0][y-1];

		for (j = 1; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			pDst_image[3 * tmp_ind] = pDst_image[3 * tmp_ind + 1] =
			pDst_image[3 * tmp_ind + 2] = (pDst_image[3 * (tmp_ind - 1)] +
			                              pG_image[3 * (tmp_ind - 1)]);
			//neboli: OUTPUT_BW[x][y] = OUTPUT_BW[x-1][y] + Grad_X[x-1][y]; 
		}
	}
}

/* --------------------------------------------------------------------------- *
 * Calibration of the output values                                            *
 * --------------------------------------------------------------------------- */
void TMOColor2Gray::calibrate(TMOImage& src_image, TMOImage& dst_image){
	long i, j, tmp_y, xmax = src_image.GetWidth(),
	     ymax = src_image.GetHeight();
	double* pSrc_image = src_image.GetData(),
	      * pDst_image=dst_image.GetData();
	double SUM_L_new = 0, SUM_L_old = 0,
	       SUM_L2_new = 0, SUM_L_newL_old = 0;
	double	A = 0, B = 0;

	//assert: src_image je v Yxy
	//assert: dst_image je v RGB a to v stupnich sedi

	for (i = 0; i < ymax; ++i) {
		tmp_y = i * xmax;
		for (j = 0; j < xmax ; ++j) {
			SUM_L_new += pDst_image[3 * (tmp_y + j)];
			SUM_L_old += pSrc_image[3 * (tmp_y + j)];
			SUM_L2_new += pDst_image[3 * (tmp_y + j)] * pDst_image[3 * (tmp_y + j)];
			SUM_L_newL_old += pDst_image[3 * (tmp_y + j)] * pSrc_image[3 * (tmp_y + j)];
		}
	}

	if (SUM_L2_new - SUM_L_new * SUM_L_new != 0)
		B = (SUM_L_newL_old - SUM_L_new * SUM_L_old) /
		    (SUM_L2_new - SUM_L_new * SUM_L_new);
	A = (SUM_L_old - B * SUM_L_new) / (xmax * ymax); 
	printf("Normalization: A+B*L == %g+%g*L\n", A, B);

	for (i = 0; i < ymax; ++i) {
		tmp_y = i * xmax;
		for (j = 0; j < xmax; ++j) {
			pDst_image[3 * (tmp_y + j)]= 0.01 * (A + B * pDst_image[3 * (tmp_y + j)]);
			pDst_image[3 * (tmp_y + j) + 1]= .01 * (A + B * pDst_image[3 * (tmp_y + j) + 1]);
			pDst_image[3 * (tmp_y + j) + 2]= .01 * (A + B * pDst_image[3 * (tmp_y + j) + 2]);

			pSrc_image[3 * (tmp_y + j)] *= 0.01;
			pSrc_image[3 * (tmp_y + j) + 1] = pSrc_image[3 * (tmp_y + j)];
			pSrc_image[3 * (tmp_y + j) + 2] = pSrc_image[3 * (tmp_y + j)];
		}
	}
}

/* --------------------------------------------------------------------------- *
 * An implementation of Poisson solver                                         *
 * --------------------------------------------------------------------------- */
int TMOColor2Gray::Transform()
{
	int i=0, j=0;
	double L_max=0.,  // max luminance 
		   L_min=0.,  // min luminance
		   L_temp=0.;   // 
	double dStonits=pSrc->GetStonits();
	pSrc->Convert(TMO_Yxy);
	pSrc->GetMinMaxAvg(&L_min, &L_max, &L_temp);
	fprintf(stdout, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * dStonits, L_max * dStonits, L_temp * dStonits);
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB, true);
	xmax=pSrc->GetWidth();
	ymax=pSrc->GetHeight();
	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();
	int max,k;
	int tmp_y,tmp_ind;
	char filename[500];
	strcpy(filename, pSrc->GetFilename());



///////////////////////////////////////////////////////////////////////////////////////
// Vypocet gradientu

	int max_square_pow2=std::max(xmax, ymax);
	double log2=log((double)max_square_pow2) / log(2.0);
	if(max_square_pow2 > pow(2.0, floor(log2))) 
		max_square_pow2=pow(2.0, ceil(log2));
	max=max_square_pow2*max_square_pow2;
	int jmax = max_square_pow2;
	++jmax;

	double** H=alloc_image1(max_square_pow2, max_square_pow2);
	multiply_image1(H, max_square_pow2, max_square_pow2, 0.0); 
	vect2D *nablaH=new vect2D[max];
	double tmp_divG;
	TMOImage G_image;
	TMOImage divG_image;

	/*if(!NEW_FORMULA){
	/////////////////////////////////////////////////////////////
	//////vypocet gradientu z luminance
	//// vypocte hodnoty luminance pro vsechny pixely a zaroven log luminance (H)
	pSrc->Convert(TMO_Yxy);
	for (i = 0; i < ymax; i++)
	{
		for (int j = 0; j < xmax; j++)
		{
			double L_w = *pSourceData++;
			double x_w = *pSourceData++;
			double y_w = *pSourceData++;

			H[i][j]=TAKE_LOG(L_w);
			//H[i][j]=L_w;
	
			*pDestinationData++ = L_w ;  
			*pDestinationData++ = x_w; //colors remain unchanged
			*pDestinationData++ = y_w; //colors remain unchanged
		}
	}	
	
	// vypocte hodnoty gradientu H (nabla H)
	// a odhadne parametr alpha jako 0.1 * Avg ||\Nabla H_k(x, y)||
	double avg_gradient=0;
	for (i = 0; i < ymax; i++)
	{
		tmp_y = i*xmax;
		for (j = 0; j < xmax; j++)
		{
			tmp_ind=j+tmp_y;
			nablaH[tmp_ind].x=(((j+1)==xmax)?0.0:(H[i][j+1]-H[i][j]));
			nablaH[tmp_ind].y=(((i+1)==ymax)?0.0:(H[i+1][j]-H[i][j]));
			avg_gradient+=sqrt(nablaH[tmp_ind].x*nablaH[tmp_ind].x + nablaH[tmp_ind].y*nablaH[tmp_ind].y);
		}
	}
	avg_gradient/=(xmax*ymax);
	free_image1(H, max_square_pow2, max_square_pow2);
	
	// vypocteme hodnoty divG
	//a taky si je ulozime G a divG -> dump:
	G_image.New(xmax, ymax); //Gx,Gy,DivG
	divG_image.New(xmax, ymax);
	double *pG_image = G_image.GetData();
	double *pdivG_image = divG_image.GetData();
	int image_tmp_y=0, image_tmp_ind=0;

	for ( i = 0; i < ymax; i++ )
	{
		tmp_y = i*xmax;
		image_tmp_y = i*jmax;
		for ( j = 0; j < xmax; j++ )
		{
			tmp_ind=j+tmp_y;
			image_tmp_ind=j+image_tmp_y;
			tmp_divG=(((j-1)<0)?(nablaH[tmp_ind].x):(nablaH[tmp_ind].x-nablaH[tmp_ind-1].x));
			tmp_divG=tmp_divG+(((i-1)<0)?(nablaH[tmp_ind].y):(nablaH[tmp_ind].y-nablaH[(i-1)*xmax+j].y));//

			pG_image[3*tmp_ind]=fabs(nablaH[tmp_ind].x); //Gx
			pG_image[3*tmp_ind+1]=fabs(nablaH[tmp_ind].y); //Gy
			pG_image[3*tmp_ind+2]=0; //

			pdivG_image[3*tmp_ind]=tmp_divG;
			pdivG_image[3*tmp_ind+1]=tmp_divG;
			pdivG_image[3*tmp_ind+2]=tmp_divG;
		}
	}

	G_image.SetFilename(filename);
	G_image.SaveWithSuffix("G");
	//G_image.SaveWithSuffix("G", TMO_RAW);	

	divG_image.SetFilename(filename);
	divG_image.SaveWithSuffix("divG");
	//divG_image.SaveWithSuffix("divG", TMO_RAW);	
	}
	else*/
	{
	/////////////////////////////////////////////////////////////
	//////vypocet gradientu z Coloroidu
	//prevod do XYZ
	READ_SPECTRUM_DATAE();
	READ_COLOROID_DATAE();
	READ_COLOR2GRAY_DATAE();
	_7_basic_fi_computation();

	double X, Y, Z;
	for (i = 0; i < ymax; i++) {
		tmp_y = i*xmax;
		for (j = 0; j < xmax; ++j) {
			tmp_ind = tmp_y + j;

			RGB709_XYZ(Cdisplay_01_Clinear(pSourceData[3 * tmp_ind]),
			           Cdisplay_01_Clinear(pSourceData[3 * tmp_ind + 1]),
			           Cdisplay_01_Clinear(pSourceData[3 * tmp_ind + 2]),
			           &X, &Y, &Z);
			//XYZ_Lab_(X, Y, Z, Xn, Yn, Zn, &pSourceData[3*tmp_ind], &pSourceData[3*tmp_ind+1], &pSourceData[3*tmp_ind+2]);
			pSourceData[3 * tmp_ind] = X;
			pSourceData[3 * tmp_ind + 1] = Y;
			pSourceData[3 * tmp_ind + 2] = Z;
		}
	}

	//// vypocte hodnoty gradientu H (nabla H)
	TMOImage threshold_image;
	threshold_image.New(xmax, ymax);
	double* threshold_data = threshold_image.GetData();
	double avg_gradient = 0;
	for (i = 0; i < ymax; ++i) {
		tmp_y = i*xmax;
		for (j = 0; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			nablaH[tmp_ind].x = ((j + 1) == xmax) ? 0. :
			                    formulaColoroid(pSourceData, i, j + 1, i, j, threshold_data);
				//(H[i][j+1]-H[i][j]));
			nablaH[tmp_ind].y = ((i + 1) == ymax) ? 0. :
			                    formulaColoroid(pSourceData, i+1,j,i,j, threshold_data);
				//(H[i+1][j]-H[i][j]));
			avg_gradient += sqrt(nablaH[tmp_ind].x * nablaH[tmp_ind].x +
			                     nablaH[tmp_ind].y * nablaH[tmp_ind].y);
		}
	}

	avg_gradient /= (xmax * ymax);
	free_image1(H, max_square_pow2, max_square_pow2);
	threshold_image.SetFilename(filename);
	threshold_image.SaveWithSuffix("_threshold");

	// vypocteme hodnoty divG
	//a taky si je ulozime G a divG -> dump:
	G_image.New(xmax, ymax); //Gx,Gy,DivG
	divG_image.New(xmax, ymax);
	double* pG_image = G_image.GetData();
	double* pdivG_image = divG_image.GetData();
	int image_tmp_y = 0,
	    image_tmp_ind = 0;

	for (i = 0; i < ymax; ++i) {
		tmp_y = i * xmax;
		image_tmp_y = i * jmax;
		for (j = 0; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			image_tmp_ind = j + image_tmp_y;
			tmp_divG = ((j - 1) < 0) ? nablaH[tmp_ind].x : (nablaH[tmp_ind].x - nablaH[tmp_ind - 1].x);
			tmp_divG = tmp_divG + (((i - 1) < 0) ? (nablaH[tmp_ind].y) : (nablaH[tmp_ind].y - nablaH[(i - 1) * xmax + j].y));

			pG_image[3 * tmp_ind] = nablaH[tmp_ind].x; //Gx
			pG_image[3 * tmp_ind + 1] = nablaH[tmp_ind].y; //Gy
			pG_image[3 * tmp_ind + 2] = 0; //

			pdivG_image[3 * tmp_ind] = tmp_divG;
			pdivG_image[3 * tmp_ind + 1] = tmp_divG;
			pdivG_image[3 * tmp_ind + 2] = tmp_divG;
		}
	}

	G_image.SetFilename(filename);
	G_image.SaveWithSuffix("G");
	//G_image.SaveWithSuffix("G", TMO_EXR);	

	divG_image.SetFilename(filename);
	divG_image.SaveWithSuffix("divG");
	//divG_image.SaveWithSuffix("divG", TMO_RAW);	

	//inconsistencyCorrection(G_image, divG_image, eps);
	correct_grad(G_image, eps);

	G_image.SetFilename(filename);
	G_image.SaveWithSuffix("G_corrected");
	divG_image.SetFilename(filename);
	divG_image.SaveWithSuffix("divG_corrected");

	GFintegration(G_image, *pDst);
	//GFintegrationOp(G_image, *pDst);

	// XXX WTF IS THIS vvvv SHIT???
	double* pDst_image = pDst->GetData();
	for (i = 0; i < pDst->GetHeight() ; ++i) {
		tmp_y = i * pDst->GetWidth();
		for (j = 0; j < pDst->GetWidth() ; ++j) {
			pDst_image[3 * (tmp_y + j)] = .01 * pDst_image[3 * (tmp_y + j)] * pDst_image[3 * (tmp_y + j)];
			//pDst_image[3*(tmp_y+j)+1]=0.01*pDst_image[3*(tmp_y+j)+1]*pDst_image[3*(tmp_y+j)+1];
			//pDst_image[3*(tmp_y+j)+2]=0.01*pDst_image[3*(tmp_y+j)+2]*pDst_image[3*(tmp_y+j)+2];
			pDst_image[3 * (tmp_y + j) + 1] = pDst_image[3 * (tmp_y + j)];
			pDst_image[3 * (tmp_y + j) + 2] = pDst_image[3 * (tmp_y + j)];
		}
	}

	//pDst->Convert(TMO_Yxy, true);
	//double L_temp2;
	//pDst->GetMinMaxAvg(&L_min, &L_max, &L_temp2);
	////fprintf(stdout, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * dStonits, L_max * dStonits, L_temp * dStonits);


	//pDst->Convert(TMO_RGB);

	//calibrate(*pSrc, *pDst);

	pSrc->Convert(TMO_RGB, true);
	//pSrc->SaveWithSuffix("-inputL", TMO_EXR);
	pSrc->SaveWithSuffix("-inputL");

	pDst->Convert(TMO_Yxy);
	pDst->GetMinMaxAvg(&L_min, &L_max, &L_temp);
	fprintf(stdout, "Min luminance: %g[cd/m^2]\nMax luminance: %g[cd/m^2]\nAvg luminance: %g[cd/m^2]\n", L_min * dStonits, L_max * dStonits, L_temp * dStonits);

	pDst->Convert(TMO_RGB);
	//pDst->CorrectGamma(gamma);

	return 0;

	}
}
