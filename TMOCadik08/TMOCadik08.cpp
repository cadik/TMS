//(c)Martin Cadik
//03--05/2007 - Prague

// XXX
//POZOR - neni dodelane:
//	double ZETA [psi_xmax][psi_ymax];
//	double PSI [psi_xmax][psi_ymax];
// -jsou "natvrdo" maximalni dimenze - pro mcdesk-test.hdr

/* --------------------------------------------------------------------------- *
 * TMOCadik08.cpp: implementation of the TMOCadik08 class.                       *
 * --------------------------------------------------------------------------- */
#include <iostream>
#include <vector>
#include <cmath>
#include "./TMOCadik08.h"
//#include "Laszlo/_common.h"
#include "common/text_loader.h"

static double Cdisplay_01_Clinear(double Cdisplay)
{
	//both are normalized to [0,1]
	//default 709 gamma = 2.2
	double Clinear;

	Clinear =  Cdisplay <= 0.081 ?  Cdisplay / 4.5 : 
	pow((Cdisplay + 0.099) / 1.099, 2.2); 

	return(Clinear);

}

//______________________________________________________________________________
TMOCadik08::TMOCadik08() :
	simd{},
	step{simd.create_command_queue()},
	exe{simd.create_program(com::text_loader{"TMOCadik08/resources/kernels/"
	                                         "color2gray.cl"}().c_str(),
	                        "-ITMOCadik08/resources/kernels")},
	wgs{simd.get_wgs()},
	dim{std::sqrt(wgs)}
{
	SetName(L"Cadik08");
	SetDescription(L"Cadik08 in gradient domain - version 0.1");

	normalize.SetName(L"normalize");
	normalize.SetDescription(L"Normalize output luminance (using 1/Lo_max) to be in <0, 1>");
	normalize.SetDefault(true);
	normalize = true;
	this->Register(normalize);

	gamma.SetName(L"gamma");
	gamma.SetDescription(L"Gamma Correction Factor [-]");
	gamma.SetDefault(2.2);
	gamma = 1.;
	gamma.SetRange(1.0e-3, 1.0e+2);
	this->Register(gamma);

	s.SetName(L"s");
	s.SetDescription(L"overprojection step");
	s.SetDefault(1.);
	s = 1.;
	s.SetRange(0., 1.99);
	this->Register(s);
	
	eps.SetName(L"eps");
	eps.SetDescription(L"epsilon for gradient correction");
	eps.SetDefault(.1);
	eps = .1;
	eps.SetRange(0., 10.);
	this->Register(eps);
}

//______________________________________________________________________________
TMOCadik08::~TMOCadik08()
{
}

//______________________________________________________________________________
int TMOCadik08::Transform()
{
	std::cout << "processing:" << std::endl;
	int i = 0, j = 0;
	pDst->Convert(TMO_RGB, true);
	long xmax = pSrc->GetWidth();
	long ymax = pSrc->GetHeight();
	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();
	int max, k;
	int tmp_y, tmp_ind;

	const std::string filename{pSrc->GetFilename()};

	int max_square_pow2 = std::max(xmax, ymax);
	double log2 = log((double) max_square_pow2) / log(2.);
	if(max_square_pow2 > pow(2., floor(log2))) 
		max_square_pow2 = pow(2., ceil(log2));
	max = max_square_pow2 * max_square_pow2;
	int jmax = max_square_pow2;
	++jmax;

	std::vector<vect2D> nablaH(max);
	TMOImage G_image;

	/////////////////////////////////////////////////////////////
	//////vypocet gradientu z Coloroidu
	//prevod do XYZ
	model.READ_SPECTRUM_DATAE();
	model.READ_COLOROID_DATAE();
	model.READ_COLOR2GRAY_DATAE();
	model._7_basic_fi_computation();

	double X, Y, Z;
	for (i = 0; i < ymax; ++i) {
		tmp_y = i * xmax;
		for (j = 0; j < xmax; ++j) {
			tmp_ind = tmp_y + j;

			model.RGB709_XYZ(Cdisplay_01_Clinear(pSourceData[3 * tmp_ind]),
			                 Cdisplay_01_Clinear(pSourceData[3 * tmp_ind + 1]),
			                 Cdisplay_01_Clinear(pSourceData[3 * tmp_ind + 2]),
			                 &X, &Y, &Z);
			pSourceData[3 * tmp_ind] = X;
			pSourceData[3 * tmp_ind + 1] = Y;
			pSourceData[3 * tmp_ind + 2] = Z;
		}
	}

	//// vypocte hodnoty gradientu H (nabla H)
	TMOImage threshold_image;
	threshold_image.New(xmax, ymax);
	double* threshold_data = threshold_image.GetData();
	//double avg_gradient = 0;
	for (i = 0; i < ymax; ++i) {
		tmp_y = i * xmax;
		for (j = 0; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			nablaH[tmp_ind].x = ((j + 1) == xmax) ? 0. :
			                    formula_coloroid(pSourceData, i, j + 1, i, j, xmax, threshold_data);
				//(H[i][j+1]-H[i][j]));
			nablaH[tmp_ind].y = ((i + 1) == ymax) ? 0. :
			                    formula_coloroid(pSourceData, i + 1, j, i, j, xmax, threshold_data);
				//(H[i+1][j]-H[i][j]));
		}
	}

	threshold_image.SetFilename(filename.c_str());
	threshold_image.SaveWithSuffix("_threshold");

	G_image.New(xmax, ymax); //Gx,Gy,DivG
	double* pG_image = G_image.GetData();
	int image_tmp_y = 0;

	for (i = 0; i < ymax; ++i) {
		tmp_y = i * xmax;
		image_tmp_y = i * jmax;
		for (j = 0; j < xmax; ++j) {
			tmp_ind = j + tmp_y;

			pG_image[3 * tmp_ind] = nablaH[tmp_ind].x; //Gx
			pG_image[3 * tmp_ind + 1] = nablaH[tmp_ind].y; //Gy
			pG_image[3 * tmp_ind + 2] = 0.; //
		}
	}

	G_image.SetFilename(filename.c_str());
	G_image.SaveWithSuffix("G");

	//inconsistencyCorrection(G_image, eps); // XXX CPU version
	correct_grad(G_image, eps); // zkonverguje, ale chova se divne...

	G_image.SetFilename(filename.c_str());
	G_image.SaveWithSuffix("G_corrected");

	//GFintegration(G_image, *pDst); // XXX CPU version
	integrate2x(G_image, *pDst);
	//trans_range(*pDst, 0., 255.);

	pDst->Convert(TMO_RGB);
	// XXX what is this??? vvvvv 
	/*double* pDst_image = pDst->GetData();
	for (i = 0; i < pDst->GetHeight() ; ++i) {
		tmp_y = i * pDst->GetWidth();
		for (j = 0; j < pDst->GetWidth() ; ++j) {
			pDst_image[3 * (tmp_y + j)] = .01 * pDst_image[3 * (tmp_y + j)] * pDst_image[3 * (tmp_y + j)];
			//pDst_image[3*(tmp_y+j)+1]=0.01*pDst_image[3*(tmp_y+j)+1]*pDst_image[3*(tmp_y+j)+1];
			//pDst_image[3*(tmp_y+j)+2]=0.01*pDst_image[3*(tmp_y+j)+2]*pDst_image[3*(tmp_y+j)+2];
			pDst_image[3 * (tmp_y + j) + 1] = pDst_image[3 * (tmp_y + j)];
			pDst_image[3 * (tmp_y + j) + 2] = pDst_image[3 * (tmp_y + j)];
		}
	}*/

	calibrate(*pSrc, *pDst);
	pDst->SetFilename(filename.c_str());
	pDst->SaveWithSuffix("_out");

	//pDst->CorrectGamma(gamma);

	return 0;
}

//______________________________________________________________________________
double TMOCadik08::formula_coloroid(const double* const data,
                                    const long y1, const long x1,
                                    const long y2, const long x2,
                                    const long xmax,
                                    double* const threshold_data)
{
	const double threshold = 1.;
	double X1, Y1, Z1,
	       X2, Y2, Z2,
	       A_hue, T, V,
	       grad_luminance,
	       dA, dT, dV;

	long tmp_ind_1 = y1 * xmax + x1,
	     tmp_ind_2 = y2 * xmax + x2;

	X1 = data[3 * tmp_ind_1];
	X2 = data[3 * tmp_ind_2];
	Y1 = data[3 * tmp_ind_1 + 1];
	Y2 = data[3 * tmp_ind_2 + 1];
	Z1 = data[3 * tmp_ind_1 + 2];
	Z2 = data[3 * tmp_ind_2 + 2];

	threshold_data[3 * tmp_ind_1] = data[3 * tmp_ind_1] * .01;
	threshold_data[3 * tmp_ind_1 + 1] = data[3 * tmp_ind_1] * .01;
	threshold_data[3 * tmp_ind_1 + 2] = data[3 * tmp_ind_1] * .01;

	//the next function defines the grad value, starting from (X1,Y1,Z1) and (X2,Y2,Z2) neighbor pixels
	//grad_LUMINANCE is measured in Coloroid luminance V
	//the integral of consitent image will be obtained in Coloroid V luminance
	//CIE Y = 0.01 * (V * V) is the formula to Y !!!!!
	//pipeline: from Y the rgb and after the gamma-RGB

	//51 45 46   blue
	//12 50 82   yellow  A, T, V
	//X1 = 25.87;  Y1 = 21.16;  Z1 = 79.7;
	//X2 = 61.81;  Y2 = 67.24;  Z2 = 23.2;

	grad_luminance = model.LUMINANCE_GRAD(X1, Y1, Z1,
	                                X2, Y2, Z2,
	                                &dA, &dT, &dV);

	//A , T , V  are the separately not important hue, saturation and luminance parts

	//CURRENTLY:    #define WEIGHT_LIGHTNESS_CHROMINANCE 3.0 (color_datae.h)
	//this weight regultes the luminance chrominance ratio
	//for higher chrominance effect it has to be increased the above default value

	return grad_luminance;
}

//______________________________________________________________________________
void TMOCadik08::correct_grad(TMOImage& g, const double eps) const
{
	const unsigned rows = g.GetHeight(),
	               cols = g.GetWidth();
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols * 3 *
	                                         sizeof(double))},
	                 err{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        rows * cols * sizeof(double))};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   g.GetData())};

	double e_max;
	const double s = this->s.GetDouble();
	do {
		e_max = 0.;

		exe["calc_error"].set_args(grad, err, rows, cols);
		status = step.ndrange_kernel(exe["calc_error"], {},
		                             {rows, cols}, {dim, dim},
		                             {status});

		for (unsigned n = 0; n < 4; ++n) {
			exe["correct_grad"].set_args(grad, err, s, rows, cols, n);
			status = step.ndrange_kernel(exe["correct_grad"], {},
			                             {rows, cols}, {dim, dim},
			                             {status});
		}

		status = reduce("reduce_absmax", err, rows * cols, e_max,
		                {status});

		std::cout << "e_max: " << e_max << std::endl;
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * 3 * sizeof(double),
	                          g.GetData(), {status});
	step.wait({status});
}

//______________________________________________________________________________
cl::event TMOCadik08::reduce(const std::string type, const cl::buffer& in,
                             const unsigned n, double& out,
                             const cl::event_list pending) const
{
	// number of work-groups needed to reduce the problem
	unsigned m = std::ceil(static_cast<float>(com::math::ceil2mul(n,
	                       2 * wgs)) / (2.f * wgs));
	// partial arrays
	const cl::buffer tmp{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        m * sizeof(double))};

	exe[type].set_args(in, cl::local_mem{2 * wgs * sizeof(cl_double)},
	                   tmp, n);
	cl::event status{step.ndrange_kernel(exe[type], {}, {m * wgs}, {wgs},
	                                     pending)};
	if (m > 1)
		status = reduce(type, tmp, m, out, {status});
	else
		status = step.read_buffer(tmp, 0, sizeof(cl_double),
		                          &out, {status});

	return status;
}

//______________________________________________________________________________
void TMOCadik08::integrate2x(TMOImage& g, TMOImage& o) const
{
	const unsigned rows = o.GetHeight(),
	               cols = o.GetWidth();
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols * 3 *
	                                         sizeof(double))},
	                 out{simd.create_buffer(CL_MEM_READ_WRITE |
	                                        CL_MEM_USE_HOST_PTR,
	                                        rows * cols * 3 * sizeof(double),
	                                        o.GetData())};
	cl::event status{step.write_buffer<false>(grad, 0, rows * cols *
	                                          3 * sizeof(double),
	                                          g.GetData())};

	// sequentially initialize the first column in +y direction.
	// this has a strict data-dependency, and can't be done in parallel
	double* g_data = g.GetData(),
	      * o_data = o.GetData();
	o_data[0] = o_data[1] = o_data[2] = 0.;
	for (size_t y = 1; y < rows; ++y) {
		const size_t i = y * cols;
		if (y > 0)
			o_data[3 * i] =
			o_data[3 * i + 1] =
			o_data[3 * i + 2] = o_data[3 * (i - cols)] +
			                    g_data[3 * (i - cols) + 1];
	}

	// for each row, perform integration in +x direction
	exe["integrate2x"].set_args(grad, out, rows, cols);
	status = step.ndrange_kernel(exe["integrate2x"], {},
	                             {rows}, {wgs}, {status});

	// XXX
	//status = step.read_buffer(out, 0, 3 * rows * cols * sizeof(double),
	//                          o.GetData(), {status});
	step.wait({status});
}

//______________________________________________________________________________
void TMOCadik08::trans_range(TMOImage& i, const double new_min,
                             const double new_max) const
{
	const unsigned rows = i.GetHeight(),
	               cols = i.GetWidth();
	const cl::buffer in{simd.create_buffer(CL_MEM_READ_WRITE,
	                                       rows * cols * 3 *
	                                       sizeof(double))};
	cl::event status{step.write_buffer(in, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   i.GetData())};

	double old_min, old_max;
	status = reduce("reduce_min", in, 3 * rows * cols, old_min, {status});
	status = reduce("reduce_max", in, 3 * rows * cols, old_max, {status});

	exe["trans_range"].set_args(in, old_min, old_max, new_min, new_max,
	                            rows * cols);
	status = step.ndrange_kernel(exe["trans_range"], {}, {rows * cols},
	                             {wgs}, {status});

	status = step.read_buffer(in, 0, rows * cols * 3 * sizeof(double),
	                          i.GetData(), {status});
}

/* --------------------------------------------------------------------------- *
 * Gradient inconsistency correction -- CPU                                    *
 * --------------------------------------------------------------------------- */
void TMOCadik08::inconsistencyCorrection(TMOImage& G_image,
                                            const double eps)
{
	long xmax = G_image.GetWidth(),
	     ymax = G_image.GetHeight();
	long tmp_y, tmp_ind, i, j;
	double* pG_image = G_image.GetData();
	double E = 0., maxE = 0.;

	do {
		maxE = 0.;
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

				E *= .25 * s;

				pG_image[3 * tmp_ind] -= E;	//Gx
				pG_image[3 * tmp_ind+1] += E;	//Gy
				if(j + 1 < xmax)
					pG_image[3 * (tmp_ind + 1) + 1] -= E;		//Gy(i+1,j)
				if(i + 1 < ymax)
					pG_image[3 * (tmp_ind + xmax)] += E;	//Gx(i,j+1)
			}
		}
		std::cout << "maxE: " << maxE << std::endl;
	} while (maxE > eps);
}

/* --------------------------------------------------------------------------- *
 * Gradient field integration -- CPU                                           *
 * --------------------------------------------------------------------------- */
void TMOCadik08::GFintegration(TMOImage& G_image, TMOImage& Dst_image)
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
void TMOCadik08::calibrate(TMOImage& src_image, TMOImage& dst_image){
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

//______________________________________________________________________________
// MAXIMUM ERROR SELECTION
/*void TMOCadik08::correct_grad_mer(TMOImage& g, const double eps)
{
	// maximum error selection
	const unsigned rows = g.GetHeight(),
	               cols = g.GetWidth();
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols * 3 *
	                                         sizeof(double))},
	                 err{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        rows * cols * sizeof(double))},
	                 dummy{simd.create_buffer(CL_MEM_READ_WRITE, 1)};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   g.GetData())};

	double e_max;
	do {
		e_max = 0.;

		exe["calc_error"].set_args(grad, err, rows, cols);
		status = step.ndrange_kernel(exe["calc_error"], {0, 0},
		                             {rows, cols}, {dim, dim},
		                             {status});

		unsigned uv = 0;
		status = reduce_maxi(err, dummy, rows * cols, 0, e_max,
		                     uv, {status});

		// correct the one with e_max; identified by uv
		exe["correct_gradi"].set_args(grad, err, s.GetDouble(),
		                              uv, rows, cols);
		status = step.ndrange_kernel(exe["correct_gradi"], {0, 0},
		                             {rows, cols}, {dim, dim}, {status});

		std::cout << "e_max: " << e_max << std::endl;
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * 3 * sizeof(double),
	                          g.GetData());
}

//______________________________________________________________________________
// REDUCE_MAX outputting e_max index also
cl::event TMOCadik08::reduce_maxi(const cl::buffer& in,
                                     const cl::buffer& inds,
                                     const unsigned n,
                                     const unsigned use_inds,
                                     double& out, unsigned& uv,
                                     const std::vector<cl::event> pending)
{
	// number of work-groups needed to reduce the problem
	unsigned m = std::ceil(static_cast<float>(com::math::ceil2mul(n,
	                       2 * wgs)) / (2.f * wgs));
	// partial maximas
	const cl::buffer maximas{simd.create_buffer(CL_MEM_READ_WRITE,
	                                            m * sizeof(double))},
	                 coords{simd.create_buffer(CL_MEM_READ_WRITE,
	                                           m * sizeof(unsigned))};

	exe["reduce_maxi"].set_args(in, inds, use_inds,
	                            cl::local_mem{2 * wgs * sizeof(cl_double)},
	                            cl::local_mem{2 * wgs * sizeof(cl_uint)},
	                            maximas, coords, n);
	cl::event status{step.ndrange_kernel(exe["reduce_maxi"], {0},
	                                     {m * wgs}, {wgs},
	                                     pending)};
	if (m > 1)
		status = reduce_max(maximas, coords, m, 1, out, uv, {status});
	else {
		status = step.read_buffer(maximas, 0, sizeof(cl_double),
		                          &out, {status});
		status = step.read_buffer(coords, 0, sizeof(cl_uint),
		                          &uv, {status});
	}

	return status;
}*/
//______________________________________________________________________________
// correct_grad version of interleaved elementary loops using local memory
/*
void TMOCadik08::correct_grad(TMOImage& g, const double eps) const
{
	// maximum error selection
	const unsigned rows = g.GetHeight(),
	               cols = g.GetWidth();
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols * 3 *
	                                         sizeof(double))},
	                 err{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        rows * cols * sizeof(double))};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   g.GetData())};

	double e_max;
	do {
		e_max = 0.;

		exe["correct_grad"].set_args(grad,
		                             cl::local_mem{2 * (dim + 1) *
		                             (dim + 1)}, err, s.GetDouble(),
		                             rows, cols);
		status = step.ndrange_kernel(exe["correct_grad"], {},
		                             {rows, cols}, {dim, dim},
		                             {status});

		status = reduce("reduce_absmax", err, rows * cols, e_max,
		                {status});

		std::cout << "e_max: " << e_max << std::endl;
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * 3 * sizeof(double),
	                          g.GetData());
}*/

