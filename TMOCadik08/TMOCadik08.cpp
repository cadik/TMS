// original authors: Laszlo Neumann,
//                   Martin Cadik

//(c)Martin Cadik
//03--05/2007 - Prague

/* --------------------------------------------------------------------------- *
 * TMOCadik08.cpp: implementation of the TMOCadik08 class.                       *
 * --------------------------------------------------------------------------- */
#include <iostream>
#include <vector>
#include <cmath>
#include "./TMOCadik08.h"
#include "common/text_loader.h"
#include "quadtree.h"
#include "morton.h"

//______________________________________________________________________________
static double Cdisplay01Clinear(double Cdisplay)
{
	//both are normalized to [0,1]
	//default 709 gamma = 2.2
	double Clinear;

	Clinear = Cdisplay <= 0.081 ?  Cdisplay / 4.5 : 
	          pow((Cdisplay + 0.099) / 1.099, 2.2); 

	return Clinear;
}

//______________________________________________________________________________
TMOCadik08::TMOCadik08() :
	simd{CL_DEVICE_TYPE_DEFAULT},
#ifdef PROFILE
        step{simd.create_command_queue(CL_QUEUE_PROFILING_ENABLE)},
#else
	step{simd.create_command_queue()},
#endif
	exe{simd.create_program(com::text_loader{"TMOCadik08/resources/kernels/"
	                                         "color2gray.cl"}().c_str(),
	                        "-ITMOCadik08/resources/kernels")},
	wgs{simd.get_wgs()},
	dim{std::sqrt(wgs)},
	maba{simd.get_maba()}
{
	model.readSpectrumDatae();
	model.readColoroidDatae();
	model.readColor2GrayDatae();
	model._7_basic_fi_computation();

	SetName(L"Cadik08");
	SetDescription(L"Color to grayscale conversion in gradient domain");

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

	type.SetName(L"type");
	type.SetDescription(L"type of gradient correction:\n"
	                      "cyc - naive cyclical correction,\n"
	                      "cb - chessboard-style correction,\n"
	                      "cbloc - chessboard-style correction,\n"
	                      "quad - hierarchical quadtree correction,\n"
	                      "cpu - sequential version");
	type.SetDefault("cbloc");
	type = "cbloc";
	this->Register(type);
}

//______________________________________________________________________________
TMOCadik08::~TMOCadik08()
{
}

//______________________________________________________________________________
int TMOCadik08::Transform()
{
	std::cerr << "processing:" << std::endl;
	pDst->Convert(TMO_RGB, true);
	long xmax = pSrc->GetWidth(),
	     ymax = pSrc->GetHeight();
	double* pSourceData = pSrc->GetData(),
	      * pDestinationData = pDst->GetData();

	const std::string filename{pSrc->GetFilename()};

	int max_square_pow2 = std::max(xmax, ymax);
	double log2 = log((double) max_square_pow2) / log(2.);
	if(max_square_pow2 > pow(2., floor(log2))) 
		max_square_pow2 = pow(2., ceil(log2));
	int max = max_square_pow2 * max_square_pow2;
	int jmax = max_square_pow2;
	++jmax;

	// vypocet gradientu z Coloroidu
	// prevod do XYZ
	double X, Y, Z;
	for (int i = 0; i < ymax; ++i) {
		int tmp_y = i * xmax;
		for (int j = 0; j < xmax; ++j) {
			int tmp_ind = tmp_y + j;

			model.rgb709Xyz(Cdisplay01Clinear(pSourceData[3 * tmp_ind]),
			                Cdisplay01Clinear(pSourceData[3 * tmp_ind + 1]),
			                Cdisplay01Clinear(pSourceData[3 * tmp_ind + 2]),
			                &X, &Y, &Z);
			pSourceData[3 * tmp_ind] = X;
			pSourceData[3 * tmp_ind + 1] = Y;
			pSourceData[3 * tmp_ind + 2] = Z;
		}
	}

	TMOImage G_image;
	G_image.New(xmax, ymax);
	double* pG_image = G_image.GetData();

	if (type.GetString() == "cyc" || type.GetString() == "cb" || type.GetString() == "cbloc" || type.GetString() != "hier") {
		std::vector<vec2d> nablaH(max);
		for (int i = 0; i < ymax; ++i) {
			int tmp_y = i * xmax;
			for (int j = 0; j < xmax; ++j) {
				int tmp_ind = j + tmp_y;
				// H[i][j + 1] - H[i][j]:
				nablaH[tmp_ind].x = (j + 1) == xmax ? 0. :
								      formulaColoroid(pSourceData, i, j + 1, i, j, xmax);
				// H[i + 1][j] - H[i][j]:
				nablaH[tmp_ind].y = (i + 1) == ymax ? 0. :
								      formulaColoroid(pSourceData, i + 1, j, i, j, xmax);
			}
		}

		std::cerr << type.GetString() << std::endl;
		if (type.GetString() == "cb")
			correctGradCb(nablaH, ymax, xmax, eps);
		else if (type.GetString() == "cbloc")
			correctGradCbLoc(nablaH, ymax, xmax, eps);
		else if (type.GetString() == "cyc")
			correctGradCyc(nablaH, ymax, xmax, eps);
		else
			inconsistencyCorrection(nablaH, ymax, xmax, eps);

		for (int i = 0; i < ymax; ++i) {
			int tmp_y = i * xmax;
			int image_tmp_y = i * jmax;
			for (int j = 0; j < xmax; ++j) {
				int tmp_ind = j + tmp_y;

				pG_image[3 * tmp_ind] = nablaH[tmp_ind].x; // Gx
				pG_image[3 * tmp_ind + 1] = nablaH[tmp_ind].y; // Gy
				pG_image[3 * tmp_ind + 2] = 0.;
			}
		}
	}
	else if (type.GetString() == "hier") {
		quadtree nablaH(std::max(xmax, ymax));
		for (int i = 0; i < ymax; ++i)
			for (int j = 0; j < xmax; ++j) {
				const morton z = morton2d_encode(j, i);
				// H[i][j + 1] - H[i][j]:
				nablaH[z].x = (j + 1) == xmax ? 0. :
								formulaColoroid(pSourceData, i, j + 1, i, j, xmax);
				// H[i + 1][j] - H[i][j]:
				nablaH[z].y = (i + 1) == ymax ? 0. :
								formulaColoroid(pSourceData, i + 1, j, i, j, xmax);
			}

		correctGradHier(nablaH, eps);

		for (int i = 0; i < ymax; ++i) {
			int tmp_y = i * xmax;
			int image_tmp_y = i * jmax;
			for (int j = 0; j < xmax; ++j) {
				int tmp_ind = j + tmp_y;

				const morton z = morton2d_encode(j, i);
				pG_image[3 * tmp_ind] = nablaH[z].x; // Gx
				pG_image[3 * tmp_ind + 1] = nablaH[z].y; // Gy
				pG_image[3 * tmp_ind + 2] = 0.;
			}
		}
	}

	G_image.SetFilename(filename.c_str());
	G_image.SaveWithSuffix("G_corrected");

	integrate2x(G_image, *pDst);
	//GFintegration(G_image, *pDst);

	pDst->Convert(TMO_RGB);

	calibrate(*pSrc, *pDst);
	pDst->SetFilename(filename.c_str());

	//pDst->CorrectGamma(gamma);

	return 0;
}

//______________________________________________________________________________
double TMOCadik08::formulaColoroid(const double* const data,
                                   const long y1, const long x1,
                                   const long y2, const long x2,
                                   const long xmax)
{
	long i = y1 * xmax + x1,
	     j = y2 * xmax + x2;

	const double X1 = data[3 * i],
	             X2 = data[3 * j],
	             Y1 = data[3 * i + 1],
	             Y2 = data[3 * j + 1],
	             Z1 = data[3 * i + 2],
	             Z2 = data[3 * j + 2];

	//the next function defines the grad value, starting from (X1,Y1,Z1) and (X2,Y2,Z2) neighbor pixels
	//grad_LUMINANCE is measured in Coloroid luminance V
	//the integral of consitent image will be obtained in Coloroid V luminance
	//CIE Y = 0.01 * (V * V) is the formula to Y !!!!!
	//pipeline: from Y the rgb and after the gamma-RGB

	//51 45 46   blue
	//12 50 82   yellow  A, T, V
	//X1 = 25.87;  Y1 = 21.16;  Z1 = 79.7;
	//X2 = 61.81;  Y2 = 67.24;  Z2 = 23.2;

	//A , T , V  are the separately not important hue, saturation and luminance parts

	//CURRENTLY:    #define WEIGHT_LIGHTNESS_CHROMINANCE 3.0 (color_datae.h)
	//this weight regultes the luminance chrominance ratio
	//for higher chrominance effect it has to be increased the above default value
	double dA, dT, dV;
	return model.luminanceGrad(X1, Y1, Z1,
	                           X2, Y2, Z2,
	                           &dA, &dT, &dV);
}

//==============================================================================
// chessboard version
//______________________________________________________________________________
void TMOCadik08::correctGradCb(std::vector<vec2d>& g, const unsigned rows,
                               const unsigned cols, const double eps) const
{
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols *
	                                         sizeof(vec2d))},
	                 err{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        rows * cols *
	                                        sizeof(double))};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   sizeof(vec2d),
	                                   g.data())};

#ifdef PROFILE
	int it = 0;
#endif
	double e_max;
	do {
		e_max = 0.;

		for (unsigned char mode = 0; mode < 3; ++mode) {
			exe["cb_correct_grad"].set_args(grad, err,
			                                s.GetDouble(),
			                                rows, cols, mode);
			status = step.ndrange_kernel(exe["cb_correct_grad"],
			                             {}, {rows, cols},
			                             {dim, dim},
			                             {status});
#ifdef PROFILE
			accumulate(status);
#endif
		}

		status = reduce("reduce", err, rows * cols, e_max,
		                {status});

		std::cerr << "e_max: " << e_max << std::endl;
#ifdef PROFILE
		++it;
#endif
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * sizeof(vec2d),
	                          g.data());
#ifdef PROFILE
	std::cerr << "PROFILE: proccessing time: " << acc / 1e9 << " [s]" << std::endl;
#endif
#ifdef PROFILE
	std::cerr << "iteration count: " << it << std::endl;
#endif
}

//______________________________________________________________________________
void TMOCadik08::correctGradCbLoc(std::vector<vec2d>& g, const unsigned rows,
                                  const unsigned cols, const double eps) const
{
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols *
	                                         sizeof(vec2d))},
	                 err{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        rows * cols *
	                                        sizeof(double))};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   sizeof(vec2d),
	                                   g.data())};

#ifdef PROFILE
	int it = 0;
#endif
	double e_max;
	do {
		e_max = 0.;

		exe["cb_correct_grad_l"].set_args(grad, cl::local_mem{(dim + 1) * (dim + 1) * sizeof(vec2d)}, err,
		                                s.GetDouble(), rows, cols);
		status = step.ndrange_kernel(exe["cb_correct_grad_l"],
		                             {}, {rows, cols},
		                             {dim, dim},
		                             {status});
#ifdef PROFILE
		accumulate(status);
#endif

		status = reduce("reduce", err, rows * cols, e_max,
		                {status});

		std::cerr << "e_max: " << e_max << std::endl;
#ifdef PROFILE
		++it;
#endif
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * sizeof(vec2d),
	                          g.data());
#ifdef PROFILE
	std::cerr << "PROFILE: proccessing time: " << acc / 1e9 << " [s]" << std::endl;
#endif
#ifdef PROFILE
	std::cerr << "iteration count: " << it << std::endl;
#endif
}

//==============================================================================
// hierarchical version
//______________________________________________________________________________
cl::event TMOCadik08::evalQuadtree(const cl::buffer& root,
                                   const unsigned height, const unsigned size,
                                   vec2d* const g,
                                   cl::event_list pending) const
{
	cl::event status;
	for (size_t i = height - 2; i; --i, pending = {status}) {
		const unsigned offseti = i > 1 ? ((4 * pow(4, i) - 1) / 3) : 5,
		               sizei = ((4 * pow(4, i + 1) - 1) / 3) - offseti,
		               offsetk = ((4 * pow(4, i + 2) - 1) / 3);
		unsigned aligni;
		const cl::buffer rooti{root.create_sub<vec2d>(0,
		                                              offseti,
		                                              offsetk - offseti,
		                                              maba,
		                                              &aligni)};

		exe["hie_avg_grad"].set_args(rooti, aligni, sizei);
		status = step.ndrange_kernel(exe["hie_avg_grad"],
		                             {}, {sizei}, {wgs},
		                             pending);
#ifdef PROFILE
		accumulate(status);
#endif
	}

	status = step.read_buffer(root, 0, size * sizeof(vec2d),
				  g, {status});

	return status;
}

//______________________________________________________________________________
void TMOCadik08::correctGradHier(quadtree& nablaH, const double eps) const
{
	const cl::buffer root{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         nablaH.size() * sizeof(vec2d))};
	cl::event status{step.write_buffer(root, 0,
	                                   nablaH.size() * sizeof(vec2d),
	                                   nablaH.data())};

#ifdef PROFILE
	int it = 0;
#endif
	double e_max;
	do {
		e_max = 0.;

		// (re-)calculate gradient average in coarser levels
		status = evalQuadtree(root, nablaH.get_height(),
		                      nablaH.size(), nablaH.data(), {status});

		const morton tmp[16]{0, 1, 2, 3, 4, 5, 6, 7, 8,
		                     9, 10, 11, 12, 13, 14, 15};
		//const morton tmp[4]{0, 1, 2, 3};
		cl::buffer parent_index{simd.create_buffer(CL_MEM_READ_WRITE |
		                                           CL_MEM_COPY_HOST_PTR,
		                                           16 * sizeof(morton),
		                                           tmp)};
		unsigned n_a_i = 64;
		for (unsigned i = 2; i < nablaH.get_height() - 1 && n_a_i; ++i) {
			const unsigned offseti = (4 * pow(4, i) - 1) / 3,
		                       sizei = ((4 * pow(4, i + 1) - 1) / 3) -
			                       offseti;
			unsigned aligni;
			const cl::buffer rooti{root.create_sub<vec2d>(0,
			                                              offseti,
			                                              sizei,
			                                              maba,
			                                              &aligni)},
			                 flagi{simd.create_buffer(CL_MEM_READ_WRITE,
			                                          n_a_i * sizeof(cl_uchar))},
			                 zsi{simd.create_buffer(CL_MEM_READ_WRITE,
			                                        n_a_i * sizeof(morton))};
			exe["hie_calc_err"].set_args(rooti, aligni, sizei, parent_index, zsi,
			                              flagi, eps, n_a_i);
			status = step.ndrange_kernel(exe["hie_calc_err"],
			                             {}, {n_a_i},
			                             {wgs}, {status});
#ifdef PROFILE
			accumulate(status);
#endif

			const cl::buffer offsets{simd.create_buffer(CL_MEM_READ_WRITE,
			                                            n_a_i * sizeof(cl_uint))};
			exe["mutate"].set_args(flagi, offsets, n_a_i);
			status = step.ndrange_kernel(exe["mutate"],
			                             {}, {n_a_i},
			                             {wgs}, {status});
			status = scan("prescan", offsets, n_a_i, {status});
			unsigned n_a_j;
			status = step.read_buffer(offsets, (n_a_i - 1) * sizeof(cl_uint), sizeof(cl_uint),
	                                          &n_a_j, {status});
			unsigned char add = 0;
			status = step.read_buffer(flagi, (n_a_i - 1) * sizeof(cl_uchar), sizeof(cl_uchar),
	                                          &add, {status});
			n_a_j += add;
			if (n_a_j) {
				// extract Morton codes of the parent indices for active nodes in the lower level
				parent_index = simd.create_buffer(CL_MEM_READ_WRITE,
				                                  n_a_j * sizeof(morton));
				exe["hie_tag_active"].set_args(flagi, offsets, zsi,
				                               parent_index, n_a_i);
#ifdef PROFILE
				accumulate(status);
#endif
				status = step.ndrange_kernel(exe["hie_tag_active"],
				                             {}, {n_a_i}, {wgs},
				                             {status});
#ifdef PROFILE
				accumulate(status);
#endif
			}
			n_a_i = 4 * n_a_j;
		}

		if (!n_a_i)
			break;

		const cl::buffer err{simd.create_buffer(CL_MEM_READ_WRITE,
		                                        n_a_i * sizeof(cl_double))};
		unsigned align0;
		const cl::buffer root0{root.create_sub<vec2d>(0,
		                                              nablaH.size() - nablaH.len(),
		                                              nablaH.len(),
		                                              maba,
		                                              &align0)};

		for (unsigned char mode = 0; mode < 3; ++mode) {
			exe["hie_correct_grad"].set_args(root0, align0, parent_index,
			                                 err, eps, s.GetDouble(),
			                                 nablaH.len(), mode, n_a_i);
			status = step.ndrange_kernel(exe["hie_correct_grad"],
			                             {}, {n_a_i}, {wgs},
			                             {status});
#ifdef PROFILE
			accumulate(status);
#endif
		}


		status = reduce("reduce", err, n_a_i, e_max,
		                {status});

		std::cerr << "e_max: " << e_max << std::endl;
#ifdef PROFILE
		++it;
#endif
	} while (e_max > eps);

	status = step.read_buffer(root, 0, nablaH.size() * sizeof(vec2d),
	                          nablaH.data(), {status});

#ifdef PROFILE
	std::cerr << "PROFILE: proccessing time: " << acc / 1e9 << " [s]" << std::endl;
#endif
#ifdef PROFILE
	std::cerr << "iteration count: " << it << std::endl;
#endif
}

//==============================================================================
// naive version
//______________________________________________________________________________
void TMOCadik08::correctGradCyc(std::vector<vec2d>& g, const unsigned rows,
                                const unsigned cols, const double eps) const
{
	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols * sizeof(vec2d))},
	                 err{simd.create_buffer(CL_MEM_READ_WRITE,
	                                        rows * cols * sizeof(double))};
	cl::event status{step.write_buffer(grad, 0, rows * cols * sizeof(vec2d),
	                                   g.data())};

#ifdef PROFILE
	int it = 0;
#endif
	double e_max;
	const double s = this->s.GetDouble();
	do {
		e_max = 0.;

		exe["cyc_calc_error"].set_args(grad, err, rows, cols);
		status = step.ndrange_kernel(exe["cyc_calc_error"], {},
		                             {rows, cols}, {dim, dim},
		                             {status});
#ifdef PROFILE
		accumulate(status);
#endif

		for (unsigned n = 0; n < 4; ++n) {
			exe["cyc_correct_grad"].set_args(grad, err, s, rows, cols, n);
			status = step.ndrange_kernel(exe["cyc_correct_grad"], {},
			                             {rows, cols}, {dim, dim},
			                             {status});
#ifdef PROFILE
			accumulate(status);
#endif
		}

		status = reduce("reduce", err, rows * cols, e_max,
		                {status});

		std::cerr << "e_max: " << e_max << std::endl;
#ifdef PROFILE
		++it;
#endif
	} while (e_max > eps);

	status = step.read_buffer(grad, 0, rows * cols * sizeof(vec2d),
	                          g.data(), {status});
	step.wait({status});

#ifdef PROFILE
	std::cerr << "PROFILE: proccessing time: " << acc / 1e9 << " [s]" << std::endl;
	std::cerr << "iteration count: " << it << std::endl;
#endif
}

#ifdef PROFILE
void TMOCadik08::accumulate(const cl::event& status) const
{
	cl_ulong start, end;
	step.wait({status});
	status.profiling_info(CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start);
	status.profiling_info(CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end);
	acc += end - start;
}
#endif

//==============================================================================
// double integration
//______________________________________________________________________________
void TMOCadik08::integrate2x(TMOImage& g, TMOImage& o) const
{
	const unsigned rows = o.GetHeight(),
	               cols = o.GetWidth();

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

	const cl::buffer grad{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         rows * cols * 3 *
	                                         sizeof(double))},
	                 out{simd.create_buffer(CL_MEM_READ_WRITE |
	                                        CL_MEM_COPY_HOST_PTR,
	                                        rows * cols * 3 * sizeof(double),
	                                        o.GetData())};
	cl::event status{step.write_buffer(grad, 0, rows * cols *
	                                   3 * sizeof(double),
	                                   g.GetData())};

	// for each row, perform integration in +x direction
	exe["integrate2x"].set_args(grad, out, rows, cols);
	status = step.ndrange_kernel(exe["integrate2x"], {},
	                             {rows}, {wgs}, {status});
#ifdef PROFILE
	accumulate(status);
#endif

	status = step.read_buffer(out, 0, rows * cols * 3 * sizeof(double),
	                          o.GetData(), {status});

	step.wait({status});
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
	double A = 0, B = 0;

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
			pDst_image[3 * (tmp_y + j)]= .01 * (A + B * pDst_image[3 * (tmp_y + j)]);
			pDst_image[3 * (tmp_y + j) + 1]= .01 * (A + B * pDst_image[3 * (tmp_y + j) + 1]);
			pDst_image[3 * (tmp_y + j) + 2]= .01 * (A + B * pDst_image[3 * (tmp_y + j) + 2]);

			pSrc_image[3 * (tmp_y + j)] *= .01;
			pSrc_image[3 * (tmp_y + j) + 1] = pSrc_image[3 * (tmp_y + j)];
			pSrc_image[3 * (tmp_y + j) + 2] = pSrc_image[3 * (tmp_y + j)];
		}
	}
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
#ifdef PROFILE
	accumulate(status);
#endif
	if (m > 1)
		status = reduce(type, tmp, m, out, {status});
	else
		status = step.read_buffer(tmp, 0, sizeof(cl_double),
		                          &out, {status});

	return status;
}

//______________________________________________________________________________
cl::event TMOCadik08::scan(const std::string type, const cl::buffer& in,
                           const unsigned n, const cl::event_list pending) const
{
	const unsigned m = std::ceil(static_cast<float>(com::math::ceil2mul(n,
	                                                2 * wgs)) / (2 * wgs));

	const cl::buffer sums{simd.create_buffer(CL_MEM_READ_WRITE,
	                                         m * sizeof(cl_uint))};

	exe[type].set_args(in, sums, cl::local_mem{2 * wgs * sizeof(cl_uint)},
	                   n);
	cl::event status{step.ndrange_kernel(exe[type], {0}, {m * wgs}, {wgs},
	                                     pending)};
#ifdef PROFILE
	accumulate(status);
#endif

	if (m > 1) {
		status = scan(type, sums, m, {status});
		exe["add"].set_args(in, sums, n);
		status = step.ndrange_kernel(exe["add"], {2 * wgs}, {n},
		                             {wgs}, {status});
#ifdef PROFILE
		accumulate(status);
#endif
	}

	return status;
}


/* --------------------------------------------------------------------------- *
 * Gradient inconsistency correction -- CPU                                    *
 * --------------------------------------------------------------------------- */
void TMOCadik08::inconsistencyCorrection(std::vector<vec2d>& grad,
                                         const long ymax, const long xmax,
                                         const double eps)
{
#ifdef PROFILE
	int it = 0;
#endif
	double e_max;
	do {
#ifdef PROFILE
		++it;
#endif
		e_max = 0.;
		for (int i = 0; i < ymax; ++i) {
			int tmp_y = i * xmax;
			for (int j = 0; j < xmax; ++j) {
				int tmp_ind = j + tmp_y;
				double e = grad[tmp_ind].x - //Gx
				           grad[tmp_ind].y + //Gy
				           ((j + 1< xmax) ? (grad[tmp_ind + 1].y) : 0.) - //Gy(i+1,j)
				           ((i + 1< ymax) ? (grad[tmp_ind + xmax].x) : 0.); //Gx(i,j+1)

				if (fabs(e) > e_max)
					e_max = fabs(e);

				e *= .25 * s;

				grad[tmp_ind].x -= e; //Gx
				grad[tmp_ind].y += e; //Gy
				if (j + 1 < xmax)
					grad[tmp_ind + 1].y -= e; //Gy(i+1,j)
				if (i + 1 < ymax)
					grad[tmp_ind + xmax].x += e; //Gx(i,j+1)
			}
		}
		std::cerr << "e_max: " << e_max << std::endl;
	} while (e_max > eps);
#ifdef PROFILE
	std::cerr << "iteration count: " << it << std::endl;
#endif
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
			//neboli: OUTPUT_BW[0][y] = OUTPUT_BW[0][y - 1] + Grad_Y[0][y - 1];

		for (j = 1; j < xmax; ++j) {
			tmp_ind = j + tmp_y;
			pDst_image[3 * tmp_ind] = pDst_image[3 * tmp_ind + 1] =
			pDst_image[3 * tmp_ind + 2] = (pDst_image[3 * (tmp_ind - 1)] +
			                              pG_image[3 * (tmp_ind - 1)]);
			//neboli: OUTPUT_BW[x][y] = OUTPUT_BW[x - 1][y] + Grad_X[x - 1][y]; 
		}
	}
}
