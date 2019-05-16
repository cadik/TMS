#include <fstream>

#include "constructQpMatrices.h"
#include "osqpOptimization.h"

std::vector<cv::Mat> optimizationWithOsqp(cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels)
{
	std::cout << "QP optimization with OSQP" << '\n';
	int w1 = detailImage.cols;
	int h1 = detailImage.rows;
	std::vector<cv::Mat> sAndT;
	cv::Mat s;
	s = cv::Mat::ones(h1, w1, CV_32F);
	cv::Mat t;
	t = cv::Mat::zeros(h1, w1, CV_32F);
	
	double r1 = 200, r2 = 500;
	int mainSize = h1*w1*2; // number of qp variables (n)
	// std::cout << "n.o. qp variables (mainSize) is: " << mainSize << '\n';

	std::vector<triplet_t> P_triplets = getHessianTriplets(detailImage, weight1, weight2, r1, r2);

	// arrays for creating sparse Hessian matrix (for OSQP it is 'P')
	std::vector<c_int> P_i, P_p;
	std::vector<c_float> P_x;
	// int res = convert_collapse(P_triplets, P_i, P_p, P_x);
	int res = convert(P_triplets, P_i, P_p, P_x);
	// std::cout << "\nrow indices:" << '\n';
	// for(auto i : P_i) std::cout << i << ' ';
	// std::cout << "\ncolumn first indices:" << '\n';
	// for(auto p : P_p) std::cout << p << ' ';
	// std::cout << "\nvalues:" << '\n';
	// for(auto x : P_x) std::cout << x << ' ';
	// std::cout << std::endl;
	c_int P_nnz = P_x.size();
	std::cout << "Hessian P_triplets size: " << P_triplets.size() << '\n';
	std::cout << "Hessian NNZ elements count: " << P_nnz << '\n';
	P_triplets.clear();

	// FIXME allocate these arrays dynamically on heap? (with using automatic memory management (unique_ptr, ...))
	c_float *q = (c_float *)c_malloc(mainSize * sizeof(c_float)); // dense array for linear part of cost function (size n)
	for (int i = 0; i < mainSize; i++) { q[i] = 1.0f; }
	c_float *l = (c_float *)c_malloc(h1*w1*3 * sizeof(c_float)); // lower bound A (for each rgb color of pixel one lower bound)
	c_float *u = (c_float *)c_malloc(h1*w1*3 * sizeof(c_float)); // upper bound A (for each rgb color of pixel one upper bound)

	// arrays for creating A matrix (constraint) as sparse
	std::vector<c_int> A_i, A_p;
	std::vector<c_float> A_x;
	std::vector<triplet_t> A_triplets;

	// create also A matrix as sparse

	int counterA = 0; // constraints counter
	for (int j = 0; j < h1; j++) {
		for (int i = 0; i < w1; i++) {
			/*
				Setting A (constraints)
			*/
			/*
				Setting constraint -B_i <= t_i + s_iD_i <= 1 - B_i
				not done at borders, because for each constraint there is only one lower/upper bound
			*/
			A_triplets.push_back({i + j*w1, counterA, (double)(detailChannels[0].at<float>(j, i)/255.0)});
			A_triplets.push_back({i + j*w1 + w1*h1, counterA, 1.0});
			l[j*w1 + i] = -(double)(baseChannels[0].at<float>(j, i)/255.0); // lower bound A
			u[j*w1 + i] = 1.0 - (double)(baseChannels[0].at<float>(j, i)/255.0); // upper bound A
			counterA++;
			A_triplets.push_back({i + j*w1, counterA, (double)(detailChannels[1].at<float>(j, i)/255.0)});
			A_triplets.push_back({i + j*w1 + w1*h1, counterA, 1.0});
			l[1*h1*w1 + j*w1 + i] = -(double)(baseChannels[1].at<float>(j, i)/255.0); // lower bound A
			u[1*h1*w1 + j*w1 + i] = 1.0 - (double)(baseChannels[1].at<float>(j, i)/255.0); // upper bound A
			counterA++;
			A_triplets.push_back({i + j*w1, counterA, (double)(detailChannels[2].at<float>(j, i)/255.0)});
			A_triplets.push_back({i + j*w1 + w1*h1, counterA, 1.0});
			l[2*h1*w1 + j*w1 + i] = -(double)(baseChannels[2].at<float>(j, i)/255.0); // lower bound A
			u[2*h1*w1 + j*w1 + i] = 1.0 - (double)(baseChannels[2].at<float>(j, i)/255.0); // upper bound A
			counterA++;
		}
	}

	res = convert(A_triplets, A_i, A_p, A_x);
	// for(auto i : A_i) std::cout << i << ' ';
	// std::cout << std::endl;
	// for(auto p : A_p) std::cout << p << ' ';
	// std::cout << std::endl;
	// for(auto x : A_x) std::cout << x << ' ';
	// std::cout << std::endl;
	c_int A_nnz = A_x.size();
	A_triplets.clear();

	// solve quadratic programming problem
	// Problem settings
	std::cout << "after creating A" << '\n';
	OSQPSettings * settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

	// Structures
	OSQPWorkspace * work;  // Workspace
	OSQPData * data;  // OSQPData

	std::cout << "before populating data" << '\n';
	// Populate data
	data = (OSQPData *)c_malloc(sizeof(OSQPData));
	data->n = mainSize;		// number of variables
	data->m = h1*w1*3;		// number of constraints (one for each color of pixel)
	data->P = csc_matrix(data->n, data->n, P_nnz, P_x.data(), P_i.data(), P_p.data());
	data->q = q;			// dense array for linear part of cost function (size n)
	data->A = csc_matrix(data->m, data->n, A_nnz, A_x.data(), A_i.data(), A_p.data());
	data->l = l;
	data->u = u;
	std::cout << "after populating data" << '\n';

// void dump_vec(c_float *v, c_int len, const char *file_name);
	dump_vec(q, mainSize, "q-vector.txt");
	dump_vec(l, h1*w1*3, "l-vector.txt");
	dump_vec(u, h1*w1*3, "u-vector.txt");
	dump_csc_matrix(data->P, "P-matrix.txt");
	dump_csc_matrix(data->A, "A-matrix.txt");
	// print_csc_matrix(data->P, "P-matrix");
	// print_csc_matrix(data->A, "A-matrix");

	std::cout << "before setting default settings" << '\n';
	// Define Solver settings as default
	osqp_set_default_settings(settings);
	// settings->max_iter = 10000;

	std::cout << "before setting workspace" << '\n';
	// Setup workspace
	work = osqp_setup(data, settings);
	// if some failure occured, return an empty vector
	if(work == OSQP_NULL)
		return sAndT;

	std::cout << "before solving" << '\n';
	// Solve Problem
	osqp_solve(work);
	std::cout << "after solving" << '\n';

	/*
		Put the results into matrices
	*/
	int lineCounter = 0;
	int parametersCounterS = 0;
	int parametersCounterT = 0;
	for (int i = 0; i < mainSize; i++) {
			/*
				s parameters are the first, t are after them
			*/
			if (lineCounter < w1*h1) {
				s.at<float>(parametersCounterS / w1, parametersCounterS % w1) = work->solution->x[parametersCounterS];
				parametersCounterS++;				 			    
			} else {				
				t.at<float>(parametersCounterT / w1, parametersCounterT % w1) = work->solution->x[w1*h1 + parametersCounterT];
				parametersCounterT++; 
			}
			lineCounter++;
	}


	// Clean workspace
	osqp_cleanup(work);
	c_free(data->A);
	c_free(data->P);
	c_free(data);
	c_free(q);
	c_free(l);
	c_free(u);
	c_free(settings);

	sAndT.push_back(s);
	sAndT.push_back(t);
	return sAndT;
}
