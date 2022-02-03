/* -------------------------------------------------------------------------- *
 * constructQpMatrices.h: constructing Hessian sparse matrix for QP           *
 * Authors: Tomas Hudziec (2019)                                              *
 *          - triplet and column_major_sort structs                           *
 *          - convert functions                                               *
 *          Pavel Sedlar (2018)                                               *
 *          - getHessianTriplets function                                     *
 * -------------------------------------------------------------------------- */
#ifndef CONSTRUCT_QP_MATRICES
#define CONSTRUCT_QP_MATRICES

#include <opencv2/opencv.hpp>
#include <algorithm> // sort
#include <vector>
#include <osqp.h>

// triplet x, y, value for non-zero elements of sparse matrix
typedef struct triplet
{
	int x;
	int y;
	double v;
	friend std::ostream &operator<<(std::ostream &stream, const triplet &t)
	{
		stream << '(' << t.x << ',' << t.y << ") = " << t.v;
		return stream;
	}
	// comparing, if indices are the same
	inline bool operator==(const triplet &other)
	{
		return x == other.x && y == other.y;
	}
} triplet_t;

struct column_major_sort
{
	inline bool operator()(const triplet_t &t1, const triplet_t &t2)
	{
		return (t1.x < t2.x) || ((t1.x == t2.x) && (t1.y < t2.y));
	}
};

int convert(std::vector<triplet_t> &triplets, std::vector<c_int> &ir, std::vector<c_int> &jc, std::vector<c_float> &val);
int convert_collapse(std::vector<triplet_t> &triplets, std::vector<c_int> &ir, std::vector<c_int> &jc, std::vector<c_float> &val);

std::vector<triplet_t> getHessianTriplets(cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, double r1, double r2);

#endif /* end of include guard: CONSTRUCT_QP_MATRICES */
