#include "constructQpMatrices.h"

// converts non-zero sparse matrix elements in form of triplets <x,y,value>
// to Compressed Column (CC) Sparse Matrix File Format:
// ir[] - row indices (its size is tripletsSize)
// val[] - values (its size is tripletsSize)
// jc[] - indices to first entry of columns (its size is number of columns + 1)
int convert(std::vector<triplet_t> &triplets, std::vector<c_int> &ir, std::vector<c_int> &jc, std::vector<c_float> &val) {
	// sort triplets by rows and columns - in column-major order
	std::sort(triplets.begin(), triplets.end(), column_major_sort());
	
	// clear vectors ir, jc, val
	ir.clear(); jc.clear(); val.clear();

	// first entry of column is always first element
	size_t tripletsSize = triplets.size();
	jc.push_back(0);
	for (std::vector<triplet_t>::size_type i = 0; i < tripletsSize - 1; i++) {
		ir.push_back(triplets[i].y);
		val.push_back(triplets[i].v);
		// push index of first item in column to jc vector
		if (triplets[i].x != triplets[i+1].x) {
			jc.push_back(c_int(i+1));
		}
	}
	// push last elements to ir and val
	ir.push_back(triplets[tripletsSize-1].y);
	val.push_back(triplets[tripletsSize-1].v);
	// push last index + 1 to jc TODO comment better
	jc.push_back(c_int(tripletsSize));
	
	return 0;
}

// converts non-zero sparse matrix elements in form of triplets <x,y,value>
// to Compressed Column (CC) Sparse Matrix File Format:
// ir[] - row indices (its size is tripletsSize)
// val[] - values (its size is tripletsSize)
// jc[] - indices to first entry of columns (its size is number of columns + 1)
// and collapses (sums) duplicates
int convert_collapse(std::vector<triplet_t> &triplets, std::vector<c_int> &ir, std::vector<c_int> &jc, std::vector<c_float> &val) {
	// sort triplets by rows and columns - in column-major order
	std::sort(triplets.begin(), triplets.end(), column_major_sort());
	
	// clear vectors ir, jc, val
	ir.clear(); jc.clear(); val.clear();

	// first entry of column is always first element
	size_t tripletsSize = triplets.size();
	int duplicates = 0;
	jc.push_back(0);
	for (std::vector<triplet_t>::size_type i = 0; i < tripletsSize - 1; i++) {
		c_float value = triplets[i].v; // TODO templating? (make sum possible to be of any type)
		// summing duplicates, they are one after other, since triplets vector is sorted
		while (triplets[i] == triplets[i+1] && i < tripletsSize - 1) {
			value += triplets[i+1].v;
			i++;
			duplicates++;
		}
		ir.push_back(triplets[i].y);
		val.push_back(value);
		// push index of first item in column to jc vector
		if (triplets[i].x != triplets[i+1].x) {
			jc.push_back(c_int(int(i)+1 - duplicates));
		}
	}
	// add last item, if it is not duplicate
	if (!(triplets[tripletsSize-1] == triplets[tripletsSize-2])) {
		// push last elements to ir and val
		ir.push_back(triplets[tripletsSize-1].y);
		val.push_back(triplets[tripletsSize-1].v);
		// push last index + 1 to jc TODO comment better
		jc.push_back(c_int(tripletsSize));
	}
	
	return 0;
}


std::vector<triplet_t> getHessianTriplets(cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, double r1, double r2)
{
	int h1 = detailImage.rows;
	int w1 = detailImage.cols;
	
	// array of triplets of nonzero elemets of sparse matrix
	std::vector<triplet_t> triplets;

	/*
		Setting s
	*/
	for (int j = 0; j < h1; j++) {
		for (int i = 0; i < w1; i++) {
			
			double sPowerTerm1 = -((double)(pow(detailImage.at<float>(j, i)/255.0, 2)));
			//double sPowerTerm2 = r1*4*weight1.at<float>(j, i);
			double sPowerTerm2 = 0.0;
			double sPowerTerm3 = 0.0;
			double sPowerTerm4 = 0.0;
			double sPowerTerm5 = 0.0;
			double sPowerTerm6 = 0.0;
			double sPowerTerm7 = 0.0;
			double sPowerTerm8 = 0.0;
			double sPowerTerm9 = 0.0;
			double sPowerTerm10 = 0.0;
			double sPowerTerm11 = 0.0;
			double sPowerTerm12 = 0.0;
			double sPowerTerm13 = 0.0;
			double sumOfPowerTerms = 0.0;

			if (i > 1) {
				sPowerTerm2 = 2*r1*(double)weight1.at<float>(j, i - 1);
				triplets.push_back({w1 * j + i, w1 * j + (i - 2), -sPowerTerm2});
				triplets.push_back({w1 * j + (i - 2), w1 * j + i, -sPowerTerm2});
			}
			if (j < (h1 - 2)) {
				sPowerTerm3 = 2*r1*(double)weight1.at<float>(j + 1, i);
				triplets.push_back({w1 * (j + 2) + i, w1 * j + i, -sPowerTerm3});
				triplets.push_back({w1 * j + i, w1 * (j + 2) + i, -sPowerTerm3});
			}

			if (i < (w1 - 2)) {
				sPowerTerm4 = 2*r1*(double)weight1.at<float>(j, i + 1);
				triplets.push_back({w1 * j + (i + 2), w1 * j + i, -sPowerTerm4});
				triplets.push_back({w1 * j + i, w1 * j + (i + 2), -sPowerTerm4});
			}

			if (j > 1) {
				sPowerTerm5 = 2*r1*(double)weight1.at<float>(j - 1, i);
				triplets.push_back({w1 * j + i, w1 * (j  - 2) + i, -sPowerTerm5});
				triplets.push_back({w1 * (j  - 2) + i, w1 * j + i, -sPowerTerm5});
			}
			
			if (i == 0) {
				sPowerTerm6 = 2*r1*(double)weight1.at<float>(j, i);
				triplets.push_back({w1 * j + i, w1 * j + i + 1, -sPowerTerm6});
				triplets.push_back({w1 * j + i + 1, w1 * j + i, -sPowerTerm6});
			}

			if (i == 1) {
				sPowerTerm7 = 2*r1*(double)weight1.at<float>(j, i - 1);
				triplets.push_back({w1 * j + i, w1 * j + i - 1, -sPowerTerm7});
				triplets.push_back({w1 * j + i - 1, w1 * j + i, -sPowerTerm7});
			}

			if (j == 0) {
				sPowerTerm8 = 2*r1*(double)weight1.at<float>(j, i);
				triplets.push_back({w1 * j + i, w1 * (j + 1) + i, -sPowerTerm8});
				triplets.push_back({w1 * (j + 1) + i, w1 * j + i, -sPowerTerm8});
			}

			if (j == 1) {
				sPowerTerm9 = 2*r1*(double)weight1.at<float>(j - 1, i);
				triplets.push_back({w1 * j + i, w1 * (j - 1) + i, -sPowerTerm9});
				triplets.push_back({w1 * (j - 1) + i, w1 * j + i, -sPowerTerm9});
			}

			if (i == w1 - 1) {
				sPowerTerm10 = 2*r1*(double)weight1.at<float>(j, i);
				triplets.push_back({w1 * j + i, w1 * j + i - 1, -sPowerTerm10});
				triplets.push_back({w1 * j + i - 1, w1 * j + i, -sPowerTerm10});
			}

			if (i == w1 - 2) {
				sPowerTerm11 = 2*r1*(double)weight1.at<float>(j, i + 1);
				triplets.push_back({w1 * j + i, w1 * j + i + 1, -sPowerTerm11});
				triplets.push_back({w1 * j + i + 1, w1 * j + i, -sPowerTerm11});
			}

			if (j == h1 - 1) {
				sPowerTerm12 = 2*r1*(double)weight1.at<float>(j, i);
				triplets.push_back({w1 * j + i, w1 * (j - 1) + i, -sPowerTerm12});
				triplets.push_back({w1 * (j - 1) + i, w1 * j + i, -sPowerTerm12});
			}

			if (j == h1 - 2) {
				sPowerTerm13 = 2*r1*(double)weight1.at<float>(j + 1, i);
				triplets.push_back({w1 * j + i, w1 * (j + 1) + i, -sPowerTerm13});
				triplets.push_back({w1 * (j + 1) + i, w1 * j + i, -sPowerTerm13});
			}

			sumOfPowerTerms = sPowerTerm1 +
							sPowerTerm2 +
							sPowerTerm3 +
							sPowerTerm4 +
							sPowerTerm5 +
							sPowerTerm6 + 
							sPowerTerm7 + 
							sPowerTerm8 + 
							sPowerTerm9 +
							sPowerTerm10 +
							sPowerTerm11 +
							sPowerTerm12 +
							sPowerTerm13; 
			
			triplets.push_back({j * w1 + i, j * w1 + i, sumOfPowerTerms});
		}
	}
	

	/*
		Setting t
	*/
	for (int j = 0; j < h1; j++) {
		for (int i = 0; i < w1; i++) {
			//double sPowerTerm01 = 4*r2*(double)weight2.at<float>(j, i);
			double sPowerTerm01 = 0.0;
			double sPowerTerm02 = 0.0;
			double sPowerTerm03 = 0.0;
			double sPowerTerm04 = 0.0;
			double sPowerTerm05 = 0.0;
			double sPowerTerm06 = 0.0;
			double sPowerTerm07 = 0.0;
			double sPowerTerm08 = 0.0;
			double sPowerTerm09 = 0.0;
			double sPowerTerm010 = 0.0;
			double sPowerTerm011 = 0.0;
			double sPowerTerm012 = 0.0;
			double sumOfPowerTerms1 = 0.0;

			if (i > 1) {
				sPowerTerm01 = 2*r2*(double)weight2.at<float>(j, i - 1); //check
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + h1 * j + (i - 2), -sPowerTerm01});
				triplets.push_back({w1*h1 + w1 * j + (i - 2), w1*h1 + w1 * j + i, -sPowerTerm01});
			}
			if (j < (h1 - 2)) {
				sPowerTerm02 = 2*r2*(double)weight2.at<float>(j + 1, i);
				triplets.push_back({w1*h1 + w1 * (j + 2) + i, w1*h1 + w1 * j + i, -sPowerTerm02});
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * (j + 2) + i, -sPowerTerm02});
			}

			if (i < (w1 - 2)) {
				sPowerTerm03 = 2*r2*(double)weight2.at<float>(j, i + 1);
				triplets.push_back({w1*h1 + w1 * j + (i + 2), w1*h1 + w1 * j + i, -sPowerTerm03});
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * j + (i + 2), -sPowerTerm03});
			}

			if (j > 1) {
				sPowerTerm04 = 2*r2*(double)weight2.at<float>(j - 1, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * (j  - 2) + i, -sPowerTerm04});
				triplets.push_back({w1*h1 + w1 * (j  - 2) + i, w1*h1 + w1 * j + i, -sPowerTerm04});
			}

			if (i == 0) {
				sPowerTerm05 = 2*r2*(double)weight2.at<float>(j, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * j + i + 1, -sPowerTerm05});
				triplets.push_back({w1*h1 + w1 * j + i + 1, w1*h1 + w1 * j + i, -sPowerTerm05});
			}

			if (i == 1) {
				sPowerTerm06 = 2*r2*(double)weight2.at<float>(j, i - 1);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * j + i - 1, -sPowerTerm06});
				triplets.push_back({w1*h1 + w1 * j + i - 1, w1*h1 + w1 * j + i, -sPowerTerm06});
			}
		
			if (j == 0) {
				sPowerTerm07 = 2*r2*(double)weight2.at<float>(j, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * (j + 1) + i, -sPowerTerm07});
				triplets.push_back({w1*h1 + w1 * (j + 1) + i, w1*h1 + w1 * j + i, -sPowerTerm07});
			}

			if (j == 1) {
				sPowerTerm08 = 2*r2*(double)weight2.at<float>(j - 1, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * (j - 1) + i, -sPowerTerm08});
				triplets.push_back({w1*h1 + w1 * (j - 1) + i, w1*h1 + w1 * j + i, -sPowerTerm08});
			}

			if (i == w1 - 1) {
				sPowerTerm09 = 2*r2*(double)weight2.at<float>(j, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * j + i - 1, -sPowerTerm09});
				triplets.push_back({w1*h1 + w1 * j + i - 1, w1*h1 + w1 * j + i, -sPowerTerm09});
			}

			if (i == w1 - 2) {
				sPowerTerm010 = 2*r2*(double)weight2.at<float>(j, i + 1);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * j + i + 1, -sPowerTerm010});
				triplets.push_back({w1*h1 + w1 * j + i + 1, w1*h1 + w1 * j + i, -sPowerTerm010});
			}

			if (j == h1 - 1) {
				sPowerTerm011 = 2*r2*(double)weight2.at<float>(j, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * (j - 1) + i, -sPowerTerm011});
				triplets.push_back({w1*h1 + w1 * (j - 1) + i, w1*h1 + w1 * j + i, -sPowerTerm011});
			}

			if (j == h1 - 2) {
				sPowerTerm012 = 2*r2*(double)weight2.at<float>(j + 1, i);
				triplets.push_back({w1*h1 + w1 * j + i, w1*h1 + w1 * (j + 1) + i, -sPowerTerm012});
				triplets.push_back({w1*h1 + w1 * (j + 1) + i, w1*h1 + w1 * j + i, -sPowerTerm012});
			}

			sumOfPowerTerms1 = sPowerTerm01 +
							sPowerTerm02 +
							sPowerTerm03 +
							sPowerTerm04 +
							sPowerTerm05 +
							sPowerTerm06 +
							sPowerTerm07 +
							sPowerTerm08 +
							sPowerTerm09 +
							sPowerTerm010 +
							sPowerTerm011 +
							sPowerTerm012;
			triplets.push_back({w1*h1 + j * w1 + i, w1*h1 + j * w1 + i, sumOfPowerTerms1});

		}
	}
	// double values; original code had it, but without it it produces little bit more vivid output
	// std::for_each(triplets.begin(), triplets.end(), [](triplet_t &t){ t.v*=2; });

	return triplets;
}
