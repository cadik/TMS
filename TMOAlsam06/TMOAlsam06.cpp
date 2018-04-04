/* --------------------------------------------------------------------------- *
 * TMOAlsam06.cpp: implementation of the Alsam A., Kosvas O.         		   *
 * Grey colour sharpening												   	   *
 * --------------------------------------------------------------------------- */

#include "TMOAlsam06.h"

#include <algorithm>
#include <math.h> 
#include <iostream>

using namespace Eigen;

TMOAlsam06::TMOAlsam06()
{
	SetName(L"Alsam06");								
	SetDescription(L"Generates greyscale image with color separation and texture enhancement.");
}

TMOAlsam06::~TMOAlsam06()
{
}

int TMOAlsam06::Transform()
{
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			

	int i;
	int pixelCount = pSrc->GetHeight() * pSrc->GetWidth();

	/**
	 *                      R: G: B:
	 *             pixel1 | R1 G1 B1 |
	 * Matrix p =  pixel2 | R2 G2 B3 |
	 *             pixel3 | R3 G3 B3 |
	 */
	MatrixXd p(pixelCount, 3);

	for (i = 0; i < pixelCount; i++) {
		p(i, 0) = *pSourceData++;
		p(i, 1) = *pSourceData++;
		p(i, 2) = *pSourceData++;
	}

	/**
	 * Matrix pc = p^T * p (correlation matrix)
	 */
	MatrixXd pc(3,3);
	pc = p.transpose() * p;

	EigenSolver<MatrixXd> es(pc);
	
	/**
	 * Matrix u = eigen vectors of matrix pc.
	 * pc = u * d * u^T
	 * Matrix d = diagonal matrix with eigen values of pc.
	 */
	MatrixXd u(3,3);
	u = es.pseudoEigenvectors();

	/**
	 * Matrix pu = u1 * u1^T
	 * projection operator of most significant vector u1
	 */
	MatrixXd pu(3,3);
	pu = u.col(0) * u.col(0).transpose();

	double g;

	for (i = 0; i < pixelCount; i++) {
		pSrc->ProgressBar(i, pixelCount);

		/**
		 * Matrix ppu = projection of each pixel onto u1
		 */
		MatrixXd ppu(1, 3);
		ppu = p.row(i) * pu;

		g = ppu.squaredNorm();

		*pDestinationData++ = g;
		*pDestinationData++ = g;
		*pDestinationData++ = g;
	}

	pSrc->ProgressBar(i, pixelCount);
	return 0;
}

