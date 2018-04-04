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

	kernelSizeParameter.SetName(L"kernel");
	kernelSizeParameter.SetDescription(L"Kernel size of Gaussian filter applied to calculate high frequency masks.");
	kernelSizeParameter.SetDefault(15);
	kernelSizeParameter=15;
	kernelSizeParameter.SetRange(3, 50);
	this->Register(kernelSizeParameter);

	standardDeviationParameter.SetName(L"deviation");
	standardDeviationParameter.SetDescription(L"Standard deviation of Gaussian filter applied to calculate high frequency masks.");
	standardDeviationParameter.SetDefault(1);
	standardDeviationParameter=1;
	standardDeviationParameter.SetRange(1, 10);
	this->Register(standardDeviationParameter);
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
	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();

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
	MatrixXd gu(width, height);

	for (i = 0; i < pixelCount; i++) {
		//pSrc->ProgressBar(i, pixelCount);

		/**
		 * Matrix ppu = projection of each pixel onto u1
		 */
		MatrixXd ppu(1, 3);
		ppu = p.row(i) * pu;

		g = ppu.squaredNorm();

		/**pDestinationData++ = g;
		*pDestinationData++ = g;
		*pDestinationData++ = g;*/
		gu(floor(i / width), i % height) = g;
	}

	//pSrc->ProgressBar(i, pixelCount);
	MatrixXd gauss(kernelSizeParameter.GetInt(), kernelSizeParameter.GetInt());
	double sum = 0;
	double mean = kernelSizeParameter.GetInt() / 2;

	for (int x = 0; x < kernelSizeParameter.GetInt(); x++) {
		for (int y = 0; y < kernelSizeParameter.GetInt(); y++) {
			//gauss(x,y) = exp(-0.5 * (pow(x, 2) + pow(y, 2)) / pow(standardDeviationParameter, 2)) / (2 * M_PI * standardDeviationParameter * standardDeviationParameter);
			gauss(x,y) = exp(-0.5 * (pow((x - mean) / standardDeviationParameter.GetInt(), 2.0) + pow((y - mean) / standardDeviationParameter.GetInt(), 2.0))) / (2 * M_PI * standardDeviationParameter.GetInt() * standardDeviationParameter.GetInt());
			sum += gauss(x,y);
		}
	}

	for (int x = 0; x < kernelSizeParameter; x++) {
		for (int y = 0; y < kernelSizeParameter; y++) {
			gauss(x,y) /= sum;
		}
	}

	return 0;
}

