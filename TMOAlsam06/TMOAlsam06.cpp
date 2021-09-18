/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       Brno 2018                                              *
*                                                                              *
*                       Grey colour sharpening                                 *
*                       Implementation of the Alsam A., Kosvas O.              *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOAlsam06.cpp
 * @brief Grey colour sharpening
 * @class TMOALsam06
 */

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

	if (kernelSizeParameter % 2 == 0) {
		kernelSizeParameter = kernelSizeParameter - 1;
	}

	MatrixXd ir(width, height);
	MatrixXd ig(width, height);
	MatrixXd ib(width, height);

	for (i = 0; i < pixelCount; i++) {
		ir(i % width, floor(i / width)) = *pSourceData++;
		ig(i % width, floor(i / width)) = *pSourceData++;
		ib(i % width, floor(i / width)) = *pSourceData++;
	}

	pSourceData = pSrc->GetData();

	double minOriginal = std::numeric_limits<double>::max();
	double maxOriginal = std::numeric_limits<double>::min();

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

		if (p(i,0) > maxOriginal) {
			maxOriginal = p(i,0);
		}

		if (p(i,0) < minOriginal) {
			minOriginal = p(i,0);
		}

		if (p(i,1) > maxOriginal) {
			maxOriginal = p(i,1);
		}

		if (p(i,1) < minOriginal) {
			minOriginal = p(i,1);
		}

		if (p(i,2) > maxOriginal) {
			maxOriginal = p(i,2);
		}

		if (p(i,2) < minOriginal) {
			minOriginal = p(i,2);
		}
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
		pSrc->ProgressBar(i, pixelCount * 3);

		/**
		 * Matrix ppu = projection of each pixel onto u1
		 */
		MatrixXd ppu(1, 3);
		ppu = p.row(i) * pu;

		g = ppu.squaredNorm();
		gu(i % width, floor(i / width)) = g;
	}

	MatrixXd gauss(kernelSizeParameter.GetInt(), kernelSizeParameter.GetInt());
	//double sum = 0;
	double mean = kernelSizeParameter.GetInt() / 2;
	int flooredMean = floor(mean);

	for (int x = 0; x < kernelSizeParameter.GetInt(); x++) {
		for (int y = 0; y < kernelSizeParameter.GetInt(); y++) {
			//gauss(x,y) = exp(-0.5 * (pow(x, 2) + pow(y, 2)) / pow(standardDeviationParameter, 2)) / (2 * M_PI * standardDeviationParameter * standardDeviationParameter);
			gauss(x,y) = exp(-0.5 * (pow((x - mean) / standardDeviationParameter.GetInt(), 2.0) + pow((y - mean) / standardDeviationParameter.GetInt(), 2.0))) / (2 * M_PI * standardDeviationParameter.GetInt() * standardDeviationParameter.GetInt());
			//sum += gauss(x,y);
		}
	}

	/*
	//normalization
	for (int x = 0; x < kernelSizeParameter; x++) {
		for (int y = 0; y < kernelSizeParameter; y++) {
			gauss(x,y) /= sum;
		}
	}
	*/

	MatrixXd gb(width, height);

	// convolution
	for (int y = 0; y < height; y++) {
		pSrc->ProgressBar(pixelCount + y * width, pixelCount * 3);
		for (int x = 0; x < width; x++) {
			double value = 0;
			int members = 0;

			for (int i = -flooredMean; i <= flooredMean; i++) {
				for (int j = -flooredMean; j <= flooredMean; j++) {
					if (x + i < 0 || y + j < 0 || x + i >= width || y + j >= height) {
						continue;
					}

					value += gu(x + i, y + j) * gauss(i + flooredMean, j + flooredMean);
					members++;
				}
			}

			gb(x,y) = value / members;
		}
	}

	/*
	* hr, hg, hb - high frequency masks for red, green and blue color channels
	*/
	MatrixXd hr(width, height);
	MatrixXd hg(width, height);
	MatrixXd hb(width, height);

	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			hr(x,y) = ir(x,y) - gb(x,y);
			hg(x,y) = ig(x,y) - gb(x,y);
			hb(x,y) = ib(x,y) - gb(x,y);
		}
	}

	double minNew = std::numeric_limits<double>::max();
	double maxNew = std::numeric_limits<double>::min();
	double imageNew[pixelCount];
	int index = 0;

	for (int y = 0; y < height; y++) {
		pSrc->ProgressBar(2 * pixelCount + y * width, pixelCount * 3);
		for (int x = 0; x < width; x++) {
			imageNew[index] = gu(x,y) + (hr(x,y) + hg(x,y) + hb(x,y)) / 3;

			if (imageNew[index] > maxNew) {
				maxNew = imageNew[index];
			}

			if (imageNew[index] < minNew) {
				minNew = imageNew[index];
			}

			index++;
		}
	}

	for (int i = 0; i < pixelCount; i++) {
		double value = ((imageNew[i] - minNew)/(maxNew - minNew)) * (maxOriginal - minOriginal) + minOriginal;

		*pDestinationData++ = value;
		*pDestinationData++ = value;
		*pDestinationData++ = value;
	}

	return 0;
}

