/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       VYF project                                            *
*                       Author: Matus Bicanovsky                               *
*                       Brno 2025                                              *
*                                                                              *
*                       Implementation of the TMOHsin11 class                  *
*             Color to grayscale transform preserving natural order of hues    *
*******************************************************************************/
/* --------------------------------------------------------------------------- *
 * TMOHsin11.cpp: implementation of the TMOHsin11 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOHsin11.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOHsin11::TMOHsin11()
{
	SetName(L"Hsin11");					
	SetDescription(L"Color to grayscale image operator based on paper Color to grayscale transform preserving natural order of hues"); 
}

TMOHsin11::~TMOHsin11()
{
}
//function implementing the first stage of the algorthm, Global Mapping 
cv::Mat TMOHsin11::globalMapping(cv::Mat& input)
{
	double betaDegrees = 344.0;             		//value of beta in degrees from the paper
	double beta = betaDegrees * (M_PI / 180.0); 	//convert to radians

	int rows = input.rows;
	int cols = input.cols;

	cv::Mat Iw(rows, cols, CV_32F);          //matrix to store result of equation 2 from paper (intesity Iw)
	cv::Mat Cw(rows, cols, CV_32F);          //matrix to store result of equation 3 from paper (chrominance Cw)

	cv::Mat input32F;
	input.convertTo(input32F, CV_32FC3);		

	cv::Mat hsvInput;
	cv::cvtColor(input32F, hsvInput, cv::COLOR_BGR2HSV_FULL); 		//convert to HSV color space
	
	//calculations of weighted intensity and chrominance
	for(int y = 0; y < rows; y++)
	{
		for(int x = 0; x < cols; x++)
		{
			cv::Vec3f bgrPix = input32F.at<cv::Vec3f>(y, x); 				//get pixel value in BGR format
			cv::Vec3f hsvPix = hsvInput.at<cv::Vec3f>(y, x); 			//get pixel value in HSV format

			float R = bgrPix[2];  
			float G = bgrPix[1];  
			float B = bgrPix[0];  

			float I = 0.5*R + (1.0/3.0)*G + (1.0/6.0)*B;  				//equation 2 from paper for weighted intensity
			Iw.at<float>(y,x) = I;

			float S = hsvPix[1];   									//Saturation value from HSV color space
			float theta = hsvPix[0] * (M_PI / 180.0); 					//convert hue to radians
			float Hw = 0.5 * std::cos(theta + beta) + 0.5;             //equation 4
			Cw.at<float>(y,x) = I * S * Hw;                            //equation 3
		}
	}

	// Lambda search section, equations 8-10 from paper

	//compute matrix T used in equation 8 and 9
	cv::Mat T(rows, cols, CV_32F);
	for (int y = 0; y < rows; ++y)
	{
		for (int x = 0; x < cols; ++x)
		{
			cv::Vec3f bgrPix = input32F.at<cv::Vec3f>(y, x); 

			float R = bgrPix[2];
			float G = bgrPix[1];
			float B = bgrPix[0];

			//check for grayscale pixel
			bool isGray = std::abs(R - G) < 1e-5 && std::abs(G - B) < 1e-5;

			if (isGray)
			{
				T.at<float>(y, x) = 0.0f;			//black out the gray pixels
			}
			else
			{
				//use same intensity equation as Iw
				float I = 0.5f * R + (1.0f / 3.0f) * G + (1.0f / 6.0f) * B;
				T.at<float>(y, x) = I;
			}
		}
	}

	double maxC, minC;
	cv::minMaxLoc(Cw, &minC, &maxC); 						//find min and max of Cw
	double lambdaMax = 0.5 / std::max(maxC, 1e-7); 			//max value of Lambda defined in paper

	double bestLambda = 0.0;
	double bestStd = -1;
	cv::Mat tmp, gx, gy, grad;

	for(double lambda = 0.0; lambda <= lambdaMax; lambda += 0.01)
	{
		tmp = T + lambda * Cw;						//equation 9
		cv::threshold(tmp, tmp, 1.0, 1.0, cv::THRESH_TRUNC);

		cv::Sobel(tmp, gx, CV_32F, 1, 0, 3);		//gradient in x direction
		cv::Sobel(tmp, gy, CV_32F, 0, 1, 3);		//gradient in y direction
		cv::magnitude(gx, gy, grad);	
		
		double minG, maxG;
		cv::minMaxLoc(grad, &minG, &maxG, nullptr, nullptr);
		
		cv::Scalar m,s;
		cv::meanStdDev(grad, m, s);				//calculate mean and std of gradient
	
		if(s[0] > bestStd)
		{
			bestStd = s[0];						//update best std
			bestLambda = lambda;				//update best lambda
		}
	}
	std::cerr << "Best lambda = " << bestLambda << std::endl;
	return Iw + bestLambda * Cw;			//return the final result of equation 7 from paper
}

//function implementing conjugate gradient method solver using Eigen library and its functions
cv::Mat TMOHsin11::conjugateGrad(cv::Mat& input)
{
	//flatten into 1D array for Eigen
	int rows = input.rows;
	int cols = input.cols;
	int size = rows * cols;
	auto idx = [cols](int y, int x) { return y*cols + x; };	


	Eigen::SparseMatrix<double> A(size, size);				//sparse matrix for Laplacian operator
	A.reserve(Eigen::VectorXi::Constant(size, 5));			//reserve space for 5 elements per row (current pixel and 4 neighbors)
	Eigen::VectorXd b(size);								//vector for the right side of the equation Ax = b

	for(int y = 0; y < rows; y++)
	{
		for(int x = 0; x < cols; x++)
		{
			int i = idx(y, x);
			double diagonal = 0.0;

			auto add = [&](int yy, int xx)
			{
				if(yy < 0 || yy >= rows || xx < 0 || xx >= cols) return;		//skip out of bounds "neighbors"
				int j = idx(yy, xx);
				A.coeffRef(i, j) += -1.0;	//off-diagonal weight
				diagonal += 1.0;
				
							
			};

			add(y-1, x);	//top neighbor
			add(y+1, x);	//bottom neighbor
			add(y, x-1);	//left neighbor
			add(y, x+1);	//right neighbor

			A.coeffRef(i, i) = diagonal;
			b[i] = input.at<double>(y, x);	//set the right side of Ax=b to the input value
		}
	}
	A.makeCompressed();	
	//conjugate gradient solver
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	cg.setTolerance(1e-10);					//set tolerance for convergence
	cg.setMaxIterations(400);				//set max iterations
	cg.compute(A);							//compute the sparse matrix
	Eigen::VectorXd x = cg.solve(b);		//solve the system of equations
	//reshape vector back to cv::Mat
	cv::Mat result(rows, cols, CV_64F);
	for(int y = 0; y < rows; y++)
	{
		for(int x2 = 0; x2 < cols; x2++)
		{
			int i = idx(y, x2);
			result.at<double>(y, x2) = x[i];	
		}
	}
	return result;	
}
//function implementing second stage of algorithm, the Local lightness optimization based on input image and grayscale result from global mapping stage
cv::Mat TMOHsin11::localLightness(cv::Mat& grayscale, cv::Mat& colorInput)
{
	cv::Mat G;
	grayscale.convertTo(G, CV_64F);

	cv::Mat input32F;
	colorInput.convertTo(input32F, CV_32FC3);	
	//convert color input to Lab color space to get L* channel
	cv::Mat lab;
	cv::cvtColor(input32F, lab, cv::COLOR_BGR2Lab);	//convert to Lab color space
	std::vector<cv::Mat> channels;
	cv::split(lab, channels);						//split into Lab into channels
	cv::Mat L_channel;
	channels[0].convertTo(L_channel, CV_64F, 1.0/100.0);	//scale values in L* channel to [0,1]

	//constants for mask window
	int kernelRadius = 12;
	int kernelSize = 25;			//neigborhood of 25x25 

	double sigma = 1.1;
	cv::Mat g1 = cv::getGaussianKernel(kernelSize, sigma, CV_64F);		//Gaussian kernel for local lightness
	cv::Mat g = g1 * g1.t();											//2D Gaussian kernel

	cv::threshold(G, G, 1.0, 1.0, cv::THRESH_TRUNC);	//truncate values to [0,1]

	//matrix for storing local contrast
	cv::Mat H(grayscale.rows, grayscale.cols, CV_64F, 0.0);
	
	//double for-loop over the neighborhood (m,n) of each pixel
	for(int m = -kernelRadius; m <= kernelRadius; m++)
	{
		for(int n = -kernelRadius; n <= kernelRadius; n++)
		{
			if(n == 0 && m == 0) continue;			//skip comparing to itself

			double weight = g.at<double>(m + kernelRadius, n + kernelRadius);	//get weight from Gaussian kernel for equation 14
			//for every pixel and current neighbour we calculate ΔG and ΔL* and then choose one with larger value
			for(int y = 0; y < grayscale.rows; y++)
			{
				int yy = y - m;
				if (yy < 0 || yy >= grayscale.rows) continue;	//skip if out of bounds
				for(int x = 0; x < grayscale.cols; x++)
				{
					int xx = x - n;
					if (xx < 0 || xx >= grayscale.cols) continue;	//skip if out of bounds	
					
					//compute local contrast from grayscale image and L* channel, equations 11 and 12
					double deltaG = G.at<double>(y, x) - G.at<double>(yy, xx);	//equation 11
					double deltaL = L_channel.at<double>(y, x) - L_channel.at<double>(yy, xx);	//equation 12

					//choose stronger local contrast
					double delta = (std::abs(deltaG) >= std::abs(deltaL)) ? deltaG : deltaL;	//equation 13
					H.at<double>(y, x) += static_cast<double>(weight * delta);	//equation 14
				}
			}
		}
	}
	H -= cv::mean(H)[0];								//subtract mean from local contrast image
	cv::Mat I = conjugateGrad(H);						//call to conjugate gradient solver
	cv::normalize(I, I, 0, 1, cv::NORM_MINMAX);			//normalize to [0,1]
	return I;
}

/* --------------------------------------------------------------------------- *
 * main function of implementation of this color-to-grayscale operator *
 * --------------------------------------------------------------------------- */
int TMOHsin11::Transform()
{

	double *pSourceData = pSrc->GetData();		
	double *pDestinationData = pDst->GetData(); 
												
	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();
	cv::Mat inputImage(height, width, CV_64FC3);

	for(int j = 0; j < height; j++)
	{
		for(int i = 0; i < width; i++)
		{
			double R = *pSourceData++;
			double G = *pSourceData++;
			double B = *pSourceData++;

			inputImage.at<cv::Vec3d>(j, i) = cv::Vec3d(B, G, R);    //storing as BGR for opencv functions
		}
	}
	
	cv::Mat globalResult = globalMapping(inputImage);	//call to global mapping function

	cv::Mat result = localLightness(globalResult, inputImage);	//call to local lightness function
	
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double val = result.at<double>(j, i);	

			// and store results to the destination image
			*pDestinationData++ = val;
			*pDestinationData++ = val;
			*pDestinationData++ = val;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());

	return 0;
}
