/* --------------------------------------------------------------------------- *
 * TMOHsin11.cpp: implementation of the TMOHsin11 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOHsin11.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOHsin11::TMOHsin11()
{
	SetName(L"Hsin11");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOHsin11::~TMOHsin11()
{
}
cv::Mat TMOHsin11::globalMapping(cv::Mat& input)
{
	double betaDegrees = 344.0;             		//value of beta in degrees from the paper
	double beta = betaDegrees * (M_PI / 180.0); 	//convert to radians

	int rows = input.rows;
	int cols = input.cols;

	cv::Mat Iw(rows, cols, CV_32F);          //matrix to store result of equation 2 from paper
	cv::Mat Cw(rows, cols, CV_32F);          //matrix to store result of equation 3 from paper

	cv::Mat input32F;
	input.convertTo(input32F, CV_32FC3);		//convert input to double precision

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

			float I = 0.5*R + (1.0/3.0)*G + (1.0/6.0)*B;  				//equation 2 from paper
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
			cv::Vec3f bgrPix = input32F.at<cv::Vec3f>(y, x);  // Assume input32F is in [0,1]

			float R = bgrPix[2];
			float G = bgrPix[1];
			float B = bgrPix[0];

			// Check for grayscale pixel: R ≈ G ≈ B
			bool isGray = std::abs(R - G) < 1e-5 && std::abs(G - B) < 1e-5;

			if (isGray)
			{
				T.at<float>(y, x) = 0.0f;
			}
			else
			{
				// Use same intensity equation as Iw
				float I = 0.5f * R + (1.0f / 3.0f) * G + (1.0f / 6.0f) * B;
				T.at<float>(y, x) = I;
			}
		}
	}

	double maxC, minC;
	cv::minMaxLoc(Cw, &minC, &maxC); 						//find min and max of Cw
	double lambdaMax = 0.5 / std::max(maxC, 1e-7); 			//defined in paper

	double bestLambda = 0.0;
	double bestStd = -1;
	cv::Mat tmp, gx, gy, grad;

	cv::minMaxLoc(Cw, &minC, &maxC);
	std::cerr << "Cw range: [" << minC << ", " << maxC << "]" << std::endl;

	for(double lambda = 0.0; lambda <= lambdaMax; lambda += 0.01)
	{
		tmp = T + lambda * Cw;	//equation 9
		cv::threshold(tmp, tmp, 1.0, 1.0, cv::THRESH_TRUNC);

		cv::Sobel(tmp, gx, CV_32F, 1, 0, 3);		//gradient in x direction
		cv::Sobel(tmp, gy, CV_32F, 0, 1, 3);		//gradient in y direction
		cv::magnitude(gx, gy, grad);	
		
		double minG, maxG;
		cv::minMaxLoc(grad, &minG, &maxG, nullptr, nullptr);
		
		cv::Scalar m,s;
		cv::meanStdDev(grad, m, s);		//calculate mean and std of gradient
	
		if(s[0] > bestStd)
		{
			bestStd = s[0];						//update best std
			bestLambda = lambda;				//update best lambda
		}
	}
	std::cerr << "Max lambda = " << lambdaMax << std::endl;	//print max lambda value
	std::cerr << "Best lambda = " << bestLambda << std::endl;	//print best lambda value
	return Iw + bestLambda * Cw;			//return the final result of equation 7 from paper
}

cv::Mat TMOHsin11::conjugateGrad(cv::Mat& input)
{
	//flatten into 1D array for Eigen
	int rows = input.rows;
	int cols = input.cols;
	int size = rows * cols;
	auto idx = [cols](int y, int x) { return y*cols + x; };	//indexing function for 1D array

	//sparse matrix for the Laplacian operator
	Eigen::SparseMatrix<double> A(size, size);	//sparse matrix for Laplacian operator
	A.reserve(Eigen::VectorXi::Constant(size, 5));	//reserve space for 5 non-zero elements per row
	Eigen::VectorXd b(size);	//vector for the right-hand side of the equation

	for(int y = 0; y < rows; y++)
	{
		for(int x = 0; x < cols; x++)
		{
			int i = idx(y, x);
			double diagonal = 0.0;

			auto add = [&](int yy, int xx)
			{
				if(yy < 0 || yy >= rows || xx < 0 || xx >= cols) return;
				int j = idx(yy, xx);
				A.coeffRef(i, j) += -1.0;	//off-diagonal weight
				diagonal += 1.0;
				
							
			};

			add(y-1, x);	//top neighbor
			add(y+1, x);	//bottom neighbor
			add(y, x-1);	//left neighbor
			add(y, x+1);	//right neighbor

			A.coeffRef(i, i) = diagonal;
			b[i] = input.at<double>(y, x);	//set the right-hand side to the input value
		}
	}
	A.makeCompressed();	
	double maxRowSum = 0.0;

	for (int i = 0; i < A.rows(); ++i)
	{
		double s = 0.0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it)
			s += it.value();

		maxRowSum = std::max(maxRowSum, std::abs(s));
	}

	std::cerr << "Max |row-sum(A)| = " << std::scientific << maxRowSum << '\n';
	int i0 = 0;
	double rowsum = 0;
	for (Eigen::SparseMatrix<double>::InnerIterator it(A,i0); it; ++it)
    	rowsum += it.value();
	std::cerr << "row sum (0,0) = " << rowsum << std::endl;
	
	//conjugate gradient solver
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	cg.setTolerance(1e-10);	//set tolerance for convergence
	cg.setMaxIterations(400);	//set max iterations
	cg.compute(A);	//compute the sparse matrix
	Eigen::VectorXd x = cg.solve(b);	//solve the system of equations
	std::cerr << "max|Ax-b| = "
          << (A * x - b).cwiseAbs().maxCoeff() << std::endl;
	//reshape vector back to cv::Mat
	cv::Mat result(rows, cols, CV_64F);
	for(int y = 0; y < rows; y++)
	{
		for(int x2 = 0; x2 < cols; x2++)
		{
			int i = idx(y, x2);
			result.at<double>(y, x2) = x[i];	//set pixel value in result matrix
		}
	}
	return result;	
}

cv::Mat TMOHsin11::localLightness(cv::Mat& grayscale, cv::Mat& colorInput)
{
	cv::Mat G;
	grayscale.convertTo(G, CV_64F);

	cv::Mat input32F;
	colorInput.convertTo(input32F, CV_32FC3);	//convert input to double precision
	//convert color input to Lab color space to get L* channel
	cv::Mat lab;
	cv::cvtColor(input32F, lab, cv::COLOR_BGR2Lab);	//convert to Lab color space
	std::vector<cv::Mat> channels;
	cv::split(lab, channels);						//split into L*, a*, b* channels
	cv::Mat L_channel;
	channels[0].convertTo(L_channel, CV_64F, 1.0/100.0);	//scale values in L* channel to [0,1]

	//constants for mask window
	int kernelRadius = 6;
	int kernelSize = 13;

	double sigma = 1.1;
	cv::Mat g1 = cv::getGaussianKernel(kernelSize, sigma, CV_64F);		//Gaussian kernel for local lightness
	cv::Mat g = g1 * g1.t();											//2D Gaussian kernel

	double sumW = cv::sum(g)[0];          // g is CV_64F
	std::cerr << "Sum G = " << sumW << '\n';   // should print 1.00000

	//cv::Mat G = grayscale.clone();	//copy grayscale image for local lightness calculation
	cv::threshold(G, G, 1.0, 1.0, cv::THRESH_TRUNC);	//truncate values to [0,1]

	//matrix for storing local contrast
	cv::Mat H(grayscale.rows, grayscale.cols, CV_64F, 0.0);
	
	int G_amount = 0;	//counter for local contrast from grayscale image
	int L_amount = 0;	//counter for local contrast from L* channel
	//double for loop over the neighborhood (m,n) of each pixel
	for(int m = -kernelRadius; m <= kernelRadius; m++)
	{
		for(int n = -kernelRadius; n <= kernelRadius; n++)
		{
			if(n == 0 && m == 0) continue;			//skip comparing to itself

			double weight = g.at<double>(m + kernelRadius, n + kernelRadius);	//get weight from Gaussian kernel for equation 14
			//for every pixel we calculate ΔG and ΔL* and then choose one with larger value
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
					if(std::abs(deltaG) >= std::abs(deltaL))
					{
						G_amount += 1;
					}
					else{
						L_amount += 1;
					}
					H.at<double>(y, x) += static_cast<double>(weight * delta);	//equation 14
				}
			}
		}
	}
	
	std::cerr << "G_amount="<<G_amount<<"  L_amount="<<L_amount<<std::endl;
	double min, max;
	

	H -= cv::mean(H)[0];	//subtract mean from local contrast image

	cv::minMaxLoc(H, &min, &max);	//find min and max of local contrast
	std::cerr << "H ="<<min<<"  max="<<max
          << "  mean="<<cv::mean(H)[0] << std::endl;

	cv::Mat H8;
	H.convertTo(H8, CV_8U, 255/(max - min), -255*min/(max - min));	//convert to 8-bit for display
	cv::imwrite("H.png", H8);	//save local contrast image for debugging
	
	
	cv::Mat I = conjugateGrad(H);	//call to conjugate gradient solver
	I -= cv::mean(I)[0];	//subtract mean from result
	cv::Mat flat; I.reshape(0,1).copyTo(flat);
	cv::sort(flat,flat,cv::SORT_ASCENDING);
    double p1 = flat.at<double>( flat.total()*0.01 );
    double p99= flat.at<double>( flat.total()*0.99 );
    I = (I-p1)/(p99-p1); 
	I.setTo(0,I<0); I.setTo(1,I>1);
	return I;
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOHsin11::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
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
	//store the globalResult
	cv::Mat globalResultU8;
	globalResult.convertTo(globalResultU8, CV_8UC3, 255.0);	//convert to 8-bit for display
	cv::imwrite("globalResult.png", globalResultU8);	//save global result 

	cv::Mat result = localLightness(globalResult, inputImage);	//call to local lightness function
	double pY, px, py;
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double val = result.at<double>(j, i);	//get pixel value from global result

			// and store results to the destination image
			*pDestinationData++ = val;
			*pDestinationData++ = val;
			*pDestinationData++ = val;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0;
}
