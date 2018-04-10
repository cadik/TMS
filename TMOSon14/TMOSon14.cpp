/* --------------------------------------------------------------------------- *
 * TMOSon14.cpp: implementation of the TMOSon14 class.   *
 * --------------------------------------------------------------------------- */
#include "TMOSon14.h"
// #include <fftw3.h>
//#include <opencv2/core/core.hpp>
// #include "L0minimization.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;
using namespace Eigen;
/* --------------------------------------------------------------------------- *
 * Constructor se<ves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
/*
	pokus pallas
*/
#include <opencv2/opencv.hpp>
#include <stdlib.h>     /* srand, rand */


// #include "glog/logging.h"
#include "pallas/simulated_annealing.h"

#include <cfloat>
/*
#include "pallas/cooling_schedule.h"
#include "pallas/history_concept.h"
#include "pallas/step_function.h"
#include "pallas/types.h"
#include "pallas/internal/state.h"
#include "pallas/internal/metropolis.h"
#include "pallas/scoped_ptr.h"*/

#include <stdio.h>
#include <time.h>

cv::Mat myOwn2DFilterOptimization(cv::Mat image, const double* x, int rows, int cols)
{  
    //std::cout << image.rows << std::endl;
   // std::cout << image.cols << std::endl;
	cv::Mat filteredImage;
	filteredImage = cv::Mat::zeros(rows, cols, CV_32F);

	for(int j = 0; j < rows; ++j)
	{
		for(int i =0; i < cols; ++i)
		{
			
            //std::cout << "ahoj" << std::endl;
			long double sum = 0;
			cv::Mat gaussFilter;
           //  std::cout << x[cols * i + j] << std::endl;
			if (round(x[cols * i + j]) != 0) {
				gaussFilter = cv::Mat::zeros(round(x[cols * i + j]) * 2 + 1, round(x[cols * i + j]) * 2 + 1, CV_32F);
				gaussFilter = cv::getGaussianKernel(
					round(x[cols * i + j]) * 2 + 1,
					round(x[cols * i + j]) / 3.0,
					CV_32F
				);
			} else {
				gaussFilter = cv::Mat::ones(1, 1, CV_32F);
			}
			
            //std::cout << x[cols * i + j] << std::endl;
			cv::Mat gauss2DFilter = gaussFilter * gaussFilter.t();
			
            int kCenter = (int)(round(x[cols * i + j]));
			int kSize = (int)(round(x[cols * i + j] * 2 + 1));      
              
			for(int m = 0; m < kSize; ++m)    
			{
				int mm = kSize - 1 - m;      

				for(int n = 0; n < kSize; ++n) 
				{
					int nn = kSize - 1 - n;  

					int jj = j + (m - kCenter);
					int ii = i + (n - kCenter);

					if( jj >= 0 && jj < rows && ii >= 0 && ii < cols ) {
						long double mut1 = image.at<float>(jj, ii);
						long double mut2 = gauss2DFilter.at<float>(mm, nn);
						sum += mut1 * mut2;
					}
				}
			}
			filteredImage.at<float>(j, i) = sum;
		}
	}

    return filteredImage;
}

cv::Mat getSecondTerm (const double* x, int rows, int cols, float ypsilon) {
	cv::Mat a;
	a = cv::Mat::zeros(cols, rows, CV_32F);

	for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
           // std::cout << x[cols * j + i] << std::endl;
		   double b1 = 0;
		   double b2 = 0;
		   if ((j - 1 >= 0) && (j + 1 < rows)) {
			   // std::cout << b1 << std::endl;    
               b2 = pow(x[cols * i + j - 1] - x[cols * i + j + 1], 2);   
		   }

		   if ((i - 1 >= 0) && (i + 1 < cols)) {
			   // std::cout << b1 << std::endl;			          
			   b1 = pow(x[cols * (i - 1) + j] - x[cols * (i + 1) + j], 2);
		   }
		   a.at<float>(i, j) = ypsilon*(b1 + b2);
        }
    }
    // std::cout << a << std::endl;
    return a;
}

class sigmaClass : public pallas::GradientCostFunction {
public:
    int Width;
    int Height;
    cv::Mat OriginalImage;
    cv::Mat SmoothenedImage;

    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {

        float y = 1e-5;
        long double sum1 = 0;
        cv::Mat prvniCast = myOwn2DFilterOptimization(SmoothenedImage, parameters, Width, Height);
        long double sum2 = 0;
        cv::Mat druhaCast = getSecondTerm(parameters, Width, Height, y);
        for (int j = 0; j < Width; j++) {
            for (int i = 0; i < Height; i++) {
                sum1 += pow(prvniCast.at<float>(j, i) - OriginalImage.at<float>(j, i), 2);
                sum2 += druhaCast.at<float>(j, i);
                // std::cout << x[Width * j + i] << std::endl;
            }
        }
        cost[0] = sum1 + sum2;

        return true;
    }

    virtual int NumParameters() const { 
        return Width * Height; 
    }
    sigmaClass(int width, int height, cv::Mat originalImage, cv::Mat smoothenedImage);
};

sigmaClass::sigmaClass(int width, int height, cv::Mat originalImage, cv::Mat smoothenedImage) {
    this->Width = width;
    this->Height = height;
    this->OriginalImage = originalImage;
    this->SmoothenedImage = smoothenedImage;
}

cv::Mat optimizeForSigma(int height, int width, cv::Mat originalImage, cv::Mat smoothenedImage) {
    // define the starting point for the optimization
    srand ( time(NULL) );
    double parameters[height*width];

    for (int i = 0; i < height*width; i++) {
        parameters[i] = rand() % 10;
    }

    // set up global optimizeForSigmar options only initialization
    // is need to accept the default options
    pallas::SimulatedAnnealing::Options options;

    // Increase the number of iterations
    // SA often requires many iterations to get
    // with 1 to 2 significant figures of the
    // optimal solution
    options.max_iterations = 50; // 50
    options.dwell_iterations = 50;
    options.max_stagnant_iterations = 50;

    // set a high initial temperature to allow
    // the SA algorithm to fully explore the
    // parameter space
    options.cooling_schedule_options.initial_temperature = 1000;

    // quit the optimization of cost gets within
    // 3 significant figures of global minimum
    options.minimum_cost = 0.001;
    // define custom step function which will bound the
    // randomized candidate solution in order to limit
    // the search and speed up convergence
    double upper_bounds [height*width];
    std::fill_n(upper_bounds, height*width, 10);
    //std::cout << upper_bounds[1564] << std::endl;
    double lower_bounds [height*width];
    std::fill_n(lower_bounds, height*width, 0);
    //std::cout << lower_bounds[26] << std::endl;

    double step_size = 1.0;

    pallas::scoped_ptr<pallas::StepFunction> step_function (new pallas::BoundedStepFunction(step_size,
                                                                                            upper_bounds,
                                                                                            lower_bounds,
                                                                                            height*width));
    options.set_step_function(step_function);

    // initialize a summary object to hold the
    // optimization details
    pallas::SimulatedAnnealing::Summary summary;
   // Summary summary;
    // create a problem from your cost function
    // TOTO TAM NARVAT
    pallas::GradientProblem problem(new sigmaClass(width, height, originalImage, smoothenedImage));

    // solve the problem and store the optimal position
    // in parameters and the optimization details in
    // the summary
    pallas::Solve(options, problem, parameters, &summary);

    cv::Mat a;
	a = cv::Mat::zeros(height, width, CV_32F);
    // std::cout << "Global minimum found at:" << std::endl;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            a.at<float>(j, i) = parameters[width * j + i];
        //std::cout << "Parameter " << i << " : " << parameters[i] << std::endl;
        }
    }

    return a;
}
/*
	pokus pallas
*/
TMOSon14::TMOSon14()
{
	SetName(L"Son14");						// TODO - Insert operator name
	SetDescription(L"Art-Photography detail enhancement");	// TODO - Insert description
	/**
	 * Mu - Parameter
	 **/
	mu.SetName(L"Mu");				// TODO - Insert parameters names
	mu.SetDescription(L"Represents rate Mu for detail maximization");	// TODO - Insert parameter descriptions
	mu.SetDefault(0.5);							// TODO - Add default values
	mu=0.5;
	mu.SetRange(0.0,1.0);				// TODO - Add acceptable range if needed
	this->Register(mu);
	
}

TMOSon14::~TMOSon14()
{
}


/* --------------------------------------------------------------------------- *
 * This overloaded function
 * -------------------------------------------------------- */
int TMOSon14::Transform()
{ 
	ofstream myfile;
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	/*
	 * Base matrix 
	 **/
    cv::Mat r;
	cv::Mat g;
	cv::Mat b;
	
	r = cv::Mat::zeros(height, width, CV_32F);
	g = cv::Mat::zeros(height, width, CV_32F);
	b = cv::Mat::zeros(height, width, CV_32F);
	
	
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_RGB);								// This is format of Y as luminance
	// pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
											// of three doubles representing
	/*
	 * Fill base matrix
	 * */
	 int j;
	for (j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
			r.at<float>(j,i) = *pSourceData++; 
			g.at<float>(j,i) = *pSourceData++;  //getting separate RGB channels
			b.at<float>(j,i) = *pSourceData++;
		}
	}
    
	pSrc->ProgressBar(j, pSrc->GetHeight());
    /*
	 * For L0 smoothing
	 * */
	std::vector<cv::Mat> array_to_merge;

    array_to_merge.push_back(b*256);
    array_to_merge.push_back(g*256);
    array_to_merge.push_back(r*256);

    cv::Mat originalImage;
    
    cv::merge(array_to_merge, originalImage);
   
	/*
	 * Base Decomposition
	 **/
	 
	/*
	 * Phase 1 - Original L0 Smoothing
	 **/
	std::cout << "Base Phase1" << std::endl;
	std::cout << "Base Phase1" << std::endl;	
	// std::vector<cv::Mat> result = minimizeL0Gradient(color);
	cv::Mat basePhase1 = minimizeL0Gradient1(originalImage);
	 /*
     * split basePhase1, experiment; 
     **/
     
	cv::Mat basePhase1Chan[3];
	cv::split(basePhase1, basePhase1Chan);
	
    (basePhase1Chan[0]).convertTo(basePhase1Chan[0], CV_32F);
    (basePhase1Chan[1]).convertTo(basePhase1Chan[1], CV_32F);
    (basePhase1Chan[2]).convertTo(basePhase1Chan[2], CV_32F);  
	/*
	 * Phase 2 - L0 smooting with adaptive lambda matrix
	 **/	
	std::cout << "Base Phase2" << std::endl;
	cv::Mat gradientFrom1stSmoothing = getGradientMagnitude(basePhase1);
	cv::Mat adaptiveLambdaMatrix1 = getAdaptiveLambdaMatrix(gradientFrom1stSmoothing, height, width);
    cv::Mat basePhase2 = minimizeL0GradientSecondFaze(originalImage, adaptiveLambdaMatrix1, height, width);
	
	
	/*
     * split basePhase2; 
     **/
     
	cv::Mat basePhase2Chan[3];
	cv::split(basePhase2, basePhase2Chan);
	
    (basePhase2Chan[0]).convertTo(basePhase2Chan[0], CV_32F);
    (basePhase2Chan[1]).convertTo(basePhase2Chan[1], CV_32F);
    (basePhase2Chan[2]).convertTo(basePhase2Chan[2], CV_32F);
    
    
    /*
     * Phase 3 -- getting final base layer
     **/
     
    cv::Mat sumOfCostsBase = getSumOfCosts(basePhase2Chan[0], basePhase2Chan[1], basePhase2Chan[2], height, width);
    /*myfile.open ("./TMOSon14/txt/sumOfCostsBase.txt");
	myfile << sumOfCostsBase;
	myfile.close();*/
	
	cv::Mat sumOfCostsOriginal = getSumOfCosts(r, g, b, height, width);
	std::cout << "Base Phase3" << std::endl;
	cv::Mat sigmaMap = optimizeForSigma(height, width, sumOfCostsOriginal, sumOfCostsBase);
	//cv::Mat sigmaMap = stochasticOptimizationForGetSigma(sumOfCostsBase/256.0, sumOfCostsOriginal, height, width, 50000);
	
	cv::Mat basePhase3R = myOwn2DFilter(r*256, sigmaMap, height, width);
	cv::Mat basePhase3G = myOwn2DFilter(g*256, sigmaMap, height, width);
	cv::Mat basePhase3B = myOwn2DFilter(b*256, sigmaMap, height, width);
	std::cout << "Base phase -- COMPLETED" << std::endl;

     
    cv::Mat detailLayerR = getDetailLayer(r*256, basePhase3R, height, width);
    cv::Mat detailLayerG = getDetailLayer(g*256, basePhase3G, height, width);
    cv::Mat detailLayerB = getDetailLayer(b*256, basePhase3B, height, width);

    cv::Mat sumOfDetail = getSumOfCosts(detailLayerR, detailLayerG, detailLayerB, height, width);
    cv::Mat sumOfBase = getSumOfCosts(basePhase3R, basePhase3G, basePhase3B, height, width);
    

    std::vector<cv::Mat> array_to_merge1;

    array_to_merge1.push_back(basePhase3R);
    array_to_merge1.push_back(basePhase3G);
    array_to_merge1.push_back(basePhase3B);

    cv::Mat baseImage;
    
    cv::merge(array_to_merge1, baseImage);
    
    cv::Mat gradientOfBaseLayer = getGradientMagnitude(baseImage);

     
    cv::Mat r1Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 200);
    cv::Mat r2Layer = getWeightsFromBaseLayer(gradientOfBaseLayer, height, width, 500);
	
	std::vector<cv::Mat> detail;
	detail.push_back((detailLayerR.clone()) / 256.0);
	detail.push_back((detailLayerG.clone()) / 256.0);
	detail.push_back((detailLayerB.clone()) / 256.0);
	
	/*std::vector<cv::Mat> ST = detailMaximalization(sumOfBase/256.0, sumOfDetail/256.0, r1Layer, r2Layer, height, width, 1, detail);	
	std::cout << "Detail maximalization -- COMPLETED" << std::endl;

	cv::Mat detailMaximizedLayerR = getDetailControl(basePhase3R, detailLayerR, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerG = getDetailControl(basePhase3G, detailLayerG, ST[0], ST[1], mu, height, width);
    cv::Mat detailMaximizedLayerB = getDetailControl(basePhase3B, detailLayerB, ST[0], ST[1], mu, height, width);*/

/*  myfile.open ("./TMOSon14/txt/sumOfCostsBase.txt");
	myfile << sumOfCostsBase;
	myfile.close();*/
	/*
	 * Function for control details enhancement of picture 
	 **/

	/*
	 * Showing picture (shows blurred picture)
	 **/
	/* for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{													// simple variables
			// and store results to the destination image
			*pDestinationData++ = ((basePhase3R).at<float>(j,i) + (detailLayerR).at<float>(j,i)) / 256.0;
			*pDestinationData++ = ((basePhase3G).at<float>(j,i) + (detailLayerG).at<float>(j,i)) / 256.0;
			*pDestinationData++ = ((basePhase3B).at<float>(j,i) + (detailLayerB).at<float>(j,i)) / 256.0;
		}
	}*/
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{													// simple variables
			// and store results to the destination image
			// *pDestinationData++ = (detailMaximizedLayerR).at<float>(j,i) / 256.0;// + (detailChan[2]).at<float>(j,i)) / 256.0;
			// *pDestinationData++ = (detailMaximizedLayerG).at<float>(j,i) / 256.0;// + (detailChan[1]).at<float>(j,i)) / 256.0;
			// *pDestinationData++ = (detailMaximizedLayerB).at<float>(j,i) / 256.0;// + (detailChan[0]).at<float>(j,i)) / 256.0;
		
			*pDestinationData++ = (basePhase3R).at<float>(j,i) / 256.0;// + (detailChan[2]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = (basePhase3G).at<float>(j,i) / 256.0;// + (detailChan[1]).at<float>(j,i)) / 256.0;
			*pDestinationData++ = (basePhase3B).at<float>(j,i) / 256.0;//
		}
	}
	pDst->Convert(TMO_RGB);
return 0;
}
