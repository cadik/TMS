// Pallas Solver
// Copyright 2015. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: ryan.latture@gmail.com (Ryan Latture)

#include "glog/logging.h"

// Each solver is defined in its own header file.
// include the solver you wish you use:
#include "pallas/simulated_annealing.h"
#include <opencv2/opencv.hpp>
// define a problem you wish to solve by inheriting
// from the pallas::GradientCostFunction interface
// and implementing the Evaluate and NumParameters methods.
#include <stdio.h>
#include <time.h>
#include <cfloat>
#include <stdlib.h>     /* srand, rand */
//#include "pokus.h"

using namespace cv;
using namespace std;
using namespace pallas;


cv::Mat myOwn2DFilter(cv::Mat image, const double* x, int rows, int cols)
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
        cv::Mat prvniCast = myOwn2DFilter(SmoothenedImage, parameters, Width, Height);
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

    pallas::SimulatedAnnealing::Options options;

    options.max_iterations = 50;
    options.dwell_iterations = 50;
    options.max_stagnant_iterations = 50;

    options.cooling_schedule_options.initial_temperature = 1000;

    options.minimum_cost = 0.001;

    double upper_bounds [height*width];
    std::fill_n(upper_bounds, height*width, 10);
    double lower_bounds [height*width];
    std::fill_n(lower_bounds, height*width, 0);

    double step_size = 1.0;

    pallas::scoped_ptr<pallas::StepFunction> step_function (new pallas::BoundedStepFunction(step_size,
                                                                                            upper_bounds,
                                                                                            lower_bounds,
                                                                                            height*width));
    options.set_step_function(step_function);

    pallas::SimulatedAnnealing::Summary summary;
    pallas::GradientProblem problem(new sigmaClass(width, height, originalImage, smoothenedImage));
    pallas::Solve(options, problem, parameters, &summary);

    cv::Mat a;
	a = cv::Mat::zeros(height, width, CV_32F);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            a.at<float>(j, i) = parameters[width * j + i];
        }
    }

    return a;
}

cv::Mat getThirdBaseLayer(cv::Mat image, cv::Mat sigmaMap, int rows, int cols)
{
	cv::Mat filteredImage;
	filteredImage = cv::Mat::zeros(rows, cols, CV_32F);
	for(int j = 0; j < rows; ++j)
	{
		for(int i =0; i < cols; ++i)
		{
			/*
			 * getting gaussian kernel
			 * */
			long double sum = 0;
			cv::Mat gaussFilter;
			if (sigmaMap.at<float>(j, i) != 0) {
				gaussFilter = cv::Mat::zeros(sigmaMap.at<float>(j, i) * 2 + 1, sigmaMap.at<float>(j, i) * 2 + 1, CV_32F);
				gaussFilter = cv::getGaussianKernel(
					sigmaMap.at<float>(j, i) * 2 + 1,
					sigmaMap.at<float>(j, i) / 3.0,
					CV_32F
				);
			} else {
				gaussFilter = cv::Mat::ones(1, 1, CV_32F);
			}
			
			
			cv::Mat gauss2DFilter = gaussFilter * gaussFilter.t();
			
			int kCenter = (int)(sigmaMap.at<float>(j, i));
			int kSize = (int)(sigmaMap.at<float>(j, i) * 2 + 1);
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

// jejich fce
class Rosenbrock : public pallas::GradientCostFunction {
public:
    virtual ~Rosenbrock() {}

    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {
        const double x = parameters[0];
        const double y = parameters[1];

        cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
        if (gradient != NULL) {
            gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
            gradient[1] = 200.0 * (y - x * x);
        }
        return true;
    }

    virtual int NumParameters() const { return 2; }
};

/*
    Detail Maximalization./
*/

class detailClass : public pallas::GradientCostFunction {
public:
    int Width;
    int Height;
    cv::Mat DetailImage;
    cv::Mat Weight1;
    cv::Mat Weight2;

    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {

        /*
            pridat omezeni pro vydeleni pro V kanal
        */
        /**
         *  ||s_iD_i||^2
         */
        long double sum1 = 0;
        for (int j = 0; j < Height; j++) {
            for (int i = 0; i < Width; i++) {
                sum1 -= pow(DetailImage.at<float>(j, i) * parameters[Width * j + i], 2);
            }
        }
        /**
         *  r1*sum(w_i||\/s_i||^2
         */
        long double sum2 = 0;

        for (int j = 0; j < Height; j++) {
            for (int i = 0; i < Width; i++) {
                long double nablaPowerS = 0;
                if (i - 1 < 0 && j - 1 < 0) {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i + 1)], 2) -
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j - 1) + i]) + 
                        pow(parameters[Width * (j - 1) + i], 2);
                } else if (i + 1 == Width && j + 1 == Height) {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i - 1)], 2) -
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j + 1) + i]) + 
                        pow(parameters[Width * (j + 1) + i], 2);
                } else if (i - 1 < 0) {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i + 1)], 2) -
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j - 1) + i]) + 
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j - 1) + i], 2) -
                        2*(parameters[Width * (j - 1) + i])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j + 1) + i], 2);
                } else if (i + 1 == Width) {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i - 1)], 2) +
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j - 1) + i]) -
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j - 1) + i], 2) -
                        2*(parameters[Width * (j - 1) + i])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j + 1) + i], 2);
                } else if (j - 1 < 0) {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i - 1)], 2) -
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * j + (i + 1)]) - 
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * j + (i + 1)], 2) +
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j + 1) + i], 2);
                } else if (j + 1 == Height) {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i - 1)], 2) -
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * j + (i + 1)]) + 
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j - 1) + i]) +
                        pow(parameters[Width * j + (i + 1)], 2) -
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j - 1) + i]) +
                        pow(parameters[Width * (j - 1) + i], 2);
                } else {
                    nablaPowerS = 
                        pow(parameters[Width * j + (i - 1)], 2) -
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * j + (i + 1)]) + 
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j - 1) + i]) -
                        2*(parameters[Width * j + (i - 1)])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * j + (i + 1)], 2) -
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j - 1) + i]) +
                        2*(parameters[Width * j + (i + 1)])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j - 1) + i], 2) +
                        2*(parameters[Width * (j - 1) + i])*(parameters[Width * (j + 1) + i]) +
                        pow(parameters[Width * (j + 1) + i], 2);
                }

                sum2 += Weight1.at<float>(j, i)*nablaPowerS;
            }
        }
        /**
         *  r1*sum(w_i||\/s_i||^2
         */
        long double sum3 = 0;

        for (int j = 0; j < Height; j++) {
            for (int i = 0; i < Width; i++) {
                long double nablaPowerT = 0;
                if (i - 1 < 0 && j - 1 < 0) {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i + 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) + 
                        pow(parameters[(Width*Height) + Width * (j - 1) + i], 2);
                } else if (i + 1 == Width && j + 1 == Height) {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i - 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) + 
                        pow(parameters[(Width*Height) + Width * (j + 1) + i], 2);
                } else if (i - 1 < 0) {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i + 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) + 
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j - 1) + i], 2) -
                        2*(parameters[(Width*Height) + Width * (j - 1) + i])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j + 1) + i], 2);
                } else if (i + 1 == Width) {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i - 1)], 2) +
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) -
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j - 1) + i], 2) -
                        2*(parameters[(Width*Height) + Width * (j - 1) + i])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j + 1) + i], 2);
                } else if (j - 1 < 0) {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i - 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * j + (i + 1)]) - 
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * j + (i + 1)], 2) +
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j + 1) + i], 2);
                } else if (j + 1 == Height) {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i - 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * j + (i + 1)]) + 
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) +
                        pow(parameters[(Width*Height) + Width * j + (i + 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j - 1) + i], 2);
                } else {
                    nablaPowerT = 
                        pow(parameters[(Width*Height) + Width * j + (i - 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * j + (i + 1)]) + 
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) -
                        2*(parameters[(Width*Height) + Width * j + (i - 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * j + (i + 1)], 2) -
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j - 1) + i]) +
                        2*(parameters[(Width*Height) + Width * j + (i + 1)])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j - 1) + i], 2) +
                        2*(parameters[(Width*Height) + Width * (j - 1) + i])*(parameters[(Width*Height) + Width * (j + 1) + i]) +
                        pow(parameters[(Width*Height) + Width * (j + 1) + i], 2);
                }

                sum3 += Weight2.at<float>(j, i)*nablaPowerT;
            }
        }
        // works
       /* if (parameters[0] > 2) {
            return false;
        }*/
        //std::cout << "dssd" << std::endl;
        cost[0] = sum1 + sum2*200 + sum3*500;// sum1 + sum2 * 200; //+ 200*sum2;// + sum2;

        return true;
    }

    virtual int NumParameters() const { 
        return Width * Height * 2; 
    }
    detailClass(int width, int height, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2);
};

detailClass::detailClass(int width, int height, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2) {
    this->Width = width;
    this->Height = height;
    this->DetailImage = detailImage;
    this->Weight1 = weight1;
    this->Weight2 = weight2;
}

std::vector<cv::Mat> optimizeForGettingSAndTparameters(int height, int width, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2) {
    google::InstallFailureSignalHandler();
    std::cout << "done" << std::endl;
    // define the starting point for the optimization
    srand ( time(NULL) );
    double parameters[height*width*2];

    for (int i = 0; i < height*width*2; i++) {
        parameters[i] = rand() % 10 - 5;// / (rand() % 10);
    }

    pallas::SimulatedAnnealing::Options options;

    options.max_iterations = 50;
    options.dwell_iterations = 50;
    options.max_stagnant_iterations = 50;

    options.cooling_schedule_options.initial_temperature = 1000;

   // options.minimum_cost = 0.1;

    double upper_bounds [height*width*2];
    std::fill_n(upper_bounds, height*width*2, 5);
    double lower_bounds [height*width*2];
    std::fill_n(lower_bounds, height*width*2, -5);

    double step_size = 0.25;
    std::cout << "done" << std::endl;
    pallas::scoped_ptr<pallas::StepFunction> step_function (new pallas::BoundedStepFunction(step_size,
                                                                                            upper_bounds,
                                                                                            lower_bounds,
                                                                                            height*width*2));
    options.set_step_function(step_function);
    std::cout << "dofdssfdne" << std::endl;
    pallas::SimulatedAnnealing::Summary summary;
    pallas::GradientProblem problem(new detailClass(width, height, detailImage, weight1, weight2));
    std::cout << "dsdffsaaaaone" << std::endl;

    pallas::Solve(options, problem, parameters, &summary);
    std::vector<cv::Mat> sAndT;
    std::cout << summary.FullReport() << std::endl;
    cv::Mat s;
	s = cv::Mat::zeros(height, width, CV_32F);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            s.at<float>(j, i) = parameters[width * j + i];
        }
    }
    sAndT.push_back(s);
    cv::Mat t;
    t = cv::Mat::zeros(height, width, CV_32F);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            t.at<float>(j, i) = parameters[(width*height) + width * j + i];
        }
    }
    sAndT.push_back(t);
    return sAndT;
}

cv::Mat detailLayer(const cv::Mat &image1, const cv::Mat &image, int rows, int cols) {
    cv::Mat detail; 
    detail = cv::Mat::zeros(rows, cols, CV_32F);  
    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            detail.at<float>(j, i) = abs(image1.at<float>(j, i) - image.at<float>(j, i));
        }
    }
    return detail;
}

cv::Mat getGradient(const cv::Mat &image, int rows, int cols) {
    cv::Mat gradient; 
    gradient = cv::Mat::zeros(rows, cols, CV_32F);  
    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            double x = 0;
            double y = 0;
            if (i > 0 && i < cols - 1) {
                x = image.at<float>(j, i - 1) - image.at<float>(j, i + 1);
            } else {
                x = image.at<float>(j, i);
            }

            if (j > 0 && j < rows - 1) {
                y = image.at<float>(j - 1, i) - image.at<float>(j + 1, i);
            } else {
                y = image.at<float>(j, i);
            }

            gradient.at<float>(j, i) = x + y;
        }
    }
    return gradient;
}
/*
 * Function for getting weights for detail maximilization
 **/
cv::Mat getWeightsFromBaseLayer(const cv::Mat &gradient, int rows, int cols, int r){
    
	double a = 0.2;
    // initialize
    cv::Mat weights; 
    weights = cv::Mat::zeros(rows, cols, CV_32F);  
    
    for (int j = 0; j < rows; j++) {
		for (int i = 0; i < cols; i++) {
			if (abs(gradient.at<float>(j, i)/a) <= 1) {
				weights.at<float>(j, i) = pow(1 - pow(abs(gradient.at<float>(j, i)/a), 3), 3);
				
				if (weights.at<float>(j, i) * r <= 2.0) {
					weights.at<float>(j, i) += 2.0/(float)r;
				}
			} else {
				weights.at<float>(j, i) = 0;
			}
		}	
	}
    return weights;
}

//constrol
cv::Mat getDetailLayer(const cv::Mat &orig, const cv::Mat &base, int rows, int cols) 
{		
		cv::Mat detailLayer;
		detailLayer = cv::Mat::zeros(rows, cols, CV_32F);

		for(int j=0; j<rows; j++){
            for(int i=0; i<cols; i++){
                detailLayer.at<float>(j, i) = orig.at<float>(j, i) - base.at<float>(j, i);
            }
        }
        
        return detailLayer;
}

cv::Mat getDetailControl(const cv::Mat &base, const cv::Mat &detail,const cv::Mat &s,const cv::Mat &t,float mu, int rows, int cols) 
{
		cv::Mat detailLayer;
		detailLayer = cv::Mat::zeros(rows, cols, CV_32F);
		(s).convertTo(s, CV_32F);
		(t).convertTo(t, CV_32F);
		for(int j=0; j<rows; j++){
            for(int i=0; i<cols; i++){
                detailLayer.at<float>(j, i) = (mu * s.at<float>(j, i) + (1 - mu)) * detail.at<float>(j, i) + base.at<float>(j, i) + mu * t.at<float>(j, i);
            }
        }
        
        return detailLayer;
}

int main(int argc, char** argv) {
    std::cout << "hello wordl" << std::endl;
    cv::Mat originalImage = cv::imread("../data/dahlia.png");
    Mat originalImage_splited[3];   //destination array
    split(originalImage,originalImage_splited);//split source
    (originalImage_splited[0]).convertTo(originalImage_splited[0], CV_32F);
    (originalImage_splited[1]).convertTo(originalImage_splited[1], CV_32F);
    (originalImage_splited[2]).convertTo(originalImage_splited[2], CV_32F);
    std::cout << "dahlia separated" << std::endl;

    cv::Mat smoothedImage = cv::imread("../data/fin.png");
    Mat smoothedImage_splited[3];   //destination array
    split(smoothedImage,smoothedImage_splited);
    (smoothedImage_splited[0]).convertTo(smoothedImage_splited[0], CV_32F);
    (smoothedImage_splited[1]).convertTo(smoothedImage_splited[1], CV_32F);
    (smoothedImage_splited[2]).convertTo(smoothedImage_splited[2], CV_32F);

    std::cout << "dahlia smoothed separated" << std::endl;
    int width = originalImage.cols;
    int height = originalImage.rows;
    cv::Mat sumOriginalImage = cv::Mat::zeros(height, width, CV_32F);
    cv::Mat sumSmoothedImage = cv::Mat::zeros(height, width, CV_32F);

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            sumOriginalImage.at<float>(j, i) = (originalImage_splited[0]).at<float>(j, i) + originalImage_splited[1].at<float>(j, i) + originalImage_splited[2].at<float>(j, i);
            sumSmoothedImage.at<float>(j, i) = (smoothedImage_splited[0]).at<float>(j, i) + smoothedImage_splited[1].at<float>(j, i) + smoothedImage_splited[2].at<float>(j, i);
           // std::cout << originalImage_splited[0].at<float>(j, i) << std::endl;
        }
    }

    cv::Mat sumOfDetail = detailLayer(sumOriginalImage, sumSmoothedImage, height, width);
    cv::Mat gradient = getGradient(sumSmoothedImage, height, width); 
    cv::Mat weight1 = getWeightsFromBaseLayer(gradient/256.0, height, width, 200);
    cv::Mat weight2 = getWeightsFromBaseLayer(gradient/256.0, height, width, 500);


    std::cout << "done" << std::endl;
    std::vector<cv::Mat> kokotina = optimizeForGettingSAndTparameters(height, width, sumOfDetail, weight1, weight2);
    std::cout << kokotina[0].at<float>(1,2) << std::endl;

    std::cout << "getDetailLayer" << std::endl;
    cv::Mat detailB = getDetailLayer(originalImage_splited[0], smoothedImage_splited[0], height, width);
    cv::Mat detailG = getDetailLayer(originalImage_splited[1], smoothedImage_splited[1], height, width);
    cv::Mat detailR = getDetailLayer(originalImage_splited[2], smoothedImage_splited[2], height, width);
    //  std::cout << detailR<< std::endl;
    std::cout << "detailMaximalization" << std::endl;
    cv::Mat detailMaximizedLayerB = getDetailControl(smoothedImage_splited[0], detailB, kokotina[0], kokotina[1], 0.5, height, width);
    cv::Mat detailMaximizedLayerG = getDetailControl(smoothedImage_splited[1], detailG, kokotina[0], kokotina[1], 0.5, height, width);
    cv::Mat detailMaximizedLayerR = getDetailControl(smoothedImage_splited[2], detailR, kokotina[0], kokotina[1], 0.5, height, width);
    
    //std::cout << detailMaximizedLayerB << std::endl;
    cv::Mat fin_img;
    vector<Mat> channels;
    channels.push_back(detailMaximizedLayerB);
    channels.push_back(detailMaximizedLayerG);
    channels.push_back(detailMaximizedLayerR);
    merge(channels, fin_img);

    cv::imwrite("../data/finDetail.png", fin_img);
    //std::cout << "summary of both channels complete" << std::endl;
    //cv::Mat sigma = optimizeForSigma(height, width, sumOriginalImage, sumSmoothedImage);
    // std::cout << sigma << std::endl;
    //std::cout << "sigma is complete" << std::endl;

   /* cv::Mat basePhase3R = getThirdBaseLayer(originalImage_splited[0], sigma, height, width);
	cv::Mat basePhase3G = getThirdBaseLayer(originalImage_splited[1], sigma, height, width);
	cv::Mat basePhase3B = getThirdBaseLayer(originalImage_splited[2], sigma, height, width);

    cv::Mat fin_img;
    vector<Mat> channels;
    channels.push_back(basePhase3R);
    channels.push_back(basePhase3G);
    channels.push_back(basePhase3B);
    merge(channels, fin_img);

    cv::imwrite("../data/fin.png", fin_img);*/
    return 0;
}