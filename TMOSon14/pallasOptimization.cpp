#include "pallasOptimization.h"

/*
    Part 1: functions for getting adaptive sigma map
*/
/*
    Function for getting first term, aka "2D filter with adaptive masks for each pixel"
*/
cv::Mat getFirstTerm(cv::Mat image, const double* x, int rows, int cols)
{  
	cv::Mat filteredImage;
	filteredImage = cv::Mat::zeros(rows, cols, CV_32F);

	for(int j = 0; j < rows; ++j)
	{
		for(int i =0; i < cols; ++i)
		{			
			long double sum = 0;
			cv::Mat gaussFilter;
            /*
                Setting up the right gaussian kernel by size and sigma value acording to "x" parameter
            */
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

/*
    Function for second term, aka getting first derivates of sigmaMap
*/
cv::Mat getSecondTerm (const double* x, int rows, int cols, float ypsilon) {
	cv::Mat a;
	a = cv::Mat::zeros(cols, rows, CV_32F);

	for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
		   double b1 = 0;
		   double b2 = 0;
		   if ((j - 1 >= 0) && (j + 1 < rows)) {
               b2 = pow(x[cols * i + j - 1] - x[cols * i + j + 1], 2);   
		   }

		   if ((i - 1 >= 0) && (i + 1 < cols)) {
			   b1 = pow(x[cols * (i - 1) + j] - x[cols * (i + 1) + j], 2);
		   }
		   a.at<float>(i, j) = ypsilon*(b1 + b2);
        }
    }
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

        /*
            Creating objective function
        */
        float y = 1e-5;
        long double sum1 = 0;
        cv::Mat firstTerm = getFirstTerm(SmoothenedImage, parameters, Width, Height);
        long double sum2 = 0;
        cv::Mat secondTerm = getSecondTerm(parameters, Width, Height, y);
        for (int j = 0; j < Width; j++) {
            for (int i = 0; i < Height; i++) {
                sum1 += pow(firstTerm.at<float>(j, i) - OriginalImage.at<float>(j, i), 2);
                sum2 += secondTerm.at<float>(j, i);
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

cv::Mat optimizeForSigma(int height, int width, cv::Mat originalImage, cv::Mat smoothenedImage, int optim1Iteration) {
    // define the starting point for the optimization
    /*
        Defining randomly starting points for getting such a best solution
    */
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
    options.max_iterations = optim1Iteration; // 50
    options.dwell_iterations = optim1Iteration;
    options.max_stagnant_iterations = optim1Iteration;

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
    /*
        Bounds set up in interval <0, 10>
    */
    double upper_bounds [height*width];
    std::fill_n(upper_bounds, height*width, 10);
    double lower_bounds [height*width];
    std::fill_n(lower_bounds, height*width, 0);

    /*
        Stepsize set up to 1 for getting integer values
    */
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
    pallas::GradientProblem problem(new sigmaClass(width, height, originalImage, smoothenedImage));

    // solve the problem and store the optimal position
    // in parameters and the optimization details in
    // the summary
    pallas::Solve(options, problem, parameters, &summary);

    /*
        Adding to cv matrix
    */
    cv::Mat a;
	a = cv::Mat::zeros(height, width, CV_32F);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            a.at<float>(j, i) = parameters[width * j + i];
        }
    }

    return a;
}

/*
    Part 2: functions for maximalizing details
*/

/*
    Class for creating objective functions
*/
class detailClass : public pallas::GradientCostFunction {
public:
    int Width;
    int Height;
    cv::Mat DetailImage;
    cv::Mat Weight1;
    cv::Mat Weight2;
    std::vector<cv::Mat> BaseChannels;
    std::vector<cv::Mat> DetailChannels;

    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {

        /*
            Constraint - NOT working in pallas
        */
        
        /*for (int j = 0; j < Height; j++) {
            for (int i = 0; i < Width; i++) {
                if (-10 > (((parameters[(Width*Height) + Width * j + i]) + BaseChannels[0].at<float>(j, i)) + (parameters[Width * j + i])*DetailChannels[0].at<float>(j, i))) {
                    return false;
                } else if ((((parameters[(Width*Height) + Width * j + i]) + BaseChannels[0].at<float>(j, i)) + (parameters[Width * j + i])*DetailChannels[0].at<float>(j, i)) > 265) {
                    return false;
                }

               if (-10 > (((parameters[(Width*Height) + Width * j + i]) + BaseChannels[1].at<float>(j, i)) + (parameters[Width * j + i])*DetailChannels[1].at<float>(j, i))) {
                    return false;
                } else if ((((parameters[(Width*Height) + Width * j + i]) + BaseChannels[1].at<float>(j, i)) + (parameters[Width * j + i])*DetailChannels[1].at<float>(j, i)) > 265) {
                    return false;
                }

                if (-10 > (((parameters[(Width*Height) + Width * j + i]) + BaseChannels[2].at<float>(j, i)) + (parameters[Width * j + i])*DetailChannels[2].at<float>(j, i))) {
                    return false;
                } else if ((((parameters[(Width*Height) + Width * j + i]) + BaseChannels[2].at<float>(j, i)) + (parameters[Width * j + i])*DetailChannels[2].at<float>(j, i)) > 265) {
                    return false;
                }
            }
        }*/
        /**
         * Setting up objective function
         * /
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
                /*
                    Top-left corner
                */
                if ((i - 1 < 0 && j - 1 < 0) || (i + 1 == Width && j + 1 == Height)) {
                    nablaPowerS = 0.0;
                /*
                    Width boundary
                */
                } else if ((i - 1 < 0)||(i + 1 == Width)) {
                    nablaPowerS = pow((parameters[Width * (j - 1) + i]) + (parameters[Width * (j + 1) + i]), 2);
                /*
                    Height boundary
                */
                } else if ((j - 1 < 0) + (j + 1 == Height)) {
                    nablaPowerS = pow((parameters[Width * j + (i - 1)]) + (parameters[Width * j + (i + 1)]), 2);
                } else {
                    nablaPowerS = pow((parameters[Width * j + (i - 1)]) + (parameters[Width * j + (i + 1)]), 2) +
                                  pow((parameters[Width * (j - 1) + i]) + (parameters[Width * (j + 1) + i]), 2);
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
                /*
                    Top-left corner
                */
                if ((i - 1 < 0 && j - 1 < 0) || (i + 1 == Width && j + 1 == Height)) {
                    nablaPowerT = 0.0;
                /*
                    Width boundary
                */
                } else if ((i - 1 < 0)||(i + 1 == Width)) {
                    nablaPowerT = pow((parameters[(Width*Height) + Width * (j - 1) + i]) + (parameters[(Width*Height) + Width * (j + 1) + i]), 2);
                /*
                    Height boundary
                */
                } else if ((j - 1 < 0) + (j + 1 == Height)) {
                    nablaPowerT = pow((parameters[(Width*Height) + Width * j + (i - 1)]) + (parameters[(Width*Height) + Width * j + (i + 1)]), 2);
                } else {
                    nablaPowerT = pow((parameters[(Width*Height) + Width * j + (i - 1)]) + (parameters[(Width*Height) + Width * j + (i + 1)]), 2) +
                                  pow((parameters[(Width*Height) + Width * (j - 1) + i]) + (parameters[(Width*Height) + Width * (j + 1) + i]), 2);
                }

                sum3 += Weight2.at<float>(j, i)*nablaPowerT;
            }
        }
        cost[0] = sum1 + sum2*200 + sum3*500;// sum1 + sum2 * 200; //+ 200*sum2;// + sum2;

        return true;
    }

    virtual int NumParameters() const { 
        return Width * Height * 2; 
    }
    detailClass(int width, int height, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels);
};

detailClass::detailClass(int width, int height, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels) {
    this->Width = width;
    this->Height = height;
    this->DetailImage = detailImage;
    this->Weight1 = weight1;
    this->Weight2 = weight2;
    this->BaseChannels = baseChannels;
    this->DetailChannels = detailChannels;
}

/*
    Detail maximalization
*/
std::vector<cv::Mat> optimizeForGettingSAndTparameters(int height, int width, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels, int optim2Iteration) {
    //google::InstallFailureSignalHandler();
    // define the starting point for the optimization as much similar to constraint
    srand ( time(NULL) );
    double parameters[height*width*2];
    std::cout << "generating random values" << std::endl;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            while (true) {
                bool condition1 = false;
                bool condition2 = false;
                bool condition3 = false;

                parameters[(width*height) + width * j + i] = (rand() % 4 - 2) / (float)(rand() % 2 + 1);
                parameters[width * j + i] = (rand() % 4 - 2) / (float)(rand() % 2 + 1); // s
                if (0 <= (((parameters[(width*height) + width * j + i]) + baseChannels[0].at<float>(j, i)) + (parameters[width * j + i])*detailChannels[0].at<float>(j, i))) {
                    condition1 = true;
                } else if (255 < (((parameters[(width*height) + width * j + i]) + baseChannels[0].at<float>(j, i)) + (parameters[width * j + i])*detailChannels[0].at<float>(j, i))) {
                    condition1 = false;
                }

                if (0 <= (((parameters[(width*height) + width * j + i]) + baseChannels[1].at<float>(j, i)) + (parameters[width * j + i])*detailChannels[1].at<float>(j, i))) {
                    condition2 = true;
                } else if (255 < (((parameters[(width*height) + width * j + i]) + baseChannels[1].at<float>(j, i)) + (parameters[width * j + i])*detailChannels[1].at<float>(j, i))) {
                    condition2 = false;
                }

                if (0 <= (((parameters[(width*height) + width * j + i]) + baseChannels[2].at<float>(j, i)) + (parameters[width * j + i])*detailChannels[2].at<float>(j, i))) {
                    condition3 = true;
                } else if (255 < (((parameters[(width*height) + width * j + i]) + baseChannels[2].at<float>(j, i)) + (parameters[width * j + i])*detailChannels[2].at<float>(j, i))) {
                    condition3 = false;
                }

                if (condition1 == true && condition2 == true && condition3 == true) {
                    // std::cout << parameters[(width*height) + width * j + i] + "" << std::endl;
                    break;
                }
            }
        }
    }

    /*for (int i = 0; i < height*width*2; i++) {
        parameters[i] = 0;
    }*/
    std::cout << "random values generated" << std::endl;
    pallas::SimulatedAnnealing::Options options;

    options.max_iterations = optim2Iteration;// 250
    options.dwell_iterations = optim2Iteration;
    options.max_stagnant_iterations = optim2Iteration;

    options.cooling_schedule_options.initial_temperature = 1000;

   // options.minimum_cost = 0.1;

    double upper_bounds [height*width*2];
    std::fill_n(upper_bounds, height*width*2, 3);
    double lower_bounds [height*width*2];
    std::fill_n(lower_bounds, height*width*2, -3);

    double step_size = 0.05;
    pallas::scoped_ptr<pallas::StepFunction> step_function (new pallas::BoundedStepFunction(step_size,
                                                                                            upper_bounds,
                                                                                            lower_bounds,
                                                                                            height*width*2));
    options.set_step_function(step_function);
    pallas::SimulatedAnnealing::Summary summary;
    pallas::GradientProblem problem(new detailClass(width, height, detailImage, weight1, weight2, baseChannels, detailChannels));

    pallas::Solve(options, problem, parameters, &summary);
    std::cout << summary.FullReport() << std::endl;
    std::vector<cv::Mat> sAndT;
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