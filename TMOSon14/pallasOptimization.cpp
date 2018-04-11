#include "pallasOptimization.h"

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
    options.max_iterations = 1; // 50
    options.dwell_iterations = 1;
    options.max_stagnant_iterations = 1;

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