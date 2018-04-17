#include "baseAndDetailDecomposition.h"

/*
 * L0 Smoothing phase1 (for getting only one image)
 **/
cv::Mat minimizeL0Gradient1(const cv::Mat &src){
    int rows = src.rows;
    int cols = src.cols;
    std::vector<cv::Mat> src_channels;
    cv::split(src, src_channels);

    int num_of_channels = src_channels.size();    
    std::vector<cv::Mat> S_channels(num_of_channels), I_channels(num_of_channels), S_U8_channels(num_of_channels);
    for(int i=0; i<num_of_channels; i++){
        src_channels[i].convertTo(I_channels[i], CV_32FC1);
        I_channels[i] *= 1./255;
        I_channels[i].copyTo(S_channels[i]);            
    }

    cv::Mat S, H, V, grad_x, grad_y;
    std::vector<cv::Mat> S_mats;
    float beta = beta0;   
    S = cv::Mat(rows, cols, CV_32FC1);
    H = cv::Mat(rows, cols, CV_32FC1);
    V = cv::Mat(rows, cols, CV_32FC1);
    grad_x = cv::Mat::zeros(rows, cols, CV_32FC1);
    grad_y = cv::Mat::zeros(rows, cols, CV_32FC1);      
    init(rows, cols);

    while(beta < beta_max){

        for(int i=0; i<num_of_channels; i++){
            optimize(S_channels[i], I_channels[i], H, V, grad_x, grad_y, beta);
        }

        beta = beta*kappa;

        for(int i=0; i<num_of_channels; i++){
            cv::convertScaleAbs(S_channels[i], S_U8_channels[i], 255.0);
        }        
        cv::merge(S_U8_channels, S);  

		S_mats.push_back(S.clone());
    }
    return S;
}

/*
 * Function for getting gradient magnitude from image 
 **/
cv::Mat getGradientMagnitude(const cv::Mat &src)
{
	int rows = src.rows;
    int cols = src.cols;
    
    std::vector<cv::Mat> src_channels;
    cv::split(src, src_channels);

    int num_of_channels = src_channels.size();    
    std::vector<cv::Mat> S_channels(num_of_channels);
    std::vector<cv::Mat> I_channels(num_of_channels);

    cv::Mat S, grad_x, grad_y, gradient;
    std::vector<cv::Mat> S_mats;    
    S = cv::Mat(rows, cols, CV_32FC1);
    grad_x = cv::Mat::zeros(rows, cols, CV_32FC1);
    grad_y = cv::Mat::zeros(rows, cols, CV_32FC1); 
    gradient = cv::Mat::zeros(rows, cols, CV_32FC1);  

    for(int i=0; i<num_of_channels; i++){
        src_channels[i].convertTo(I_channels[i], CV_32FC1);
        I_channels[i] *= 1./255;
        I_channels[i].copyTo(S_channels[i]);
        computeGradient(S_channels[i], grad_x, grad_y);
        for(int j=0; j<rows; j++){
            for(int i=0; i<cols; i++){
                float gx = grad_x.at<float>(j, i);
                float gy = grad_y.at<float>(j, i);
                gradient.at<float>(j, i) += sqrt(pow(gx, 2) + pow(gy, 2));
            }
        }
    }
    return gradient;
}

/*
 * get gradient from one layer
 **/
 cv::Mat getGradientMagnitudeFromOneLayer(const cv::Mat &src)
{
	int rows = src.rows;
    int cols = src.cols;
    

    cv::Mat S, grad_x, grad_y, gradient;   
    S = cv::Mat(rows, cols, CV_32FC1);
    grad_x = cv::Mat::zeros(rows, cols, CV_32FC1);
    grad_y = cv::Mat::zeros(rows, cols, CV_32FC1); 
    gradient = cv::Mat::zeros(rows, cols, CV_32FC1);  


	computeGradient(src, grad_x, grad_y);
	for(int j=0; j<rows; j++){
		for(int i=0; i<cols; i++){
			float gx = grad_x.at<float>(j, i);
			float gy = grad_y.at<float>(j, i);
			gradient.at<float>(j, i) += sqrt(pow(gx, 2) + pow(gy, 2));
		}
	}
    return gradient;
}

/*
 * Function for getting adaptive lambda matrix from gradient magnitude 
 **/
cv::Mat getAdaptiveLambdaMatrix(const cv::Mat &gradient, 
								int rows, 
								int cols) 
{
    cv::Mat adaptiveLambdaMatrix;
    adaptiveLambdaMatrix = cv::Mat(rows, cols, CV_32FC1);

    double a = 0.2;
    double sigma = 0.1;
    double epsilon = 0.00001;
    

    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            /*
                Bisquare function
            */
            double tmp = 0.0;
            double u = gradient.at<float>(j, i) - a;
            if ( u < -sigma) {
                tmp = 1/3;
            } else if ((-sigma <= u) && (u < 0.0)) {
                tmp = ((pow(u, 2))/(pow(sigma, 2))) - ((pow(u, 4))/(pow(sigma, 4))) + ((pow(u, 6))/(3*(pow(sigma, 6))));
            } else if (u >= 0.0) {
                tmp = 0;
            }
            /*
				Final count of adaptiveLambdaMatrix
            */
            adaptiveLambdaMatrix.at<float>(j, i) = 3*(0.01 - epsilon)*tmp + epsilon;
        }
    } 
    return adaptiveLambdaMatrix;
}

/*
 * Function for optimizing image with adaptive lambda matrix
 **/
void optimizeWithAdaptiveLambdaMatrix(cv::Mat &S, 
              const cv::Mat &I, 
              cv::Mat &H, 
              cv::Mat &V, 
              cv::Mat &grad_x,
              cv::Mat &grad_y,
              cv::Mat &lambdaMatrix,
              float &beta)
{
    int rows = S.rows;
    int cols = S.cols;

    computeGradient(S, grad_x, grad_y);

    for(int j=0; j<rows; j++){
        for(int i=0; i<cols; i++){
            float gx = grad_x.at<float>(j, i);
            float gy = grad_y.at<float>(j, i);
            float val = gx*gx + gy*gy;        

            if(val < lambdaMatrix.at<float>(j, i)/beta){
                H.at<float>(j, i) = V.at<float>(j, i) = 0;
            }
            else{          
                H.at<float>(j, i) = gx;
                V.at<float>(j, i) = gy;
            }      
        }            
    }  

    computeS(S, I, H, V, beta);
}

/*
 * Function for minimizing image with adaptive Lambda matrix
 **/
cv::Mat minimizeL0GradientSecondFaze(const cv::Mat &src, cv::Mat lambdaMatrix1, int rows, int cols){
    std::vector<cv::Mat> src_channels;
    cv::split(src, src_channels);

    int num_of_channels = src_channels.size();    
    std::vector<cv::Mat> S_channels(num_of_channels), I_channels(num_of_channels), S_U8_channels(num_of_channels);
    for(int i=0; i<num_of_channels; i++){
        src_channels[i].convertTo(I_channels[i], CV_32FC1);
        I_channels[i] *= 1./255;
        I_channels[i].copyTo(S_channels[i]);            
    }

    cv::Mat S, H, V, grad_x, grad_y;
    std::vector<cv::Mat> S_mats;
    float beta = beta0;
    S = cv::Mat(rows, cols, CV_32FC1);
    H = cv::Mat(rows, cols, CV_32FC1);
    V = cv::Mat(rows, cols, CV_32FC1);
    grad_x = cv::Mat::zeros(rows, cols, CV_32FC1);
    grad_y = cv::Mat::zeros(rows, cols, CV_32FC1);   

    init(rows, cols);

    do {

        for(int i=0; i<num_of_channels; i++){
            optimizeWithAdaptiveLambdaMatrix(S_channels[i], I_channels[i], H, V, grad_x, grad_y, lambdaMatrix1, beta);
        }

        beta = beta*kappa;

        for(int i=0; i<num_of_channels; i++){
            cv::convertScaleAbs(S_channels[i], S_U8_channels[i], 255.0);
        }
        cv::merge(S_U8_channels, S); 
        if (beta >= beta_max) {
            S_mats.push_back(S.clone()); 
        }
    } while (beta < beta_max);
    
    return S;
}

/*
 * Function for getting detail layer from base layer
 **/
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
/*
 * Function for detail maximilization 
 **/
 
std::vector<cv::Mat> detailMaximalization(const cv::Mat &base, const cv::Mat &detailSum, const cv::Mat &r1Weight, const cv::Mat &r2Weight,  int rows, int cols, int counter, std::vector<cv::Mat> detail) 
{
	srand (time(NULL));
    std::vector<cv::Mat> S_mats;
    
    long double tmpSum = 1e45;
    
    int r1 = 200;
    int r2 = 500;
    
    cv::Mat S;
	cv::Mat T;
    
	cv::Mat	tmpS = cv::Mat::zeros(rows, cols, CV_32FC1);
	cv::Mat	tmpT = cv::Mat::zeros(rows, cols, CV_32FC1);
    for (int a = 0; a < counter; a++) {
		/*
		 * Getting random layer
		 * */
		long double sumR1 = 0;
		long double sumR2 = 0;
		S = cv::Mat::zeros(rows, cols, CV_32FC1);
		T = cv::Mat::zeros(rows, cols, CV_32FC1);
		long double sum = 0;
		for(int j = 0; j < rows; j++){
			for(int i = 0; i < cols; i++){
				/*
				 * Constraint condition 
				 **/
				while(true) {
					double aS = rand() % 10000;
					double bS = rand() % 10000;
					double cS = rand() % 1000;
					
					S.at<float>(j, i) = ((double)aS - (double)bS) / (double)cS;
					
					double aT = rand() % 10000;
					double bT = rand() % 10000;
					double cT = rand() % 1000;
					
					T.at<float>(j, i) = ((double)aT - (double)bT) / (double)cT;
					
					bool cond1 = false;
					bool cond2 = false;
					if ((base.at<float>(j,i) + T.at<float>(j, i) + S.at<float>(j, i) * (detail[0]).at<float>(j, i)) >= 0 && (base.at<float>(j,i) + T.at<float>(j, i) + S.at<float>(j, i) * (detail[0]).at<float>(j, i)) <= 1) {
						cond1 = true;
					}
					
					if ((cond1)&&(base.at<float>(j,i) + T.at<float>(j, i) + S.at<float>(j, i) * (detail[1]).at<float>(j, i)) >= 0 && (base.at<float>(j,i) + T.at<float>(j, i) + S.at<float>(j, i) * (detail[1]).at<float>(j, i)) <= 1) {
						cond2 = true;
					}
					
					if ((cond1)&&(cond2)&&(base.at<float>(j,i) + T.at<float>(j, i) + S.at<float>(j, i) * (detail[2]).at<float>(j, i)) >= 0 && (base.at<float>(j,i) + T.at<float>(j, i) + S.at<float>(j, i) * (detail[2]).at<float>(j, i)) <= 1) {
						break;
					}
					
					// break;
				}				
			}
		}
		/*
		 * getting gradients of T and S
		 * */
		cv::Mat gradientT = getGradientMagnitudeFromOneLayer(T);
		cv::Mat gradientS = getGradientMagnitudeFromOneLayer(S);
		
		for(int jj = 0; jj < rows; jj++){
			for(int ii = 0; ii < cols; ii++){
				sum -= pow(S.at<float>(jj, ii) * detailSum.at<float>(jj, ii), 2);
				sumR1 += r1Weight.at<float>(jj, ii) * pow(gradientS.at<float>(jj, ii), 2);
				sumR2 += r2Weight.at<float>(jj, ii) * pow(gradientT.at<float>(jj, ii), 2);
			}
		}
		
		sumR1 *= r1;
		sumR2 *= r2;
		
		if ((sum + sumR1 + sumR2) < tmpSum) {
			tmpSum = sum + sumR1 + sumR2;
			tmpS = S;
			tmpT = T;
		}
	}
    S_mats.push_back(tmpS.clone());
    S_mats.push_back(tmpT.clone());
    
    return S_mats;
}

/*
 * Function for detail control
 **/
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

/*
 * Function for finding parameters from objective function for getting sigma map
 **/
cv::Mat stochasticOptimizationForGetSigma(cv::Mat base, cv::Mat original, int rows, int cols, int counter) 
{
		float y = 1e-5;
		srand (time(NULL));
		
		cv::Mat randomLayer;		
		randomLayer = cv::Mat::zeros(rows, cols, CV_32F);
		
		cv::Mat tmpLayer;		
		tmpLayer = cv::Mat::zeros(rows, cols, CV_32F);
		
		long double tmpSum = 1e30;
		
		for (int a = 0; a < counter; a++) {
			cv::Mat blurredImage;		
			blurredImage = cv::Mat::zeros(rows, cols, CV_32F);
			/*
			 * Getting random layer with discrete values
			 **/
			for(int j=0; j<rows; j++){
				for(int i=0; i<cols; i++){
					randomLayer.at<float>(j, i) = rand() % 10;
				}
			}
			
			/*
			 * Getting blurred version of image
			 **/
			cv::Mat rSmoothed = myOwn2DFilter(base, randomLayer, rows, cols);
			
			/*
			 * Getting derivations
			 **/
			
			cv::Mat randomLayerMagnitude = getGradientMagnitudeFromOneLayer(randomLayer);
			
			/*
			 * Minimalization
			 * */
			long double sum = 0.0;
			 
			for(int j=0; j<rows; j++){
				for(int i=0; i<cols; i++){
					sum += pow(rSmoothed.at<float>(j, i) - original.at<float>(j, i), 2) + y * pow(randomLayerMagnitude.at<float>(j, i), 2);
				}
			}
			
			if (sum < tmpSum) {
				tmpSum = sum;
				tmpLayer = randomLayer;
			}
			
		}
			
        return tmpLayer;
}

/*
 * Function for getting sum of costs for image optimizing
 **/
cv::Mat getSumOfCosts(cv::Mat r, cv::Mat g, cv::Mat b, int rows, int cols) 
{
	cv::Mat sum;
	sum = cv::Mat::zeros(rows,cols, CV_32F);
	
	for (int j = 0; j < rows; ++j) 
	{
		for (int i = 0; i < cols; ++i) 
		{
			sum.at<float>(j, i) = r.at<float>(j, i) + g.at<float>(j, i) + b.at<float>(j, i); 
		}
	}
	
	return sum;
}

/*
 * Convolution 2D filter with different sigma kernel for each pixel
 **/
cv::Mat myOwn2DFilter(cv::Mat image, cv::Mat sigmaMap, int rows, int cols)
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

