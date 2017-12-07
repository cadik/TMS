#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <Eigen/Sparse>
#include <math.h>       /* sqrt */
#include "L0minimization.h"

/*redef*/ //std::string input_file, out_dir;
// optimization params
/*redef*/ // float lambda;
/*redef*/ // float beta0;
/*redef*/ // float beta_max;
/*redef*/ // float kappa;
// buffers for solving linear system
/*redef*/ // Eigen::SparseMatrix<float> A0, E;
/*redef*/ // Eigen::SparseMatrix<float> GX, GY;
/*redef*/ // Eigen::VectorXf S_vec, I_vec, H_vec, V_vec;

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
    // Compute Gradient
    computeGradient(S, grad_x, grad_y);
    // Computing h, v
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
    // Computing S
    computeS(S, I, H, V, beta);
}

cv::Mat getAdaptiveLambdaMatrix(const cv::Mat &gradientFromFirstSmoothing, int rows, int cols) {

    cv::Mat adaptiveLambdaMatrix;
    adaptiveLambdaMatrix = cv::Mat(rows, cols, CV_32FC1);

    double a = 0.2;
    double sigma = 0.1;
    double epsilon = 0.00001;
    

    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            /*
                Mezivypocet, bisquare funkce
            */
            double tmp = 0.0;
            double u = gradientFromFirstSmoothing.at<float>(j, i) - a;
            // std::cout << u << std::endl;
            if ( u < -sigma) {
                tmp = 1/3;
            } else if ((-sigma <= u) && (u < 0.0)) {
                tmp = ((pow(u, 2))/(pow(sigma, 2))) - ((pow(u, 4))/(pow(sigma, 4))) + ((pow(u, 6))/(3*(pow(sigma, 6))));
            } else if (u >= 0.0) {
                tmp = 0;
            }
           // std::cout << tmp << std::endl;
            /*
            Konecny vypocet adaptivni lambdy
            */
            // std::cout << 3*(lambda - epsilon)*tmp + epsilon << std::endl;
            adaptiveLambdaMatrix.at<float>(j, i) = 3*(0.1 - epsilon)*tmp + epsilon;
        }
    } 
    // std::cout << adaptiveLambdaMatrix << std::endl;
    return adaptiveLambdaMatrix;

}

std::vector<cv::Mat> minimizeL0GradientSecondFaze(const cv::Mat &src, cv::Mat lambdaMatrix1, int rows, int cols){
    std::vector<cv::Mat> src_channels;
    cv::split(src, src_channels);

    int num_of_channels = src_channels.size();    
    std::vector<cv::Mat> S_channels(num_of_channels), I_channels(num_of_channels), S_U8_channels(num_of_channels);
    for(int i=0; i<num_of_channels; i++){
        src_channels[i].convertTo(I_channels[i], CV_32FC1);
        I_channels[i] *= 1./255;
        I_channels[i].copyTo(S_channels[i]);            
    }
    // initialize
    cv::Mat S, H, V, grad_x, grad_y;
    std::vector<cv::Mat> S_mats;
    float beta = beta0;
    int count = 0;    
    S = cv::Mat(rows, cols, CV_32FC1);
    H = cv::Mat(rows, cols, CV_32FC1);
    V = cv::Mat(rows, cols, CV_32FC1);
    grad_x = cv::Mat::zeros(rows, cols, CV_32FC1);
    grad_y = cv::Mat::zeros(rows, cols, CV_32FC1);   
    /***********************************/
    init(rows, cols);
    // main loop
    do {
        // minimize L0 gradient
        for(int i=0; i<num_of_channels; i++){
            optimizeWithAdaptiveLambdaMatrix(S_channels[i], I_channels[i], H, V, grad_x, grad_y, lambdaMatrix1, beta);
        }
        // Update param
        beta = beta*kappa;
        std::cout << "iteration #" << count++ << " beta: " << beta << std::endl;

        for(int i=0; i<num_of_channels; i++){
            cv::convertScaleAbs(S_channels[i], S_U8_channels[i], 255.0);
        }
        cv::merge(S_U8_channels, S); 
        if (beta >= beta_max) {
            S_mats.push_back(S.clone()); 
        }
    } while (beta < beta_max);
    std::cout << "complete" << std::endl;
    return S_mats;
}

cv::Mat getGradientFromFirstSmoothing(const cv::Mat &src, int rows, int cols){
    std::vector<cv::Mat> src_channels;
    cv::split(src, src_channels);

    int num_of_channels = src_channels.size();    
    std::vector<cv::Mat> S_channels(num_of_channels);
    std::vector<cv::Mat> I_channels(num_of_channels);
    // initialize
    cv::Mat S, grad_x, grad_y, gradient;
    std::vector<cv::Mat> S_mats;    
    S = cv::Mat(rows, cols, CV_32FC1);
    grad_x = cv::Mat::zeros(rows, cols, CV_32FC1);
    grad_y = cv::Mat::zeros(rows, cols, CV_32FC1); 
    gradient = cv::Mat::zeros(rows, cols, CV_32FC1);  
    /***********************************/
    init(rows, cols);
    // main loop
    // getG
    for(int i=0; i<num_of_channels; i++){
        // Compute Gradient
        src_channels[i].convertTo(I_channels[i], CV_32FC1);
        I_channels[i] *= 1./255;
        I_channels[i].copyTo(S_channels[i]);
        computeGradient(S_channels[i], grad_x, grad_y);
        // Computing h, v
        // std::cout << src_channels[i] << std::endl;
        for(int j=0; j<rows; j++){
            for(int i=0; i<cols; i++){
                float gx = grad_x.at<float>(j, i);
                float gy = grad_y.at<float>(j, i);
                gradient.at<float>(j, i) += sqrt(pow(gx, 2) + pow(gy, 2));
            }
        }
    }
    // std::cout << gradient << std::endl;
    return gradient;
}

int main(){
    // parse user input
    lambda = 0.01;
    beta0 = 2*lambda;
    beta_max = 10000;
    kappa = 2.0;
    // read input image
    /*
    tady si budu menit input file a out dir
    */
    input_file = "./data/obrazky/road_big.jpg";
    out_dir = "./data/obrazky/";
    cv::Mat img = cv::imread(input_file, 1);
    if(img.empty()){
        std::cout << "can't read input image " << std::endl;
        return -1;
    }
    int rows = img.rows;
    int cols = img.cols;
    /*
    vytvorim labda matici pro zacatek
    */
    std::stringstream ss;
    
   /* cv::Mat lambdaMatrix1;
    lambdaMatrix1 = cv::Mat(rows, cols, CV_32FC1);
    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            lambdaMatrix1.at<float>(j, i) = lambda;
        }
    } 
    // L0 gradient minimization, PRVNI PRUCHOD
    std::vector<cv::Mat> results = minimizeL0Gradient(img, lambdaMatrix1, rows, cols);
    
    std::string first_l0_smoothing_picture = "road_big_smooth.png";
    ss << out_dir << first_l0_smoothing_picture;
    cv::imwrite(ss.str(), results[0]);*/
    
    cv::Mat img2 = cv::imread("./data/obrazky/road_big_smooth.png", 1);
    if(img2.empty()){
        std::cout << "can't read input image " << std::endl;
        return -1;
    }
    // vytvoreni lambda mapy + gradientu
    cv::Mat gradientFrom1stSmoothing = getGradientFromFirstSmoothing(img2, rows, cols);
    // std::cout << gradientFrom1stSmoothing << std::endl;
    // std::cout << gradientFrom1stSmoothing << std::endl;
    /*//std::cout << "bdssd" << std::endl;
    cv::FileStorage file("some_name.txt", cv::FileStorage::WRITE);
    //std::cout << "bdskbsfd" << std::endl;
    // Write to file!
    file << "neco" << gradientFrom1stSmoothing;
    cv::Mat adaptiveLambdaMatrix1 = getAdaptiveLambdaMatrix(gradientFrom1stSmoothing, rows, cols);
    // L0 gradient minimization, DRUHY PRUCHOD
    // std::cout << adaptiveLambdaMatrix1 << std::endl;
    std::vector<cv::Mat> results1 = minimizeL0Gradient(img, adaptiveLambdaMatrix1, rows, cols);
    std::string second_l0_smoothing_picture = "road_big_smooth_.png";
    ss << out_dir << second_l0_smoothing_picture;
    cv::imwrite(ss.str(), results1[0]);*/
    
    
    return 0;
}
