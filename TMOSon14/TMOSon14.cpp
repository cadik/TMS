/* --------------------------------------------------------------------------- *
 * TMOSon14.cpp: implementation of the TMOSon14 class.   *
 * --------------------------------------------------------------------------- */
#include "TMOSon14.h"
#include <fftw3.h>
//#include <opencv2/core/core.hpp>


#include <opencv2/objdetect/objdetect.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
// #include "L0minimization.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;
/* --------------------------------------------------------------------------- *
 * Constructor se<ves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOSon14::TMOSon14()
{
	SetName(L"Son14");						// TODO - Insert operator name
	SetDescription(L"Art-Photography detail enhancement");	// TODO - Insert description
	/**
	 * Kappa - Parameter
	 * */
	kappa.SetName(L"Kappa");				// TODO - Insert parameters names
	kappa.SetDescription(L"Represents rate Kappa for modified l0 smoothing");	// TODO - Insert parameter descriptions
	kappa.SetDefault(2);							// TODO - Add default values
	kappa=2.;
	kappa.SetRange(-10.0,10.0);				// TODO - Add acceptable range if needed
	this->Register(kappa);
	
}

TMOSon14::~TMOSon14()
{
}

/**
 * From L0-smoothing variables and fuctions
 */

// optimization params

// buffers for solving linear system
Eigen::SparseMatrix<float> A0, E;
Eigen::SparseMatrix<float> GX, GY;
Eigen::VectorXf S_vec, I_vec, H_vec, V_vec;

float lambda = 0.01;
float beta0 = 2*lambda;
float beta_max = 10000;
float kappa = 2.0;
bool exact = true;
int iter_max = 1000;
// functions
void buildGradientMatrix(Eigen::SparseMatrix<float> &G, 
                         const int rows,
                         const int cols,
                         const std::vector<std::pair<int, float> > x_indices, 
                         const std::vector<std::pair<int, float> > y_indices
                         )
{
    int num_of_variables = rows*cols;
    std::vector<Eigen::Triplet<float> > coeffcients;
    bool compute_x = x_indices.empty() ? false : true;
    bool compute_y = y_indices.empty() ? false : true;

    G = Eigen::SparseMatrix<float>(num_of_variables, num_of_variables);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            int index = i*cols+j;
            int n_index, num_indices;
            if(compute_x){
                num_indices = x_indices.size();
                for(int k=0; k<num_indices; k++){
                    n_index = index + x_indices[k].first;                    
                    if(n_index >= num_of_variables){
                        continue;
                    }
                    coeffcients.push_back(Eigen::Triplet<float>(index, n_index, x_indices[k].second));
                }
            }
            if(compute_y){
                num_indices = y_indices.size();
                for(int k=0; k<num_indices; k++){                    
                    n_index = (i+y_indices[k].first)*cols+j;
                    if(n_index >= num_of_variables){
                        continue;
                    }
                    coeffcients.push_back(Eigen::Triplet<float>(index, n_index, y_indices[k].second));
                }                    
            }
        }
    }
    G.setFromTriplets(coeffcients.begin(), coeffcients.end());
}

void constructSparseIdentityMatrix(Eigen::SparseMatrix<float> &mat, const int &num_of_variables){
    mat = Eigen::SparseMatrix<float>(num_of_variables, num_of_variables);
    std::vector<Eigen::Triplet<float> > coeffcients;
    for(int i=0; i<num_of_variables; i++){
        coeffcients.push_back(Eigen::Triplet<float>(i, i, 1.0f));        
    }
    mat.setFromTriplets(coeffcients.begin(), coeffcients.end());
}

void init(const int &rows, const int &cols)
{
    int num_of_variables = rows*cols;

    // build gradient matrix
    std::vector<std::pair<int, float> > indices;
    indices.push_back(std::pair<int, float>(0, 1.0f));
    indices.push_back(std::pair<int, float>(1, -1.0f));
    buildGradientMatrix(GX, rows, cols, indices, std::vector<std::pair<int, float> >());
    buildGradientMatrix(GY, rows, cols, std::vector<std::pair<int, float> >(), indices);

    A0 = (GX.transpose()*GX+GY.transpose()*GY);

    constructSparseIdentityMatrix(E, num_of_variables);    

    S_vec = Eigen::VectorXf::Zero(rows*cols);    
    I_vec = Eigen::VectorXf::Zero(rows*cols);        
    H_vec = Eigen::VectorXf::Zero(rows*cols);    
    V_vec = Eigen::VectorXf::Zero(rows*cols);        
}

void vec2CvMat(const Eigen::VectorXf &vec, cv::Mat &mat, const int &rows, const int &cols){
    for(int i=0; i<rows; i++){
        float *ptr = reinterpret_cast<float*>(mat.data+mat.step*i);
        for(int j=0; j<cols; j++){
            *ptr = vec[i*cols+j];
            ++ptr;
        }
    }   
}
void cvMat2Vec(const cv::Mat &mat, Eigen::VectorXf &vec){
    int rows = mat.rows;
    int cols = mat.cols;

    for(int i=0; i<rows; i++){
        float *ptr = reinterpret_cast<float*>(mat.data+mat.step*i);
        for(int j=0; j<cols; j++){
            vec[i*cols+j] = *ptr;
            ++ptr;
        }
    }    
}

void computeGradient(const cv::Mat &mat, cv::Mat &grad_x, cv::Mat &grad_y){
    int rows = mat.rows;
    int cols = mat.cols;  

    for(int i=0; i<rows-1; i++){ 
        float *ptr = reinterpret_cast<float*>(mat.data+mat.step*i);
        float *n_ptr = reinterpret_cast<float*>(mat.data+mat.step*(i+1));
        float *gx_ptr = reinterpret_cast<float*>(grad_x.data+grad_x.step*i);
        float *gy_ptr = reinterpret_cast<float*>(grad_y.data+grad_y.step*i);
        for(int j=0; j<cols-1; j++){        
            *gx_ptr = *ptr - *(ptr+1);
            *gy_ptr = *ptr - *n_ptr;
            ++ptr;
            ++n_ptr;
            ++gx_ptr;
            ++gy_ptr;
        }
    }
}

void computeS(cv::Mat &S, 
              const cv::Mat &I,
              const cv::Mat &H,
              const cv::Mat &V,
              const float &beta)
{
    int rows = S.rows;
    int cols = S.cols;

    //boost::timer t;

    cvMat2Vec(I, I_vec);
    cvMat2Vec(H, H_vec);
    cvMat2Vec(V, V_vec);

    //std::cout << "\t\t mat2vec " << t.elapsed() << " sec" << std::endl;    
    //t.restart();

    // build linear system As=b
    Eigen::SparseMatrix<float> A = beta*A0 + E;
    Eigen::VectorXf b = I_vec + beta*(GX.transpose()*H_vec+GY.transpose()*V_vec);

    // solve linear system
    if(exact){
        Eigen::SimplicialLLT<Eigen::SparseMatrix<float> > solver;
        solver.compute(A);
        if(solver.info()!=Eigen::Success) {
            std::cout << "decomposition failed" << std::endl;
        }    
        S_vec = solver.solve(b);
    }
    else{
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
        S_vec = solver.compute(A).solve(b);        
    }
    //std::cout << "\t\t solve linear system " << t.elapsed() << " sec" << std::endl;
    //t.restart();    

    // update S
    vec2CvMat(S_vec, S, rows, cols);
    //std::cout << "\t\t vec2mat " << t.elapsed() << " sec" << std::endl;        
}

void optimize(cv::Mat &S, 
              const cv::Mat &I, 
              cv::Mat &H, 
              cv::Mat &V, 
              cv::Mat &grad_x,
              cv::Mat &grad_y,
              float &beta)
{
    int rows = S.rows;
    int cols = S.cols;

    //boost::timer t;

    // Compute Gradient
    computeGradient(S, grad_x, grad_y);
    // std::cout << "\t compute gradient " << t.elapsed() << " sec" << std::endl;
    // t.restart();

    // Computing h, v
    for(int j=0; j<rows; j++){
        for(int i=0; i<cols; i++){
            float gx = grad_x.at<float>(j, i);
            float gy = grad_y.at<float>(j, i);
            float val = gx*gx + gy*gy;        

            if(val < lambda/beta){
                H.at<float>(j, i) = V.at<float>(j, i) = 0;
            }
            else{          
                H.at<float>(j, i) = gx;
                V.at<float>(j, i) = gy;
            }      
        }            
    }
    // std::cout << "\t compute h, v " << t.elapsed() << " sec" << std::endl;    
    // t.restart();    

    // Computing S
    computeS(S, I, H, V, beta);

    //std::cout << "\t compute S " << t.elapsed() << " sec" << std::endl;    
}

std::vector<cv::Mat> minimizeL0Gradient(const cv::Mat &src){
	ofstream myfile;

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
    init(rows, cols);
	myfile.open ("exampssle.txt");
	myfile << "kokot";
	myfile.close();
    // main loop
    while(beta < beta_max){
		myfile.open ("exampssle.txt");
		myfile << beta;
		myfile.close();
        //boost::timer t;
        // minimize L0 gradient
        for(int i=0; i<num_of_channels; i++){
            optimize(S_channels[i], I_channels[i], H, V, grad_x, grad_y, beta);
        }
        // Update param
        beta = beta*kappa;
        std::cout << "iteration #" << count++ << " beta: " << beta << std::endl;

        for(int i=0; i<num_of_channels; i++){
            cv::convertScaleAbs(S_channels[i], S_U8_channels[i], 255.0);
        }        
        cv::merge(S_U8_channels, S);        
        S_mats.push_back(S.clone());
        if(count >= iter_max){
            break;
        }
        //std::cout << "iteration: " << t.elapsed() << " sec" << std::endl;
    }
    return S_mats;
}
/*
	end of L0-smothing
*/
/*
    My functions
*/
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
/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOSon14::Transform()
{
	ofstream myfile;

	myfile.open ("exampsssle.txt");
	myfile << "1Writing this to a file.\n";
	myfile.close();
	
	// Initialy images are in RGB format, but you can 
	// convert it into other format
	//pSrc->Convert(TMO_Yxy);								// This is format of Y as luminance
	//pDst->Convert(TMO_Yxy);								// x, y as color information

	double* pSourceData = pSrc->GetData();				// You can work at low level data
	double* pDestinationData = pDst->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components
	// double r, g, b;
    
	cv::Mat r;
	cv::Mat g;
	cv::Mat b;
	
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	
	r = cv::Mat::zeros (height, width, CV_64FC1);
	g = cv::Mat::zeros (height, width, CV_64FC1);
	b = cv::Mat::zeros (height, width, CV_64FC1);
	// cv::Mat pSourceData1 = pSrc->GetData();
	/**
	 * Detail/Base decomposition 
	 */
	/**
	 * Base decomposition, fase 1 (original L0-Smoothing) 
	 */
	// std::cout 
	// cv::Mat img = cv::imread(pSrc);
	// std::vector<cv::Mat> results = minimizeL0Gradient(pSrc);
    // std::vector<cv::Mat> results = minimizeL0Gradient();
	for (int j = 0; j < height; j++)
	{
		pSrc->ProgressBar(j, height);
		for (int i = 0; i < width; i++)
		{
			r.at<double>(j, i) = ceil(255 * *pSourceData++); 
			g.at<double>(j, i) = ceil(255 * *pSourceData++);  //getting separate RGB channels
			b.at<double>(j, i) = ceil(255 * *pSourceData++);
		}
	}
	myfile.open ("exampsssle.txt");
	myfile << "2Writing this to a file.\n";
	myfile.close();
	// std::cout << r << std::endl;
	std::vector<cv::Mat> array_to_merge;

    array_to_merge.push_back(b);
    array_to_merge.push_back(g);
    array_to_merge.push_back(r);

    cv::Mat color;

    cv::merge(array_to_merge, color);
    myfile.open ("color.txt");
	myfile << color;
	myfile.close();
	std::vector<cv::Mat> results = minimizeL0Gradient(color);

	cv::Mat src_channels[3];
    cv::split(results[(int)results.size() - 1], src_channels);
	// int j;
    myfile.open ("color2.txt");
	myfile << results[(int)results.size() - 1];
	myfile.close();

    cv::Mat pgm_double;
    pgm_double = (results[(int)results.size() - 1]);


    cv::Mat chan[3];
    cv::split(pgm_double, chan);

    myfile.open ("color3.txt");
	myfile << pgm_double;
	myfile.close();

    myfile.open ("color31.txt");
	myfile << chan[0];
	myfile.close();

    myfile.open ("color32.txt");
	myfile << chan[1];
	myfile.close();

    myfile.open ("color33.txt");
	myfile << chan[2];
	myfile.close();

	for (int j = 0; j < height; j++)
	{
		// int i;	// You can provide progress bar
		for (int i = 0; i < width; i++)
		{
            *pDestinationData++ = chan[0].at<float>(j, i) * 1.0/255.0;
			*pDestinationData++ = chan[1].at<float>(j, i) * 1.0/255.0;
			*pDestinationData++ = chan[2].at<float>(j, i) * 1.0/255.0;
		}
	}

    
	// pDst->Convert(TMO_RGB);
	return 0;
}

