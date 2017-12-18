/****************************************************************************************
 * 																						*		
 * This is a C++ implementation of "Image Smoothing via L0 Gradient Minimization", 		*
 * Li Xu, Cewu Lu, Yi Xu, Jiaya Jia, SIGGRAPH ASIA 2011.								*
 * **************************************************************************************
 * AUTHORS OF CODE: https://github.com/daikiyamanaka									*
 * EDIT: values, functions and libraries including 										*
 * boost" commented because of not neccesary using (better compiling)					*
 * GITHUB LINK: https://github.com/daikiyamanaka/L0-gradient-smoothing					*
 * 																						*	
 ****************************************************************************************/
#include "L0minimization.h"

float lambda = 0.01;
float beta0 = lambda * 2;
float beta_max = 10000;
float kappa = 2.0;
bool exact = false;
int iter_max = 1000;

Eigen::SparseMatrix<float> A0, E;
Eigen::SparseMatrix<float> GX, GY;
Eigen::VectorXf S_vec, I_vec, H_vec, V_vec;

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

    // main loop
    while(beta < beta_max){
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
