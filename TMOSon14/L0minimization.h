/**
 * Zde uz jsou zahrnuty vsechny fce z originalniho l0minimalization
 */
#include <fstream>
#include <opencv2/opencv.hpp>
#include <Eigen/Sparse>

extern float lambda;
extern float beta0;
extern float beta_max;
extern float kappa;
extern bool exact;
extern int iter_max;

extern Eigen::SparseMatrix<float> A0, E;
extern Eigen::SparseMatrix<float> GX, GY;
extern Eigen::VectorXf S_vec, I_vec, H_vec, V_vec;

/**
 * functions
 */
extern void buildGradientMatrix(Eigen::SparseMatrix<float> &G, 
                         const int rows,
                         const int cols,
                         const std::vector<std::pair<int, float> > x_indices, 
                         const std::vector<std::pair<int, float> > y_indices
                         );
                         
extern void constructSparseIdentityMatrix(Eigen::SparseMatrix<float> &mat, const int &num_of_variables);

extern void init(const int &rows, const int &cols);

extern void vec2CvMat(const Eigen::VectorXf &vec, cv::Mat &mat, const int &rows, const int &cols);

extern void cvMat2Vec(const cv::Mat &mat, Eigen::VectorXf &vec);

extern void computeGradient(const cv::Mat &mat, cv::Mat &grad_x, cv::Mat &grad_y);

extern void computeS(cv::Mat &S, 
              const cv::Mat &I,
              const cv::Mat &H,
              const cv::Mat &V,
              const float &beta);

extern void optimize(cv::Mat &S, 
              const cv::Mat &I, 
              cv::Mat &H, 
              cv::Mat &V, 
              cv::Mat &grad_x,
              cv::Mat &grad_y,
              float &beta);

extern std::vector<cv::Mat> minimizeL0Gradient(const cv::Mat &src);

                        

