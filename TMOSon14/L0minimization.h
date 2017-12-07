/**
 * Zde uz jsou zahrnuty vsechny fce z originalniho l0minimalization
 */
#include <fstream>
#include <opencv2/opencv.hpp>
// #include <boost/filesystem.hpp>
// #include <boost/program_options.hpp>
// #include <boost/timer.hpp>
#include <Eigen/Sparse>

// namespace po = boost::program_options;
// namespace fs = boost::filesystem;

// std::string input_file, out_dir, config_file;
float lambda;
float beta0;
float beta_max;
float kappa;
bool exact;
int iter_max;
// buffers for solving linear system
/*Eigen::SparseMatrix<float> A0, E;
Eigen::SparseMatrix<float> GX, GY;
Eigen::VectorXf S_vec, I_vec, H_vec, V_vec;*/

void parseCommandLine(int argc, char** argv);

void parseConfigFile(const std::string &config_filename);

void saveConfigFile(const std::string filename);

void buildGradientMatrix(Eigen::SparseMatrix<float> &G, 
                         const int rows,
                         const int cols,
                         const std::vector<std::pair<int, float> > x_indices, 
                         const std::vector<std::pair<int, float> > y_indices
                         );

void constructSparseIdentityMatrix(Eigen::SparseMatrix<float> &mat, const int &num_of_variables);

void init(const int &rows, const int &cols);

void vec2CvMat(const Eigen::VectorXf &vec, cv::Mat &mat, const int &rows, const int &cols);

void cvMat2Vec(const cv::Mat &mat, Eigen::VectorXf &vec);

void computeGradient(const cv::Mat &mat, cv::Mat &grad_x, cv::Mat &grad_y);

void computeS(cv::Mat &S, 
              const cv::Mat &I,
              const cv::Mat &H,
              const cv::Mat &V,
              const float &beta);

void optimize(cv::Mat &S, 
              const cv::Mat &I, 
              cv::Mat &H, 
              cv::Mat &V, 
              cv::Mat &grad_x,
              cv::Mat &grad_y,
              float &beta);

std::vector<cv::Mat> minimizeL0Gradient(const cv::Mat &src);