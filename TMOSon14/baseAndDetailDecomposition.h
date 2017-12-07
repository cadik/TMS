/**
 * Zde uz jsou zahrnuty vsechny fce z originalniho l0minimalization
 */
#include <fstream>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer.hpp>
#include <Eigen/Sparse>

#include "L0minimization.h";

namespace po = boost::program_options;
namespace fs = boost::filesystem;

std::string input_file, out_dir, config_file;

// optimization params
float lambda;
float beta0;
float beta_max;
float kappa;
bool exact;
int iter_max = 1000;

// buffers for solving linear system
Eigen::SparseMatrix<float> A0, E;
Eigen::SparseMatrix<float> GX, GY;
Eigen::VectorXf S_vec, I_vec, H_vec, V_vec;

void optimizeWithAdaptiveLambdaMatrix(cv::Mat &S, 
    const cv::Mat &I, 
    cv::Mat &H, 
    cv::Mat &V, 
    cv::Mat &grad_x,
    cv::Mat &grad_y,
    cv::Mat &lambdaMatrix,
    float &beta);

cv::Mat getAdaptiveLambdaMatrix(const cv::Mat &gradientFromFirstSmoothing, int rows, int cols);

std::vector<cv::Mat> minimizeL0GradientSecondFaze(const cv::Mat &src, cv::Mat lambdaMatrix1, int rows, int cols);

cv::Mat getGradientFromFirstSmoothing(const cv::Mat &src, int rows, int cols);