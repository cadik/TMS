#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <Eigen/Sparse>

extern float beta_max;
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
                                const std::vector<std::pair<int, float>> x_indices,
                                const std::vector<std::pair<int, float>> y_indices);

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
                     const cv::Mat &U,
                     cv::Mat &H,
                     cv::Mat &V,
                     cv::Mat &orig_grad_x,
                     cv::Mat &orig_grad_y,
                     cv::Mat &grad_x,
                     cv::Mat &grad_y,
                     float &beta,
                     float &eta,
                     float &lambda);

// extern std::vector<cv::Mat> minimizeL0Gradient(const cv::Mat &src);
// extern cv::Mat minimizeL0Gradient(const cv::Mat &src);
extern cv::Mat minimizeL0Gradient(const cv::Mat &src, float eta = 4.0, float lambda = 0.01, float kappa = 2.0);

extern cv::Mat calcNeighbourhoodVariance(const cv::Mat &I, int r);
