/**
 * Zde uz jsou zahrnuty vsechny fce z originalniho l0minimalization
 */
#include <fstream>
#include <opencv2/opencv.hpp>
#include <Eigen/Sparse>
#include "L0minimization.h"

extern cv::Mat minimizeL0Gradient1(const cv::Mat &src);

extern cv::Mat getAdaptiveLambdaMatrix(const cv::Mat &gradient, 
									   int rows, 
									   int cols);

extern void optimizeWithAdaptiveLambdaMatrix(cv::Mat &S, 
										     const cv::Mat &I, 
										     cv::Mat &H, 
										     cv::Mat &V, 
										     cv::Mat &grad_x,
										     cv::Mat &grad_y,
										     cv::Mat &lambdaMatrix,
										     float &beta);

extern cv::Mat minimizeL0GradientSecondFaze(const cv::Mat &src, 
											cv::Mat lambdaMatrix1, 
											int rows, 
											int cols);

extern cv::Mat getGradientMagnitude(const cv::Mat &src);

extern cv::Mat getDetailLayer(const cv::Mat &orig, 
							  const cv::Mat &base, 
							  int rows, 
							  int cols);
							  
extern cv::Mat getWeightsFromBaseLayer(const cv::Mat &gradient, 
									   int rows, 
									   int cols,
									   int r);
							  
extern cv::Mat getDetail(const cv::Mat &base, 
						 const cv::Mat &detail,
					     const cv::Mat &s,
						 const cv::Mat &t,
						 float mu, 
						 int rows, 
						 int cols);
						 
extern cv::Mat stochasticOptimizationForGetSigma(cv::Mat base, cv::Mat original, int rows, int cols, int counter);

extern cv::Mat getSumOfCosts(cv::Mat r, cv::Mat g, cv::Mat b, int rows, int cols); 

extern cv::Mat myOwn2DFilter(cv::Mat image, cv::Mat sigmaMap, int rows, int cols);
