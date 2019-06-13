/* -------------------------------------------------------------------------- *
 * baseAndDetailDecomposition.h: image decomposition to base & detail layers  *
 * Author: Pavel Sedlar (2018)                                                *
 *          - initial code based on L0 minimization                           *
 * Editor: Tomas Hudziec (2019)                                               *
 *          - commented out unused code                                       *
 * -------------------------------------------------------------------------- */
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
										     float &beta, 									
										     cv::Mat &lambdaMatrix);

extern cv::Mat minimizeL0GradientSecondFaze(const cv::Mat &src, 
											cv::Mat lambdaMatrix1, 
											int rows1, 
											int cols1);

extern cv::Mat getGradientMagnitude(const cv::Mat &src);

// extern cv::Mat getGradientMagnitudeFromOneLayer(const cv::Mat &src);
							  
extern cv::Mat getWeightsFromBaseLayer(const cv::Mat &gradient, 
									   int rows, 
									   int cols,
									   int r);
									   
// extern std::vector<cv::Mat> detailMaximalization(const cv::Mat &base, const cv::Mat &detailSum, const cv::Mat &r1Weight, const cv::Mat &r2Weight, int rows, int cols, int counter, std::vector<cv::Mat> detail); 
							  
extern cv::Mat getDetailControl(const cv::Mat &base, 
						 const cv::Mat &detail,
					     const cv::Mat &s,
						 const cv::Mat &t,
						 float mu, 
						 int rows, 
						 int cols);
						 
// extern cv::Mat stochasticOptimizationForGetSigma(cv::Mat base, cv::Mat original, int rows, int cols, int counter);

// extern cv::Mat myOwn2DFilter(cv::Mat image, cv::Mat sigmaMap, int rows, int cols);
