/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                      A Spatial Post-Processing Algorithm                     *
 *                         for Images of Night Scenes                           *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2025                                 *
 *                                                                              *
 *******************************************************************************/
#include "TMO.h"

class TMOThompson02 : public TMO
{
public:
	/**
	 * @brief Constructor
	 */
	TMOThompson02();

	/**
	 * @brief Destructor
	 */
	virtual ~TMOThompson02();

	/**
	 * @brief Function to reduce luminance based on the average luminance and
	 * maximum luminance
	 * @param Y luminance
	 * @param YLogAvg average luminance
	 * @param Ymax maximum luminance
	 * @return double: reduced luminance
	 */
	double luminanceReduction(double Y, double YLogAvg, double Ymax);

	/**
	 * @brief Function to get scotopic luminance matrix
	 * @return cv::Mat: scotopic luminance matrix
	 */
	cv::Mat getScotopicLuminanceMat();

	/**
	 * @brief Function to get mesopic factor
	 * @param L: luminance
	 * @return double: mesopic factor
	 */
	double getMesopicFactor(double L);

	/**
	 * @brief Function to convert TMO data to matrix
	 * @return cv::Mat: matrix
	 */
	cv::Mat TMO2mat();

	/**
	 * @brief Function to convert matrix to TMO data
	 * @param input: matrix
	 */
	void mat2TMO(cv::Mat &input);

	/**
	 * @brief Function to apply Gaussian blur
	 * @param input: input matrix
	 * @param sigma: sigma value
	 * @return cv::Mat: blurred matrix
	 */
	cv::Mat applyGaussianBlur(cv::Mat &input, double sigma);

	/**
	 * @brief Function to apply night filter
	 * @param input: input matrix
	 * @return cv::Mat: filtered matrix
	 */
	cv::Mat applyNightFilter(cv::Mat &input);

	/**
	 * @brief Function to add Gaussian noise
	 * @param input: input matrix
	 * @return cv::Mat: noisy matrix
	 */
	cv::Mat addGaussianNoise(cv::Mat &input);

	/**
	 * @brief Function to apply the tone mapping operator
	 * @return int: 0 = success, 1 = error
	 */
	virtual int Transform();

protected:
	// parameters
	TMODouble sigmaBlur = 2.0;	   // sigma value for Gaussian blur
	TMODouble gammaEdge = 1.25;	   // gamma value for edge enhancement
	TMODouble sigmaNoise = 0.0125; // sigma value for Gaussian noise
	TMODouble rcf = 1.0;		   // range compression factor
	TMODouble mf = 0.0;			   // mesopic factor

	// constants
	// RGB value for bluish grey
	std::vector<double> bluishGreyRGB = {0.92, 0.97, 1.27};
};