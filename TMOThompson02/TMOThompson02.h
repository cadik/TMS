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
	// structure to hold a matrix and a double
	typedef struct matAndDouble
	{
		cv::Mat mat;
		double d;
	} matAndDouble;

	/**
	 * @brief Constructor
	 */
	TMOThompson02();

	/**
	 * @brief Destructor
	 */
	virtual ~TMOThompson02();

	/**
	 * @brief Function to get luminance matrix
	 * @return matAndDouble: structure with luminance matrix and maximum luminance
	 */
	matAndDouble getLuminanceMat();

	/**
	 * @brief Function to get displayable luminance
	 * @param L: luminance
	 * @param maxLuminance: maximum luminance
	 * @return double: displayable luminance
	 */
	double getDisplayableLuminance(double L, double maxLuminance);

	/**
	 * @brief Function to map luminance
	 * @param luminanceMat: luminance matrix
	 * @param maxLuminance: maximum luminance
	 */
	void mapLuminance(cv::Mat &luminanceMat, double maxLuminance);

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
	TMODouble sigmaBlur = 1.6;	   // sigma value for Gaussian blur
	TMODouble gammaEdge = 1.25;	   // gamma value for edge enhancement
	TMODouble sigmaNoise = 0.0125; // sigma value for Gaussian noise

	// constants
	// RGB value for bluish grey
	std::vector<double> bluishGreyRGB = {1.05, 0.97, 1.27};
};
