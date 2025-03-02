/*******************************************************************************
 *                                                                              *
 *                         Brno University of Technology                        *
 *                       Faculty of Information Technology                      *
 *                                                                              *
 *                   Fast Local Tone Mapping, Summed-Area Tables                *
 *                          and Mesopic Vision Simulation                       *
 * 																			    *
 *                                 Bachelor thesis                              *
 *             Author: Jan Findra [xfindr01 AT stud.fit.vutbr.cz]               *
 *                                    Brno 2024                                 *
 *                                                                              *
 *******************************************************************************/

#include "TMO.h"

class TMOSlomp12 : public TMO
{
public:
	/**
	 * @brief Constructor
	 */
	TMOSlomp12();

	/**
	 * @brief Destructor
	 */
	virtual ~TMOSlomp12();

	/**
	 * @brief Function to convert the image to logarithmic luminance matrix
	 * @return cv::Mat: logarithmic luminance matrix
	 */
	cv::Mat TMOImageToLogLuminanceMat();

	/**
	 * @brief Function to log the luminance image
	 * @param luminanceMat pointer to the luminance matrix
	 */
	void logLuminanceImage(cv::Mat *luminanceMat);

	/**
	 * @brief Function to perform full mipmap
	 * @param mat pointer to the matrix
	 * @return double: average value
	 */
	double fullMipmap(cv::Mat *mat);

	/**
	 * @brief Function to scale the luminance matrix
	 * @param luminanceMat pointer to the luminance matrix
	 * @param keyValue key value
	 */
	void scaleLuminance(cv::Mat *luminanceMat, double keyValue);

	/**
	 * @brief Function to scale the luminance image
	 * @param luminanceMat pointer to the luminance matrix
	 */
	void scaledLuminanceImage(cv::Mat *luminanceMat);

	/**
	 * @brief Function to compute the box filter
	 * @param SAT pointer to the summed-area table
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param s scale
	 * @return double: box filter value
	 */
	double boxFilter(cv::Mat *SAT, int x, int y, int s);

	/**
	 * @brief Function to get the normalized difference
	 * @param conv0 first convolution
	 * @param conv1 second convolution
	 * @param s scale
	 * @return double: normalized difference
	 */
	double getNormalizedDifference(double conv0, double conv1, int s);

	/**
	 * @brief Function to get the maximum scale
	 * @param SAT pointer to the summed-area table
	 * @param x x coordinate
	 * @param y y coordinate
	 * @param averageValue average value which was subtracted from the luminance matrix
	 * @return int: maximum scale
	 */
	int getMaxScale(cv::Mat *SAT, int x, int y, double averageValue);

	/**
	 * @brief Function to get the red response value
	 * @param illuminance illuminance
	 * @return double: red response value
	 */
	double redResponseValue(double illuminance);

	/**
	 * @brief Function to compute the arithmetic luminance average
	 * @return double: arithmetic luminance average
	 */
	double arithLuminanceAverage();

	/**
	 * @brief Function to apply the tone mapping operator
	 * @return int: 0 = success, 1 = error
	 */
	virtual int Transform();

protected:
	// parameters
	TMOBool local;	 // local or global luminance compression
	TMOBool mesopic; // mesopic vision color change
	TMOBool varying; // spatially-varying or spatially-uniform mesopic vision reproduction operator

	// constants
	double alpha = 0.72;
	double phi = 8.;
	double epsilon = 0.025;
};
