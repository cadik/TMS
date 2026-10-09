/* --------------------------------------------------------------------------- *
 * TMOZhangWan24.h: Color-to-gray image conversion using salient colors and    *
 * radial basis functions                                                      *
 *                                                                             *
 * Paper: L. Zhang, Y. Wan: Color-to-gray image conversion using salient       *
 *        colors and radial basis functions, 2024                              *
 * Author of the implementation: Maria Novakova                                *
 * --------------------------------------------------------------------------- */

#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace cv;

class TMOZhangWan24 : public TMO
{
public:
	TMOZhangWan24();
	virtual ~TMOZhangWan24();
	virtual int Transform();

private:
	// the three steps of the method
	void quantizeColors(cv::Mat &image, int &k, int max_k, float theta_0, float theta_1); // Sec. 3.1
	void ordering(int method);                                                            // Sec. 3.2
	cv::Mat createGrayScale(cv::Mat image, float sigma, bool edgeFix);                    // Sec. 3.3

	// part 1 helpers
	cv::Vec3f computePrincipalDirection(const cv::Mat &image);
	void expandCentroids(cv::Vec3f c_0, int &k, cv::Mat img, std::vector<cv::Vec3f> *centers);
	float euclideanDistance(const cv::Vec3f &color1, const cv::Vec3f &color2);
	std::vector<std::vector<std::pair<cv::Point, cv::Vec3f>>> clusterImage(cv::Mat image, std::vector<cv::Vec3f> centroids);
	void actualizeCenters(std::vector<cv::Vec3f> *centers, std::vector<std::vector<std::pair<cv::Point, cv::Vec3f>>> clusters);
	void MSE(std::vector<cv::Vec3f> centers, std::vector<std::vector<std::pair<cv::Point, cv::Vec3f>>> clusters, std::vector<float> *individualMSE);
	float MSE_k(cv::Mat image, std::vector<cv::Vec3f> centers, std::vector<std::vector<std::pair<cv::Point, cv::Vec3f>>> clusters);
	float g(cv::Mat image, std::vector<std::pair<cv::Point, cv::Vec3f>> cluster);
	float MSEG_k(cv::Mat image, std::vector<cv::Vec3f> centers, std::vector<std::vector<std::pair<cv::Point, cv::Vec3f>>> clusters);
	float Entropy(cv::Mat img);
	bool classification(float E);
	float M_k(float MSE_k, float MSEG_k);

	// part 2 helpers
	float rgb2grayOfCenter(const cv::Vec3f &labColor);
	float weightedEuclidean(const cv::Vec3f &color1, const cv::Vec3f &color2);

	// part 3 helpers
	float gaussianKernel(const cv::Vec3f &color1, const cv::Vec3f &color2, float sigma);
	float laplaceKernel(const cv::Vec3f &color1, const cv::Vec3f &color2, float sigma);
	float clamp(float value);
	float getGreyValue(cv::Vec3f img_color, std::vector<float> a, float sigma);
	cv::Mat repairEdgePixels(const cv::Mat &image, const cv::Mat &F); // extension, NOT part of the paper

	// state between the steps
	std::vector<cv::Vec3f> centers;                                     // salient colors (normalised Lab)
	std::vector<std::vector<std::pair<cv::Point, cv::Vec3f>>> clusters; // pixels of each salient color
	std::vector<float> grayvalues;                                      // ordered gray value of each salient color
	cv::Mat rgbImage;                                                   // input as BGR float in [0,1]

protected:
	TMOInt maxK;
	TMODouble sigma;
	TMOInt orderingMethod;
	TMOBool edgeRepair;
};
