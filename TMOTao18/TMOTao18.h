#include "TMO.h"
#include "TMOv.h"
#include <cmath>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

class TMOTao18 : public TMOv
{
public:
	TMOTao18();
	virtual ~TMOTao18();
	virtual int Transform();
	virtual int TransformVideo();

	void computeProximity(cv::Mat& currentFrame, cv::Mat& previousFrame, float& deltaL, float& deltaC, cv::Mat& mask);
	float computeEntropy(const cv::Mat& hist);
	std::vector<int> classify(const std::vector<cv::Vec2f>& proximityValues);
	double gaussianLikelihood(const cv::Vec2d& ui, const cv::Vec2d& mean, const cv::Mat& cov, double prior);
	double logLikelihood(const std::vector<cv::Vec2f>& dataPoints, const std::vector<int>& clusters, const std::vector<cv::Vec2f>& means, const std::vector<cv::Mat>& covariances, const std::vector<double>& priors);
	bool hasConverged(double L, double prevL, double threshold);
	cv::Mat applyHPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray, double phi);
	cv::Mat applyMPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray);
	cv::Mat applyLPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray, double beta);

protected:
	TMODouble dParameter;
};
