#include "TMO.h"
#include "TMOv.h"
#include <cmath>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <algorithm>
#include <vector>

class TMOTao18 : public TMOv
{
public:
	TMOTao18();
	virtual ~TMOTao18();
	virtual int Transform();
	virtual int TransformVideo();

	void computeProximity(const cv::Mat& currentFrame, const cv::Mat& previousFrame, double& deltaL, double& deltaC);
	double computeEntropy(const cv::Mat& hist);
	int classify(double deltaL, double deltaC);
	double gaussianLikelihood(const cv::Vec2d& ui, const cv::Vec2d& mean, const cv::Mat& cov, double prior);
	double logLikelihood(const std::vector<cv::Vec2d>& dataPoints, const std::vector<int>& clusters, const std::vector<cv::Vec2d>& means, const std::vector<cv::Mat>& covariances, const std::vector<double>& priors);
	bool hasConverged(double L, double prevL, double threshold);
	cv::Mat applyHPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray, double phi);

protected:
	TMODouble dParameter;
};
