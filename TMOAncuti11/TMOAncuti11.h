#include "TMO.h"
#include "TMOv.h"
#include <opencv2/opencv.hpp>
#include <opencv2/saliency.hpp>


struct scaleDifference{
	cv::Mat diff;         //difference map
	int c;                //center scale
	int s;                //surround scale
};


class TMOAncuti11 : public TMOv
{
public:
	TMOAncuti11();
	virtual ~TMOAncuti11();
	virtual int Transform();
	virtual int TransformVideo();

	

	cv::Mat decolorization(cv::Mat &input, double eta, double phi);
	cv::Mat computeSaliencyMap(cv::Mat &input, bool color);
	double offsetAngleSelection(cv::Mat &input);
	double averageSaliencyInRegion(cv::Mat& saliencyMap, cv::Point center, int radius);
	std::vector<cv::Point> findSalientPoint(cv::Mat &saliencyMap, int topN, int radius);
	std::vector<cv::Mat> buildPyramid(cv::Mat& input, int levels);
	cv::Mat intensityChannel(cv::Mat& input);
	std::vector<scaleDifference> mapsDifference(std::vector<cv::Mat>& input);
	cv::Mat sumNormalized(std::vector<scaleDifference>& input, cv::Size size);
	cv::Mat ittiNormalize(cv::Mat& input);
	void computeColorChannels(cv::Mat& input, cv::Mat& Rp, cv::Mat& Gp, cv::Mat& Bp, cv::Mat& Yp);
	cv::Mat gaborFilter(cv::Mat& input, float angle);

protected:
	TMODouble dParameter;
};
