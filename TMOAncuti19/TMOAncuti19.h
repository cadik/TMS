#include "TMO.h"

#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>

class TMOAncuti19 : public TMO
{
public:
	TMOAncuti19();
	virtual ~TMOAncuti19();
	virtual int Transform();

protected:
	bool LoadICAModel(const std::string &modelFile);
    void DecomposeInput(TMOImage* src, cv::Mat &red, cv::Mat &green, cv::Mat &blue);
    void ComputeWeightMaps(const cv::Mat &red, const cv::Mat &green, const cv::Mat &blue,
                           cv::Mat &normW_R, cv::Mat &normW_G, cv::Mat &normW_B);
	cv::Mat ComputeLocalWeightMap(const cv::Mat &channel, int patchSize = 5);
	double ComputeEntropy(const cv::Mat &patch);
	cv::Mat ComputeGlobalWeightMap(const cv::Mat &channel, int patchSize);
	double ComputeGlobalWeightAtPixel(const cv::Mat &channel, int r, int c, int patchSize);
    void BuildPyramids(const cv::Mat &channel, std::vector<cv::Mat> &lapPyr, int numLevels = 4);
    void BuildWeightPyramids(const cv::Mat &weight, std::vector<cv::Mat> &gaussPyr, int numLevels = 4);
    cv::Mat ReconstructFromPyramid(const std::vector<cv::Mat> &fusedPyr);
	cv::Mat NormalizeImage(const cv::Mat &input);
    void WriteOutput(TMOImage* dst, const cv::Mat &result);
	
	TMODouble alphaParameter;
	cv::Mat mIcaBasis;
	cv::Mat mMu;
	cv::Mat mSigma;
};
