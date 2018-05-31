#include <qpOASES.hpp>
#include <cmath>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

extern std::vector<cv::Mat> optimizationWithOases(int height, int width, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels);