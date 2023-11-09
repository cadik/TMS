#include "TMO.h"

class TMOQueiroz06 : public TMO
{
public:
	TMOQueiroz06();
	virtual ~TMOQueiroz06();
	virtual int Transform();
	cv::Mat srcToMat(TMOImage *image);
	cv::Mat resizeToSquare(const cv::Mat &input);
	void PrintImage(const cv::Mat input, std::string name);
	void HaarWavelet(cv::Mat &src,cv::Mat &dst,int level);
	void InvHaarWavelet(cv::Mat &src,cv::Mat &dst,int level);

protected:
	TMODouble dParameter;
	TMOInt dVariant;
	TMOInt dSavePartialResults;
};
