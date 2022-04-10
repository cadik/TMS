#include "TMO.h"
#include "ColorspaceConverter.h"
#include "MatUtil.h"

#define CHANNELSCNT 3
#define NC 1
#define IC 2
#define RF 3

#define EDGEAWARESMOOTH 0
#define DETAILENHC 1
#define STYLIZATION 2
#define PENCILSKETCH 3

class TMOGastal11 : public TMO
{
public:
	TMOGastal11();
	virtual ~TMOGastal11();
	virtual int Transform();

protected:
    TMOInt numIter;
    TMOInt filterType;
    TMOInt filterAppl;
    TMODouble sigma_s;
    TMODouble sigma_r;
    TMODouble shadeFactor;
    TMOBool pencilColor;

private:
    ColorspaceConverter converter;
    MatUtil matUtil;

    Mat TMOImage2Mat(TMOImage* pSrc);

    void domainTransform (Mat srcMat,
                          Mat &dHdx,
                          Mat &dVdy,
                          Mat &ct_H,
                          Mat &ct_V,
                          double sigma_s,
                          double sigma_r,
                          uint8_t filterOper);

    void boxFilter (Mat &srcMat,
                    Mat &domainPosition,
                    Mat &lowerPos,
                    Mat &upperPos,
                    Mat &lowerIdx,
                    Mat &upperIdx,
                    double boxRadius);

    void NCfilter (Mat srcMat, Mat &outMat, Mat &domainPosition, double radius);

    void ICfilter (Mat srcMat, Mat &outMat, Mat &domainPosition, double radius);

    void Rfilter(Mat srcMat, Mat &outMat, Mat &domainPosition, double sigma_H);

    void filterOperationNC(Mat srcMat,
                           Mat &resMat,
                           double sigma_s,
                           double sigma_r,
                           uint8_t numIter);

    void filterOperationIC(Mat srcMat,
                           Mat &resMat,
                           double sigma_s,
                           double sigma_r,
                           uint8_t numIter);

    void filterOperationRF(Mat srcMat,
                           Mat &resMat,
                           double sigma_s,
                           double sigma_r,
                           uint8_t numIter);

    void filterOperation(Mat srcMat,
                         Mat &resMat,
                         double sigma_s,
                         double sigma_r,
                         uint8_t filterType,
                         uint8_t numIter);

    void edgeAwareSmoothing(Mat srcMat,
                            Mat &resMat,
                            double sigma_s,
                            double sigma_r,
                            uint8_t filterType,
                            uint8_t numIter);

    void detailEnhancement(Mat srcMat,
                           Mat &resMat,
                           double sigma_s,
                           double sigma_r,
                           uint8_t filterType,
                           uint8_t numIter);

    void findMagnitude(Mat srcMat, Mat &mag);

    void stylization(Mat srcMat,
                     Mat &resMat,
                     double sigma_s,
                     double sigma_r,
                     uint8_t filterType,
                     uint8_t numIter);

    void pencilSketch(Mat srcMat,
                      Mat &sketch,
                      Mat &colorSketch,
                      double sigma_s,
                      double sigma_r,
                      double shadeFactor,
                      bool color,
                      uint8_t numIter);

    void pencilDrawing(Mat srcMat,
                       Mat &resMat,
                       double sigma_s,
                       double sigma_r,
                       double shadeFactor,
                       bool color,
                       uint8_t numIter);
};
