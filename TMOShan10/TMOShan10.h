#include <opencv2/opencv.hpp>
#include <Eigen/Sparse>
#include "TMO.h"
#include "MatUtil.h"

using namespace std;
using namespace cv;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;

class TMOShan10 : public TMO
{
public:
	TMOShan10();
	virtual ~TMOShan10();
	virtual int Transform();

protected:
   TMODouble beta1;
   TMODouble beta2;
   TMODouble beta3;
   TMODouble sSat;
   TMOInt winSize;
   TMOInt levelNum;

private:
   MatUtil matUtil;

   /*
    * Create guidance map.
    */
   Mat generateGuidanceMap(Mat lumo,
                           int winSize,
                           double beta1,
                           double beta2,
                           const double beta3,
                           const double kappa);

   /*
    * Returns linear coefficients defined in Eq. (17) in the appendix
    * of the article.
    */
   Mat getLinearCoefficients(Mat lum,
                             Mat lumo,
                             Mat map,                             
                             uint8_t multiGridFilt,
                             int winSize,
                             double epsilon);

   /*
    * Upsample by linear coefficients defined in Eq. (17) in the appendix
    * of the article.
    */
   Mat upsampleByLinearCoefficients(Mat lum,
                                    Mat lumo,
                                    Mat bLumo,
                                    Mat map,                                    
                                    uint8_t multiGridFilt,
                                    int winSize,
                                    double epsilon);

   /*
    * Creates matrices S and B defined in the article.
    */
   std::pair<Mat, SpMat> getLaplacian(Mat lumo,
                                      Mat map,
                                      uint8_t multiGridFilter,
                                      int winSize,
                                      double epsilon);

   /*
    * Method has two important parts:
    * - creating the matrices S and B
    * - solving the linear system
    */
   Mat solveLinearSystemFire(Mat lumo,
                             Mat constsMap,
                             Mat constsValue,
                             Mat map,
                             uint8_t multiGridFilter,
                             int winSize,
                             double epsilon,
                             double lambda);

   /*
    * A multigrid architecture that solves the large linear system.
    */
   Mat solveLinearSystem(Mat lumo,
                         Mat constsMap,
                         Mat constsValue,
                         Mat map,                         
                         uint8_t levelNum,
                         int winSize,
                         double epsilon,
                         double lambda);
};
