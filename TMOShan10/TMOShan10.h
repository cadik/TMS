#include <opencv2/opencv.hpp>
#include "TMO.h"

using namespace std;
using namespace cv;

class TMOShan10 : public TMO
{
public:
	TMOShan10();
	virtual ~TMOShan10();
	virtual int Transform();

protected:
   TMODouble beta1;
   TMODouble beta2;
   TMODouble sSat;
   TMOInt winSize;

private:
    Mat generateGuidanceMap(Mat lumo,
                            int winSize,
                            double beta1,
                            double beta2,
                            const double beta3,
                            const double kappa);
};
