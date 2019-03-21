#include "TMO.h"

#include <opencv2/core/core.hpp> // Basic OpenCV structures (cv::Mat, Scalar)
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"

#include <vector>

#include "Cluster.h"

class TMOCheryl11 : public TMO
{
public:
    TMOCheryl11();
    virtual ~TMOCheryl11();
    virtual int Transform();

protected:
    // TMO code
    TMODouble dParameter;
    TMOInt iClusterCount; // TODO Heuristic to estimate number of colors.
    
    // Cheryl11 code
    cv::Mat inputImg;
    std::vector<Cheryl11::Cluster> clusters;
    
    void clusterize(bool showClusteredImg);
};
