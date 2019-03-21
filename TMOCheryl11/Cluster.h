#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>

#include <opencv2/core/core.hpp> // Basic OpenCV structures (cv::Mat, Scalar)

namespace Cheryl11 {

    class Pixel {
    public:
        Pixel() { }
        Pixel(cv::Mat color, cv::Mat position)
        {
            this->color = color;
            this->position = position;
        }
        ~Pixel();
        
    protected:
        cv::Mat color;
        cv::Mat position;
    };
    
    class Cluster {
    public:
        //Cluster() {}
        Cluster(int rows, int cols, int color_dimension = 3);
        ~Cluster();

        void addPixel(cv::Mat color, cv::Mat position);
        
        cv::Mat getAverageCoordinates();
        int getCoordinatesSize()
        {
            return coordinates.rows * coordinates.cols;
        }
        cv::Mat getAverageColor();
        cv::Mat getClusterImage()
        {
            return clusterImage;
        }
        
        void setColorCenter(cv::Mat center) { colorCenter = center; }
    protected:
        cv::Mat colorCenter; // Same as averageColor
        cv::Mat clusterImage;
        
        cv::Mat colors;
        cv::Mat coordinates;
        
        std::vector<Pixel> pixels;
    };
   
} // namespace Cheryl11

#endif /* CLUSTER_H */
