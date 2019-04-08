#ifndef CHERYL11_CLUSTER_H
#define CHERYL11_CLUSTER_H

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
        void makeRegionMask();
        bool isPixelOwner(int row, int col);
        bool isPixelMasked(int row, int col);
        void makeAverageCoordinates();
        
        cv::Mat getAverageCoordinates();
        int getCoordinatesSize()
        {
            return coordinates.rows * coordinates.cols;
        }
        cv::Mat getAverageColor();
        double getMappedColor();
        cv::Mat getClusterImage()
        {
            return clusterImage;
        }

        void setColorCenter(cv::Mat center);
        
        float nearestClusterPathLenght;
        
    protected:
        cv::Mat colorCenter; // Same as averageColor
        double mappedColorCenter; // 'm' in equation
        
        cv::Mat clusterImage;
        cv::Mat regionMask;
        
        cv::Mat colors; // List of colors
        cv::Mat coordinates; // List of coordinates for colors
        std::vector<Pixel> pixels; // The same as above but completed
        
        cv::Mat averageCoordinates; // For path between clusters
    };
   
} // namespace Cheryl11

#endif /* CHERYL11_CLUSTER_H */
