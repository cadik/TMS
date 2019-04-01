/* --------------------------------------------------------------------------- *
 * Cluster.cpp: implementation of the Cluster class.   *
 * --------------------------------------------------------------------------- */
#include "Cluster.h"

#include <iostream>

#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"

using namespace Cheryl11;

Cluster::Cluster(int rows, int cols, int color_dimension)
{
    colorCenter = cv::Mat::zeros(1, color_dimension, CV_32F);
    clusterImage = cv::Mat::zeros(rows, cols, CV_32FC3);
    
    colors = cv::Mat(0, color_dimension, CV_32F);
    coordinates = cv::Mat(0, 2, CV_32F);
}

Cluster::~Cluster()
{
}

void Cluster::addPixel(cv::Mat color, cv::Mat position)
{
    colors.push_back(color);
    coordinates.push_back(position);
    
    clusterImage.at<cv::Vec3f>(position.at<float>(0, 0), position.at<float>(0, 1)) = color;
}

void Cluster::makeRegionMask()
{
    cv::extractChannel(clusterImage, regionMask, 0);
    cv::threshold(regionMask, regionMask, 0.0, 255.0, cv::THRESH_BINARY);
    regionMask.convertTo(regionMask, CV_8U);
    
    cv::Laplacian(regionMask, regionMask, CV_8U, 1, 1, 0); // Delated region by 1px
    // Maybe cv::Canny() could be more interesting
    // Or make blur
}

bool Cluster::isPixelOwner(int row, int col)
{
    if (clusterImage.at<float>(row, col) > 0)
    {
        return true;
    }
    return false;
}

bool Cluster::isPixelMasked(int row, int col)
{
    if (regionMask.at<float>(row, col) > 0)
    {
        return true;
    }
    return false;
}

void Cluster::makeAverageCoordinates()
{
    cv::reduce(coordinates, averageCoordinates, 0, cv::REDUCE_AVG, CV_32F);
}

cv::Mat Cluster::getAverageCoordinates()
{
    return averageCoordinates;
}

cv::Mat Cluster::getAverageColor()
{
    return colorCenter;
}
