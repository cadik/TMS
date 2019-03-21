/* --------------------------------------------------------------------------- *
 * Cluster.cpp: implementation of the Cluster class.   *
 * --------------------------------------------------------------------------- */

#include "Cluster.h"

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

cv::Mat Cluster::getAverageCoordinates()
{
    cv::Mat averageCoordinates;
    cv::reduce(coordinates, averageCoordinates, 0, cv::REDUCE_AVG, CV_32F);
    return averageCoordinates;
}

cv::Mat Cluster::getAverageColor()
{
    return colorCenter;
}
