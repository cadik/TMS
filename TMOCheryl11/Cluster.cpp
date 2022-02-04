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

    nearestClusterPathLenght = 0.0f;
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
    if (clusterImage.at<cv::Vec3f>(row, col)[0] > 0)
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

void Cluster::makeCovarianceMatrix()
{
    cv::calcCovarMatrix(colors, covMatrix, mu, CV_COVAR_NORMAL | CV_COVAR_ROWS); // TODO 'mu' can be changed by 'colorCenter' but 'COVAR_USE_AVG' flag must be set.

    covMatrix = covMatrix / (colors.rows - 1);

    invCovMatrix = covMatrix.inv();

    /* 
    std::cerr << "cov: " << std::endl;
    std::cerr << covMatrix << std::endl;
.
    std::cerr << "cov inverse: " << std::endl;
    std::cerr << invCovMatrix << std::endl;

    std::cerr << "mu: " << std::endl;
    std::cerr << mu << std::endl;
     * */
}

double Cluster::getWeight(cv::Mat u)
{
    u = u.reshape(1, 1);
    u.convertTo(u, CV_64F);

    cv::Mat U = u - mu;
    cv::Mat M = invCovMatrix * u.t();
    cv::Mat ret = 1 / (U * M);

    return ret.at<double>();
}

cv::Mat Cluster::getAverageCoordinates()
{
    return averageCoordinates;
}

double Cluster::getMappedColor()
{
    return mappedColorCenter;
}

cv::Mat Cluster::getAverageColor()
{
    return colorCenter;
}

void Cluster::setColorCenter(cv::Mat center)
{
    colorCenter = center;

    cv::Mat mapped_color = cv::Mat(1, 1, CV_32FC3);
    mapped_color.at<cv::Vec3f>(0, 0) = colorCenter;
    cv::cvtColor(mapped_color, mapped_color, cv::COLOR_Luv2BGR);
    cv::cvtColor(mapped_color, mapped_color, cv::COLOR_BGR2GRAY);
    mappedColorCenter = mapped_color.at<float>(0, 0);
}
