/* --------------------------------------------------------------------------- *
 * TMOCheryl11.cpp: implementation of the TMOCheryl11 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOCheryl11.h"

#include <opencv2/core/core.hpp> // Basic OpenCV structures (cv::Mat, Scalar)
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"

#include <iostream>

using namespace std;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOCheryl11::TMOCheryl11()
{
    SetName(L"Cheryl11");						// TODO - Insert operator name
    SetDescription(L"Cluster-Based Color Space Optimizations");	// TODO - Insert description

    dParameter.SetName(L"Luminance");				// TODO - Insert parameters names
    dParameter.SetDescription(L"Value of multiplayer for Y of Yxy color scheme.");	// TODO - Insert parameter descriptions
    dParameter.SetDefault(1);							// TODO - Add default values
    dParameter=1.;
    dParameter.SetRange(-1000.0,1000.0);				// TODO - Add acceptable range if needed
    this->Register(dParameter);
}

TMOCheryl11::~TMOCheryl11()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOCheryl11::Transform()
{
    // Source image is stored in local parameter pSrc
    // Destination image is in pDst

    // Initialy images are in RGB format, but you can
    // convert it into other format
    pSrc->Convert(TMO_LAB);
    pDst->Convert(TMO_LAB);

    double* pSourceData = pSrc->GetData();				// You can work at low level data
    double* pDestinationData = pDst->GetData();			// Data are stored in form of array
                                                        // of three doubles representing
                                                        // three colour components
    double pL, pA, pB;

    int j = 0;
    for (j = 0; j < pSrc->GetHeight(); j++)
    {
        pSrc->ProgressBar(j, pSrc->GetHeight());	// You can provide progress bar
        for (int i = 0; i < pSrc->GetWidth(); i++)
        {
            pL = *pSourceData++;
            pA = *pSourceData++;
            pB = *pSourceData++;

            // Here you can use your transform
            // expressions and techniques...
            pL *= dParameter;							// Parameters can be used like
            pA = 00.0;
            pB = 00.0;

            // and store results to the destination image
            *pDestinationData++ = pL;
            *pDestinationData++ = pA;
            *pDestinationData++ = pB;
        }
    }
    
    pSrc->ProgressBar(j, pSrc->GetHeight());
    pDst->Convert(TMO_RGB);

    clusterize(true);
    
    return 0;
}

void TMOCheryl11::clusterize(bool showClusteredImg = false)
{
    const int max_clusters = 6;
    cv::Mat resultColorTable(max_clusters, 1, CV_8UC3);
    resultColorTable.at<char>(0, 0) = 0;
    resultColorTable.at<char>(0, 1) = 0;
    resultColorTable.at<char>(0, 2) = 255;
    
    resultColorTable.at<char>(1, 0) = 0;
    resultColorTable.at<char>(1, 1) = 255;
    resultColorTable.at<char>(1, 2) = 0;
    
    resultColorTable.at<char>(2, 0) = 255;
    resultColorTable.at<char>(2, 1) = 0;
    resultColorTable.at<char>(2, 2) = 0;
    
    resultColorTable.at<char>(3, 0) = 0;
    resultColorTable.at<char>(3, 1) = 0;
    resultColorTable.at<char>(3, 2) = 255/2;
    
    resultColorTable.at<char>(4, 0) = 0;
    resultColorTable.at<char>(4, 1) = 255/2;
    resultColorTable.at<char>(4, 2) = 0;
    
    resultColorTable.at<char>(5, 0) = 255/2;
    resultColorTable.at<char>(5, 1) = 0;
    resultColorTable.at<char>(5, 2) = 0;

    cv::Mat points(pSrc->GetHeight(), pSrc->GetWidth(), CV_64FC3, pSrc->GetData());
    points.convertTo(points, CV_32FC3); // Converting for cv::kmeans assert.
    points = points.reshape(1, points.total()); // TODO Must be reshaped? Try make better first loading image to cv::Mat.

    const int clusterCount = 7; // TODO Heuristic to estimate number of colors.
    vector<int> labels;
    cv::TermCriteria criteria(CV_TERMCRIT_ITER, 10, 1.0);
    cv::Mat1f centers;

    double compactness = cv::kmeans(points, clusterCount, labels, criteria, 3, cv::KMEANS_PP_CENTERS, centers);
    cerr << "Compactness: " << compactness << endl; // TODO more compactness better results? (very good: 1.05613e+07)
    cerr << "Centers: " << centers << endl;

    cv::Mat imgResult(pSrc->GetHeight(), pSrc->GetWidth(), CV_8UC3);
    imgResult = imgResult.reshape(1, imgResult.total());
    
    for (int i = 0; i < pSrc->GetHeight() * pSrc->GetWidth(); i++)
    {
        // By cluster center color
        imgResult.at<char>(i, 0) = centers(labels[i], 0);
        imgResult.at<char>(i, 1) = centers(labels[i], 1);
        imgResult.at<char>(i, 2) = centers(labels[i], 2);
        /**/

/*
        // By special table:
        imgResult.at<char>(i, 0) = resultColorTable.at<char>(labels[i], 0);
        imgResult.at<char>(i, 1) = resultColorTable.at<char>(labels[i], 1);
        imgResult.at<char>(i, 2) = resultColorTable.at<char>(labels[i], 2);
        /**/
    }

    imgResult = imgResult.reshape(3, pSrc->GetHeight());
    imgResult.convertTo(imgResult, CV_8U);


    if (showClusteredImg) {
        imshow("clusters", imgResult);
        cv::waitKey();
    }
}