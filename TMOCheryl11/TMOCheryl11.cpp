/* --------------------------------------------------------------------------- *
 * TMOCheryl11.cpp: implementation of the TMOCheryl11 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOCheryl11.h"

#include <opencv2/core/core.hpp> // Basic OpenCV structures (cv::Mat, Scalar)
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"

#include <iostream>
#include <vector>

#include "Cluster.h"

using namespace Cheryl11;
using namespace std;

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOCheryl11::TMOCheryl11()
{
    SetName(L"Cheryl11");
    SetDescription(L"Cluster-Based Color Space Optimizations");

    dParameter.SetName(L"Luminance");				// TODO - Insert parameters names
    dParameter.SetDescription(L"Value of multiplayer for Y of Yxy color scheme.");	// TODO - Insert parameter descriptions
    dParameter.SetDefault(1);							// TODO - Add default values
    dParameter = 1.0;
    dParameter.SetRange(-1000.0, 1000.0);				// TODO - Add acceptable range if needed
    this->Register(dParameter);
    
    iClusterCount.SetName(L"Clustering");
    iClusterCount.SetDescription(L"Number of clustering input image.");
    iClusterCount.SetDefault(7);
    iClusterCount = 7;
    iClusterCount.SetRange(0, 256);
    this->Register(iClusterCount);
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

/*
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
*/

    // Solution:
    inputImg = cv::imread(pSrc->GetFilename()); // default makes: CV_8UC3, BGR, 3 channels, values 0 to 255
    cv::imshow("inputImg", inputImg);
    inputImg.convertTo(inputImg, CV_32FC3);
    inputImg *= 1./255;
    cv::cvtColor(inputImg, inputImg, cv::COLOR_BGR2Luv);
    cv::imshow("inputImg - Luv", inputImg); // CV_32FC3, LUV, 3 channels
    cerr << "Size: " << inputImg.rows << "x" << inputImg.cols <<endl;
    
    clusterize(true);
    
    return 0;
}

void TMOCheryl11::clusterize(bool showClusteredImg = false)
{
    for (int i = 0; i < iClusterCount.GetInt(); i++) {
        clusters.push_back(Cluster(inputImg.rows, inputImg.cols));
    }
    
    // Points to k-means
    cv::Mat points = inputImg;
    points = points.reshape(1, points.total()); // TODO Must be reshaped? Try make better first loading image to cv::Mat.

    vector<int> labels; // Labels of results
    cv::TermCriteria criteria(CV_TERMCRIT_ITER, 10, 1.0);
    cv::Mat1f centers; // results of k-means

    // k-means for LAB color clusters
    double compactness = cv::kmeans(points, iClusterCount.GetInt(), labels, criteria, 3, cv::KMEANS_PP_CENTERS, centers);
    cerr << "Compactness: " << compactness << endl; // TODO more compactness better results? (very good: 1.05613e+07)
    cerr << "Color centers: " << centers << endl;

    cv::Mat tmpCoordinates(1, 2, CV_32F);

    cv::Mat imgResult(inputImg.rows, inputImg.cols, CV_32FC3);
    imgResult = imgResult.reshape(1, imgResult.total());
    
    for (int i = 0; i < inputImg.rows * inputImg.cols; i++)
    {
        // Colorize clusters by center color
        imgResult.at<float>(i, 0) = centers(labels[i], 0);
        imgResult.at<float>(i, 1) = centers(labels[i], 1);
        imgResult.at<float>(i, 2) = centers(labels[i], 2);

        int row = i / inputImg.cols;
        int col = i - (row * inputImg.cols);
        tmpCoordinates.at<float>(0, 0) = row;
        tmpCoordinates.at<float>(0, 1) = col;
        clusters.at(labels[i]).addPixel(points.row(i), tmpCoordinates);
    }

    imgResult = imgResult.reshape(3, inputImg.rows);

    for (int i = 0; i < iClusterCount; i++)
    {
        clusters.at(i).makeAverageCoordinates();
        
        cerr << "Cluster spatial mean: ";
        cv::Mat average = clusters.at(i).getAverageCoordinates();
        cv::Point2f d(average.at<float>(0, 1), average.at<float>(0, 0));
        cerr << d << endl;

        cv::circle(imgResult, d, 20, cv::Scalar(255, 0, 0), 1, cv::LineTypes::LINE_AA);
    }
    
    makeGraph();
    
    const vector<Graph::Edge> &edges = graph.getEdges();
    for (int i = 0; i < edges.size(); i++)
    {
        cv::Mat average = clusters.at(edges.at(i).c0).getAverageCoordinates();
        cv::Point2f p1(average.at<float>(0, 1), average.at<float>(0, 0));
        average = clusters.at(edges.at(i).c1).getAverageCoordinates();
        cv::Point2f p2(average.at<float>(0, 1), average.at<float>(0, 0));

        cv::line(imgResult, p1, p2, cv::Scalar(255, 0, 0), 1, cv::LineTypes::LINE_AA);
    }
    
    if (showClusteredImg) {
        cv::cvtColor(imgResult, imgResult, cv::COLOR_Luv2BGR);
        //imgResult *= 255;
        //imgResult.convertTo(imgResult, CV_8UC3);
        imshow("clusters", imgResult);
        
        cv::Mat test = clusters[0].getClusterImage(); // 13 for sun
        cv::cvtColor(test, test, cv::COLOR_Luv2BGR);
        test *= 255;
        test.convertTo(test, CV_8UC3);
        imshow("test", test);
        
        cv::cvtColor(inputImg, inputImg, cv::COLOR_Luv2BGR);
        inputImg *= 255;
        inputImg.convertTo(inputImg, CV_8UC3);
        cv::imshow("input_3", inputImg); // CV_8UC3, BGR, 3 channels, values 0 to 255
        
        cv::waitKey();
    }
}

void TMOCheryl11::makeGraph()
{
    for (int i = 0; i < iClusterCount.GetInt(); i++)
    {
        clusters[i].makeRegionMask();
        
        int histogramValuesCounter = 0;
        vector<int> histogram(iClusterCount.GetInt(), 0); // This histogram count neighbour pixels per clusters
        
        for (int r = 0; r < inputImg.rows; r++) {
            for (int c = 0; c < inputImg.cols; c++)
            {
                if (clusters[i].isPixelMasked(r, c))
                {
                    for (int j = 0; j < iClusterCount.GetInt(); j++)
                    {
                        if (clusters[j].isPixelOwner(r, c))
                        {
                            histogramValuesCounter += 1;
                            histogram.at(j) += 1;
                        }
                    }
                }
            }
        }
        
        // Histogram results and edge making
        float threshold = histogramValuesCounter / float(iClusterCount.GetInt());
        for (int h = 0; h < iClusterCount.GetInt(); h++)
        {
            //cerr << "[" << i << "] hist: " << histogram.at(h) << endl;
            if (threshold < histogram.at(h))
            {
                graph.addEdge(i, h);
            }
        }
        //cerr << "Histogram counter: " << histogramValuesCounter << " threshold: " << threshold << endl;
    }
}