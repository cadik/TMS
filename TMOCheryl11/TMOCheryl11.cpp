/* --------------------------------------------------------------------------- *
 * TMOCheryl11.cpp: implementation of the TMOCheryl11 class.                   *
 * --------------------------------------------------------------------------- */

#include "TMOCheryl11.h"

#include <opencv2/core/core.hpp> // Basic OpenCV structures (cv::Mat, Scalar)
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"

#include "optim.hpp"

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

    iClusterCount.SetName(L"Clusters");
    iClusterCount.SetDescription(L"Number of clustering input image. <1, 256>");
    iClusterCount.SetDefault(7);
    iClusterCount = 7;
    iClusterCount.SetRange(1, 256);
    this->Register(iClusterCount);
    
    dContrast_k.SetName(L"Contrast");
    dContrast_k.SetDescription(L"Contrast enhacement. <0.0, 1.0>");
    dContrast_k.SetDefault(0.2);
    dContrast_k = 0.2;
    dContrast_k.SetRange(0.0, 1.0);
    this->Register(dContrast_k);
    
    dRegularization_w.SetName(L"Regularization");
    dRegularization_w.SetDescription(L"Increase to preserve the realistic nature of the scene. <0.0, 1.0>");
    dRegularization_w.SetDefault(0.8);
    dRegularization_w = 0.8;
    dRegularization_w.SetRange(0.0, 1.0);
    this->Register(dRegularization_w);
    
    dJND.SetName(L"JND");
    dJND.SetDescription(L"Just Noticeable Difference. <0.0, 1.0>");
    dJND.SetDefault(0.001);
    dJND = 0.001;
    dJND.SetRange(0.0, 1.0);
    this->Register(dJND);
}

TMOCheryl11::~TMOCheryl11()
{
}

double opt_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);

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
    
    cv::Mat inputGrey;
    cv::cvtColor(inputImg, inputGrey, cv::COLOR_BGR2GRAY);
    cv::imshow("inputGrey", inputGrey);

    inputImg.convertTo(inputImg, CV_32FC3);
    inputImg *= 1./255;
    cv::cvtColor(inputImg, inputImg, cv::COLOR_BGR2Luv);
    cv::imshow("inputImg - Luv", inputImg); // CV_32FC3, LUV, 3 channels
    cerr << "Size: " << inputImg.rows << "x" << inputImg.cols <<endl;
    
    cv::Mat mapped_result = clusterize(true);
    
    // initial values:
    OptimData optim_data;
    optim_data.clusters = &clusters;
    optim_data.graph = &graph;
    optim_data.k = dContrast_k;
    optim_data.w = dRegularization_w;
    
    arma::vec x = arma::zeros(clusters.size(), 1) + 0.5; // Init at 0.5 -> is it necessary?? TODO -> 0.0 sometimes makes negative numbers in results

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    bool success = optim::de(x, opt_fn, &optim_data);
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (success) {
        std::cerr << "de: Optimization completed successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        std::cerr << "de: Optimization completed unsuccessfully." << std::endl;
    }
    arma::cerr << "\nde: solution:\n" << x << arma::endl;

    inputGrey.convertTo(inputGrey, CV_32F);
    inputGrey *= 1./255;

    cv::Mat img_result = cv::Mat::zeros(inputImg.rows, inputImg.cols, CV_32FC3);
    cv::Mat img_result2 = inputGrey.clone();
    for (int r = 0; r < inputImg.rows; r++)
    {
        for (int c = 0; c < inputImg.cols; c++)
        {
            for (int j = 0; j < clusters.size(); j++)
            {
                img_result2.at<float>(r, c) += (x.at(j) - mapped_result.at<float>(j)) * clusters[j].getWeight(cv::Mat(inputImg.at<cv::Vec3f>(r, c)));
                
                if (clusters[j].isPixelOwner(r, c))
                {
                    float optimized_color = (float)abs(x.at(j) - 1); // TODO: results are not between 0.0 and 1.0 -> maybe it is ok for blending...
                    img_result.at<cv::Vec3f>(r, c)[0] = optimized_color;
                    img_result.at<cv::Vec3f>(r, c)[1] = optimized_color;
                    img_result.at<cv::Vec3f>(r, c)[2] = optimized_color;
                }
            }
            //img_result2.at<float>(r, c) = abs(img_result2.at<float>(r, c));
        }
    }
    cv::imshow("optimized", img_result);
    cv::imshow("optimized_2", img_result2);
    
    cv::cvtColor(inputImg, inputImg, cv::COLOR_Luv2BGR);
    inputImg *= 255;
    inputImg.convertTo(inputImg, CV_8UC3);
    cv::imshow("inputImg - back", inputImg); // CV_8UC3, BGR, 3 channels, values 0 to 255

    cv::waitKey();
    return 0;
}

double opt_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    OptimData *optim_data = (OptimData*) opt_data;
    std::vector<Cheryl11::Cluster> *clusters = (std::vector<Cheryl11::Cluster>*)optim_data->clusters;
    Cheryl11::Graph *graph = (Cheryl11::Graph*)optim_data->graph;
    const vector<Graph::Edge> &edges = graph->getEdges();

    // equation (2)
    double Et = 0.0;
    for (int i = 0; i < graph->getEdgesCount(); i++)
    {
        double Tau_ij = edges.at(i).lenght;
        
        double x_j = vals_inp(edges.at(i).c1);
        double x_i = vals_inp(edges.at(i).c0);
        
        double m_ij = abs(clusters->at(edges.at(i).c0).getMappedColor() - clusters->at(edges.at(i).c1).getMappedColor());
        double k = optim_data->k; // 0.1 - 1.0 ... controls the amount of contrast enhancement
        double mapped_L2 = m_ij;
        double psi_ij = edges.at(i).psi;
        double o_ij = edges.at(i).colorL2;
        double a_ij = k * (abs(psi_ij * (o_ij - mapped_L2)) / graph->getScaleFactor());
        double t_ij = (m_ij + a_ij); // m_ij is the same as M_ij -> equation is shorter...

        // equation (2)
        Et += Tau_ij * pow((x_j - x_i) - t_ij, 2);
    }

    // equation (4)
    double w = optim_data->w; // increase value to preserve the realistic nature of the scene

    double Em = 0.0;
    for (int i = 0; i < clusters->size(); i++)
    {
        double Tau_i = 1 - exp(-1 * (clusters->at(i).nearestClusterPathLenght + 0.1));
        double x_i = vals_inp(i);
        double m_i = clusters->at(i).getMappedColor();
        
        // equation (4)
        Em += Tau_i * pow(x_i - m_i, 2);
    }
    
    // equation (1)
    double obj_val = Et + w * Em;

    return obj_val;
}

cv::Mat TMOCheryl11::clusterize(bool showClusteredImg = false)
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
    // Important! - k-means are not deterministic
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
        
        clusters.at(i).makeCovarianceMatrix();
    }
    
    // Set color center and nearest cluster
    for (int i = 0; i < iClusterCount; i++)
    {
        // Color center
        cv::Mat tmpColor(1, 3, CV_32F);
        tmpColor.at<float>(0, 0) = centers(i, 0);
        tmpColor.at<float>(0, 1) = centers(i, 1);
        tmpColor.at<float>(0, 2) = centers(i, 2);
        clusters.at(i).setColorCenter(tmpColor);
        
        // Nearest Cluster
        cv::Mat average = clusters.at(i).getAverageCoordinates();
        float tmpX = average.at<float>(0, 1);
        float tmpY = average.at<float>(0, 0);
        float tmpLenghtPath;
        for (int j = 0; j < iClusterCount; j++)
        {
            if (j != i)
            {
                tmpLenghtPath = sqrt(pow(tmpX - clusters.at(j).getAverageCoordinates().at<float>(0, 1), 2) + pow(tmpY - clusters.at(j).getAverageCoordinates().at<float>(0, 0), 2));
                if (tmpLenghtPath < clusters.at(i).nearestClusterPathLenght || clusters.at(i).nearestClusterPathLenght == 0)
                {
                    clusters.at(i).nearestClusterPathLenght = tmpLenghtPath;
                }
            }
        }
        //cerr << "Cluster nearest ("<<i<<"): " << clusters.at(i).nearestClusterPathLenght << endl;
    }
    
    makeGraph();
    
    const vector<Graph::Edge> &edges = graph.getEdges();
    double max_scale_factor = 0.0;
    for (int i = 0; i < edges.size(); i++)
    {
        cv::Mat average = clusters.at(edges.at(i).c0).getAverageCoordinates();
        cv::Point2f p1(average.at<float>(0, 1), average.at<float>(0, 0));
        average = clusters.at(edges.at(i).c1).getAverageCoordinates();
        cv::Point2f p2(average.at<float>(0, 1), average.at<float>(0, 0));

        cv::line(imgResult, p1, p2, cv::Scalar(255, 0, 0), 1, cv::LineTypes::LINE_AA);
        
        //cerr << edges.at(i).lenght << " " << edges.at(i).c0 << "," << edges.at(i).c1 << endl;
        
        double m_ij = abs(clusters.at(edges.at(i).c0).getMappedColor() - clusters.at(edges.at(i).c1).getMappedColor());
        double kappa = 80.0;
        double c = dJND; // JND try to change it between 0.0 and 1.0
        double mapped_L2 = m_ij;
        double psi_ij = 1.0 / (1.0 + exp(-kappa * (mapped_L2 - c)));
        graph.setPsi(i, psi_ij);
        double o_ij = edges.at(i).colorL2;
        if (abs(psi_ij * (o_ij - mapped_L2)) > max_scale_factor) {
            max_scale_factor = abs(psi_ij * (o_ij - mapped_L2));
        }
    }
    graph.setScaleFactor(max_scale_factor);
    
    cv::cvtColor(imgResult, imgResult, cv::COLOR_Luv2BGR);
    //imgResult *= 255;
    //imgResult.convertTo(imgResult, CV_8UC3);
    
    cv::cvtColor(imgResult, imgResult, cv::COLOR_BGR2GRAY);
    
    if (showClusteredImg) {
        imshow("clusters", imgResult);
        
        cv::Mat test = clusters[0].getClusterImage(); // 13 for sun
        cv::cvtColor(test, test, cv::COLOR_Luv2BGR);
        test *= 255;
        test.convertTo(test, CV_8UC3);
        imshow("test", test);
        
        cv::waitKey();
    }
    
    //return imgResult;
    
    cv::Mat c = centers.reshape(3, iClusterCount);
    cv::cvtColor(c, c, cv::COLOR_Luv2BGR);
    cv::cvtColor(c, c, cv::COLOR_BGR2GRAY);
    return c;
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
                float lenght = sqrt(pow(
                                 clusters.at(i).getAverageCoordinates().at<float>(0, 0) -
                                 clusters.at(h).getAverageCoordinates().at<float>(0, 0),
                                2)
                             +
                             pow(
                                 clusters.at(i).getAverageCoordinates().at<float>(0, 1) -
                                 clusters.at(h).getAverageCoordinates().at<float>(0, 1),
                                2)
                            );// edge lenght
                float colorL2 = sqrt(pow(
                                 clusters.at(i).getAverageColor().at<float>(0, 0) -
                                 clusters.at(h).getAverageColor().at<float>(0, 0),
                                2)
                             +
                             pow(
                                 clusters.at(i).getAverageColor().at<float>(0, 1) -
                                 clusters.at(h).getAverageColor().at<float>(0, 1),
                                2)
                            +
                             pow(
                                 clusters.at(i).getAverageColor().at<float>(0, 2) -
                                 clusters.at(h).getAverageColor().at<float>(0, 2),
                                2)
                            );// edge lenght
                graph.addEdge(i, h, lenght, colorL2);
            }
        }
        //cerr << "Histogram counter: " << histogramValuesCounter << " threshold: " << threshold << endl;
    }
}