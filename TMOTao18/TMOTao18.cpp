/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       Diploma thesis                                         *
*                       Author: Matus Bicanovsky                               *
*                       Brno 2024                                              *
*                                                                              *
*                       Implementation of the TMOTao18 class                   *
*       Video Decolorization Using Visual Proximity Coherence Optimization     *
*******************************************************************************/
/* --------------------------------------------------------------------------- *
 * TMOTao18.cpp: implementation of the TMOTao18 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOTao18.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOTao18::TMOTao18()
{
	SetName(L"Tao18");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

TMOTao18::~TMOTao18()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOTao18::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pY, px, py;

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			pY *= dParameter; // Parameters can be used like
							  // simple variables

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}

float TMOTao18::computeEntropy(const cv::Mat& hist)
{
	float entropy = 0.0;
	float total = cv::sum(hist)[0];
	for (int i = 0; i < hist.rows; i++)
	{
		float p = hist.at<float>(i) / total;
		if (p != 0)
		{
			entropy -= p * log(p);
		}
	}
	return entropy;
}

void TMOTao18::computeProximity(const cv::Mat& currentFrame, const cv::Mat& previousFrame, float& deltaL, float& deltaC)
{
	cv::Mat currLab, prevLab;
	//convert to Lab color space
	cv::cvtColor(currentFrame, currLab, cv::COLOR_BGR2Lab);
	cv::cvtColor(previousFrame, prevLab, cv::COLOR_BGR2Lab);

	//calculate difference between current and previous frame
	cv::Mat diffLab = cv::abs(currLab - prevLab);

	//calculate the histogram of the difference
	std::vector<cv::Mat> diffLabChannels(3);
	cv::normalize(diffLab, diffLab, 0.0, 255.0, cv::NORM_MINMAX, CV_32F);
	cv::split(diffLab, diffLabChannels);
	cv::Mat histL, histA, histB;
	int histSize = 256;
	float range[] = {0,256};
	const float* histRange = {range};
	cv::calcHist(&diffLabChannels[0], 1, 0, cv::Mat(), histL, 1, &histSize, &histRange);
	cv::calcHist(&diffLabChannels[1], 1, 0, cv::Mat(), histA, 1, &histSize, &histRange);
	cv::calcHist(&diffLabChannels[2], 1, 0, cv::Mat(), histB, 1, &histSize, &histRange);

	//calculate the entropy of the histograms to get decolorization proximity values
	deltaL = std::abs(computeEntropy(histL));
	float deltaA = computeEntropy(histA);
	float deltaB = computeEntropy(histB);

	deltaC = sqrt((deltaA * deltaA + deltaB * deltaB) / 2.0);
}

cv::Mat TMOTao18::applyHPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray, double phi)
{
	//convert rgb frames to Lab color space
	cv::Mat currLab, prevLab;
	cv::cvtColor(currentFrame, currLab, cv::COLOR_BGR2Lab);
	cv::cvtColor(previousFrame, prevLab, cv::COLOR_BGR2Lab);

	//split lab channels
	std::vector<cv::Mat> currLabChannels(3), prevLabChannels(3);
	cv::split(currLab, currLabChannels);
	cv::split(prevLab, prevLabChannels);

	//calculate optical flow between current and previous frame
	cv::Mat flow;
	cv::calcOpticalFlowFarneback(prevLabChannels[0], currLabChannels[0], flow, 0.5, 3, 15, 3, 5, 1.2, 0);

	//calculate differential refinement frame D_i
	cv::Mat differentialRefinement = cv::Mat::zeros(currentFrame.size(), CV_64F);
	for(int y = 0; y < currentFrame.rows; y++)
	{
		for(int x = 0; x < currentFrame.cols; x++)
		{
			cv::Point2f flowAt = flow.at<cv::Point2f>(y, x);
			int newX = cv::borderInterpolate(x + flowAt.x, currentFrame.cols, cv::BORDER_REFLECT_101);
            int newY = cv::borderInterpolate(y + flowAt.y, currentFrame.rows, cv::BORDER_REFLECT_101); 
			//calculate chrominance diff
			double chromaDiff = 0.0;
			for(int k = 1; k < 3; k++){  //only a and b channels
				chromaDiff += currLabChannels[k].at<double>(newY, newX) - prevLabChannels[k].at<double>(y, x);
			}
			chromaDiff /= 2.0;
			//calculate the luminance difference
            double luminanceDiff = currLabChannels[0].at<double>(y, x) - prevLabChannels[0].at<double>(y, x);

            //calculate the differential refinement value
            differentialRefinement.at<double>(y, x) = phi * chromaDiff + (1 - phi) * luminanceDiff;
		}
	}
}

std::vector<int> TMOTao18::classify(const std::vector<cv::Vec2f>& proximityValues)
{
	//set number of clusters
	const int k = 3;


	//convert the data points to a cv::Mat of type CV_32F
	fprintf(stderr, "size of proximityValues: %d\n", proximityValues.size());
    cv::Mat dataMat(proximityValues.size(), 2, CV_32F);
    for (size_t i = 0; i < proximityValues.size(); i++) {
        dataMat.at<float>(static_cast<int>(i), 0) = proximityValues[i][0];
		dataMat.at<float>(static_cast<int>(i), 1) = proximityValues[i][1];
		//fprintf(stderr, "Matrices values: %f %f at %d\n", dataMat.at<float>(static_cast<int>(i), 0), dataMat.at<float>(static_cast<int>(i), 1), i);
    }
	fprintf(stderr, "Data matrix created\n");
	// Debug prints to check the matrix dimensions and type
	fprintf(stderr, "Data matrix size: %d x %d\n", dataMat.rows, dataMat.cols);
	fprintf(stderr, "Data matrix type: %d\n", dataMat.type());

	//initialize cluster assignemnts with k-means++
	std::vector<int> clusters(proximityValues.size());
	cv::TermCriteria criteria(cv::TermCriteria::Type(1+2), 10, 1.0);
    cv::kmeans(dataMat, k, clusters, criteria, 3, cv::KMEANS_PP_CENTERS);
	fprintf(stderr, "K-means clustering done\n");
	//print what is stored in clusters
	//for(int i = 0; i < clusters.size(); i++)
	//{
	//	fprintf(stderr, "Cluster %d: %d\n", i, clusters[i]);
	//}
	//calculate mean, covariance and prior probability for each cluster
	std::vector<cv::Vec2f> means(k);
	std::vector<cv::Mat> covariances(k);
	std::vector<double> priors(k);

	for (int j = 0; j < k; j++) {
        cv::Mat clusterData(0, 2, CV_32F);
        for (size_t i = 0; i < proximityValues.size(); i++) {
            if (clusters[i] == j) {
                cv::Mat row = (cv::Mat_<float>(1, 2) << proximityValues[i][0], proximityValues[i][1]);
                clusterData.push_back(row);
            }
        }
		fprintf(stderr, "Cluster data size: %d x %d\n", clusterData.rows, clusterData.cols);
        if (!clusterData.empty()) {
            cv::Mat meanMat, covarMat;
            cv::calcCovarMatrix(clusterData, covarMat, meanMat, cv::COVAR_NORMAL | cv::COVAR_ROWS);
            covarMat = covarMat / (clusterData.rows - 1); // Normalize the covariance matrix
            means[j] = cv::Vec2f(meanMat.at<double>(0, 0), meanMat.at<double>(0, 1));
            covariances[j] = covarMat;
            priors[j] = static_cast<double>(clusterData.rows) / proximityValues.size();
        } else {
            covariances[j] = cv::Mat::eye(2, 2, CV_64F); // Identity matrix if no data points
            means[j] = cv::Vec2f(0, 0);
            priors[j] = 0;
        }
    }
	fprintf(stderr, "Cluster statistics calculated\n");
	// Debug prints to check the initial means, covariances, and priors
    for (int j = 0; j < k; j++) {
		fprintf(stderr, "Initial mean for cluster %d: %lf %lf\n", j, means[j][0], means[j][1]);
		fprintf(stderr, "Initial covariance for cluster %d: %lf %lf %lf %lf\n", j, covariances[j].at<double>(0, 0), covariances[j].at<double>(0, 1), covariances[j].at<double>(1, 0), covariances[j].at<double>(1, 1));
		fprintf(stderr, "Initial prior for cluster %d: %lf\n", j, priors[j]);
    }

	//set penalty parameters
	const double lambda = 1.0;
	const double xi = 1.0;

	bool convergance = false;
	//iteratively update the cluster assignments
	int iter = 0;
	double L = 0.0;
	double prevL = -std::numeric_limits<double>::infinity();
	std::vector<int> newClusters = clusters;
    while (true) {
		int consecutiveCount = 0;
		int lastCluster = newClusters[0];
        // Reassign each data point based on GMM probability
        
        for (size_t i = 0; i < proximityValues.size(); i++) {
            cv::Vec2d ui(proximityValues[i][0], proximityValues[i][1]);
            std::vector<double> clusterProbabilities(k);

            for (int j = 0; j < k; j++) {
                if (newClusters[i] == j) {
        			if (lastCluster == j) {
            			consecutiveCount++;
        			} else {
            			consecutiveCount = 1;
        			}
        			lastCluster = j;
				}
    			double penalty = (consecutiveCount > 5) ? lambda / (1 + std::exp(-xi * (consecutiveCount - 5))) : 0.0;
                clusterProbabilities[j] = gaussianLikelihood(ui, means[j], covariances[j], priors[j]) - penalty;
            }

            newClusters[i] = std::distance(clusterProbabilities.begin(), std::max_element(clusterProbabilities.begin(), clusterProbabilities.end()));
        }

        // Update cluster statistics (mean, covariance, and prior)
        for (int j = 0; j < k; j++) {
            cv::Mat clusterData(0, 2, CV_32F);
            for (size_t i = 0; i < proximityValues.size(); i++) {
                if (newClusters[i] == j) {
                    cv::Mat row = (cv::Mat_<float>(1, 2) << proximityValues[i][0], proximityValues[i][1]);
                    clusterData.push_back(row);
                }
            }
            if (!clusterData.empty()) {
                cv::Mat meanMat;
                cv::calcCovarMatrix(clusterData, covariances[j], meanMat, cv::COVAR_NORMAL | cv::COVAR_ROWS);
                means[j] = cv::Vec2f(meanMat.at<float>(0, 0), meanMat.at<float>(0, 1));
                priors[j] = static_cast<double>(clusterData.rows) / proximityValues.size();
            } else {
                covariances[j] = cv::Mat::eye(2, 2, CV_64F); // Identity matrix if no data points
                means[j] = cv::Vec2f(0, 0);
                priors[j] = 0;
            }
        }

        // Calculate the log-likelihood to check for convergence
        L = logLikelihood(proximityValues, newClusters, means, covariances, priors);
		fprintf(stderr, "Iteration %d prevL %lf currentL %lf diff %lf\n", iter, prevL, L, std::abs(L - prevL));
        // Check for convergence
        if (hasConverged(L, prevL, 1e-3)) {
            break;
        }
		
        prevL = L;
        clusters = newClusters;
		iter++;
    }

    // Return the cluster assignments for all frames
    return clusters;

}

double TMOTao18::gaussianLikelihood(const cv::Vec2d& ui, const cv::Vec2d& mean, const cv::Mat& cov, double prior)
{
    cv::Vec2d diff = ui - mean;
    cv::Mat diffMat = cv::Mat(diff).reshape(1);               //convert diff to a 2x1 matrix
    cv::Mat covInv = cov.inv();
	
    cv::Mat tmpMat = diffMat.t() * covInv * diffMat;
    double tmp = tmpMat.at<double>(0, 0);
    double likelihood = std::exp(-0.5 * tmp) / std::sqrt(std::pow(2 * CV_PI, diffMat.rows) * cv::determinant(cov));
    return likelihood;
}

double TMOTao18::logLikelihood(const std::vector<cv::Vec2f>& dataPoints, const std::vector<int>& clusters, const std::vector<cv::Vec2f>& means, const std::vector<cv::Mat>& covariances, const std::vector<double>& priors)
{
    double L = 0.0;
    for (size_t i = 0; i < dataPoints.size(); i++) {
        int j = clusters[i];
        cv::Vec2d dataPoint(dataPoints[i][0], dataPoints[i][1]);
        L += std::log(priors[j]) * gaussianLikelihood(dataPoint, means[j], covariances[j], priors[j]);
    }
    return L;
}

bool TMOTao18::hasConverged(double L, double prevL, double threshold)
{
    return std::abs(L - prevL) < threshold;
}

int TMOTao18::TransformVideo()
{
	int width = vSrc->GetWidth();
	int height = vSrc->GetHeight();
	cv::VideoCapture vid = vSrc->getVideoCaptureObject();
	cv::Mat currentFrame, previousFrame;
	std::vector<cv::Vec2f> proximityValues;
	int tmp = 0;
	float minDeltaL, minDeltaC, maxDeltaL, maxDeltaC;
	minDeltaC = 100.0;
	minDeltaL = 100.0;
	maxDeltaC = 0.0;
	maxDeltaL = 0.0;
	for(int cnt = 0; cnt < vSrc->GetTotalNumberOfFrames(); cnt++)
	{
		vSrc->GetMatVideoFrame(vid, cnt, currentFrame);
		if(cnt == 0)
		{

		}
		else
		{
			float deltaL, deltaC;
			computeProximity(currentFrame, previousFrame, deltaL, deltaC);
			proximityValues.push_back(cv::Vec2f(static_cast<float>(deltaL), static_cast<float>(deltaC)));
			fprintf(stderr, "Frame %d processed deltaL %f deltaC %f\n", cnt, proximityValues[tmp][0], proximityValues[tmp][1]);
			if(proximityValues[tmp][0] < minDeltaL)
			{
				minDeltaL = proximityValues[tmp][0];
			}
			if(proximityValues[tmp][0] > maxDeltaL)
			{
				maxDeltaL = proximityValues[tmp][0];
			}
			if(proximityValues[tmp][1] < minDeltaC)
			{
				minDeltaC = proximityValues[tmp][1];
			}
			if(proximityValues[tmp][1] > maxDeltaC)
			{
				maxDeltaC = proximityValues[tmp][1];
			}
			tmp++;
		}
		//fprintf(stderr, "Frame %d processed\n", cnt);
		previousFrame = currentFrame.clone();
	}
	fprintf(stderr, "Min deltaL %f Max deltaL %f Min deltaC %f Max deltaC %f\n", minDeltaL, maxDeltaL, minDeltaC, maxDeltaC);
	std::vector<int> classifications = classify(proximityValues);
	std::vector<int> cluster0, cluster1, cluster2;
	for(int i = 0; i < classifications.size(); i++)
	{
		if(classifications[i] == 0)
		{
			cluster0.push_back(i);
		}
		else if(classifications[i] == 1)
		{
			cluster1.push_back(i);
		}
		else
		{
			cluster2.push_back(i);
		}
	}
	fprintf(stderr, "Cluster 0 size: %d\n", cluster0.size());
	fprintf(stderr, "Cluster 1 size: %d\n", cluster1.size());
	fprintf(stderr, "Cluster 2 size: %d\n", cluster2.size());
}
