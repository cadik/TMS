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

double TMOTao18::computeEntropy(const cv::Mat& hist)
{
	double entropy = 0.0;
	for (int i = 0; i < hist.rows; i++)
	{
		float p = hist.at<float>(i);
		if (p > 0)
		{
			entropy -= p * log(p);
		}
	}
	return entropy;
}

void TMOTao18::computeProximity(const cv::Mat& currentFrame, const cv::Mat& previousFrame, double& deltaL, double& deltaC)
{
	cv::Mat currLab, prevLab;
	//convert to Lab color space
	cv::cvtColor(currentFrame, currLab, cv::COLOR_BGR2Lab);
	cv::cvtColor(previousFrame, prevLab, cv::COLOR_BGR2Lab);

	//calculate difference between current and previous frame
	cv::Mat diffLab = cv::abs(currLab - prevLab);

	//calculate the histogram of the difference
	std::vector<cv::Mat> diffLabChannels(3);
	cv::split(diffLab, diffLabChannels);
	cv::Mat histL, histA, histB;
	int histSize = 256;
	float range[] = {0,256};
	const float* histRange = {range};
	cv::calcHist(&diffLabChannels[0], 1, 0, cv::Mat(), histL, 1, &histSize, &histRange);
	cv::calcHist(&diffLabChannels[1], 1, 0, cv::Mat(), histA, 1, &histSize, &histRange);
	cv::calcHist(&diffLabChannels[2], 1, 0, cv::Mat(), histB, 1, &histSize, &histRange);

	//calculate the entropy of the histograms to get decolorization proximity values
	deltaL = computeEntropy(histL);
	double deltaA = computeEntropy(histA);
	double deltaB = computeEntropy(histB);

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

int TMOTao18::classify(double deltaL, double deltaC)
{
	//set number of clusters
	const int k = 3;

	//initialize data points
	std::vector<cv::Vec2d> dataPoints;
	dataPoints.push_back(cv::Vec2d(deltaL, deltaC));
	//initialize cluster assignemnts with k-means++
	std::vector<int> clusters;
	cv::kmeans(dataPoints, k, clusters, cv::TermCriteria(), 3, cv::KMEANS_PP_CENTERS);

	//calculate mean, covariance and prior probability for each cluster
	std::vector<cv::Vec2d> means(k);
	std::vector<cv::Mat> covariances(k);
	std::vector<double> priors(k);

	 for (int j = 0; j < k; j++) {
        std::vector<cv::Vec2d> clusterData;
        for (size_t i = 0; i < dataPoints.size(); i++) {
            if (clusters[i] == j) {
                clusterData.push_back(dataPoints[i]);
            }
        }
        cv::Mat clusterMat(clusterData.size(), 2, CV_64F, &clusterData[0]);
        cv::calcCovarMatrix(clusterMat, covariances[j], means[j], cv::COVAR_NORMAL | cv::COVAR_ROWS);
        priors[j] = static_cast<double>(clusterData.size()) / dataPoints.size();
    }

	//set penalty parameters
	const double lambda = 1.0;
	const double xi = 1.0;

	bool convergance = false;
	//iteratively update the cluster assignments
	double prevL = -std::numeric_limits<double>::infinity();
    while (true) {
        //assign each data point based on GMM probability
        std::vector<int> newClusters(dataPoints.size());
        for (size_t i = 0; i < dataPoints.size(); i++) {
            cv::Vec2d ui = dataPoints[i];
            std::vector<double> clusterProbabilities(k);

            for (int j = 0; j < k; j++) {
                int mj = std::count(newClusters.begin(), newClusters.end(), j);
				double penalty = 0.0;
				if(mj > 5){
					penalty = lambda / (1 + std::exp(-xi * (mj - 5)));
				}
                clusterProbabilities[j] = gaussianLikelihood(ui, means[j], covariances[j], priors[j]) - penalty;
            }

            newClusters[i] = std::distance(clusterProbabilities.begin(), std::max_element(clusterProbabilities.begin(), clusterProbabilities.end()));
        }

        // Update cluster statistics (mean, covariance, and prior)
        for (int j = 0; j < k; j++) {
            std::vector<cv::Vec2d> clusterData;
            for (size_t i = 0; i < dataPoints.size(); i++) {
                if (newClusters[i] == j) {
                    clusterData.push_back(dataPoints[i]);
                }
            }
            cv::Mat clusterMat(clusterData.size(), 2, CV_64F, &clusterData[0]);
            cv::calcCovarMatrix(clusterMat, covariances[j], means[j], cv::COVAR_NORMAL | cv::COVAR_ROWS);
            priors[j] = static_cast<double>(clusterData.size()) / dataPoints.size();
        }

        //calculate the log-likelihood to check for convergence
        double L = logLikelihood(dataPoints, newClusters, means, covariances, priors);

        //check for convergence
        if (hasConverged(L, prevL, 1e-6)) {
            break;
        }
        prevL = L;
        clusters = newClusters;
    }

    // Return the cluster assignment for the current frame
    return clusters.back();

}

double TMOTao18::gaussianLikelihood(const cv::Vec2d& ui, const cv::Vec2d& mean, const cv::Mat& cov, double prior)
{
    cv::Vec2d diff = ui - mean;
    cv::Mat diffMat = cv::Mat(diff).reshape(1);               //convert diff to a 2x1 matrix
    cv::Mat covInv = cov.inv();
    double tmp = diffMat.t().dot(covInv * diffMat);
    double likelihood = std::exp(-0.5 * tmp) / std::sqrt(std::pow(2 * CV_PI, diffMat.rows) * cv::determinant(cov));
    return likelihood * prior;
}

double TMOTao18::logLikelihood(const std::vector<cv::Vec2d>& dataPoints, const std::vector<int>& clusters, const std::vector<cv::Vec2d>& means, const std::vector<cv::Mat>& covariances, const std::vector<double>& priors)
{
    double L = 0.0;
    for (size_t i = 0; i < dataPoints.size(); i++) {
        int j = clusters[i];
        L += std::log(gaussianLikelihood(dataPoints[i], means[j], covariances[j], priors[j]));
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
	for(int cnt = 0; cnt < vSrc->GetTotalNumberOfFrames(); cnt++)
	{
		vSrc->GetMatVideoFrame(vid, cnt, currentFrame);
		if(cnt == 0)
		{

		}
		else
		{
			double deltaL, deltaC;
			computeProximity(currentFrame, previousFrame, deltaL, deltaC);
			int decolorStrat = classify(deltaL, deltaC);
			
		}
		previousFrame = currentFrame.clone();
	}
}
