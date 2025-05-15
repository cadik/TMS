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
	SetName(L"Tao18");					 
	SetDescription(L"Color to grayscale operator for images and video, method from paper Video Decolorization Using Visual Proximity Coherence Optimization");
}

TMOTao18::~TMOTao18()
{
}
//since this operator is being used for videos i did not implement the processing of images, it would be done same as the first frame of video (so LPD with previous frame being black)
//but since i didnt get to work LPD correctly, i did not finish this part, this function is generated automatically when creating new operator plugin in the TMS
int TMOTao18::Transform()
{
	pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	pDst->Convert(TMO_Yxy); // x, y as color information

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
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


//function to compute the entropy of a histogram as mentioned in paper section Decolorization Proximity after equation 1
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
//function for computing proximity values as described in section Decolorization proximity and Algorithm 1 description
void TMOTao18::computeProximity(cv::Mat& currentFrame, cv::Mat& previousFrame, float& deltaL, float& deltaC, cv::Mat& mask)
{
	cv::Mat currLab, prevLab;
	//convert to Lab color space
	cv::cvtColor(currentFrame, currLab, cv::COLOR_BGR2Lab);
	cv::cvtColor(previousFrame, prevLab, cv::COLOR_BGR2Lab);

	//calculate difference between current and previous frame
	cv::Mat diffLab = cv::abs(currLab - prevLab);

	//calculate the histogram of the differences
	std::vector<cv::Mat> diffLabChannels(3);
	cv::normalize(diffLab, diffLab, 0.0, 255.0, cv::NORM_MINMAX, CV_32F);
	cv::split(diffLab, diffLabChannels);
	cv::Mat histL, histA, histB;
	int histSize = 256;
	float range[] = {0,256};
	const float* histRange = {range};
	
	//calculate the histograms of the channels
	cv::calcHist(&diffLabChannels[0], 1, 0, mask, histL, 1, &histSize, &histRange);
	cv::calcHist(&diffLabChannels[1], 1, 0, mask, histA, 1, &histSize, &histRange);
	cv::calcHist(&diffLabChannels[2], 1, 0, mask, histB, 1, &histSize, &histRange);
	


	//calculate the entropy of the histograms to get decolorization proximity values
	deltaL = std::abs(computeEntropy(histL));
	float deltaA = computeEntropy(histA);
	float deltaB = computeEntropy(histB);

	deltaC = sqrt((deltaA * deltaA + deltaB * deltaB) / 2.0);         //calculating deltaC according to equation in the paper in Algorithm 1
}

//function implementing Local Proximity Decolorization strategy according to paper
//however the implementation is not working properly in terms of minimizing Energy function, tried many different approaches none worked and minimized the energy function correctly
//also in the paper they did not mention how to update and select the weights for initial grayscale frame so basic equation is used
cv::Mat TMOTao18::applyLPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray, double beta)
{
	//convert rgb frames to lab
	cv::Mat currLab, prevLab, prevGrayD;
	cv::cvtColor(currentFrame, currLab, cv::COLOR_BGR2Lab);
	cv::cvtColor(previousFrame, prevLab, cv::COLOR_BGR2Lab);
	std::vector<cv::Mat> prevLabChannels(3);
	std::vector<cv::Mat> currLabChannels(3);
	cv::split(prevLab, prevLabChannels);
	cv::split(currLab, currLabChannels);
	//normalizing the channels so when calculating the spatial term we dont have significantly differnt values from grayscale values
	prevLabChannels[0] = prevLabChannels[0] / 100.0;
	currLabChannels[0] = currLabChannels[0] / 100.0;

	//initialize grayscale frame 
	cv::Mat currGray(currentFrame.size(), CV_32F);
	std::vector<cv::Mat> bgrChannels(3);
	cv::split(currentFrame, bgrChannels);
	currGray = 0.299 * bgrChannels[2] + 0.587 * bgrChannels[1] + 0.114 * bgrChannels[0];
	
	//iterative optimization to minimize energy function
	cv::Mat G;   //start with initial grayscale frame
	G = currGray.clone();
	//no further information in the paper about parameters for the conjugate gradient method so i was testing different values, none seemed to have inpact on the result
	const double learningRate = 1e-5; 
	const int maxIterations = 50;     //maximum number of iterations
	const double epsilon = 1e-3;      //convergence threshold
	cv::Mat dir = cv::Mat::zeros(G.size(), CV_32F);         //direction matrix for conjugate gradient method
	cv::Mat prevGrad = cv::Mat::zeros(G.size(), CV_32F);
	float totalEnergy = 0.0;
	float oldTotalEnergy = 0.0;
	for(int i = 0; i < maxIterations; i++)
	{
		cv::Mat gradient = cv::Mat::zeros(G.size(), CV_32F);
		totalEnergy = 0.0;
		for(int y = 0; y < G.rows; y++){
			for(int x = 0; x < G.cols; x++){
				//spatial grad , equation 6 in paper
				float eta_diff, diff;
				if(x+1 < G.cols-1)
				{
					eta_diff = currLabChannels[0].at<float>(y,x) - currLabChannels[0].at<float>(y, x+1);
					diff = G.at<float>(y, x) - G.at<float>(y, x+1) - eta_diff;
				}
				else{
					eta_diff = currLabChannels[0].at<float>(y,x) - currLabChannels[0].at<float>(y, x-1);
					diff = G.at<float>(y, x) - G.at<float>(y, x-1) - eta_diff;
				}
				float spatialGrad = 2 * diff;
				//temporal grad, equation 7 in paper
				float lum_diff = currLabChannels[0].at<float>(y,x) - prevLabChannels[0].at<float>(y,x);
				float diff2 = G.at<float>(y, x) - previousGray.at<float>(y, x) - lum_diff;
				float temporalGrad = 2 * diff2;
				//combining the spatial and temporal energy functions to get total energy, equation 8
				totalEnergy += ((1.0 - beta)*(spatialGrad)) + (beta * (temporalGrad));
				gradient.at<float>(y, x) = ((1.0 - beta) * spatialGrad) + (beta * temporalGrad);
			}
		}
		//some form of conjugate gradient method, not working properly
		if(i == 0){
			dir = -gradient;
		}
		else{
			double betCGM = cv::sum(gradient.mul(gradient))[0] / cv::sum(prevGrad.mul(prevGrad))[0];    //calculating beta for CGM
			dir = -gradient + betCGM * dir; 																//update direction									
		}
		G += learningRate * dir;				//update grayscale frame
		prevGrad = gradient.clone();
		float energyDiff = abs(totalEnergy - oldTotalEnergy);
		//check for convergence
		if(energyDiff < epsilon){
			break;
		}
		oldTotalEnergy = totalEnergy;
	}
	cv::Mat result = G.clone();
	return result;
}
//function implementing Motion Proximity Decolorization strategy according to paper
cv::Mat TMOTao18::applyMPD(const cv::Mat& currentFrame, const cv::Mat& previousFrame, const cv::Mat& previousGray)
{
	//convert input into CIELab
	cv::Mat currLab, prevLab, prevLabGray;
    cv::cvtColor(currentFrame, currLab, cv::COLOR_BGR2Lab);
    cv::cvtColor(previousFrame, prevLab, cv::COLOR_BGR2Lab);
	//split into channels
	std::vector<cv::Mat> currLabChannels(3), prevLabChannels(3);
	cv::split(currLab, currLabChannels);
	cv::split(prevLab, prevLabChannels);
	//calculate optical flow from L channel of current and previous frame
	cv::Mat flow;
    cv::calcOpticalFlowFarneback(prevLabChannels[0], currLabChannels[0], flow, 0.5, 3, 15, 3, 5, 1.2, 0);   //calculating optical flow between the current and previous frame

	//normalizing the channels
	currLabChannels[0] = currLabChannels[0] / 100.0;
	prevLabChannels[0] = prevLabChannels[0] / 100.0;
	cv::Mat curr_L = currLabChannels[0];
	cv::Mat prev_L = prevLabChannels[0];
	//initialize Ci with the L channel of current frame
    cv::Mat Ci = curr_L.clone();
	Ci.convertTo(Ci, CV_32F);
	//values of sigma_p and sigma_T were not mentioned in the paper, so these are testing values
	const float sigma_p = 0.2;            //noise variance in frame transition
    const float sigma_T = 0.5;            //noise variance in frame transition
	//parameters for gradient ascent
    const float alpha = 1e-5;              //step length for gradient ascent 
    const int maxIterations = 20;          //maximum number of iterations
	for (int iter = 0; iter < maxIterations; iter++) 
	{
        cv::Mat Ci_new = Ci.clone();

        for (int y = 0; y < currentFrame.rows; y++) {
            for (int x = 0; x < currentFrame.cols; x++) {
                cv::Point2f flowAt = flow.at<cv::Point2f>(y, x);
				
                int newX = cv::borderInterpolate(x - flowAt.x, currentFrame.cols, cv::BORDER_REFLECT_101);
                int newY = cv::borderInterpolate(y - flowAt.y, currentFrame.rows, cv::BORDER_REFLECT_101);

                //calculate Mi using the L2-norm, equation 4 
				float currLVal = curr_L.at<float>(y, x);
				float prevLVal = prev_L.at<float>(newY, newX);
				float Mi = std::sqrt(
					std::pow(currLabChannels[0].at<float>(y, x) - prevLabChannels[0].at<float>(newY, newX), 2) +
					std::pow(currLabChannels[1].at<float>(y, x) - prevLabChannels[1].at<float>(newY, newX), 2) +
					std::pow(currLabChannels[2].at<float>(y, x) - prevLabChannels[2].at<float>(newY, newX), 2)
				);
				if(Mi == 0.0){
					Mi += 0.02f;          //avoid division by zero
				}
                //calculate the gradient of P(Ci), this is probably wrong interpretation of equation 5, since results are not as expected
				float d1 = Ci.at<float>(y, x) - previousGray.at<float>(y, x);
				float d2 = Ci.at<float>(y, x) - previousGray.at<float>(newY, newX);
				float gradVal = 2.0f * d1 / (sigma_p * sigma_p) +
								2.0f * d2 / (sigma_T * sigma_T * Mi);
				
                //update Ci using gradient ascent
                Ci_new.at<float>(y, x) = Ci.at<float>(y, x) + alpha * gradVal;
            }
        }
		
        //check for convergence
        if (cv::norm(Ci_new - Ci) < 1e-3) {
            break;
        }

        Ci = Ci_new;
    }
	cv::Mat result;
	double weight = 0.5;
	cv::addWeighted(previousGray, weight, Ci, 1 - weight, 0, result);    //doing weighted sum of prev grayscale and coherence refinement frame, as descibed in paper
    return result;

}
//function implementing High Proximity Decolorization strategy according to paper
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
	//normalizing the channels
	prevLabChannels[0] = prevLabChannels[0] / 100.0;
	currLabChannels[0] = currLabChannels[0] / 100.0;
	currLabChannels[1] = (currLabChannels[1] + 128.0f) / 255.0f;  //a channel (range [-128, 127])
    currLabChannels[2] = (currLabChannels[2] + 128.0f) / 255.0f;  //b channel (range [-128, 127])
	prevLabChannels[1] = (prevLabChannels[1] + 128.0f) / 255.0f;  
    prevLabChannels[2] = (prevLabChannels[2] + 128.0f) / 255.0f; 

	//calculate optical flow between current and previous frame
	cv::Mat flow;
	cv::calcOpticalFlowFarneback(prevLabChannels[0], currLabChannels[0], flow, 0.5, 3, 15, 3, 5, 1.2, 0);

	//calculate differential refinement frame D_i, as in eqaution 3 
	cv::Mat differentialRefinement = cv::Mat::zeros(currentFrame.size(), CV_32F);
	for(int y = 0; y < currentFrame.rows; y++)
	{
		for(int x = 0; x < currentFrame.cols; x++)
		{
			cv::Point2f flowAt = flow.at<cv::Point2f>(y, x);
			int newX = cv::borderInterpolate(x + flowAt.x, currentFrame.cols, cv::BORDER_REFLECT_101);
            int newY = cv::borderInterpolate(y + flowAt.y, currentFrame.rows, cv::BORDER_REFLECT_101); 
			//calculate chrominance diff
			float chromaDiff = 0.0;
			for(int k = 1; k < 3; k++){  //only a and b channels
				chromaDiff += currLabChannels[k].at<float>(newY, newX) - prevLabChannels[k].at<float>(y, x);
			}
			//calculate the luminance difference
            float luminanceDiff = currLabChannels[0].at<float>(y, x) - prevLabChannels[0].at<float>(y, x);

            //calculate the differential refinement value
			//parameter phi is used to control balance between chromatic and luminance information, not set in paper, so it was set to 0.5
            differentialRefinement.at<float>(y, x) = phi * chromaDiff + (1 - phi) * luminanceDiff;
		}
	}
	cv::Mat result;
	cv::add(previousGray, differentialRefinement, result);
	return result;
}
//function for classification of frames for decolorization process
//described in section DC-GMM CLassifier, and Algorithm 2 description
//since i was unable to get the decolorization strategies to work, i was unable to test classifier properly
std::vector<int> TMOTao18::classify(const std::vector<cv::Vec2f>& proximityValues)
{
	//set number of clusters
	const int k = 3;
	//convert the data points to a cv::Mat of type CV_32F
    cv::Mat dataMat(proximityValues.size(), 2, CV_32F);
    for (size_t i = 0; i < proximityValues.size(); i++) {
        dataMat.at<float>(static_cast<int>(i), 0) = proximityValues[i][0];
		dataMat.at<float>(static_cast<int>(i), 1) = proximityValues[i][1];
    }
	
	//initialize cluster assignemnts with k-means
	std::vector<int> clusters(proximityValues.size());
	cv::TermCriteria criteria(cv::TermCriteria::Type(1+2), 10, 1.0);
    cv::kmeans(dataMat, k, clusters, criteria, 3, cv::KMEANS_PP_CENTERS);
	
	//calculate initial mean, covariance and prior probability for each cluster
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
        if (!clusterData.empty()) {
            cv::Mat meanMat, covarMat;
            cv::calcCovarMatrix(clusterData, covarMat, meanMat, cv::COVAR_NORMAL | cv::COVAR_ROWS);    //calculate the covariance matrix
            covarMat = covarMat / (clusterData.rows - 1); 											   //normalize the covariance matrix
            means[j] = cv::Vec2f(meanMat.at<double>(0, 0), meanMat.at<double>(0, 1));                  //calculate the mean
            covariances[j] = covarMat;
            priors[j] = static_cast<double>(clusterData.rows) / proximityValues.size();                //calculate the prior probability
        } else {
            covariances[j] = cv::Mat::eye(2, 2, CV_64F); // covar matrix if no data points
            means[j] = cv::Vec2f(0, 0);
            priors[j] = 0;
        }
    }
	
	//set penalty parameters, the values are not mentioned in the paper, so didnt know what values to use
	//but since decolorization strategies are not working, i didnt tweak and test these values to find optimal ones
	const double lambda = 1.0;
	const double xi = 1.0;

	bool convergance = false;
	//iteratively update the cluster assignments
	int iter = 0;
	double L = 0.0;
	double prevL = -std::numeric_limits<double>::infinity();
	std::vector<int> newClusters = clusters;
	std::vector<int> consecutiveCount(k, 0);
	int potentialCount = 0;
	
    while (true) {
		std::fill(consecutiveCount.begin(), consecutiveCount.end(), 0);
		int lastCluster = newClusters[0];
		
        // Reassign each data point based on GMM probability
        for (size_t i = 0; i < proximityValues.size(); i++) {
            cv::Vec2d ui(proximityValues[i][0], proximityValues[i][1]);
            std::vector<double> clusterProbabilities(k);
			int lastAssigned = (i>0) ? newClusters[i-1] : -1;
            for (int j = 0; j < k; j++) {
				potentialCount = 1;
                if(j == lastAssigned)
				{
					potentialCount = consecutiveCount[j] + 1;
				}
    			double penalty = (potentialCount > 5) ? lambda / (1.0 + std::exp(-xi * (potentialCount - 5))) : 0.0;         //calculate penalty based on the number of consecutive frames in the same cluster
                clusterProbabilities[j] = gaussianLikelihood(ui, means[j], covariances[j], priors[j]) - penalty;               //  the probability of the data point belonging to the cluster
            }

            newClusters[i] = std::distance(clusterProbabilities.begin(), std::max_element(clusterProbabilities.begin(), clusterProbabilities.end()));      //assign the data point to the cluster with the highest probability
			//counting amount of consecutive assignments to same cluster, to know whether to use penalty or not
			if(newClusters[i] == lastAssigned)
			{
				consecutiveCount[newClusters[i]]++;
			}
			else
			{
				consecutiveCount[newClusters[i]] = 1;
			}
		}
        //update cluster statistics (mean, covariance, and prior)
        for (int j = 0; j < k; j++) {
            cv::Mat clusterData(0, 2, CV_32F);
            for (size_t i = 0; i < proximityValues.size(); i++) {
                if (newClusters[i] == j) {
                    cv::Mat row = (cv::Mat_<float>(1, 2) << proximityValues[i][0], proximityValues[i][1]);       //convert the data point to a 1x2 matrix
                    clusterData.push_back(row);
                }
            }
            if (!clusterData.empty()) {
				cv::Mat meanMat, covarMat;
            	cv::calcCovarMatrix(clusterData, covarMat, meanMat, cv::COVAR_NORMAL | cv::COVAR_ROWS);    //calculate the covariance matrix
            	covarMat = covarMat / (clusterData.rows - 1); 											   //normalize the covariance matrix
            	means[j] = cv::Vec2f(meanMat.at<double>(0, 0), meanMat.at<double>(0, 1));                  //calculate the mean
            	covariances[j] = covarMat;
            	priors[j] = static_cast<double>(clusterData.rows) / proximityValues.size();
			} else {
                covariances[j] = cv::Mat::eye(2, 2, CV_64F); 			//identity matrix if no data points
                means[j] = cv::Vec2f(0, 0);
                priors[j] = 0;
            }
        }

        //calculate the log-likelihood to check for convergence
        L = logLikelihood(proximityValues, newClusters, means, covariances, priors);
        //check for convergence
        if (hasConverged(L, prevL, 1e-3)) {
            break;
        }
		
        prevL = L;                  //update the previous log-likelihood
        clusters = newClusters;     //update the cluster assignments
		iter++;
    }

    // Return the cluster assignments for all frames
    return clusters;

}
// function to calculate the likelihood of a data point given a Gaussian distribution
// in Algorithm 2 description, line 12 calculation
double TMOTao18::gaussianLikelihood(const cv::Vec2d& ui, const cv::Vec2d& mean, const cv::Mat& cov, double prior)
{
    cv::Vec2d diff = ui - mean;								  // calculate the difference between the data point and the mean
    cv::Mat diffMat = cv::Mat(diff).reshape(1);               //convert diff to a 2x1 matrix
    cv::Mat covInv = cov.inv();
	
    cv::Mat tmpMat = diffMat.t() * covInv * diffMat;
    double tmp = tmpMat.at<double>(0, 0);
    double likelihood = std::exp(-0.5 * tmp) / std::sqrt(std::pow(2 * CV_PI, diffMat.rows) * cv::determinant(cov));		// calculate the likelihood
    return likelihood;
}
// function to calculate the log-likelihood of the data points given the cluster assignments
// in Algorithm 2 description, line 24 calculation
double TMOTao18::logLikelihood(const std::vector<cv::Vec2f>& dataPoints, const std::vector<int>& clusters, const std::vector<cv::Vec2f>& means, const std::vector<cv::Mat>& covariances, const std::vector<double>& priors)
{
    double L = 0.0;
    for (size_t i = 0; i < dataPoints.size(); i++) {
        int j = clusters[i];
        cv::Vec2d dataPoint(dataPoints[i][0], dataPoints[i][1]);                       // convert the data point to a Vec2d
        L += std::log(priors[j] * gaussianLikelihood(dataPoint, means[j], covariances[j], priors[j]));       // calculate the log-likelihood
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
	//get the forground and background for the video frames
	cv::Ptr<cv::BackgroundSubtractor> bgModel = cv::createBackgroundSubtractorKNN(500, 400.0, false);
	int warmupFrames = 30;
	for(int i = 0; i < std::min(warmupFrames, vSrc->GetTotalNumberOfFrames()); i++)
	{
		//train the background model on 30 first frames
		cv::Mat frame, dummyMask;
		vSrc->GetMatVideoFrame(vid, i, frame);
		bgModel->apply(frame, dummyMask, -1);
	}
	std::vector<cv::Mat> channels;
	cv::Mat currentFrame, previousFrame, previousGray;
	std::vector<cv::Vec2f> proximityValuesFG, proximityValuesBG;


	float minDeltaL, minDeltaC, maxDeltaL, maxDeltaC;
	minDeltaC = 100.0;
	minDeltaL = 100.0;
	maxDeltaC = 0.0;
	maxDeltaL = 0.0;
	for(int cnt = 0; cnt < vSrc->GetTotalNumberOfFrames(); cnt++)
	{
		vSrc->GetMatVideoFrame(vid, cnt, currentFrame);
		if(cnt == 0) //we skip first frame since it doesnt have previous frame
		{

		}
		//compute the proximity values for pairs of current frame and its previous frame, for forground and background separately
		else
		{
			cv::Mat fgMask, bgMask;
			bgModel->apply(currentFrame, fgMask);
			cv::threshold(fgMask, bgMask, 1, 255, cv::THRESH_BINARY_INV);
			float deltaLf, deltaCf, deltaLb, deltaCb;
			computeProximity(currentFrame, previousFrame, deltaLf, deltaCf, fgMask);
			computeProximity(currentFrame, previousFrame, deltaLb, deltaCb, bgMask);

			proximityValuesFG.push_back(cv::Vec2f(static_cast<float>(deltaLf), static_cast<float>(deltaCf)));
			proximityValuesBG.push_back(cv::Vec2f(static_cast<float>(deltaLb), static_cast<float>(deltaCb)));
			
		}
		previousFrame = currentFrame.clone();
	}
	fprintf(stderr, "Proximity values calculated\n");
	//call classifier on computed proximity values to assign decolorization strategy for each frame for background and forground proximity values separately
	std::vector<int> classificationsFG = classify(proximityValuesFG);
	std::vector<int> classificationsBG = classify(proximityValuesBG);
	std::vector<int> classifications;
	for(int i = 0; i < classificationsFG.size(); i++)
	{	
		//strict approach from paper, if classification of FG and BG are same
		if(classificationsFG[i] == classificationsBG[i])
		{
			classifications.push_back(classificationsFG[i]);
		}
		else
		{
			//if classification for FG and BG are different, take the minimum of the two
			int tmp = std::min(classificationsFG[i], classificationsBG[i]);
			classifications.push_back(tmp);
		}
	}
	// helping print to see how many frames are classified as LPD, MPD, HPD
	int count0 = std::count(classifications.begin(), classifications.end(), 0);
	int count1 = std::count(classifications.begin(), classifications.end(), 1);
	int count2 = std::count(classifications.begin(), classifications.end(), 2);
	fprintf(stderr, "\nLPD: %d MPD: %d HPD: %d\n", count0, count1, count2);
	previousFrame = cv::Mat::zeros(height, width, CV_32FC3);
	previousGray = cv::Mat::zeros(height, width, CV_32F);
	cv::Mat result, normResult;
	for(int i = 0; i < vSrc->GetTotalNumberOfFrames(); i++)
	{
		vSrc->GetMatVideoFrame(vid, i, currentFrame);
		if(i == 0)
		{
			//for first frame we use LPD and we dont have classification for it because it doesnt have previous frame, and its previous grayscale frame is all black, authors did not mention how they deal with previous grayscale frame for first frame
			result = applyLPD(currentFrame, previousFrame, previousGray, 0.5);
		}
		else
		{
			//apply the decolorization strategy based on the classification
			if(classifications[i-1] == 0)
			{
				result = applyLPD(currentFrame, previousFrame, previousGray, 0.5);
			}
			else if(classifications[i-1] == 1)
			{
				result = applyMPD(currentFrame, previousFrame, previousGray);
			}
			else
			{
				result = applyHPD(currentFrame, previousFrame, previousGray, 0.5);
			}
		}
		
		normResult = result.clone();
		//clamp values to [0, 1]
		cv::min(normResult, 1.0, normResult);   
		cv::max(normResult, 0.0, normResult);   
		channels.clear();
		channels.push_back(normResult);
		channels.push_back(normResult);
		channels.push_back(normResult);
		cv::Mat finalResult;
		cv::merge(channels, finalResult);
		//storing the final frame
		vDst->setMatFrame(vDst->getVideoWriterObject(), finalResult);
		fprintf(stderr, "\rFrames %d/%d decolorized", i, vSrc->GetTotalNumberOfFrames());
		fflush(stdout);
		previousFrame = currentFrame.clone();
		previousGray = normResult.clone();
	}
	fprintf(stderr, "\n");
	return 0;
}


