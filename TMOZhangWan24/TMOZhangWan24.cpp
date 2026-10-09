/* --------------------------------------------------------------------------- *
 * TMOZhangWan24.cpp: Color-to-gray image conversion using salient colors and  *
 * radial basis functions                                                      *
 *                                                                             *
 * Paper: L. Zhang, Y. Wan: Color-to-gray image conversion using salient       *
 *        colors and radial basis functions, 2024                              *
 * Author of the implementation: Maria Novakova                                *
 *                                                                             *
 * Where the paper contradicts itself or leaves values out, this               *
 * implementation uses (see the comments in the code):                         *
 *  - RGB centroids with CIELab distance (captions of Fig. 4 and Fig. 5)       *
 *  - all values normalised to [0,1], Lab as L/100, (a+128)/255, (b+128)/255   *
 *  - squared gray error in MSEG, Eq.(5) (Fig. 7 magnitudes)                   *
 *  - absolute relative change in Eq.(4), log2 in Eq.(8)                       *
 *  - Eq.(16) solved as a k x k linear system for the weights a_j              *
 *  - sigma of the Laplace kernel is a parameter (not given in the paper)      *
 * The edgeRepair parameter is an extension, NOT part of the paper.            *
 * --------------------------------------------------------------------------- */

#include "TMOZhangWan24.h"

// CIELab of a BGR image, normalised to [0,1] like the paper normalises all values to [0,1]:
// L/100, (a+128)/255, (b+128)/255 (this is OpenCV's 8-bit Lab encoding divided by 255)
static Mat normalizedLab(const Mat &bgr8) {
    Mat bgr, lab;
    bgr8.convertTo(bgr, CV_32F, 1/255.0);
    cvtColor(bgr, lab, COLOR_BGR2Lab);                       // L in [0,100], a and b around [-128,127]
    lab = lab.reshape(1, static_cast<int>(lab.total()));
    lab.col(0) *= 1/100.0;
    lab.col(1) = (lab.col(1) + 128.0) / 255.0;
    lab.col(2) = (lab.col(2) + 128.0) / 255.0;
    return lab.reshape(3, bgr8.rows);
}

// normalised Lab (see normalizedLab) of one colour given as BGR in [0,1]
static Vec3f labOfBGR(const Vec3f &bgr) {
    Mat px(1, 1, CV_32FC3, Scalar(bgr[0], bgr[1], bgr[2])), lab;
    cvtColor(px, lab, COLOR_BGR2Lab);
    Vec3f t = lab.at<Vec3f>(0, 0);
    return Vec3f(t[0] / 100.0f, (t[1] + 128.0f) / 255.0f, (t[2] + 128.0f) / 255.0f);
}

// RGB values (BGR order, [0,1]) of the pixels of one cluster as an N x 1 three-channel matrix,
// used for the split of Eq.(1), which is done in RGB (paper, Fig. 4 and Fig. 5)
static Mat clusterPixels(const Mat &rgbImage, const vector<pair<Point, Vec3f>> &cluster) {
    Mat data(static_cast<int>(cluster.size()), 1, CV_32FC3);
    for (size_t i = 0; i < cluster.size(); i++) data.at<Vec3f>(static_cast<int>(i)) = rgbImage.at<Vec3f>(cluster[i].first);
    return data;
}

/* --------------------------------------------------------------------------- *
 * Constructor: description of the technique and its parameters                *
 * --------------------------------------------------------------------------- */
TMOZhangWan24::TMOZhangWan24()
{
	SetName(L"ZhangWan24");
	SetDescription(L"Color-to-gray image conversion using salient colors and radial basis functions (Zhang, Wan 2024)");

	maxK.SetName(L"maxK");
	maxK.SetDescription(L"Maximum number of salient (quantized) colors; <2, 40>");
	maxK.SetDefault(40);
	maxK = 40;
	maxK.SetRange(2, 40);
	this->Register(maxK);

	sigma.SetName(L"sigma");
	sigma.SetDescription(L"Sigma of the Laplace kernel (Eq. 17), on the [0,1] scale of the normalised CIELab values; not given in the paper");
	sigma.SetDefault(0.1);
	sigma = 0.1;
	sigma.SetRange(0.001, 10.0);
	this->Register(sigma);

	orderingMethod.SetName(L"ordering");
	orderingMethod.SetDescription(L"Ordering of the salient colors: 1 = rgb2gray ordering (Eq. 11), 2 = distance ordering (Algorithm 2)");
	orderingMethod.SetDefault(2);
	orderingMethod = 2;
	orderingMethod.SetRange(1, 2);
	this->Register(orderingMethod);

	edgeRepair.SetName(L"edgeRepair");
	edgeRepair.SetDescription(L"Repair the gray values of blended (anti-aliased) edge pixels; extension, NOT part of the paper");
	edgeRepair.SetDefault(false);
	edgeRepair = false;
	this->Register(edgeRepair);
}

TMOZhangWan24::~TMOZhangWan24()
{
}

/* --------------------------------------------------------------------------- *
 * The conversion: Sec. 3.1 (salient colors), Sec. 3.2 (ordering) and          *
 * Sec. 3.3 (RBF interpolation)                                                *
 * --------------------------------------------------------------------------- */
int TMOZhangWan24::Transform()
{
	pSrc->Convert(TMO_RGB);
	pDst->Convert(TMO_RGB);

	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();
	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();

	// the method works on 8-bit colors (delta = 1/255 in Eq. 1, 8-bit histogram for the entropy);
	// HDR input is scaled to [0,1] by its maximum first
	double maxValue = 0.0;
	for (int i = 0; i < width * height * 3; i++)
		maxValue = max(maxValue, pSourceData[i]);
	double scale = (maxValue > 1.0) ? 1.0 / maxValue : 1.0;

	Mat image(height, width, CV_8UC3); // BGR, as used by the OpenCV conversions below
	for (int y = 0; y < height; y++)
	{
		pSrc->ProgressBar(y, height);
		for (int x = 0; x < width; x++)
		{
			double *p = pSourceData + 3 * (y * width + x);
			image.at<Vec3b>(y, x) = Vec3b(saturate_cast<uchar>(p[2] * scale * 255.0),
										   saturate_cast<uchar>(p[1] * scale * 255.0),
										   saturate_cast<uchar>(p[0] * scale * 255.0));
		}
	}

	int k = 1;
	const float theta_0 = 0.0004f;  // paper value, synthetic images (values normalised to [0,1])
	const float theta_1 = 0.00065f; // paper value, natural images (values normalised to [0,1])

	quantizeColors(image, k, maxK, theta_0, theta_1);                             // Sec. 3.1
	ordering(orderingMethod);                                                     // Sec. 3.2
	Mat F = createGrayScale(image, static_cast<float>(sigma.GetDouble()), edgeRepair); // Sec. 3.3

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			double g = F.at<float>(y, x);
			double *d = pDestinationData + 3 * (y * width + x);
			d[0] = d[1] = d[2] = g;
		}
	}
	pSrc->ProgressBar(height, height);

	return 0;
}

/************ PART 1 ************/

/**
 * @brief First step of the conversion, determine salient colors 
 * 
 * @param image input image
 * @param k current number of salient colors
 * @param max_k maximum number for quantizied colors
 * @param theta_0 threshold value for MSE for synthetic images
 * @param theta_1 threshold value for Mk for natural images
 */
void TMOZhangWan24::quantizeColors(Mat &image, int &k, int max_k, float theta_0, float theta_1) {
    
    cv::Mat imagefloat, imageLab, grayImage, grayImage8;

    // image to grayscale
    image.convertTo(imagefloat, CV_32F, 1/255.0);
    cvtColor(imagefloat, grayImage, COLOR_BGR2GRAY);
    cvtColor(image, grayImage8, COLOR_BGR2GRAY); // 8-bit gray for the entropy histogram

    // image to CIE colors space, normalised to [0,1]
    imageLab = normalizedLab(image);

    // colours are quantised in RGB with the perceptual distance in CIELab (paper, Fig. 4b and 5b "our method"):
    // centroids are mean RGB colours, stored as their Lab values so that all distances are CIELab distances (Eq. 2)
    this->rgbImage = imagefloat;

    // determine type of the image (synthetic vs natural)
    float E = Entropy(grayImage8);
    bool classif = classification(E);

    // determine first centroid color, it is mean of all values in the image
    Scalar meanScalar = mean(imagefloat);
    float mse_k_prev, m_k_prev;

    Vec3f c_0 = labOfBGR(Vec3f(static_cast<float>(meanScalar[0]), 
                    static_cast<float>(meanScalar[1]), 
                    static_cast<float>(meanScalar[2])));

    vector<Vec3f> centers;
    centers.push_back(c_0);

    // create initial clusters
    vector<vector<pair<Point, Vec3f>>> clusters = clusterImage(imageLab, centers);

    float prevMSE_k = numeric_limits<float>::max();
    float mse_k;

    // adjust clusters and recompute centers
    do {
        mse_k = MSE_k(imageLab, centers, clusters);

        // STEP 5
        if (mse_k == 0 || ((abs(mse_k - prevMSE_k)) / prevMSE_k) < pow(10, -6)) break;

        clusters = clusterImage(imageLab, centers);

        actualizeCenters(&centers, clusters);

        prevMSE_k = mse_k;

    } while (true); // STEP 6

    float mseg_k = MSEG_k(grayImage, centers, clusters);

    float m_k = M_k(mse_k, mseg_k);

    mse_k_prev = mse_k;
    m_k_prev = m_k;

    vector<float> individualMSE;

    // save the best colors
    vector<Vec3f> best_centers;
    vector<vector<pair<Point, Vec3f>>> best_clusters;

    int position_mse;
    
    // create new centroids from selected centroid (in this case c_0)
    // (take the stored centre: after the updates above it is no longer bit-identical to c_0,
    //  and expandCentroids() finds the centre to replace by exact comparison)
    c_0 = centers[0];
    expandCentroids(c_0, k, clusterPixels(rgbImage, clusters[0]), &centers);

    // iterate till maximum number condition is not met
    while (k <= max_k) {
        // STEP 2
        clusters = clusterImage(imageLab, centers);

        // STEP 3
        actualizeCenters(&centers, clusters);

        // STEP 4
        prevMSE_k = numeric_limits<float>::max();
        do {

            mse_k = MSE_k(imageLab, centers, clusters);

            // STEP 5
            if (mse_k == 0 || ((abs(mse_k - prevMSE_k)) / prevMSE_k) < pow(10, -6)) break;

            clusters = clusterImage(imageLab, centers);

            actualizeCenters(&centers, clusters);

            prevMSE_k = mse_k;

        } while (true); // STEP 6

        // STEP 7
        best_centers = centers;
        best_clusters = clusters;

        mseg_k = MSEG_k(grayImage, centers, clusters);

        m_k = M_k(mse_k, mseg_k);

        // select centroid for new expansion, it is the one with biggest MSE
        // (MSE of the converged clusters of this iteration only)
        individualMSE.clear();
        MSE(centers, clusters, &individualMSE);

        float minMSE = -1.0f;
        position_mse = 0;

        for (size_t i = 0; i < individualMSE.size(); i++){
            if (minMSE < individualMSE[i]){
                minMSE = individualMSE[i];
                position_mse = i;
            }
        }

        c_0 = centers[position_mse];

        // principal direction of the pixels that belong to c_0 (paper, Sec. 3.1)
        expandCentroids(c_0, k, clusterPixels(rgbImage, clusters[position_mse]), &centers);

        clusters = clusterImage(imageLab, centers);

        if(classif){ // synthetic
            if(mse_k <= theta_0 && mse_k_prev > theta_0){
                break;
            }
        }else{ // natural
            if(m_k <= theta_1 && m_k_prev > theta_1){
                break;
            }
        }

        m_k_prev = m_k;
        mse_k_prev = mse_k;
    }

    this->centers = best_centers;
    this->clusters = best_clusters;
}

/************ PART 2 ************/

/**
 * @brief Ordering grayscale colors
 * 
 * @param method 1 = order by rgb2gray value of the quantized colors (paper Sec. 3.2.1, Eq. 11)
 *               2 = order by weighted distance from the basic color (paper Sec. 3.2.2, Algorithm 2)
 */
void TMOZhangWan24::ordering(int method){

    vector<vector<pair<Point, Vec3f>>> clusters = this->clusters;
    vector<Vec3f> centers = this->centers;

    float min_distance = -1.0f;
    float distance;
    int i0, i1 = 0, i2 = 0, k = centers.size();

    vector<float> grey(k); 

    // select two colors with the biggest distance (needed by method 2 only)
    for (int i = 0; i < k - 1; i++){
        Vec3f color1 = centers[i];
        for (int j = i+1; j < k; j++){
            Vec3f color2 = centers[j];
            distance = weightedEuclidean(color1, color2);

            if (min_distance < distance){
                i1 = i;
                i2 = j;
                min_distance = distance;
            }
        }
    }

    // select basic color from two colors based on lower L component
    i0 = centers[i1][0] < centers[i2][0] ? i1 : i2;

    vector<pair<int, float>> storage;

    Vec3f basic_color = centers[i0];

    // method 1 (Eq. 11): sort key is the rgb2gray value of the quantizied color
    // method 2 (Eq. 13): sort key is the distance between basic color and the quantizied color
    for (int i = 0; i < k; i++){
        Vec3f color = centers[i];
        distance = (method == 1) ? rgb2grayOfCenter(color) : weightedEuclidean(color, basic_color);
        storage.push_back(make_pair(i,distance));
    }

    // colors colors by distance
    sort(storage.begin(), storage.end(), 
        [](const pair<int, float>& a, const pair<int, float>& b) {
            if (a.second != b.second) {
                return a.second < b.second; 
            }
            return a.first < b.first;
        }
    );

    // assign evenly spaced grey colors to quantizied colors based on the sorted keys (Eq. 11 / Eq. 14)
    for (int m = 1; m <= k; m++){
        int index = storage[m-1].first;
        grey[index] = static_cast<float>(m - 1) / (k - 1);
    }

    this->grayvalues = grey;
}

/************ PART 3 ************/

/**
 * @brief Create grayscale image based on ordered gray values and radial basis function
 * 
 * @param image input image
 * @param sigma scaling parameter for Laplace Kernel
 * @param edgeFix apply the edge repair (extension, not part of the paper)
 * @return Mat gray values in [0,1] (CV_32F)
 */
Mat TMOZhangWan24::createGrayScale(Mat image, float sigma, bool edgeFix){
    vector<Vec3f> centers = this->centers;
    int k = centers.size();

    vector<float> gray = this->grayvalues; 
    vector<float> a(k);

    // determine weights a_j by solving the k x k linear system of Eq.(16):
    //   g_i = sum_j a_j * phi(x_i, x_j),  i = 1..k   <=>   Phi * a = g
    Mat Phi(k, k, CV_64F), G(k, 1, CV_64F), A;
    for (int i = 0; i < k; i++){
        G.at<double>(i) = gray[i];
        for (int j = 0; j < k; j++){
            Phi.at<double>(i, j) = laplaceKernel(centers[i], centers[j], sigma);
        }
    }
    if (!solve(Phi, G, A, DECOMP_LU)) solve(Phi, G, A, DECOMP_SVD);
    for (int i = 0; i < k; i++) a[i] = static_cast<float>(A.at<double>(i));

    float f_x;

    Mat imageLab = normalizedLab(image);

    // assign gray value to each pixel (Eq. 16-18)
    Mat F(image.rows, image.cols, CV_32F);
    for (int x = 0; x < image.rows; x++){
        for (int y = 0; y < image.cols; y++){
            Vec3f img_color = imageLab.at<Vec3f>(x,y);

            f_x = getGreyValue(img_color, a, sigma);

            F.at<float>(x, y) = clamp(f_x);
        }
    }

    // optional extension, NOT part of the paper
    if (edgeFix) F = repairEdgePixels(image, F);

    return F;
}

/************ PART 1 - HELPER FUNCTIONS ************/

/**
 * @brief Compute PCA direction from the image
 * 
 * @param image input image
 * @return Vec3f direction vector
 */
Vec3f TMOZhangWan24::computePrincipalDirection(const Mat& image) {
    Mat data = image.reshape(1, image.total()); // reshape to a single row per pixel
    data.convertTo(data, CV_32F);

    PCA pca(data, Mat(), PCA::DATA_AS_ROW, 1); // keep only the first principal component
    Vec3f D_PCA;
    for (int i = 0; i < 3; ++i) {
        D_PCA[i] = pca.eigenvectors.at<float>(0, i);
    }
    return D_PCA;   
}

/**
 * @brief Expanding centroids with the given formula, new centroids are computed from selected one with the highest MSE
 * 
 * @param c_0 coordinates of centroid from which two new are computed
 * @param k current number of quantizied colors
 * @param img RGB pixels that belong to c_0, used for the split colour and the PCA direction
 * @param centers centroids (quantizied colors)
 */
void TMOZhangWan24::expandCentroids(Vec3f c_0, int &k, Mat img, vector<Vec3f> *centers) {
    // delta = 1/255, the minimum increment of an 8-bit value normalised to [0,1] (paper, Eq. 1)
    const float delta = 1.0f / 255.0f;

    // compute PCA (in RGB, as the split is done in RGB)
    Vec3f D_pca = computePrincipalDirection(img);

    // RGB colour of c_0 = mean RGB colour of its pixels
    Scalar m = mean(img);
    Vec3f c_0_rgb(static_cast<float>(m[0]), static_cast<float>(m[1]), static_cast<float>(m[2]));

    // compute new centroids / colors (Eq. 1 in RGB), stored as Lab like all centroids
    Vec3f N_c1 = labOfBGR(c_0_rgb + delta * D_pca);
    Vec3f N_c2 = labOfBGR(c_0_rgb - delta * D_pca);

    // remove old color
    auto it = find(centers->begin(), centers->end(), c_0);
    if (it != centers->end()) {
        centers->erase(it);
    }

    // add two new colors
    centers->push_back(N_c1);
    centers->push_back(N_c2);

    // increment current quantizied colors
    k++;
}

/**
 * @brief Compute Euclidean Distance
 * 
 * @param color1 color value
 * @param color2 color value (most of the time quantizied color)
 * @return float distance
 */
float TMOZhangWan24::euclideanDistance(const Vec3f& color1, const Vec3f& color2) {
    float dL = color1[0] - color2[0];
    float da = color1[1] - color2[1];
    float db = color1[2] - color2[2];
    return sqrt(dL * dL + da * da + db * db);
}

/**
 * @brief Create clusters
 * 
 * @param image input image
 * @param centroids centroid colors
 * @return vector<vector<pair<Point, Vec3f>>> clusters of input image
 */
vector<vector<pair<Point, Vec3f>>> TMOZhangWan24::clusterImage(Mat image, vector<Vec3f> centroids) {
    vector<vector<pair<Point, Vec3f>>> clusters(centroids.size());

    for (int x = 0; x < image.rows; x++) {
        for (int y = 0; y < image.cols; y++) {
            Vec3f pixel = image.at<Vec3f>(x,y);
            float minDistance = numeric_limits<float>::max();
            int assignedCluster = 0;

            for (size_t i = 0; i < centroids.size(); i++) {
                float distance = (euclideanDistance(pixel, centroids[i]));

                if (distance < minDistance) {
                    minDistance = distance;
                    assignedCluster = i;
                }
            }

            clusters[assignedCluster].push_back(make_pair(Point(y,x),pixel));
        }
    }

    return clusters;
}

/**
 * @brief Update centers based on current colors of pixels assigned to corresponding cluster
 *        (new center = mean RGB colour of the cluster, stored as its Lab value)
 * 
 * @param centers centroids of the clusters
 * @param clusters clusters of image
 */
void TMOZhangWan24::actualizeCenters(vector<Vec3f> *centers, vector<vector<pair<Point, Vec3f>>> clusters) {
    if (centers->empty()) {
        cerr << "Centers vector is empty." << endl;
        return;
    }

    for (size_t i = 0; i < clusters.size(); ++i) {
        Vec3d meanVal(0.0, 0.0, 0.0);
        size_t clusterSize = clusters[i].size();

        if (clusterSize == 0) {
            cerr << "Cluster " << i << " is empty. Skipping mean calculation." << endl;
            continue;
        }

        for (const auto& pixel : clusters[i]) {
            meanVal += Vec3d(rgbImage.at<Vec3f>(pixel.first));
        }

        Vec3f newCenter = labOfBGR(Vec3f(meanVal / static_cast<double>(clusterSize)));

        (*centers)[i] = newCenter; 
    }
}

/**
 * @brief Compute MSE for individual clusters (for synthetic and natural)
 * 
 * @param centers centroids of the clusters
 * @param clusters clusters of image
 * @param individualMSE vector which holds MSE value for each cluster
 */
void TMOZhangWan24::MSE(vector<Vec3f> centers, vector<vector<pair<Point, Vec3f>>> clusters, vector<float> *individualMSE) {
    
    for (size_t i = 0; i < centers.size(); i++) {
        float tmp = 0;
        Vec3f color_cluster = centers[i];
        for (size_t j = 0; j < clusters[i].size(); j++) {
            Vec3f color = clusters[i][j].second;
            float diff = euclideanDistance(color, color_cluster);
            tmp += diff*diff;
        }
        if (clusters[i].size() > 0) {
            individualMSE->push_back(tmp / clusters[i].size());
        } else {
            individualMSE->push_back(0);
        }
    }
}

/**
 * @brief Compute MSE for whole image (for synthetic and natural)
 * 
 * @param image input image
 * @param centers centroids of the clusters
 * @param clusters clusters of image
 * @return float MSE value
 */
float TMOZhangWan24::MSE_k(Mat image, vector<Vec3f> centers, vector<vector<pair<Point, Vec3f>>> clusters) {
    float mse = 0.0;

    for (size_t i = 0; i < centers.size(); i++) {
        Vec3f color_cluster = centers[i];

        for (const auto& cluster_pixel : clusters[i]) {
            Point coords = cluster_pixel.first;
            Vec3f color = image.at<Vec3f>(coords);
            float diff = euclideanDistance(color, color_cluster);
            mse += diff*diff;
        }
    }

    return mse / static_cast<float>(image.total());
}

/**
 * @brief Compute average gray color of the image's clusters
 * 
 * @param image grayscale image
 * @param cluster clusters of image
 * @return float average color of cluster
 */
float TMOZhangWan24::g(Mat image, vector<pair<Point, Vec3f>> cluster){
    float Qx = 0.0;
    for (const auto& cluster_pixel : cluster) {
            Point coords = cluster_pixel.first;
            float color = image.at<float>(coords);
            Qx += color;
        }
    return Qx / cluster.size();
}

/**
 * @brief Compute MSE for grayscaled image
 * 
 * @param image grayscale image
 * @param centers centroids of the clusters
 * @param clusters clusters of image
 * @return float MSE value
 */
float TMOZhangWan24::MSEG_k(Mat image, vector<Vec3f> centers, vector<vector<pair<Point, Vec3f>>> clusters){
    float mse = 0.0;

    for (size_t i = 0; i < centers.size(); i++) {

        float gQx = g(image, clusters[i]);

        for (const auto& cluster_pixel : clusters[i]) {
            Point coords = cluster_pixel.first;
            float color = image.at<float>(coords);

            // squared difference between color in the grayscale image and computed average graycolor in the cluster
            // (Eq. 5 prints |.|, but MSEG is a mean SQUARE error and Fig. 7 shows it on the scale of MSE)
            float diff = color - gQx;
            mse += diff * diff;
        }
    }

    return mse / static_cast<float>(image.total());
}

/**
 * @brief Compute Entropy value, used for determinings type of the image
 * 
 * @param img input image
 * @return float entropy value
 */
float TMOZhangWan24::Entropy(Mat img){

    vector<int> histogram(256, 0);

    for (int i = 0; i < img.rows; i++) {
        for (int j = 0; j < img.cols; j++) {
            histogram[img.at<uchar>(i, j)]++;
        }
    }

    int total_pixels = img.rows * img.cols;

    vector<float> p(256);

    // compute probabilities of gray color in the image
    for (int i = 0; i < 256; i++) {
        p[i] = (float)histogram[i] / total_pixels;
    }

    float E = 0;
    for (size_t i = 0; i < p.size(); i++){
        if (p[i] > 0) E += p[i]*log2(p[i]); // bits (max 8 for 8-bit image), 0*log(0) is taken as 0
    }
    
    return -E;
}

/**
 * @brief Determine type of the image
 * 
 * @param E Entropy value
 * @return true return value for synthetic image
 * @return false return value for natural image
 */
bool TMOZhangWan24::classification(float E){
    return(E <= 4 ? true : false); // true for synthetic, false for natural
}

/**
 * @brief Compute combined MSE score (for natural)
 * 
 * @param MSE_k MSE of quantizied colors
 * @param MSEG_k MSE of grayscale image
 * @return float MSE value
 */
float TMOZhangWan24::M_k(float MSE_k, float MSEG_k) {
    return sqrt(MSE_k * MSEG_k);
}

/************ PART 2 - HELPER FUNCTIONS ************/

/**
 * @brief rgb2gray value of a quantizied color (Eq. 11), in [0,1]
 * 
 * @param labColor quantizied color as stored in centers (normalised Lab: L/100, (a+128)/255, (b+128)/255)
 * @return float gray value
 */
float TMOZhangWan24::rgb2grayOfCenter(const Vec3f& labColor){
    // convert back to RGB through true Lab values (L in [0,100], a and b centered on 0)
    Mat lab(1, 1, CV_32FC3, Scalar(labColor[0] * 100.0f, labColor[1] * 255.0f - 128.0f, labColor[2] * 255.0f - 128.0f)), bgr;
    cvtColor(lab, bgr, COLOR_Lab2BGR);
    Vec3f c = bgr.at<Vec3f>(0, 0);
    return 0.2989f * c[2] + 0.5870f * c[1] + 0.1140f * c[0]; // weights of the paper
}

/**
 * @brief Compute weighted Euclidean Distance
 * 
 * @param color1 color value 
 * @param color2 color value
 * @return float distance
 */
float TMOZhangWan24::weightedEuclidean(const Vec3f& color1, const Vec3f& color2){
    float dL = color1[0] - color2[0];
    float da = color1[1] - color2[1];
    float db = color1[2] - color2[2];
    return sqrt(dL * dL * 0.6 + da * da * 0.3 + db * db * 0.1);
}

/************ PART 3 - HELPER FUNCTIONS ************/

/**
 * @brief EXTENSION, NOT PART OF THE PAPER: repair the gray value of blended (anti-aliased) edge pixels.
 *
 * The method maps colours to gray values without looking at the position of a pixel. A pixel on the
 * border of two regions often has a blended colour c_p = (1-t)*c_A + t*c_B of the colours on both sides,
 * and such a colour can be mapped to a gray value outside the range of the two sides, which shows up
 * as thin dark or white lines along region borders.
 *
 * For every pixel, pairs of pixels on opposite sides (distance 1 and 2, horizontal, vertical and both
 * diagonals) are tested. If the pixel's RGB colour lies on the segment between the two colours
 * (residual <= 0.04, 0 <= t <= 1) and the two colours really differ (distance >= 0.1), the pair with
 * the largest colour difference is used. If the pixel's gray value lies outside the gray range of that
 * pair (by more than 4/255), or differs from the same blend of the two gray values, (1-t)*g_A + t*g_B,
 * by more than 8/255, it is replaced by that blend. Pixels with 4 or more neighbours of nearly the same
 * colour (difference < 0.02) belong to a flat region and are skipped: a region colour can happen to lie
 * between two other colours, which would otherwise be mistaken for a blend where several regions meet.
 * All other pixels are left unchanged. (Correcting
 * only out-of-range pixels left the in-range ones of the same border unchanged, which gave stair-stepped
 * borders.)
 *
 * @param image input image (BGR, 8-bit)
 * @param F gray values in [0,1] from Eq.(18)
 * @return Mat repaired gray values in [0,1]
 */
Mat TMOZhangWan24::repairEdgePixels(const Mat &image, const Mat &F){
    const float maxResidual = 0.04f, minEdge = 0.1f, eps = 4.0f / 255.0f, maxBlendError = 8.0f / 255.0f;
    const float sameColour = 0.01f, newRegionColour = 0.02f;
    const int dirs[4][2] = {{0, 1}, {1, 0}, {1, 1}, {1, -1}};
    const int R = image.rows, C = image.cols;
    auto inside = [&](int x, int y){ return x >= 0 && y >= 0 && x < R && y < C; };

    Mat rgb;
    image.convertTo(rgb, CV_32F, 1/255.0);
    Mat out = F.clone();
    int changed = 0;

    // flat = pixel with 4 or more near-identical neighbours (inside a region, not on a thin blended border)
    Mat flat = Mat::zeros(R, C, CV_8U);
    for (int x = 0; x < R; x++){
        for (int y = 0; y < C; y++){
            Vec3f c = rgb.at<Vec3f>(x, y);
            int same = 0;
            for (int dx = -1; dx <= 1; dx++){
                for (int dy = -1; dy <= 1; dy++){
                    if ((dx == 0 && dy == 0) || !inside(x + dx, y + dy)) continue;
                    Vec3f dn = rgb.at<Vec3f>(x + dx, y + dy) - c;
                    if (dn.dot(dn) < sameColour * sameColour) same++;
                }
            }
            if (same >= 4) flat.at<uchar>(x, y) = 1;
        }
    }
    // member = flat, or the same colour as a flat neighbour. The gray is a function of the colour,
    // so such a pixel already has the gray of its region and must not be changed.
    Mat member = flat.clone();
    for (int x = 0; x < R; x++){
        for (int y = 0; y < C; y++){
            if (member.at<uchar>(x, y)) continue;
            Vec3f c = rgb.at<Vec3f>(x, y);
            for (int dx = -1; dx <= 1; dx++){
                for (int dy = -1; dy <= 1; dy++){
                    if (!inside(x + dx, y + dy) || !flat.at<uchar>(x + dx, y + dy)) continue;
                    Vec3f dn = rgb.at<Vec3f>(x + dx, y + dy) - c;
                    if (dn.dot(dn) < sameColour * sameColour) member.at<uchar>(x, y) = 1;
                }
            }
        }
    }

    // ---- 1) thin border between two regions: blend of the two colours on opposite sides ----
    for (int x = 0; x < R; x++){
        for (int y = 0; y < C; y++){
            if (member.at<uchar>(x, y)) continue;
            Vec3f c = rgb.at<Vec3f>(x, y);
            float bestEdge = -1.0f, bestT = 0.0f, gA = 0.0f, gB = 0.0f;

            for (const auto &d : dirs){
                for (int s = 1; s <= 2; s++){
                    int xa = x - s * d[0], ya = y - s * d[1], xb = x + s * d[0], yb = y + s * d[1];
                    if (!inside(xa, ya) || !inside(xb, yb)) continue;

                    Vec3f cA = rgb.at<Vec3f>(xa, ya), cB = rgb.at<Vec3f>(xb, yb), v = cB - cA;
                    float len2 = v.dot(v);
                    if (len2 < minEdge * minEdge) continue;              // no real edge between A and B

                    float t = (c - cA).dot(v) / len2;
                    if (t < 0.0f || t > 1.0f) continue;                  // not between the two colours
                    Vec3f r = c - (cA + t * v);
                    if (r.dot(r) > maxResidual * maxResidual) continue;  // not a blend of the two colours

                    if (len2 > bestEdge){
                        bestEdge = len2; bestT = t;
                        gA = F.at<float>(xa, ya); gB = F.at<float>(xb, yb);
                    }
                }
            }

            if (bestEdge < 0.0f) continue;
            float g = F.at<float>(x, y);
            float blend = (1.0f - bestT) * gA + bestT * gB;
            if (g < min(gA, gB) - eps || g > max(gA, gB) + eps || fabs(g - blend) > maxBlendError){
                out.at<float>(x, y) = blend;
                changed++;
            }
        }
    }

    // ---- 2) junctions: where 3 or more regions meet, a pixel blends up to 3 colours. Fit its colour as a
    // convex combination of 1-3 colours of the flat regions within radius 2 and use the same weights on
    // their grays. Only these region grays are used, so a wrong border pixel cannot spread.
    int junction = 0;
    for (int x = 0; x < R; x++){
        for (int y = 0; y < C; y++){
            if (member.at<uchar>(x, y)) continue;
            vector<Vec3f> col; vector<float> gr;
            for (int dx = -2; dx <= 2; dx++){
                for (int dy = -2; dy <= 2; dy++){
                    if (!inside(x + dx, y + dy) || !flat.at<uchar>(x + dx, y + dy)) continue;
                    Vec3f cn = rgb.at<Vec3f>(x + dx, y + dy);
                    bool known = false;
                    for (const auto &k : col){ Vec3f dn = cn - k; if (dn.dot(dn) < newRegionColour * newRegionColour) known = true; }
                    if (!known){ col.push_back(cn); gr.push_back(F.at<float>(x + dx, y + dy)); }
                }
            }
            if (col.empty()) continue;

            Vec3f c = rgb.at<Vec3f>(x, y);
            const int n = col.size();
            float bestRes = 1e9f, pred = 0.0f;
            auto take = [&](float res, float p){ if (res < bestRes - 1e-4f){ bestRes = res; pred = p; } };
            for (int i = 0; i < n; i++)                                   // one region colour
                take(norm(c - col[i]), gr[i]);
            for (int i = 0; i < n; i++){                                  // blend of two
                for (int j = i + 1; j < n; j++){
                    Vec3f v = col[j] - col[i]; float l = v.dot(v);
                    if (l < 1e-8f) continue;
                    float t = (c - col[i]).dot(v) / l;
                    if (t < 0.0f || t > 1.0f) continue;
                    take(norm(c - col[i] - t * v), (1.0f - t) * gr[i] + t * gr[j]);
                }
            }
            for (int i = 0; i < n; i++){                                  // blend of three
                for (int j = i + 1; j < n; j++){
                    for (int k = j + 1; k < n; k++){
                        Vec3f u = col[j] - col[i], v = col[k] - col[i], q = c - col[i];
                        float uu = u.dot(u), uv = u.dot(v), vv = v.dot(v), det = uu * vv - uv * uv;
                        if (fabs(det) < 1e-10f) continue;
                        float a = (vv * u.dot(q) - uv * v.dot(q)) / det, b = (uu * v.dot(q) - uv * u.dot(q)) / det;
                        if (a < 0.0f || b < 0.0f || a + b > 1.0f) continue;
                        take(norm(q - a * u - b * v), (1.0f - a - b) * gr[i] + a * gr[j] + b * gr[k]);
                    }
                }
            }
            if (bestRes > maxResidual) continue;                          // colour is not a blend of the regions
            if (fabs(out.at<float>(x, y) - pred) > maxBlendError){
                if (out.at<float>(x, y) == F.at<float>(x, y)) changed++;
                out.at<float>(x, y) = pred;
                junction++;
            }
        }
    }

    cerr << "edge repair (extension): " << changed << " pixels changed ("
         << 100.0 * changed / image.total() << " %), " << junction << " of them at junctions" << endl;
    return out;
}

/**
 * @brief RBF function - Gaussian
 * 
 * @param color1 color value
 * @param color2 color value
 * @param sigma standard deviation
 * @return float 
 */
float TMOZhangWan24::gaussianKernel(const Vec3f& color1, const Vec3f& color2, float sigma){
    float euclid = euclideanDistance(color1, color2);
    return exp(- euclid * euclid / (2*sigma*sigma)); // Eq.(15) is printed without the square
}

/**
 * @brief RBF function - Laplace
 * 
 * @param color1 color value
 * @param color2 color value
 * @param sigma standard deviation
 * @return float 
 */
float TMOZhangWan24::laplaceKernel(const Vec3f& color1, const Vec3f& color2, float sigma){
    float euclid = euclideanDistance(color1, color2);
    return exp(- euclid / sigma);
}

/**
 * @brief Clamp color to [0,1]
 * 
 * @param value color value
 * @return float normalized value
 */
float TMOZhangWan24::clamp(float value) {
    if (value < 0.0f) return 0.0f;
    if (value > 1.0f) return 1.0f;
    return value;
}

/**
 * @brief Compute grayvalue with usage of RBF (this implementation uses laplace)
 * 
 * @param img_color color value
 * @param a weights
 * @param sigma standard deviation
 * @return float gray value
 */
float TMOZhangWan24::getGreyValue(Vec3f img_color, vector<float> a, float sigma){
    vector<Vec3f> centers = this->centers;

    float f_x = 0.0;
    for (size_t i = 0; i < centers.size(); i++){
        Vec3f color1 = centers[i]; 
        f_x += a[i] * laplaceKernel(img_color, color1, sigma);
    }
    return f_x;
}
