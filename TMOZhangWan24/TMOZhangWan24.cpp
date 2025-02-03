/*
* Author of the code: Mária Nováková (xnovak2w@stud.fit.vutbr.cz)
* Author of the paper: Zhang L., Wan Y.
* Date of implementation: 28.04.2024
* Description: Implementation of color-to-gray image conversion using salient colors and radial basis functions.
*/

#include "TMOZhangWan24.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOZhangWan24::TMOZhangWan24()
{
	SetName(L"ZhangWan24");
	SetDescription(L"Color-to-gray image conversion using salient colors and radial basis functions");

	max_k.SetName(L"max_k");
	max_k.SetDescription(L"Maximum number of quantized colors (clusters)");
	max_k.SetDefault(20.0);
	max_k = 20.0;
	max_k.SetRange(20.0, 40.0);
	this->Register(max_k);

   sigma.SetName(L"sigma");
	sigma.SetDescription(L"Controls the spread of the Laplace kernel's influence");
	sigma.SetDefault(5.0);
	sigma = 5.0;
	sigma.SetRange(0.0, 100.0);
	this->Register(sigma);
}

TMOZhangWan24::~TMOZhangWan24()
{
}

/////////////////////////////////// HELPERS PART 1. ///////////////////////////////////

// compute PCA, keep only the first eigenvector
Vec3d computePrincipalDirection(const Mat& image) {
   Mat data = image.reshape(1, image.total());

   PCA pca(data, Mat(), PCA::DATA_AS_ROW, 1); 
   Vec3d D_PCA;
   for (int i = 0; i < 3; ++i) {
      D_PCA[i] = pca.eigenvectors.at<double>(0, i);
   }
   return D_PCA;   
}

// expand centroids based on the principal direction and current c0 [FORMULA 3.1 (1)]
void expandCentroids(Vec3d c_0, int &k, Mat img, vector<Vec3d> *centers) {
   Vec3d sigma(1.0/255.0, 1.0/255.0, 1.0/255.0);

   Vec3d D_pca = computePrincipalDirection(img);

   Vec3d N_c1 = c_0 + sigma * D_pca;
   Vec3d N_c2 = c_0 - sigma * D_pca;

   auto it = find(centers->begin(), centers->end(), c_0);
   if (it != centers->end()) {
      centers->erase(it);
   }

   centers->push_back(N_c1);
   centers->push_back(N_c2);
   k++;
}

// compute entropy which detemrined whether the image is natural or synthetic [FORMULA 3.1 (8, 9)]
double Entropy(Mat img){

   vector<int> histogram(256, 0);

   for (int i = 0; i < img.rows; i++) {
      for (int j = 0; j < img.cols; j++) {
         histogram[int(img.at<Vec3d>(i, j)[0]*255)]++;
      }
   }

   int total_pixels = img.rows * img.cols;

   vector<double> p(256, 0);

   for (int i = 0; i < 256; i++) {
      p[i] = (double)histogram[i] / total_pixels;
   }

   double E = 0;
   for (int i = 0; i < 256; i++) {
      if (p[i] > 0.0) {    // skip zero
         E += p[i] * log(p[i]);
      }
   }
   return -E;
}

// [FORMULA 3.1 (2)]
double euclideanDistance(const Vec3d& color1, const Vec3d& color2) {
   double dL = color1[0] - color2[0];
   double da = color1[1] - color2[1];
   double db = color1[2] - color2[2];
   return sqrt(dL * dL + da * da + db * db);
}

// cluster pixels to the nearest centroid based on minimal euclidean distance
vector<vector<pair<Point, Vec3d>>> clusterImage(Mat image, vector<Vec3d> centroids) {
   vector<vector<pair<Point, Vec3d>>> clusters(centroids.size());

   for (int x = 0; x < image.rows; x++) {
      for (int y = 0; y < image.cols; y++) {
         Vec3d pixel = image.at<Vec3d>(x,y);
         double minDistance = numeric_limits<double>::max();
         int assignedCluster;

         for (size_t i = 0; i < centroids.size(); i++) {
               double distance = (euclideanDistance(pixel, centroids[i]));

               if (distance < minDistance) {
                  minDistance = distance;
                  assignedCluster = i;
               }
         }

         clusters[assignedCluster].push_back(make_pair(Point(y, x),pixel));
      }
   }



   return clusters;
}

// update centers based on average color in the cluster
void updateCenters(vector<Vec3d> *centers, vector<vector<pair<Point, Vec3d>>> clusters) {
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
         meanVal += pixel.second;
      }

      Vec3d newCenter = meanVal / static_cast<double>(clusterSize);

      (*centers)[i] = newCenter; 
   }
}

// compute MSE for k quantized colors [FORMULA 3.1 (3)]
double MSE_k(Mat image, vector<Vec3d> centers, vector<vector<pair<Point, Vec3d>>> clusters) {
   double mse = 0.0;

   for (size_t i = 0; i < centers.size(); i++) {
      Vec3d color_cluster = centers[i];

      for (const auto& cluster_pixel : clusters[i]) {
         Point coords = cluster_pixel.first;
         Vec3d color = image.at<Vec3d>(coords);
         double diff = euclideanDistance(color, color_cluster);
         mse += pow(diff, 2);
      }
   }

   return mse / static_cast<double>(image.total());
}

// [FORMULA 3.1 (6)]
double g(Mat image, vector<pair<Point, Vec3d>> cluster){
   double Qx = 0.0;
   for (const auto& cluster_pixel : cluster) {
         Point coords = cluster_pixel.first;
         double color = image.at<Vec3d>(coords)[0];
         Qx += color;

      }

   return Qx / cluster.size();
}

// compute MSEG_k for the grayscale image [FORMULA 3.1 (5)]
double MSEG_k(Mat image, vector<Vec3d> centers, vector<vector<pair<Point, Vec3d>>> clusters){
   double mse = 0.0;

   for (size_t i = 0; i < centers.size(); i++) {

      double gQx = g(image, clusters[i]);

      for (const auto& cluster_pixel : clusters[i]) {
         Point coords = cluster_pixel.first;
         double color = image.at<Vec3d>(coords)[0];

         double diff = abs(color - gQx);
         mse += diff;
      }
   }

   return mse / static_cast<double>(image.total());
}

// compute individual MSE for each cluster [FORMULA 3.1 (5)]
void MSE(vector<Vec3d> centers, vector<vector<pair<Point, Vec3d>>> clusters, vector<double> *individualMSE) {
   
   for (size_t i = 0; i < centers.size(); i++) {
      double tmp = 0;
      Vec3d color_cluster = centers[i];
      for (size_t j = 0; j < clusters[i].size(); j++) {
         Vec3d color = clusters[i][j].second;
         double diff = euclideanDistance(color, color_cluster);
         tmp += pow(diff, 2);
      }
      if (clusters[i].size() > 0) {
         individualMSE->push_back(tmp / clusters[i].size());
      } else {
         individualMSE->push_back(0);
      }
   }
}

// [FORMULA 3.1 (7)]
double M_k(double MSE_k, double MSEG_k) {
   return sqrt(MSE_k * MSEG_k);
}

/////////////////////////////////// MAIN CODE PART 1. ///////////////////////////////////

void TMOZhangWan24::quantizeColors(Mat imageLab, Mat image, Mat grayImage, int &k, double theta_0, double theta_1) {

   int position_mse;

   Scalar meanScalar = mean(imageLab);

   // compute centroid 0, it is an average of all pixels in the image
   Vec3d c_0(static_cast<double>(meanScalar[0]), 
                  static_cast<double>(meanScalar[1]), 
                  static_cast<double>(meanScalar[2]));

   vector<Vec3d> centers;
   vector<Vec3d> best_centers;
   centers.push_back(c_0);
   vector<vector<pair<Point, Vec3d>>> clusters;

   // compute new centroids Nc1,Nc2 from c_0
   expandCentroids(c_0, k, imageLab, &centers);
   vector<vector<pair<Point, Vec3d>>> best_clusters;

   vector<double> individualMSE;

   double mse_k_prev = numeric_limits<double>::max(), mse_k;
   double mseg_k, m_k_prev = numeric_limits<double>::max(), m_k;

   // compute Entropy, determine whether the image is natural or synthetic
   double E = Entropy(grayImage);

   cerr << E << endl;

   while(k <= max_k){ // outer loop, expanding clusters until max_k is reached

      cerr << "Number of clusters: " << k << endl;

      double mse_k_innerloop_prev = numeric_limits<double>::max();
      do { // inner loop, adjust the clusters, compute new centroids, iterate until satisfactory MSE_k is reached

         // assign each pixel to the nearest cluster
         clusters = clusterImage(imageLab, centers);

         // update the centroids based on average color in the cluster
         updateCenters(&centers, clusters);

         // compute MSE for the current k
         mse_k = MSE_k(imageLab, centers, clusters);

         // [FORMULA 3.1 (4)]
         if (((abs(mse_k - mse_k_innerloop_prev)) / mse_k_innerloop_prev) < pow(10, -6)) break;

         mse_k_innerloop_prev = mse_k;

      }while(true);

      // store current centroids and clusters as the best ones
      best_centers = centers;
      best_clusters = clusters;

      cerr << "MSE_K  " << mse_k << " MSE_K_PREV  " << mse_k_prev << endl;

      // outer loop break condition, break if the satisfactory MSE_k/M_k is reached [FORMULA 3.1 (10)]
      if(E > 4){
         mseg_k = MSEG_k(grayImage, centers, clusters);
         m_k = M_k(mse_k, mseg_k);
         if(m_k <= theta_1 && m_k_prev > theta_1){
               break;
         }
      }else{
         if(mse_k <= theta_0 && mse_k_prev > theta_0){
               break;
         }
      }

      // compute individual MSE for each cluster, find the one with the biggest MSE
      MSE(centers, clusters, &individualMSE);
      double minMSE = numeric_limits<double>::min();
      for (size_t i = 0; i < individualMSE.size(); i++){
         if (minMSE < individualMSE[i]){
               minMSE = individualMSE[i];
               position_mse = i;
         }
      }

      // select the cluster with the biggest MSE, set its centroid as c_0
      c_0 = centers[position_mse];

      cerr << "Selected color " << c_0 << " with MSE " << individualMSE[position_mse] << endl;

      // extract pixels that belong to the cluster with the biggest MSE
      const auto &theCluster = clusters[position_mse];
      Mat img_forPCA((int)theCluster.size(), 1, CV_64FC3);
      for (int i = 0; i < (int)theCluster.size(); i++){
         img_forPCA.at<Vec3d>(i, 0) = theCluster[i].second;
      }

      // expand the centroids, add new centroids Nc1, Nc2 based on pixels in the cluster
      expandCentroids(c_0, k, img_forPCA, &centers);

      individualMSE = vector<double>();
      position_mse = -1;

      m_k_prev = m_k;
      mse_k_prev = mse_k;

   } 

   this->centers = best_centers;
   this->clusters = best_clusters;
}



/////////////////////////////////// HELPERS PART 2. ///////////////////////////////////

// [FORMULA 3.2 (12)]
double weightedEuclidean(const Vec3d& color1, const Vec3d& color2){
   double dL = color1[0] - color2[0];
   double da = color1[1] - color2[1];
   double db = color1[2] - color2[2];
   return sqrt(dL * dL * 0.6 + da * da * 0.3 + db * db * 0.1);
}

/////////////////////////////////// MAIN CODE PART 2. ///////////////////////////////////


void TMOZhangWan24::ordering(){

   vector<vector<pair<Point, Vec3d>>> clusters = this->clusters;
   vector<Vec3d> centers = this->centers;

   double min_distance = numeric_limits<double>::min();
   double distance;
   int i0, i1, i2, k = centers.size();

   vector<double> grey(k); 

   for (int i = 0; i < k - 1; i++){
      Vec3d color1 = centers[i];
      for (int j = i+1; j < k; j++){
         Vec3d color2 = centers[j];
         distance = weightedEuclidean(color1, color2);
         // [FORMULA 3.2 (13)]
         if (min_distance < distance){
               i1 = i;
               i2 = j;
               min_distance = distance;
         }
      }
   }

   i0 = centers[i1][0] < centers[i2][0] ? i1 : i2;

   vector<pair<int, double>> storage;

   Vec3d basic_color = centers[i0];

   for (int i = 0; i < k; i++){
      Vec3d color = centers[i];
      distance = weightedEuclidean(color, basic_color);
      storage.push_back(make_pair(i,distance));
   }

   sort(storage.begin(), storage.end(), 
      [](const pair<int, double>& a, const pair<int, double>& b) {
         if (a.second != b.second) {
               return a.second < b.second; 
         }
         return a.first < b.first;
      }
   );

   for (int m = 1; m <= k; m++){
      int index = storage[m-1].first;
      // [FORMULA 3.2 (14)]
      grey[index] = static_cast<double>(m - 1) / (k - 1);
   }

   this->grayvalues = grey;
   
}

/////////////////////////////////// HELPERS PART 3. ///////////////////////////////////

// [FORMULA 3.3 (17)]
double laplaceKernel(const Vec3d& color1, const Vec3d& color2, double sigma){
   double euclid = euclideanDistance(color1, color2);
   return exp(- euclid / sigma);
}

// [FORMULA 3.3 (18)]
double clamp_value(double value) {
   if (value < 0.0f) return 0.0f;
   if (value > 1.0f) return 1.0f;
   return value;
}

double TMOZhangWan24::getGreyValue(Vec3d img_color, vector<double> a){
   vector<Vec3d> centers = this->centers;

   double f_x = 0.0;
   for (size_t i = 0; i < centers.size(); i++){
      Vec3d color1 = centers[i]; 
      // [FORMULA 3.3 (16)]

      f_x += a[i] * laplaceKernel(img_color, color1, sigma);
   }
   return f_x;
}


/////////////////////////////////// MAIN PART 3. ///////////////////////////////////

Mat TMOZhangWan24::createGrayScale(Mat imageLab, Mat image){
   vector<Vec3d> centers = this->centers;
   int k = centers.size();

   vector<double> gray = this->grayvalues; 
   vector<double> a(k); 

   double sum;

   for (int i = 0; i < k; i++){ 
      Vec3d color1 = centers[i]; 
      double greycolor1 = gray[i];
      sum = 0.0;
      for (int j = 0; j < k; j++){
               Vec3d color2 = centers[j]; 
               sum += laplaceKernel(color1, color2, sigma);             
      }
      a[i] = greycolor1 / sum;
   }

   double f_x;

   Mat grayImage = Mat::zeros(imageLab.size(), CV_64FC3);   

   for (int x = 0; x < imageLab.rows; x++){
      for (int y = 0; y < imageLab.cols; y++){
         Vec3d img_color = imageLab.at<Vec3d>(x,y);

         f_x = getGreyValue(img_color, a);
         double value = clamp_value(f_x);

         grayImage.at<Vec3d>(x, y) = Vec3d(value, value, value);
      }
   }

   return grayImage;
   
}

// Convert RGB to XYZ
void rgbToXyz(double r, double g, double b, double &x, double &y, double &z) {
   r = (r > 0.04045) ? pow((r + 0.055) / 1.055, 2.4) : (r / 12.92);
   g = (g > 0.04045) ? pow((g + 0.055) / 1.055, 2.4) : (g / 12.92);
   b = (b > 0.04045) ? pow((b + 0.055) / 1.055, 2.4) : (b / 12.92);

    r *= 100.0;
    g *= 100.0;
    b *= 100.0;

    x = r * 0.4124 + g * 0.3576 + b * 0.1805;
    y = r * 0.2126 + g * 0.7152 + b * 0.0722;
    z = r * 0.0193 + g * 0.1192 + b * 0.9505;
}

// Convert XYZ to LAB
void xyzToLab(double x, double y, double z, double &L, double &a, double &b) {
   double refX = 95.047;
   double refY = 100.000;
   double refZ = 108.883;

   x /= refX;
   y /= refY;
   z /= refZ;

    x = (x > 0.008856) ? pow(x, 1.0 / 3.0) : (7.787 * x + 16.0 / 116.0);
    y = (y > 0.008856) ? pow(y, 1.0 / 3.0) : (7.787 * y + 16.0 / 116.0);
    z = (z > 0.008856) ? pow(z, 1.0 / 3.0) : (7.787 * z + 16.0 / 116.0);

    L = (116.0 * y) - 16.0;
    a = 500.0 * (x - y);
    b = 200.0 * (y - z);
}

// Convert RGB to LAB
void rgbImageToLabImage(double* rgbData, double* labData, int width, int height) {
   for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
         int index = (i * width + j) * 3;
         double r = rgbData[index];
         double g = rgbData[index + 1];
         double b = rgbData[index + 2];

         double x, y, z, LL, aa, bb;
         rgbToXyz(r, g, b, x, y, z);
         xyzToLab(x, y, z, LL, aa, bb);

         labData[index] = LL;
         labData[index + 1] = aa;
         labData[index + 2] = bb;
      }
   }
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOZhangWan24::Transform()
{
   pSrc->Convert(TMO_RGB);
	double *pSourceData = pSrc->GetData();

   int width = pSrc->GetWidth();
   int height = pSrc->GetHeight();

   double *sourceDataRGB = new double[width * height * 3];
   double *sourceDataBW = new double[width * height * 3];
   double *sourceDataLab = new double[width * height * 3];

   memcpy(sourceDataRGB, pSourceData, width * height * 3 * sizeof(double));
   memcpy(sourceDataBW, pSourceData, width * height * 3 * sizeof(double));

	//pSrc->Convert(TMO_LAB); // This is format of Y as luminance
	//pDst->Convert(TMO_LAB); // x, y as color information

   //double *pLabData = pSrc->GetData();

   rgbImageToLabImage(sourceDataRGB, sourceDataLab, width, height);

   Mat   image(height, width, CV_64FC3, sourceDataRGB), 
         imageLab(height, width, CV_64FC3, sourceDataLab), 
         imageGray = Mat::zeros(image.size(), CV_64FC3); 

   cerr << image.at<Vec3d>(0, 0) << endl;
   cerr << imageLab.at<Vec3d>(0, 0) << endl;

   for(int x = 0; x < height; x++){
      for(int y = 0; y < width; y++){
         Vec3d color = image.at<Vec3d>(x, y); 
         // similar to opencv rgb2gray
         double grayValue = color[0] * 0.299 + color[1] * 0.587 + color[2] * 0.114;
         imageGray.at<Vec3d>(x, y) = Vec3d(grayValue, grayValue, grayValue);
      } 
   }  

   int k = 1;
   double theta_0 = 0.01;  
   double theta_1 = 0.5; 

   // 3.1 section
   quantizeColors(imageLab, image, imageGray, k, theta_0, theta_1);
   // 3.2 section
   ordering();
   // 3.3 section
   Mat output_greyimage = createGrayScale(imageLab, image);

	double *pDestinationData = pDst->GetData();
	memcpy(pDestinationData, output_greyimage.data, width*height*sizeof(double)*3);

	return 0;
}
