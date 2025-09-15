/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*          Opponent Center-Surround Contrast: Colour to Grey Conversion        *
*                                                                              *
*             Author: Peter Zdravecký [xzdrav00 AT stud.fit.vutbr.cz]          *
*                                    Brno 2024                                 *
*                                                                              *
*******************************************************************************/

#include "TMOAlsam20.h"

TMOAlsam20::TMOAlsam20()
{
	SetName(L"Alsam20");
	SetDescription(L"Opponent Center-Surround Contrast: Colour to Grey Conversion");

	pord.SetName(L"ord");
	pord.SetDescription(L"Order of the matrix norm");
	pord.SetDefault(3);
	pord = 3;
	this->Register(pord);

	pw1.SetName(L"w1");
	pw1.SetDescription(L"Linear weight for the first channel");
	pw1.SetDefault(1.0);
	pw1 = 1.0;
	this->Register(pw1);

	pw2.SetName(L"w2");
	pw2.SetDescription(L"Linear weight for the second channel");
	pw2.SetDefault(1.0);
	pw2 = 1.0;
	this->Register(pw2);

	pfile_kernel.SetName(L"kernels");
	pfile_kernel.SetDescription(L"Path to file with kernels \n"
		"		Expected data format: Dc_x, Dc_y, Ds_x, Ds_y\n"
		"		Example:\n"
		"		  1  1  25  25\n"
		"		  3  3  25  25\n");
	this->Register(pfile_kernel);
}

/**
 * @brief Load filter kernels from a file
 * 
 * @param file_path path to the file with kernels
 * @return pair<vector<Mat>, vector<Mat>> Dc and Ds kernels
 */
pair<vector<Mat>, vector<Mat>> TMOAlsam20::load_kernels(string file_path) {
	vector<Mat> Dc_list;
	vector<Mat> Ds_list;
	std::ifstream file(file_path);
	if (!file.is_open()) {
		// Default kernels from the paper
		Mat Dc = Mat::ones(1, 1, CV_64F);
		Mat Ds = Mat::ones(25, 25, CV_64F);
		Dc = Dc / 1;
		Ds = Ds / 625;
		Dc_list.push_back(Dc);
		Ds_list.push_back(Ds);

		Dc = Mat::ones(3, 3, CV_64F);
		Ds = Mat::ones(25, 25, CV_64F);
		Dc = Dc / 9;
		Ds = Ds / 625;
		for (int i = 0; i < 4; ++i) {
			Dc_list.push_back(Dc);
			Ds_list.push_back(Ds);
		}
	
		return make_pair(Dc_list, Ds_list);
	}

	string line;
	while (getline(file, line)) {
		istringstream iss(line);
		int a, b, c, d;
		if (!(iss >> a >> b >> c >> d)) { 
			fprintf(stderr, "Error: Invalid kernel format, expected 4 integers on each line\n");
			break;
		}
		Mat Dc = Mat::ones(a, b, CV_64F);
		Mat Ds = Mat::ones(c, d, CV_64F);
		Dc = Dc / (a * b);
		Ds = Ds / (c * d);
		Dc_list.push_back(Dc);
		Ds_list.push_back(Ds);
	}
	
	return make_pair(Dc_list, Ds_list);
}

TMOAlsam20::~TMOAlsam20()
{
}

/**
 * @brief Calculate the norm of a matrix
 * 
 * @param diff matrix
 * @param ord order of the norm
 * @return double norm
 */
double TMOAlsam20::calculate_norm(Mat diff, int ord = 3)
{
	Mat norm;
	pow(diff, ord, norm);
	norm = abs(norm);
	double s = sum(norm)[0];
	s = pow(s, 1.0 / ord);
	
	return s;
}

/**
 * @brief Calculate the red-green contrast
 * 
 * @param I1 red channel
 * @param I2 green channel
 * @param Dc center kernel
 * @param Ds surround kernel
 * @param ord order of the norm
 * @return contrast
 */
double TMOAlsam20::red_green_contrast(Mat I1, Mat I2, Mat Dc, Mat Ds, int ord = 3) {
	double delta_rg;
	Mat conv1, conv2;
	filter2D(I1, conv1, -1, Dc, Point(-1, -1), 0, BORDER_CONSTANT);
	filter2D(I2, conv2, -1, Ds, Point(-1, -1), 0, BORDER_CONSTANT);
	delta_rg = calculate_norm(conv1 - conv2, ord);
	return delta_rg;
}

/**
 * @brief Calculate the green-red contrast
 * 
 * @param I1 red channel
 * @param I2 green channel
 * @param Dc center kernel
 * @param Ds surround kernel
 * @param ord order of the norm
 * @return contrast
 */
double TMOAlsam20::green_red_contrast(Mat I1, Mat I2, Mat Dc, Mat Ds, int ord = 3) {
	double delta_gr;
	Mat conv1, conv2;
	filter2D(I2, conv1, -1, Dc, Point(-1, -1), 0, BORDER_CONSTANT);
	filter2D(I1, conv2, -1, Ds, Point(-1, -1), 0, BORDER_CONSTANT);
	delta_gr = calculate_norm(conv1 - conv2, ord);
	return delta_gr;
}

/**
 * @brief Calculate the blue-yellow contrast
 * 
 * @param I1 red channel
 * @param I2 green channel
 * @param I3 blue channel
 * @param Dc center kernel
 * @param Ds surround kernel
 * @param w1 weight for the first channel
 * @param w2 weight for the second channel
 * @param ord order of the norm
 * @return contrast
 */
double TMOAlsam20::blue_yellow_contrast(Mat I1, Mat I2, Mat I3, Mat Dc, Mat Ds, float w1 = 1, float w2 = 1, int ord = 3) {
	Mat contrast;
	Mat weighted_sum;
	double delta_by;
	filter2D(I3, contrast, -1, Dc, Point(-1, -1), 0, BORDER_CONSTANT);
	weighted_sum = (w1 * I1 + w2 * I2) / (w1 + w2);
	filter2D(weighted_sum, weighted_sum, -1, Ds, Point(-1, -1), 0, BORDER_CONSTANT);
	contrast = contrast - weighted_sum;
	delta_by = calculate_norm(max(contrast, Mat::zeros(contrast.size(), contrast.type())), ord);
	return delta_by;
}

/**
 * @brief Calculate the weights for the channels
 * 
 * @param I1 red channel
 * @param I2 green channel
 * @param I3 blue channel
 * @param Dc center kernel
 * @param Ds surround kernel
 * @param ord order of the norm
 * @param w1 weight for the first channel
 * @param w2 weight for the second channel
 * @return weights
 */
vector<double> TMOAlsam20::calculate_weights(Mat I1, Mat I2, Mat I3, Mat Dc, Mat Ds, int ord = 3, float w1 = 1, float w2 = 1) {
    double delta_rg = red_green_contrast(I1, I2, Dc, Ds, ord);
    double delta_gr = green_red_contrast(I1, I2, Dc, Ds, ord);
    double delta_by = blue_yellow_contrast(I1, I2, I3, Dc, Ds, w1, w2, ord);

    double delta_sum = delta_rg + delta_gr + delta_by;
    double alpha1 = delta_rg / delta_sum;
    double alpha2 = delta_gr / delta_sum;
    double alpha3 = delta_by / delta_sum;

    vector<double> alphas;
    alphas.push_back(alpha1);
    alphas.push_back(alpha2);
    alphas.push_back(alpha3);
    return alphas;
}

/**
 * @brief Convert the image to greyscale
 * 
 * @param I0 red channel
 * @param I1 green channel
 * @param I2 blue channel
 * @param alpha weights
 * @return greyscale image
 */
Mat TMOAlsam20::grayscale_image(Mat I0, Mat I1, Mat I2, vector<double> alpha) {
    Mat gs;
	gs = I0.mul(alpha[0]) + I1.mul(alpha[1]) + I2.mul(alpha[2]);	
    return gs;
}

/**
 * @brief Iteratively calculate the greyscale image
 * 
 * @param I0 red channel
 * @param I1 green channel
 * @param I2 blue channel
 * @param Dc_list list of center kernels
 * @param Ds_list list of surround kernels
 * @param ord order of the norm
 * @param w1 weight for the first channel
 * @param w2 weight for the second channel
 * @return greyscale image
 */
Mat TMOAlsam20::iterative_greyscale_image(Mat I0, Mat I1, Mat I2, vector<Mat> Dc_list, 
										vector<Mat> Ds_list, int ord = 3, double w1 = 1, double w2 = 1) {
    vector<Mat> g_list;
    for (int i = 0; i < Dc_list.size(); ++i) {
        vector<double> alpha = calculate_weights(I0, I1, I2, Dc_list[i], Ds_list[i], ord, w1, w2);
        Mat g = grayscale_image(I0, I1, I2, alpha);	
        g_list.push_back(g);
    }
	
	Mat g_mean = Mat::zeros(I0.size(), I0.type());
	for (int i = 0; i < g_list.size(); ++i) {
		g_mean += g_list[i];
	}
	g_mean = g_mean / g_list.size();
	return g_mean;
}

int TMOAlsam20::Transform()
{
	pSrc->Convert(TMO_RGB); 
	pDst->Convert(TMO_RGB);

	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();

	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();
	int N = width * height;

	Mat I(height, width, CV_64FC3, pSourceData);
	Mat I0, I1, I2;
	Mat channels[3];
	split(I, channels);
	I0 = channels[0];
	I1 = channels[1];
	I2 = channels[2];

	int ord = pord.GetInt();
	double w1 = pw1.GetDouble();
	double w2 = pw2.GetDouble();
	string file_kernel = pfile_kernel.GetString();

	vector<Mat> Dc_list, Ds_list;
	pair<vector<Mat>, vector<Mat>> kernels = load_kernels(file_kernel);
	Dc_list = kernels.first;
	Ds_list = kernels.second;

	if (Dc_list.size() != Ds_list.size()) {
		fprintf(stderr, "Error: Number of Dc and Ds kernels must be the same\n");
		return -1;
	}

	Mat g = iterative_greyscale_image(I0, I1, I2, Dc_list, Ds_list, ord, w1, w2);
		
	cv::normalize(g, g, 0, 1, cv::NORM_MINMAX, g.type());
	for (int i = 0; i < N; i++) 
	{
		pDestinationData[i * 3] = g.at<double>(i / width, i % width);
		pDestinationData[i * 3 + 1] = g.at<double>(i / width, i % width);
		pDestinationData[i * 3 + 2] = g.at<double>(i / width, i % width);
	}

	return 0;
}
