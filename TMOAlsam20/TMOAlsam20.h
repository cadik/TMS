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
#include "TMO.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace cv;
using namespace std;


class TMOAlsam20 : public TMO
{
public:
	TMOAlsam20();
	virtual ~TMOAlsam20();
	virtual int Transform();
	double calculate_norm(Mat diff, int ord);
	double red_green_contrast(Mat I1, Mat I2, Mat Dc, Mat Ds, int ord);
	double green_red_contrast(Mat I1, Mat I2, Mat Dc, Mat Ds, int ord);
	double blue_yellow_contrast(Mat I1, Mat I2, Mat I3, Mat Dc, Mat Ds, float w1, float w2 , int ord);
	vector<double> calculate_weights(Mat I1, Mat I2, Mat I3, Mat Dc, Mat Ds, int ord, float w1, float w2);
	Mat grayscale_image(Mat I0, Mat I1, Mat I2, vector<double> alpha);
	Mat iterative_greyscale_image(Mat I0, Mat I1, Mat I2, vector<Mat> Dc_list, vector<Mat> Ds_list,  int ord, double w1, double w2);
	pair<vector<Mat>, vector<Mat>> load_kernels(string file_kernel);

protected:
	TMOInt pord;
	TMODouble pw1;
	TMODouble pw2;
	TMOString pfile_kernel;
};
