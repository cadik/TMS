/* --------------------------------------------------------------------------- *
 * TMOEisemann04.cpp: implementation of the TMOEisemann04 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOEisemann04.h"



/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOEisemann04::TMOEisemann04()
{
	SetName(L"Eisemann04");
	SetDescription(L"Flash/No-flash");

	flashImagePathParameter.SetName(L"Flash");
	flashImagePathParameter.SetDescription(L"Input flash image path");
	flashImagePathParameter.SetDefault("");
	this->Register(flashImagePathParameter);
}

TMOEisemann04::~TMOEisemann04()
{
}

Mat ToMat64C3(TMOImage* image){
	return Mat(image->GetHeight(), image->GetWidth(), CV_64FC3, image->GetData());
}

Mat computeLargeScale(Mat intensity){
	double min, max;
	minMaxLoc(intensity, &min, &max);
	Mat output;
	Mat intensity32;
	intensity.convertTo(intensity32, CV_8U, 255., 0);
	bilateralFilter(intensity32, output, 9, 10,10);
	output.convertTo(output, CV_64F, 1./255., 0);
	return output;
}

void saveImage(string name, Mat image){
	Mat output;
	image.convertTo(output, CV_8UC3, 255.,0);
	cvtColor(output, output, COLOR_RGB2BGR);
	imwrite(name, output);
}

Mat computeDetail(Mat intensity, Mat largeScale){
	Mat detail = Mat();
	divide(intensity, largeScale, detail);
	return detail;
}

Mat computeColor(Mat image, Mat intensity){
	Mat color;
	Mat channels[3];
	vector<Mat> channelsV;
	split(image, channels);
	for(int i = 0; i < 3;i++){
		divide(channels[i], intensity, channels[i]);
		channelsV.push_back(channels[i]);
	}
	merge(channelsV, color);
	return color;
}

Mat computeIntensity(Mat img, Mat* img2 = NULL){
	Mat channels[3]; split(img, channels);
	Mat sums = channels[0];
	add(sums, channels[1], sums);
	add(sums, channels[1], sums);
	Mat sums3; merge(vector<Mat>{sums,sums,sums}, sums3);
	
	Mat temp = img2 == NULL ? img.clone() : img2->clone();

	divide(temp, sums3, temp);
	multiply(temp, img, temp);

	split(temp, channels);
	Mat intensity = channels[0];
	add(intensity, channels[1], intensity);
	add(intensity, channels[2], intensity);

	normalize(intensity, intensity, 0., 1., NORM_MINMAX, CV_64F);
	return intensity;
}

double sumMaxValues(Mat channel){
	double min, max, sum;
	minMaxLoc(channel, &min, &max);

	for(int j = 0; j < channel.size().width; j++)
		for(int i = 0; i < channel.size().height; i++)
			if(channel.at<double>(i,j) == max)
				sum += max;
	return sum;
}

double getWhiteBallanceChannelWeight(Mat channelF, Mat channelNF){
	return pow(0.4 * sumMaxValues(channelF) + 0.6 * sumMaxValues(channelNF), 0.2);
}

Mat whiteBalanceCorr(Mat image, Mat imageF, Mat imageNF){
	Mat channels [3]; split(image, channels);
	Mat channelsF [3]; split(imageF, channelsF);
	Mat channelsNF [3]; split(imageNF, channelsNF);
	for(int i = 0; i<3; i++){
		double weight = getWhiteBallanceChannelWeight(channelsF[i], channelsNF[i]);
		multiply(channels[i], weight, channels[i]);
	}
	Mat corrected = Mat();
	merge(vector<Mat>{channels[0], channels[1], channels[2]}, corrected);
	return corrected;
}

Mat toColor(Mat gray){
	Mat color;
	merge(vector<Mat>{gray, gray, gray}, color);
	return color;
}

int TMOEisemann04::Transform()
{
	if(flashImagePathParameter.GetString() == ""){
		cerr << "No flash image provided" << endl;
		return 1;
	}

	flashImage = TMOImage();
	flashImage.Open(flashImagePathParameter.GetString().c_str());

	Mat imageF = ToMat64C3(&flashImage);
	Mat imageNF = ToMat64C3(pSrc);

	saveImage("imageF.jpg", imageF);
	saveImage("imageNF.jpg", imageNF);

	cerr << "LOADED images\n";


	int width = pSrc->GetWidth();
	int height = pSrc->GetHeight();

	Mat intensityF =  computeIntensity(imageF);
	Mat intensityNF =  computeIntensity(imageNF, &imageF);
	Mat intensityNFLarge = computeIntensity(imageNF);

	saveImage("intensityF.jpg", toColor(intensityF));
	saveImage("intensityNF.jpg", toColor(intensityNF));
	saveImage("intensityNFLarge.jpg", toColor(intensityNFLarge));
	cerr << "COMPUTED Intensities\n";

	Mat largeScaleF = computeLargeScale(intensityF);
	Mat largeScaleNF = computeLargeScale(intensityNFLarge);

	saveImage("largeScaleF.jpg", toColor(largeScaleF));
	saveImage("largeScaleNF.jpg", toColor(largeScaleNF));
	cerr << "COMPUTED Large scale\n";

	Mat colorF = computeColor(imageF, intensityF);
	Mat colorNF = computeColor(imageNF, intensityNF);

	saveImage("colorF.jpg", colorF);
	saveImage("colorNF.jpg", colorNF);	
	cerr << "COMPUTED Color\n";

	Mat detailF = computeDetail(intensityF, largeScaleF);
	Mat detailNF = computeDetail(intensityNFLarge, largeScaleNF);

	saveImage("detailF.jpg", toColor(detailF));
	saveImage("detailNF.jpg", toColor(detailNF));
	cerr << "COMPUTED Detail\n";

	Mat resultGray, result;
	multiply(detailF, largeScaleNF, resultGray);
	result = toColor(resultGray);
	multiply(result, colorF, result);
	result = whiteBalanceCorr(result, imageF, imageNF);
	saveImage("result.jpg", result);


	double* pDestinationData = pDst->GetData();
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			Vec3d resultPixel = result.at<Vec3d>(j, i);
			for(int c = 0; c < 3; c++){
				*pDestinationData++ = resultPixel[c];
			}
		}
	}
	pDst->Convert(TMO_RGB);
	return 0;
}
