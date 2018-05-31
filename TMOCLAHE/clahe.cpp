#include "clahe.h"

#include <iostream>
#include <fstream>

using namespace std;

cv::Mat histogramEqualization(cv::Mat matrix, int height, int width,  int gridRegions, double cl) {
	double histogram[256] = {0};
	double newhistogram[256] = {0};
	double cmphistogram[256] = {0};

	cv::Mat newImage;
	newImage = cv::Mat::zeros(height, width, CV_64F);	

	/*
		Getting scale, + subValue
	*/
	double maxValue1 = 0;
	double minValue1 = 1E90;

	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			if (matrix.at<float>(j, i) < minValue1) {
				minValue1 = matrix.at<float>(j, i);
			}

			if (matrix.at<float>(j, i) > maxValue1) {
				maxValue1 = matrix.at<float>(j, i);
			}
		}
	}
	// normalized histogram, sirka binu
	double subValue1 = (maxValue1 - minValue1) / 256 ;

	/*
		Computing regions
	*/
	
	int widthNumberRegions = width / gridRegions;

	if (width % gridRegions != 0) {
		widthNumberRegions += 1;
	}
	int heightNumberRegions = height / gridRegions;
	if (height % gridRegions != 0) {
		heightNumberRegions += 1;
	}
	int counter;
	cv::Mat regionHistograms = cv::Mat::zeros(256, widthNumberRegions*heightNumberRegions, CV_32F);
	for (int j = 0; j < heightNumberRegions; j++) {	
		for (int i = 0; i < widthNumberRegions; i++) {
			/*
				Also with corners and edges
			*/
			int widthRegion = gridRegions;
			int heightRegion = gridRegions;

			if (i == widthNumberRegions - 1 && j == heightNumberRegions - 1) {
				if (height % gridRegions != 0) {
					heightRegion = height % gridRegions;
				}

				if (width % gridRegions != 0) {
					widthRegion = width % gridRegions;
				}
			} else if (i == widthNumberRegions - 1) {
				if (width % gridRegions != 0) {
					widthRegion = width % gridRegions;
				}
			}  else if (j == heightNumberRegions - 1) {
				if (height % gridRegions != 0) {
					heightRegion = height % gridRegions;
				}
			}
			int regionHistogram[256] = {0};
			double regionNewHistogram[256] = {0};
			double regionCmphistogram[256] = {0};

			/*
				Getting histogram
			*/

			int counter = 0;
			for (int jjjj = 0; jjjj < heightRegion; jjjj++) {
				for (int iiii = 0; iiii < widthRegion; iiii++)  {
					if (matrix.at<float>((j * gridRegions) + jjjj, (i * gridRegions) + iiii) == maxValue1) {
						regionHistogram[255]++;
					} else {	
						regionHistogram[(int)round(matrix.at<float>((j * gridRegions) + jjjj, (i * gridRegions) + iiii) / subValue1)]++;
					}					
				}
			}
			/*
				Clipping histogram
			*/

			double regionClippedHistogram[256] = {0};
			/*for (int o = 0; o < 256; o++) {
				regionClippedHistogram[o] = regionHistogram[o];
			}*/
			//if (cl > 0) {
				/*
					Getting average value of counts of histogram bins
				*/
				int maxRegionValue = 0;
				for (int i = 0; i < 256; i++) {
					if (regionHistogram[i] > maxRegionValue) {
						maxRegionValue = regionHistogram[i];
					}
					regionClippedHistogram[i] = regionHistogram[i];
				}
				int maxBinValueCL = (round)(maxRegionValue*cl);

				/*
					Getting number of bins to layout uniformly
				*/
				int toLayoutBins = 0;
				for (int i0 = 0; i0 < 256; i0++) {
					if (regionHistogram[i0] > maxBinValueCL) {
						toLayoutBins += regionHistogram[i0] - maxBinValueCL;
						regionClippedHistogram[i0] -= regionHistogram[i0] - maxBinValueCL;					
					}					
				}
				/*
					Uniformly coresponding bins
				*/
				int toLayoutForBin = toLayoutBins / 256;
				int restBins = toLayoutBins - toLayoutForBin * 256;
				for (int i1 = 0; i1 < 256; i1++) {
					regionClippedHistogram[i1] += toLayoutForBin;
				}

				
				for (int i2 = 0; i2 < restBins; i2++) {
					regionClippedHistogram[i2]++;
				}

			/*
				Cumulate propability
			*/
			for (int a = 0; a < 256; a++) {
				if (a == 0) {
					regionCmphistogram[a] = regionClippedHistogram[a]/(heightRegion*widthRegion);
				} else {
					regionCmphistogram[a] = regionCmphistogram[a - 1] + regionClippedHistogram[a]/(heightRegion*widthRegion);
				}
			}

			/*
				Creating new histogram
			*/
			for (int a = 0; a < 256; a++) {
				regionHistograms.at<float>(a, j * widthNumberRegions + i) = regionCmphistogram[a]*(maxValue1 - minValue1) + minValue1;
				
			}
			/*
				Save to image
			*/
			for (int b = 0; b < heightRegion; b++) {
				for (int a = 0; a < widthRegion; a++)  {
					if (matrix.at<float>((j * gridRegions) + b, (i * gridRegions) + a) == maxValue1) {
						newImage.at<float>((j * gridRegions) + b, (i * gridRegions) + a) = regionHistograms.at<float>(255, j * gridRegions + i);
					}else {
						newImage.at<float>((j * gridRegions) + b, (i * gridRegions) + a) = regionHistograms.at<float>((int)round(matrix.at<float>((j * gridRegions) + b, (i * gridRegions) + a) / subValue1), j * gridRegions + i);
					}
				}
			}		
		}
	}

	/*
		Bilinear interpolation
	*/
	int startI = gridRegions/2;
	int startJ = gridRegions/2;
	int endI = width - (width % gridRegions)/2;
	if (width % gridRegions == 0) {
		endI = width - gridRegions/2;
	}
	
	int endJ = height - (height % gridRegions)/2; 
	if (height % gridRegions == 0) {
		endJ = height - gridRegions/2;
	}
	for (int j = startJ + 1; j < endJ; j++) {
		for (int i = startI + 1; i < endI; i++) {
			double coordinatesXLess = ((i - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesXGreater = coordinatesXLess + gridRegions;
			if (coordinatesXGreater >= width) {
				coordinatesXGreater = coordinatesXLess + (width - coordinatesXLess) / 2;
			}
			double coordinatesYLess = ((j - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesYGreater = coordinatesYLess + gridRegions;
			if (coordinatesYGreater >= height) {
				coordinatesYGreater = coordinatesYLess + (height - coordinatesYLess) / 2;
			}
		

			int actualIRegion = coordinatesXLess / (gridRegions);
			int actualJRegion = coordinatesYLess / (gridRegions);
			double sizeX = i - coordinatesXLess;
			double sizeXMinusOne = (coordinatesXGreater - coordinatesXLess) - sizeX;
			double sizeY = j - coordinatesYLess;
			double sizeYMinusOne = (coordinatesYGreater - coordinatesYLess) - sizeY;
			sizeX = sizeX / (coordinatesXGreater - coordinatesXLess);
			sizeXMinusOne = sizeXMinusOne / (coordinatesXGreater - coordinatesXLess);
			sizeY = sizeY / (coordinatesYGreater - coordinatesYLess);
			sizeYMinusOne = sizeYMinusOne / (coordinatesYGreater - coordinatesYLess);

			double tmp1 = 0;
			double tmp2 = 0;
			double tmp3 = 0;
			double tmp4 = 0;
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp1 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion);
			} else {
				tmp1 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion);
			}

			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp2 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion + 1);
			} else {
				tmp2 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion + 1);
			}

			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp3 = regionHistograms.at<float>(255, widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			} else {
				tmp3 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			}
			
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp4 = regionHistograms.at<float>(255, widthNumberRegions*(actualJRegion + 1) + actualIRegion + 1);
			} else {
				tmp4 =regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*(actualJRegion + 1) + actualIRegion + 1);
			}
			double dump1 = sizeYMinusOne * (sizeXMinusOne * tmp1 + sizeX * tmp2);
			double dump2 = sizeY * (sizeXMinusOne * tmp3 + sizeX * tmp4);
		
			newImage.at<float>(j, i) = abs(dump1 + dump2);
		}
	}

	/*
		Linear interpolation top-edge
	*/

	for (int j = 0; j <= startJ; j++) {
		for (int i = startI; i <= endI; i++) {
			double coordinatesXLess = ((i - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesXGreater = coordinatesXLess + gridRegions;
			if (coordinatesXGreater >= width) {
				coordinatesXGreater = coordinatesXLess + (width - coordinatesXLess) / 2;
			}

			double coordinatesYLess = ((j - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesYGreater = coordinatesYLess + gridRegions;
			if (coordinatesYGreater >= height) {
				coordinatesYGreater = coordinatesYLess + (height - coordinatesYLess) / 2;
			}
			double sizeX = i - coordinatesXLess;
			double sizeXMinusOne = (coordinatesXGreater - coordinatesXLess) - sizeX;
			sizeX = sizeX / (coordinatesXGreater - coordinatesXLess);
			sizeXMinusOne = sizeXMinusOne / (coordinatesXGreater - coordinatesXLess);
			int actualIRegion = coordinatesXLess / (gridRegions);
			double sizeY = j - coordinatesYLess;
			double sizeYMinusOne = (coordinatesYGreater - coordinatesYLess) - sizeY;
			int actualJRegion = coordinatesYLess / (gridRegions);

			double tmp1 = 0;
			double tmp2 = 0;
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp1 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion);
			} else {
				tmp1 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion);
			}

			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp2 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion + 1);
			} else {
				tmp2 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion + 1);
			}
			newImage.at<float>(j, i) = abs(sizeXMinusOne * tmp1 + sizeX * tmp2);
		}
	}

	/*
		Linear interpolation bottom-edge
	*/
	for (int j = endJ; j < height; j++) {
		for (int i = startI; i <= endI; i++) {
			double coordinatesXLess = ((i - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesXGreater = coordinatesXLess + gridRegions;
			if (coordinatesXGreater >= width) {
				coordinatesXGreater = coordinatesXLess + (width - coordinatesXLess) / 2;
			}

			double coordinatesYLess = ((j - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesYGreater = coordinatesYLess + gridRegions;
			if (coordinatesYGreater >= height) {
				coordinatesYGreater = coordinatesYLess + (height - coordinatesYLess) / 2;
			}
			double sizeX = i - coordinatesXLess;
			double sizeXMinusOne = (coordinatesXGreater - coordinatesXLess) - sizeX;
			sizeX = sizeX / (coordinatesXGreater - coordinatesXLess);
			sizeXMinusOne = sizeXMinusOne / (coordinatesXGreater - coordinatesXLess);
			int actualIRegion = coordinatesXLess / (gridRegions);
			double sizeY = j - coordinatesYLess;
			double sizeYMinusOne = (coordinatesYGreater - coordinatesYLess) - sizeY;
			int actualJRegion = coordinatesYLess / (gridRegions);

			double tmp1 = 0;
			double tmp2 = 0;
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp1 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion);
			} else {
				tmp1 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion+ actualIRegion);
			}

			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp2 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion + 1);
			} else {
				tmp2 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion + 1);
			}
			newImage.at<float>(j, i) = abs(sizeXMinusOne * tmp1 + sizeX * tmp2);
		}
	}
	/*
		Linear interpolation left-edge
	*/
	for (int j = startJ; j <= endJ; j++) {
		for (int i = 0; i <= startI; i++) {
			double coordinatesXLess = ((i - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesXGreater = coordinatesXLess + gridRegions;
			if (coordinatesXGreater >= width) {
				coordinatesXGreater = coordinatesXLess + (width - coordinatesXLess) / 2;
			}

			double coordinatesYLess = ((j - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesYGreater = coordinatesYLess + gridRegions;
			if (coordinatesYGreater >= height) {
				coordinatesYGreater = coordinatesYLess + (height - coordinatesYLess) / 2;
			}
			double sizeX = i - coordinatesXLess;
			double sizeXMinusOne = (coordinatesXGreater - coordinatesXLess) - sizeX;
			sizeX = sizeX / (coordinatesXGreater - coordinatesXLess);
			sizeXMinusOne = sizeXMinusOne / (coordinatesXGreater - coordinatesXLess);
			int actualIRegion = coordinatesXLess / (gridRegions);
			double sizeY = j - coordinatesYLess;
			double sizeYMinusOne = (coordinatesYGreater - coordinatesYLess) - sizeY;
			sizeY = sizeY / (coordinatesYGreater - coordinatesYLess);
			sizeYMinusOne = sizeYMinusOne / (coordinatesYGreater - coordinatesYLess);
			int actualJRegion = coordinatesYLess / (gridRegions);

			double tmp1 = 0;
			double tmp2 = 0;
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp1 = regionHistograms.at<float>(255, widthNumberRegions*actualJRegion + actualIRegion);
			} else {
				tmp1 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion);
			}
			
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp2 = regionHistograms.at<float>(255, widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			} else {
				tmp2 =regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			}
			newImage.at<float>(j, i) = abs(sizeYMinusOne * tmp1 + sizeY * tmp2);
		}
	}

	/*
		Linear interpolation right-edge
	*/
	for (int j = startJ; j <= endJ; j++) {
		for (int i = endI; i < width; i++) {
			double coordinatesXLess = ((i - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesXGreater = coordinatesXLess + gridRegions;
			if (coordinatesXGreater >= width) {
				coordinatesXGreater = coordinatesXLess + (width - coordinatesXLess) / 2;
			}

			double coordinatesYLess = ((j - gridRegions/2)/gridRegions) * gridRegions + gridRegions/2;
			double coordinatesYGreater = coordinatesYLess + gridRegions;
			if (coordinatesYGreater >= height) {
				coordinatesYGreater = coordinatesYLess + (height - coordinatesYLess) / 2;
			}
			double sizeX = i - coordinatesXLess;
			double sizeXMinusOne = (coordinatesXGreater - coordinatesXLess) - sizeX;
			sizeX = sizeX / (coordinatesXGreater - coordinatesXLess);
			sizeXMinusOne = sizeXMinusOne / (coordinatesXGreater - coordinatesXLess);
			int actualIRegion = coordinatesXLess / (gridRegions);
			double sizeY = j - coordinatesYLess;
			double sizeYMinusOne = (coordinatesYGreater - coordinatesYLess) - sizeY;
			sizeY = sizeY / (coordinatesYGreater - coordinatesYLess);
			sizeYMinusOne = sizeYMinusOne / (coordinatesYGreater - coordinatesYLess);
			int actualJRegion = coordinatesYLess / (gridRegions);
			double tmp1 = 0;
			double tmp2 = 0;
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp1 = regionHistograms.at<float>(255, widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			} else {
				tmp1 = regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*actualJRegion + actualIRegion);
			}
			
			if ((int)round(matrix.at<float>(j, i) / subValue1) > 255) {
				tmp2 = regionHistograms.at<float>(255, widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			} else {
				tmp2 =regionHistograms.at<float>((int)round(matrix.at<float>(j, i) / subValue1), widthNumberRegions*(actualJRegion + 1) + actualIRegion);
			}
			newImage.at<float>(j, i) = abs(sizeYMinusOne * tmp1 + sizeY * tmp2);
		}
	}

	ofstream myfile;
	myfile.open ("example.txt");
	myfile << newImage;
	myfile.close();
    return newImage;
}

