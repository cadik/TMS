#include "clahe.h"

cv::Mat histogramEqualization(cv::Mat matrix, int height, int width,  int gridRegions, double cl) {
	double histogram[256] = {0};
	double newhistogram[256] = {0};
	double cmphistogram[256] = {0};

	cv::Mat newImage;
	newImage = cv::Mat::zeros(height, width, CV_32F);	

	/*
		Computing regions
	*/
	
	int widthNumberRegions = width / gridRegions + 1;
	int heightNumberRegions = height / gridRegions + 1;
	std::vector<double> regionSubValues;
	std::vector<double*> regionNewHistograms;
	for (int j = 0; j < heightNumberRegions; j++) {	
		for (int i = 0; i < widthNumberRegions; i++) {
			/*
				Also with corners and edges
			*/
			int widthRegion = gridRegions;
			int heightRegion = gridRegions;

			if (i == widthNumberRegions - 1 && j == heightNumberRegions - 1) {
				heightRegion = height % gridRegions;
				widthRegion = width % gridRegions;
			} else if (i == widthNumberRegions - 1) {
				widthRegion = width % gridRegions;
			}  else if (j == heightNumberRegions - 1) {
				heightRegion = height % gridRegions;
			}
			double regionHistogram[256] = {0};
			double regionNewHistogram[256] = {0};
			double regionCmphistogram[256] = {0};
			double minValue = 1E90;
			double maxValue = 0;

			/*
				Getting max and min value
			*/
			for (int jj = 0; jj < widthRegion; jj++) {
				for (int ii = 0; ii < heightRegion; ii++) {
					if (matrix.at<float>((j * gridRegions) + jj, (i * gridRegions) + ii) < minValue) {
						minValue = matrix.at<float>((j * gridRegions) + jj, (i * gridRegions) + ii);
					}

					if (matrix.at<float>((j * gridRegions) + jj, (i * gridRegions) + ii) > maxValue) {
						maxValue = matrix.at<float>((j * gridRegions) + jj, (i * gridRegions) + ii);
					}
				}
			}			
			/*
				Computing with boundaries
			*/
			/*
				Getting histogram
			*/
			double subValue = (maxValue - minValue) / 256.0;
			regionSubValues.push_back(subValue);

			int counter = 0;
			for (int jjjj = 0; jjjj < heightRegion; jjjj++) {
				for (int iiii = 0; iiii < widthRegion; iiii++)  {
					if (matrix.at<float>((j * gridRegions) + jjjj, (i * gridRegions) + iiii) == maxValue) {
						regionHistogram[255]++;
					} else {	
						regionHistogram[(int)floor(matrix.at<float>((j * gridRegions) + jjjj, (i * gridRegions) + iiii) / subValue)]++;
					}					
				}
			}
			/*
				Clipping histogram
			*/

			double regionClippedHistogram[256] = {0};
			double averageValue = 0;
			if (cl > 0) {
				/*
					Getting average value of counts of histogram bins
				*/
				for (int i = 0; i < 256; i++) {
					averageValue += regionHistogram[i];
				}
				averageValue = (int)(averageValue/256.0);

				int maxBinValueCL = (floor)(averageValue*cl);

				/*
					Getting number of bins to layout uniformly
				*/
				int toLayoutBins = 0;
				for (int i = 0; i < 256; i++) {
					if (regionHistogram[i] > maxBinValueCL) {
						toLayoutBins+= regionHistogram[i] - maxBinValueCL;
						
					}
					regionClippedHistogram[i] += regionHistogram[i];
				}

				/*
					Uniformly coresponding bins
				*/
				for (int i = 0; i < 256; i++) {
					regionClippedHistogram[i] += toLayoutBins/256;
					if (regionClippedHistogram[i] > maxBinValueCL) {
						regionClippedHistogram[i] = regionClippedHistogram[i] - maxBinValueCL;
					}
				}
			} else {
				for (int i = 0; i < 256; i++) {
					regionClippedHistogram[i] = regionHistogram[i];
				}
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
				regionNewHistogram[a] = regionCmphistogram[a]*(maxValue - minValue) + minValue;
			}
			regionNewHistograms.push_back(regionNewHistogram);
			/*
				Save to image
			*/
			/*
				Save to image
			*/
			for (int b = 0; b < heightRegion; b++) {
				for (int a = 0; a < widthRegion; a++)  {
					if (matrix.at<float>((j * gridRegions) + b, (i * gridRegions) + a) == maxValue) {
						newImage.at<float>((j * gridRegions) + b, (i * gridRegions) + a) = regionNewHistogram[255];
					}else {
						newImage.at<float>((j * gridRegions) + b, (i * gridRegions) + a) = regionNewHistogram[(int)floor(matrix.at<float>((j * gridRegions) + b, (i * gridRegions) + a) / subValue)];
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
	int endJ = height - (height % gridRegions)/2; 

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
			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion]) > 255) {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][255];			
			} else {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion])];
			}

			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1]) > 255) {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][255];
			} else {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1])];
			}

			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*(actualJRegion + 1) + actualIRegion]) > 255) {
				tmp3 = regionNewHistograms[widthNumberRegions*(actualJRegion + 1) + actualIRegion][255];
			} else {
				tmp3 = regionNewHistograms[widthNumberRegions*(actualJRegion + 1) + actualIRegion][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*(actualJRegion + 1) + actualIRegion])];
			}
			
			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*(actualJRegion + 1) + actualIRegion + 1]) > 255) {
				tmp4 = regionNewHistograms[widthNumberRegions*(actualJRegion + 1) + actualIRegion + 1][255];
			} else {
				tmp4 = regionNewHistograms[widthNumberRegions*(actualJRegion + 1) + actualIRegion + 1][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*(actualJRegion + 1) + actualIRegion + 1])];
			}
			double dump1 = sizeYMinusOne * (sizeXMinusOne * tmp1 + sizeX * tmp2);
			double dump2 = sizeY * (sizeXMinusOne * tmp3 + sizeX * tmp4);
		
			newImage.at<float>(j, i) = dump1 + dump2;
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
			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion]) > 255) {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][255];			
			} else {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion])];
			}

			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1]) > 255) {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][255];
			} else {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1])];
			}
			newImage.at<float>(j, i) = sizeX * tmp1 + sizeXMinusOne * tmp2;
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
			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion]) > 255) {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][255];			
			} else {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion])];
			}

			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1]) > 255) {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][255];
			} else {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1])];
			}
			newImage.at<float>(j, i) = sizeX * tmp1 + sizeXMinusOne * tmp2;
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
			int actualJRegion = coordinatesYLess / (gridRegions);

			double tmp1 = 0;
			double tmp2 = 0;
			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion]) > 255) {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][255];			
			} else {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion])];
			}

			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1]) > 255) {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][255];
			} else {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1])];
			}
			newImage.at<float>(j, i) = sizeX * tmp1 + sizeXMinusOne * tmp2;
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
			int actualJRegion = coordinatesYLess / (gridRegions);

			double tmp1 = 0;
			double tmp2 = 0;
			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion]) > 255) {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][255];			
			} else {
				tmp1 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion])];
			}

			if ((int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1]) > 255) {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][255];
			} else {
				tmp2 = regionNewHistograms[widthNumberRegions*actualJRegion + actualIRegion + 1][(int)floor(matrix.at<float>(j, i) / regionSubValues[widthNumberRegions*actualJRegion + actualIRegion + 1])];
			}
			newImage.at<float>(j, i) = sizeX * tmp1 + sizeXMinusOne * tmp2;
		}
	}
    return newImage;
}
