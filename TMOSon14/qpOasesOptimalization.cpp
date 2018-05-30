

#include <fstream>

#include "qpOasesOptimization.h"

/*
	Ponechány zakomentované i hodnoty vzniklé po derivaci
*/
std::vector<cv::Mat> optimizationWithOases(int height, int width, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels)
{
	
	USING_NAMESPACE_QPOASES
	/*
		Deklarace proměnných
	*/
	int w1 = width;
	int h1 = height;	
	std::vector<cv::Mat> sAndT;
    cv::Mat s;
	s = cv::Mat::ones(height, width, CV_32F);
	cv::Mat t;
    t = cv::Mat::zeros(height, width, CV_32F);

	/*
		Optimization
	*/
	
	int mainSize = h1*w1*2; // počet parametrů
  /*
		Creating D matrix
	*/
	cv::Mat DMatrix = cv::Mat::zeros(mainSize, mainSize, CV_32F);


	/*
		Setting s (podle té mapy co mám, aby jsem věděl co s  čím komunikuje v případě gradientů, jedu po prvku)
		-- Prvky jsou už pronásobeny dvěma, ono je to jedno, jestli se znásobí potom celá matice nebo samotné prvky
	*/
 for (int j = 0; j < h1; j++) {
		for (int i = 0; i < w1; i++) {
			
			double sPowerTerm1 = -((double)(pow(detailImage.at<float>(j, i)/255.0, 2))); //zkontrolovat
			//double sPowerTerm2 = 200*4*weight1.at<float>(j, i);
			double sPowerTerm2 = 0.0;
			double sPowerTerm3 = 0.0;
			double sPowerTerm4 = 0.0;
			double sPowerTerm5 = 0.0;
			double sPowerTerm6 = 0.0;
			double sPowerTerm7 = 0.0;
			double sPowerTerm8 = 0.0;
			double sPowerTerm9 = 0.0;
			double sPowerTerm10 = 0.0;
			double sPowerTerm11 = 0.0;
			double sPowerTerm12 = 0.0;
			double sPowerTerm13 = 0.0;
			double sumOfPowerTerms = 0.0;


			/* if ((i > 0 && i < width - 1) && (j > 0 && j < height - 1)) {
				sPowerTerm3 = -200*weight1.at<float>(j, i + 1);
				sPowerTerm4 = -200*weight1.at<float>(j, i - 1);
				sPowerTerm5 = -200*weight1.at<float>(j + 1, i);
				sPowerTerm6 = -200*weight1.at<float>(j - 1, i);*/
				/*sPowerTerm8 = -200*weight1.at<float>(j - 1, i + 1);
				sPowerTerm9 = -200*weight1.at<float>(j - 1, i - 1);
				sPowerTerm10 = -200*weight1.at<float>(j + 1, i + 1);
				sPowerTerm11 = -200*weight1.at<float>(j - 1, i - 1);*/
			//}

			/*
				Corners
			*/

			/*if (i == 0 && j == 0) {
				sPowerTerm3 = -200*2*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j, i + 1);
				sPowerTerm5 = -200*weight1.at<float>(j + 1, i);
			}

			if (i == 0 && j == height - 1) {
				sPowerTerm3 = -200*2*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j, i + 1);
				sPowerTerm5 = -200*weight1.at<float>(j - 1, i);
			}

			if (i == width - 1 && j == 0) {
				sPowerTerm3 = -200*2*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j, i - 1);
				sPowerTerm5 = -200*weight1.at<float>(j + 1, i);
			}

			if (i == width - 1 && j == height - 1) {
				sPowerTerm3 = -200*2*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j, i - 1);
				sPowerTerm5 = -200*weight1.at<float>(j - 1, i);
			}*/

			/*
				Boundaries
			*/

			/*if (i == 0 && (j > 0 && j < height - 1)) {
				sPowerTerm3 = -200*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j + 1, i);
				sPowerTerm5 = -200*weight1.at<float>(j - 1, i);
				sPowerTerm6 = -200*weight1.at<float>(j, i + 1);
			}

			if (i == width - 1 && (j > 0 && j < height - 1))  {
				sPowerTerm3 = -200*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j + 1, i);
				sPowerTerm5 = -200*weight1.at<float>(j - 1, i);
				sPowerTerm6 = -200*weight1.at<float>(j, i - 1);
			}

			if ((i > 0 && i < width - 1) && j == 0) {
				sPowerTerm3 = -200*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j, i + 1);
				sPowerTerm5 = -200*weight1.at<float>(j, i - 1);
				sPowerTerm6 = -200*weight1.at<float>(j + 1, i);
			}

			if ((i > 0 && i < width - 1) && j == height - 1) {
				sPowerTerm3 = -200*weight1.at<float>(j, i);
				sPowerTerm4 = -200*weight1.at<float>(j, i + 1);
				sPowerTerm5 = -200*weight1.at<float>(j, i - 1);
				sPowerTerm6 = -200*weight1.at<float>(j - 1, i);
			}*/
			

			/*
				Puvodni
			*/
			if (i > 1) {
				sPowerTerm2 = 2*200*(double)weight1.at<float>(j, i - 1);
				DMatrix.at<float>(w1 * j + i, w1 * j + (i - 2)) = -sPowerTerm2;
				DMatrix.at<float>(w1 * j + (i - 2), w1 * j + i) = -sPowerTerm2;
			}
			if (j < (h1 - 2)) {
				sPowerTerm3 = 2*200*(double)weight1.at<float>(j + 1, i);
				DMatrix.at<float>(w1 * (j + 2) + i, w1 * j + i) = -sPowerTerm3;
				DMatrix.at<float>(w1 * j + i, w1 * (j + 2) + i) = -sPowerTerm3;
			}

			if (i < (w1 - 2)) {
				sPowerTerm4 = 2*200*(double)weight1.at<float>(j, i + 1);
				DMatrix.at<float>(w1 * j + (i + 2), w1 * j + i) = -sPowerTerm4;
				DMatrix.at<float>(w1 * j + i, w1 * j + (i + 2)) = -sPowerTerm4;
			}

			if (j > 1) {
				sPowerTerm5 = 2*200*(double)weight1.at<float>(j - 1, i);
				DMatrix.at<float>(w1 * j + i, w1 * (j  - 2) + i) = -sPowerTerm5;
				DMatrix.at<float>(w1 * (j  - 2) + i, w1 * j + i) = -sPowerTerm5;
			}
			
			if (i == 0) {
				sPowerTerm6 = 2*200*(double)weight1.at<float>(j, i);
				DMatrix.at<float>(w1 * j + i, w1 * j + i + 1) = -sPowerTerm6;
				DMatrix.at<float>(w1 * j + i + 1, w1 * j + i) = -sPowerTerm6;
			}

			if (i == 1) {
				sPowerTerm7 = 2*200*(double)weight1.at<float>(j, i - 1);
				DMatrix.at<float>(w1 * j + i, w1 * j + i - 1) = -sPowerTerm7;
				DMatrix.at<float>(w1 * j + i - 1, w1 * j + i) = -sPowerTerm7;
			}

			if (j == 0) {
				sPowerTerm8 = 2*200*(double)weight1.at<float>(j, i);
				DMatrix.at<float>(w1 * j + i, w1 * (j + 1) + i) = -sPowerTerm8;
				DMatrix.at<float>(w1 * (j + 1) + i, w1 * j + i) = -sPowerTerm8;
			}

			if (j == 1) {
				sPowerTerm9 = 2*200*(double)weight1.at<float>(j - 1, i);
				DMatrix.at<float>(w1 * j + i, w1 * (j - 1) + i) = -sPowerTerm9;
				DMatrix.at<float>(w1 * (j - 1) + i, w1 * j + i) = -sPowerTerm9;
			}

			if (i == w1 - 1) {
				sPowerTerm10 = 2*200*(double)weight1.at<float>(j, i);
				DMatrix.at<float>(w1 * j + i, w1 * j + i - 1) = -sPowerTerm10;
				DMatrix.at<float>(w1 * j + i - 1, w1 * j + i) = -sPowerTerm10;
			}

			if (i == w1 - 2) {
				sPowerTerm11 = 2*200*(double)weight1.at<float>(j, i + 1);
				DMatrix.at<float>(w1 * j + i, w1 * j + i + 1) = -sPowerTerm11;
				DMatrix.at<float>(w1 * j + i + 1, w1 * j + i) = -sPowerTerm11;
			}

			if (j == h1 - 1) {
				sPowerTerm12 = 2*200*(double)weight1.at<float>(j, i);
				DMatrix.at<float>(w1 * j + i, w1 * (j - 1) + i) = -sPowerTerm12;
				DMatrix.at<float>(w1 * (j - 1) + i, w1 * j + i) = -sPowerTerm12;
			}

			if (j == h1 - 2) {
				sPowerTerm13 = 2*200*(double)weight1.at<float>(j + 1, i);
				DMatrix.at<float>(w1 * j + i, w1 * (j + 1) + i) = -sPowerTerm13;
				DMatrix.at<float>(w1 * (j + 1) + i, w1 * j + i) = -sPowerTerm13;
			}

			sumOfPowerTerms = sPowerTerm1 +
							sPowerTerm2 +
							sPowerTerm3 +
							sPowerTerm4 +
							sPowerTerm5 +
							sPowerTerm6 + 
							sPowerTerm7 + 
							sPowerTerm8 + 
							sPowerTerm9 +
							sPowerTerm10 +
							sPowerTerm11 +
							sPowerTerm12 +
							sPowerTerm13; 
			
			DMatrix.at<float>(j * w1 + i, j * w1 + i) = sumOfPowerTerms;
		}
	}
	

	/*
		Setting t
	*/

	for (int j = 0; j < h1; j++) {
		for (int i = 0; i < w1; i++) {
			//double sPowerTerm01 = 4*500*(double)weight2.at<float>(j, i);
			double sPowerTerm01 = 0.0;
			double sPowerTerm02 = 0.0;
			double sPowerTerm03 = 0.0;
			double sPowerTerm04 = 0.0;
			double sPowerTerm05 = 0.0;
			double sPowerTerm06 = 0.0;
			double sPowerTerm07 = 0.0;
			double sPowerTerm08 = 0.0;
			double sPowerTerm09 = 0.0;
			double sPowerTerm010 = 0.0;
			double sPowerTerm011 = 0.0;
			double sPowerTerm012 = 0.0;
			double sumOfPowerTerms1 = 0.0;

			/*if ((i > 0 && i < width - 1) && (j > 0 && j < height - 1)) {
				sPowerTerm02 = -500*weight2.at<float>(j, i + 1);
				sPowerTerm03 = -500*weight2.at<float>(j, i - 1);
				sPowerTerm04 = -500*weight2.at<float>(j + 1, i);
				sPowerTerm05 = -500*weight2.at<float>(j - 1, i);*/

				/*sPowerTerm06 = -500*weight2.at<float>(j + 1, i + 1);
				sPowerTerm07 = -500*weight2.at<float>(j + 1, i - 1);
				sPowerTerm08 = -500*weight2.at<float>(j - 1, i + 1);
				sPowerTerm09 = -500*weight2.at<float>(j - 1, i - 1);*/
			//}

			/*
				Corners
			*/

			/*if (i == 0 && j == 0) { 
				sPowerTerm02 = -500*2*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j, i + 1);
				sPowerTerm04 = -500*weight2.at<float>(j + 1, i);
			}

			if (i == 0 && j == height - 1) { 
				sPowerTerm02 = -500*2*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j, i + 1);
				sPowerTerm04 = -500*weight2.at<float>(j - 1, i);
			}

			if (i == width - 1 && j == 0) { 
				sPowerTerm02 = -500*2*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j, i - 1);
				sPowerTerm04 = -500*weight2.at<float>(j + 1, i);
			}

			if (i == width - 1 && j == height - 1) {
				sPowerTerm02 = -500*2*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j, i - 1);
				sPowerTerm04 = -500*weight2.at<float>(j - 1, i);
			}*/
			
			/*
				Boundaries
			*/

			/*if (i == 0 && (j > 0 && j < height - 1)) {
				sPowerTerm02 = -500*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j + 1, i);
				sPowerTerm04 = -500*weight2.at<float>(j - 1, i);
				sPowerTerm05 = -500*weight2.at<float>(j, i + 1);
			}

			if (i == width - 1 && (j > 0 && j < height - 1))  {
				sPowerTerm02 = -500*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j + 1, i);
				sPowerTerm04 = -500*weight2.at<float>(j - 1, i);
				sPowerTerm05 = -500*weight2.at<float>(j, i - 1);
			}

			if ((i > 0 && i < width - 1) && j == 0) {
				sPowerTerm02 = -500*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j, i + 1);
				sPowerTerm04 = -500*weight2.at<float>(j, i - 1);
				sPowerTerm05 = -500*weight2.at<float>(j + 1, i);
			}

			if ((i > 0 && i < width - 1) && j == height - 1) {
				sPowerTerm02 = -500*weight2.at<float>(j, i);
				sPowerTerm03 = -500*weight2.at<float>(j, i + 1);
				sPowerTerm04 = -500*weight2.at<float>(j, i - 1);
				sPowerTerm05 = -500*weight2.at<float>(j - 1, i);
			}*/

			/*
				Puvodni
			*/
			if (i > 1) {
				sPowerTerm01 = 2*500*(double)weight2.at<float>(j, i - 1); //check
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + h1 * j + (i - 2)) = -sPowerTerm01;
				DMatrix.at<float>(w1*h1 + w1 * j + (i - 2), w1*h1 + w1 * j + i) = -sPowerTerm01;
			}
			if (j < (h1 - 2)) {
				sPowerTerm02 = 2*500*(double)weight2.at<float>(j + 1, i);
				DMatrix.at<float>(w1*h1 + w1 * (j + 2) + i, w1*h1 + w1 * j + i) = -sPowerTerm02;
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * (j + 2) + i) = -sPowerTerm02;
			}

			if (i < (w1 - 2)) {
				sPowerTerm03 = 2*500*(double)weight2.at<float>(j, i + 1);
				DMatrix.at<float>(w1*h1 + w1 * j + (i + 2), w1*h1 + w1 * j + i) = -sPowerTerm03;
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * j + (i + 2)) = -sPowerTerm03;
			}

			if (j > 1) {
				sPowerTerm04 = 2*500*(double)weight2.at<float>(j - 1, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * (j  - 2) + i) = -sPowerTerm04;
				DMatrix.at<float>(w1*h1 + w1 * (j  - 2) + i, w1*h1 + w1 * j + i) = -sPowerTerm04;
			}

			if (i == 0) {
				sPowerTerm05 = 2*500*(double)weight2.at<float>(j, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * j + i + 1) = -sPowerTerm05;
				DMatrix.at<float>(w1*h1 + w1 * j + i + 1, w1*h1 + w1 * j + i) = -sPowerTerm05;
			}

			if (i == 1) {
				sPowerTerm06 = 2*500*(double)weight2.at<float>(j, i - 1);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * j + i - 1) = -sPowerTerm06;
				DMatrix.at<float>(w1*h1 + w1 * j + i - 1, w1*h1 + w1 * j + i) = -sPowerTerm06;
			}
		
			if (j == 0) {
				sPowerTerm07 = 2*500*(double)weight2.at<float>(j, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * (j + 1) + i) = -sPowerTerm07;
				DMatrix.at<float>(w1*h1 + w1 * (j + 1) + i, w1*h1 + w1 * j + i) = -sPowerTerm07;
			}

			if (j == 1) {
				sPowerTerm08 = 2*500*(double)weight2.at<float>(j - 1, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * (j - 1) + i) = -sPowerTerm08;
				DMatrix.at<float>(w1*h1 + w1 * (j - 1) + i, w1*h1 + w1 * j + i) = -sPowerTerm08;
			}

			if (i == w1 - 1) {
				sPowerTerm09 = 2*500*(double)weight2.at<float>(j, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * j + i - 1) = -sPowerTerm09;
				DMatrix.at<float>(w1*h1 + w1 * j + i - 1, w1*h1 + w1 * j + i) = -sPowerTerm09;
			}

			if (i == w1 - 2) {
				sPowerTerm010 = 2*500*(double)weight2.at<float>(j, i + 1);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * j + i + 1) = -sPowerTerm010;
				DMatrix.at<float>(w1*h1 + w1 * j + i + 1, w1*h1 + w1 * j + i) = -sPowerTerm010;
			}

			if (j == h1 - 1) {
				sPowerTerm011 = 2*500*(double)weight2.at<float>(j, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * (j - 1) + i) = -sPowerTerm011;
				DMatrix.at<float>(w1*h1 + w1 * (j - 1) + i, w1*h1 + w1 * j + i) = -sPowerTerm011;
			}

			if (j == h1 - 2) {
				sPowerTerm012 = 2*500*(double)weight2.at<float>(j + 1, i);
				DMatrix.at<float>(w1*h1 + w1 * j + i, w1*h1 + w1 * (j + 1) + i) = -sPowerTerm012;
				DMatrix.at<float>(w1*h1 + w1 * (j + 1) + i, w1*h1 + w1 * j + i) = -sPowerTerm012;
			}

			sumOfPowerTerms1 = sPowerTerm01 +
							sPowerTerm02 +
							sPowerTerm03 +
							sPowerTerm04 +
							sPowerTerm05 +
							sPowerTerm06 +
							sPowerTerm07 +
							sPowerTerm08 +
							sPowerTerm09 +
							sPowerTerm010 +
							sPowerTerm011 +
							sPowerTerm012;
			DMatrix.at<float>(w1*h1 + j * w1 + i, w1*h1 + j * w1 + i) = sumOfPowerTerms1;
							
		}
	}
	
	real_t *H = new real_t[mainSize*mainSize];
	real_t *A = new real_t[h1*w1*mainSize*3]; // potřebujeme pro každý pixel jednu podmínku
	real_t *g = new real_t[mainSize]; // to je ponecháno nule, nemá žádný vliv
	real_t *lb = new real_t[mainSize]; // dolní hranice (počet hodnot jako počet pixelů)
	real_t *ub = new real_t[mainSize]; // horní hranice (počet hodnot jako počet pixelů)
	real_t *lbA = new real_t[h1*w1*3]; // dolní hranice A (počet pixelů, pro každý pixel jedna dolní hranice)
	real_t *ubA = new real_t[h1*w1*3]; // horní hranice A (počet pixelů, pro každý pixel jedna horní hranice)

	// DMatrix = 2*DMatrix; // to už je v prvcích uděláno
	/*
		Vkládám do Hessovy matice
	*/
	DMatrix = DMatrix*2;
	for (int j = 0; j < mainSize; j++) {
		for (int i = 0; i < mainSize; i++) {
			H[j * mainSize + i] =  DMatrix.at<float>(j, i);
		}
	}
	DMatrix.release();
	/*
		Nastavuji dolní + horní hranice + nulové hodnoty g
	*/
	for (int a = 0; a < mainSize; a++) {
		g[a] = 0;		

		if (a < mainSize/2) {
			lb[a] = -5000;
			ub[a] = 5000;
		} else {
			lb[a] = -2000;
			ub[a] = 2000;
		}
		
	}

	/*
		Přednastavené A-čkovské hranice (pro sichr)
	*/
	for (int i = 0; i < h1*w1*mainSize; i++) {
		A[i] = 0.0;
	}
	
  	int counterA = 0; // počitadlo podmínek, protože se chci posouvat ob jeden řádek, tak tam je mainSiye*counterA
  	for (int j = 0; j < h1; j++) {
		for (int i = 0; i < w1; i++) {
			/*
				Nastavení A (podmínka)
			*/
			/*
				Nastavení odpovídá -B_i <= t_i + s_iD_i <= 1 - B_i
				posunuji se vždy o řádek, proto mainSize*counterA
				u hranic se to nedělá kvůli tomu, že pro každou podmínku je jenom jedna hranice dolní/horní
			*/
			A[mainSize*counterA + j * w1 + i] = (double)(detailChannels[0].at<float>(j, i)/255.0);
			A[mainSize*counterA + w1*h1 + j * w1 + i] = 1.0;
			lbA[j*w1 + i] = -(double)(baseChannels[0].at<float>(j, i)/255.0); // dolní hranice A
			ubA[j*w1 + i] = 1.0 - (double)(baseChannels[0].at<float>(j, i)/255.0); // horní hranice A
			counterA++;
			A[mainSize*counterA + j * w1 + i] = (double)(detailChannels[1].at<float>(j, i)/255.0);
			A[mainSize*counterA + w1*h1 + j * w1 + i] = 1.0;
			lbA[1*h1*w1 + + j*w1 + i] = -(double)(baseChannels[1].at<float>(j, i)/255.0); // dolní hranice A
			ubA[1*h1*w1 + + j*w1 + i] = 1.0 - (double)(baseChannels[1].at<float>(j, i)/255.0); // horní hranice A
			counterA++;
			A[mainSize*counterA + j * w1 + i] = (double)(detailChannels[2].at<float>(j, i)/255.0);
			A[mainSize*counterA + w1*h1 + j * w1 + i] = 1.0;
			lbA[2*h1*w1 + j*w1 + i] = -(double)(baseChannels[2].at<float>(j, i)/255.0); // dolní hranice A
			ubA[2*h1*w1 + j*w1 + i] = 1.0 - (double)(baseChannels[2].at<float>(j, i)/255.0); // horní hranice A
			counterA++;
		}
	}


	QProblem example( mainSize, h1*w1*3); // počet parametrů, počet podmínek

	Options options;
	example.setOptions( options );
	int_t nWSR = 50000;
	example.init(H,g,A,lb, ub,lbA,ubA, nWSR);


	
	real_t *xOpt = new real_t[mainSize];
	example.getPrimalSolution( xOpt ); 

	/*
		Vložení do vytvořených matic
	*/
	int lineCounter = 0;
	int parametersCounterS = 0;
	int parametersCounterT = 0;
	for (int i = 0; i < mainSize; i++) {
			/*
				Zjištění, jestli je to parametr s či t, s parametry jsou jako první
			*/
			if (lineCounter < w1*h1) {
				s.at<float>(parametersCounterS / w1, parametersCounterS % w1) = xOpt[parametersCounterS];
				parametersCounterS++;				 			    
			} else {				
				t.at<float>(parametersCounterT / w1, parametersCounterT % w1) = xOpt[w1*h1 + parametersCounterT];
				parametersCounterT++; 
			}
			lineCounter++;
	}

	/*
		Uvolním paměť
	*/

	delete H;
	delete xOpt;
	delete A;
	delete lb;
	delete ub;
	delete g;
	delete lbA;
	delete ubA;
	std::ofstream outfile;
	outfile.open("s.txt");//std::ios_base::app
	outfile << s;
	outfile.close();

	std::ofstream outfile1;
	outfile1.open("t.txt");//std::ios_base::app
	outfile1 << t;
	outfile1.close();
	sAndT.push_back(s);
	sAndT.push_back(t);
	return sAndT;
}