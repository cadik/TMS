#include "cgalOptimization.h"

std::vector<cv::Mat> optimizeForGettingSAndTparametersWithCgal(int height, int width, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels) {
  	
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;
		
		std::vector<cv::Mat> sAndT;
    cv::Mat s;
		s = cv::Mat::ones(height, width, CV_32F);
		cv::Mat t;
    t = cv::Mat::zeros(height, width, CV_32F);

	/*
		Optimization
	*/
	std::cout << "Getting parameters completed" << std::endl;
	Program qp (CGAL::SMALLER, false, 0, false, 0);

  /*
		Creating D matrix
	*/
	cv::Mat DMatrix = cv::Mat::zeros(width*height*2, width*height*2, CV_32F);
	/*
    Adding condition, upper lower bounds etc
  */
  int counter = 0;
  for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			/*
				Lower Upper Bounds for s parameters, (not neccesary search in whole space)
			*/
		    qp.set_l(width*j + i, true, 10);
			/* 
				Conditions for counter, 0 <= (B_i + t_i) + D_i*s_i
				also -t_i - s_i * D_i <= B_i
			*/
			qp.set_a(width*height + width * j + i, counter, -1); //ti
			qp.set_a(width * j + i, counter, -(double)(detailChannels[0].at<float>(j, i)/255.0));
			qp.set_b(counter, (double)(baseChannels[0].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(width*height + width * j + i, counter, -1); //ti
			qp.set_a(width * j + i, counter, -(double)(detailChannels[1].at<float>(j, i)/255.0));
			qp.set_b(counter, (double)(baseChannels[1].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(width*height + width * j + i, counter, -1); //ti
			qp.set_a(width * j + i, counter, -(double)(detailChannels[2].at<float>(j, i)/255.0));
			qp.set_b(counter, (double)(baseChannels[2].at<float>(j, i)/255.0));
			counter++;

			/*
				Conditions for counter, (B_i + t_i) + D_i*s_i <= 1
				also t_i + D_i*s_i <= 1 - B_i
			*/
		    qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, (double)(detailChannels[0].at<float>(j, i)/255.0));
			qp.set_b(counter, 1 - (double)(baseChannels[0].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, (double)(detailChannels[1].at<float>(j, i)/255.0));
			qp.set_b(counter, 1 - (double)(baseChannels[1].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, (double)(detailChannels[2].at<float>(j, i)/255.0));
			qp.set_b(counter, 1 - (double)(baseChannels[2].at<float>(j, i)/255.0));
			counter++;
		}
	}

	/*
		Setting s
	*/
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			double sPowerTerm1 = -2*(pow((double)(detailImage.at<float>(j, i)/765.0), 2));
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

			if (i > 1) {
				sPowerTerm2 = 2*200*(double)weight1.at<float>(j, i - 1);
				DMatrix.at<float>(width * j + i, width * j + (i - 2)) = -sPowerTerm2;
				DMatrix.at<float>(width * j + (i - 2), width * j + i) = -sPowerTerm2;
			}
			if (j < (height - 2)) {
				sPowerTerm3 = 2*200*(double)weight1.at<float>(j + 1, i);
				DMatrix.at<float>(width * (j + 2) + i, width * j + i) = -sPowerTerm3;
				DMatrix.at<float>(width * j + i, width * (j + 2) + i) = -sPowerTerm3;
			}

			if (i < (width - 2)) {
				sPowerTerm4 = 2*200*(double)weight1.at<float>(j, i + 1);
				DMatrix.at<float>(width * j + (i + 2), width * j + i) = -sPowerTerm4;
				DMatrix.at<float>(width * j + i, width * j + (i + 2)) = -sPowerTerm4;
			}

			if (j > 1) {
				sPowerTerm5 = 2*200*(double)weight1.at<float>(j - 1, i);
				DMatrix.at<float>(width * j + i, width * (j - 2) + i) = -sPowerTerm5;
				DMatrix.at<float>(width * (j - 2) + i, width * j + i) = -sPowerTerm5;
			}
			
			if (i == 0) {
				sPowerTerm6 = (double)(2*200*(double)weight1.at<float>(j, i));
				DMatrix.at<float>(width * j + i, width * j + i + 1) = -sPowerTerm6;
				DMatrix.at<float>(width * j + i + 1, width * j + i) = -sPowerTerm6;
			}

			if (i == 1) {
				sPowerTerm7 = (double)(2*200*(double)weight1.at<float>(j, i - 1));
				DMatrix.at<float>(width * j + i, width * j + i - 1) = -sPowerTerm7;
				DMatrix.at<float>(width * j + i - 1, width * j + i) = -sPowerTerm7;
			}

			if (j == 0) {
				sPowerTerm8 = (double)(2*200*(double)weight1.at<float>(j, i));
				DMatrix.at<float>(width * j + i, width * (j + 1) + i) = -sPowerTerm8;
				DMatrix.at<float>(width * (j + 1) + i, width * j + i) = -sPowerTerm8;
			}

			if (j == 1) {
				sPowerTerm9 = (double)(2*200*(double)weight1.at<float>(j - 1, i));
				DMatrix.at<float>(width * j + i, width * (j - 1) + i) = -sPowerTerm9;
				DMatrix.at<float>(width * (j - 1) + i, width * j + i) = -sPowerTerm9;
			}

			if (i == width - 1) {
				sPowerTerm10 = (double)(2*200*(double)weight1.at<float>(j, i));
				DMatrix.at<float>(width * j + i, width * j + i - 1) = -sPowerTerm10;
				DMatrix.at<float>(width * j + i - 1, width * j + i) = -sPowerTerm10;
			}

			if (i == width - 2) {
				sPowerTerm11 = (double)(2*200*(double)weight1.at<float>(j, i + 1));
				DMatrix.at<float>(width * j + i, width * j + i + 1) = -sPowerTerm11;
				DMatrix.at<float>(width * j + i + 1, width * j + i) = -sPowerTerm11;
			}

			if (j == height - 1) {
				sPowerTerm12 = (double)(2*200*(double)weight1.at<float>(j, i));
				DMatrix.at<float>(width * j + i, width * (j - 1) + i) = -sPowerTerm12;
				DMatrix.at<float>(width * (j - 1) + i, width * j + i) = -sPowerTerm12;
			}

			if (j == height - 2) {
				sPowerTerm13 = (double)(2*200*(double)weight1.at<float>(j + 1, i));
				DMatrix.at<float>(width * j + i, width * (j + 1) + i) = -sPowerTerm13;
				DMatrix.at<float>(width * (j + 1) + i, width * j + i) = -sPowerTerm13;
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
			
			DMatrix.at<float>(j * width + i, j * width + i) = sumOfPowerTerms;
		}
	}
	

	/*
		Setting t
	*/

	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
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

			if (i > 1) {
				sPowerTerm01 = 2*500*(double)weight2.at<float>(j, i - 1); //check
				DMatrix.at<float>(width*height + width * j + i, width*height + height * j + (i - 2)) = -sPowerTerm01;
				DMatrix.at<float>(width*height + width * j + (i - 2), width*height + width * j + i) = -sPowerTerm01;
			}
			if (j < (height - 2)) {
				sPowerTerm02 = 2*500*(double)weight2.at<float>(j + 1, i);
				DMatrix.at<float>(width*height + width * (j + 2) + i, width*height + width * j + i) = -sPowerTerm02;
				DMatrix.at<float>(width*height + width * j + i, width*height + width * (j + 2) + i) = -sPowerTerm02;
			}

			if (i < (width - 2)) {
				sPowerTerm03 = 2*500*(double)weight2.at<float>(j, i + 1);
				DMatrix.at<float>(width*height + width * j + (i + 2), width*height + width * j + i) = -sPowerTerm03;
				DMatrix.at<float>(width*height + width * j + i, width*height + width * j + (i + 2)) = -sPowerTerm03;
			}

			if (j > 1) {
				sPowerTerm04 = 2*500*(double)weight2.at<float>(j - 1, i);
				DMatrix.at<float>(width*height + width * j + i, width*height + width * (j - 2) + i) = -sPowerTerm04;
				DMatrix.at<float>(width*height + width * (j - 2) + i, width*height + width * j + i) = -sPowerTerm04;
			}

			if (i == 0) {
				sPowerTerm05 = (double)(2*500*(double)weight2.at<float>(j, i));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * j + i + 1) = -sPowerTerm05;
				DMatrix.at<float>(width*height + width * j + i + 1, width*height + width * j + i) = -sPowerTerm05;
			}

			if (i == 1) {
				sPowerTerm06 = (double)(2*500*(double)weight2.at<float>(j, i - 1));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * j + i - 1) = -sPowerTerm06;
				DMatrix.at<float>(width*height + width * j + i - 1, width*height + width * j + i) = -sPowerTerm06;
			}
		
			if (j == 0) {
				sPowerTerm07 = (double)(2*500*(double)weight2.at<float>(j, i));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * (j + 1) + i) = -sPowerTerm07;
				DMatrix.at<float>(width*height + width * (j + 1) + i, width*height + width * j + i) = -sPowerTerm07;
			}

			if (j == 1) {
				sPowerTerm08 = (double)(2*500*(double)weight2.at<float>(j - 1, i));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * (j - 1) + i) = -sPowerTerm08;
				DMatrix.at<float>(width*height + width * (j - 1) + i, width*height + width * j + i) = -sPowerTerm08;
			}

			if (i == width - 1) {
				sPowerTerm09 = (double)(2*500*(double)weight2.at<float>(j, i));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * j + i - 1) = -sPowerTerm09;
				DMatrix.at<float>(width*height + width * j + i - 1, width*height + width * j + i) = -sPowerTerm09;
			}

			if (i == width - 2) {
				sPowerTerm010 = (double)(2*500*(double)weight2.at<float>(j, i + 1));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * j + i + 1) = -sPowerTerm010;
				DMatrix.at<float>(width*height + width * j + i + 1, width*height + width * j + i) = -sPowerTerm010;
			}

			if (j == height - 1) {
				sPowerTerm011 = (double)(2*500*(double)weight2.at<float>(j, i));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * (j - 1) + i) = -sPowerTerm011;
				DMatrix.at<float>(width*height + width * (j - 1) + i, width*height + width * j + i) = -sPowerTerm011;
			}

			if (j == height - 2) {
				sPowerTerm012 = (double)(2*500*(double)weight2.at<float>(j + 1, i));
				DMatrix.at<float>(width*height + width * j + i, width*height + width * (j + 1) + i) = -sPowerTerm012;
				DMatrix.at<float>(width*height + width * (j + 1) + i, width*height + width * j + i) = -sPowerTerm012;
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
			DMatrix.at<float>(width*height +j * width + i, width*height + j * width + i) = sumOfPowerTerms1;
		}
	}

	cv::Mat transposedDMatrix;
	transpose(DMatrix, transposedDMatrix);

	cv::Mat newDMatrix = DMatrix * transposedDMatrix;
	/*
			Getting solution
	*/

	for (int j = 0; j < width*height*2; j++) {
		for (int i = 0; i < width*height*2; i++) {
			/*
				setting D
			*/
			if (i == j && newDMatrix.at<float>(j, i) > 0) {
				qp.set_d(j, i, newDMatrix.at<float>(j, i));
			} else {
				if (newDMatrix.at<float>(j, i) > 0 && j > i) {
					qp.set_d(j, i, newDMatrix.at<float>(j, i));
				}
			}
		}
	}
	/*
		Deleting necessary matrices from memory
	*/

	newDMatrix.release();
	DMatrix.release();

	transposedDMatrix.release();
	
    Solution solution = CGAL::solve_quadratic_program(qp, ET());
    Solution::Variable_value_iterator it = solution.variable_values_begin();
	Solution::Variable_value_iterator end = solution.variable_values_end();
	int lineCounter = 0;
	int parametersCounterS = 0;
	int parametersCounterT = 0;
	for (; it != end; ++it) {
			if (lineCounter < width*height) {
				s.at<float>(parametersCounterS / width, parametersCounterS % width) = CGAL::to_double(*it);
				parametersCounterS++;				 			    
			} else {				
				t.at<float>(parametersCounterT / width, parametersCounterT % width) = CGAL::to_double(*it);
				parametersCounterT++; 
			}
			lineCounter++;
	}

	// Checking because of not founded details
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			if (s.at<float>(j, i) == 10 && t.at<float>(j, i) == 0) {
				s.at<float>(j, i) = 1;
				t.at<float>(j, i) = 0;
			}
		}
	}

	// std::cout << t << std::endl;
	sAndT.push_back(t);
  	sAndT.push_back(s);
  	
	return sAndT;
}