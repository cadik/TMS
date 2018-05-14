#include "cgalOptimization.h"

void tryIt() {

typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;
	Program qp (CGAL::SMALLER, false, 0, false, 0);
	qp.set_d(0, 0, 2); qp.set_d (1, 1, 0.08); // !!specify 2D!!    x^2 + 4 y^2
  qp.set_c(1, -32);                                       // -32y
  qp.set_c0(64);

  Solution s = CGAL::solve_quadratic_program(qp, ET());
  assert (s.solves_quadratic_program(qp));

  /*
    Creating tmpDetailFile
  */
  ofstream fileOutputDetailMaximalization;
  fileOutputDetailMaximalization.open ("outputDetailMaximalization.txt");
  fileOutputDetailMaximalization << s;
  fileOutputDetailMaximalization.close();

  std::ifstream infile("outputDetailMaximalization.txt");
  std::string line;
  
  /*
    Getting variable values
  */
  int lineCounter = 0;
  while (std::getline(infile, line)) {
      lineCounter++;
      std::istringstream iss(line);
      if (lineCounter > 3) {
        bool isSlash = false;
        string nominator = "";
        string denominator = "";
        for (int i = 5; i < line.length(); i++) {
          if (isSlash == false) {
            if (line[i] == '/') {
              isSlash = true;
            } else {
              nominator += line[i];
            } 
          } else {
            denominator += line[i];
          }
        }
        isSlash = false;
        double number = stold(nominator);
        double number1 = stold(denominator);
        std::cout << number / number1;
        std::cout << std::endl;
      }
  }
  /*
    Removing file
  */
  if(remove("outputDetailMaximalization.txt" ) != 0 )
    std::cout <<  "Error deleting file" << std::endl;
  else
    std::cout << "File successfully deleted" << std::endl;
}

std::vector<cv::Mat> optimizeForGettingSAndTparametersWithCgal(int height, int width, cv::Mat detailImage, cv::Mat weight1, cv::Mat weight2, std::vector<cv::Mat> baseChannels, std::vector<cv::Mat> detailChannels) {
  	
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;
		
		std::vector<cv::Mat> sAndT;
    cv::Mat s;
		s = cv::Mat::zeros(height, width, CV_32F);
		cv::Mat t;
    t = cv::Mat::zeros(height, width, CV_32F);

	/*
		Optimization
	*/
	std::cout << "Getting parameters completed" << std::endl;
	Program qp (CGAL::SMALLER, false, 0, false, 0);
	
 //ykouska

/*	qp.set_a(0, 0,  1); qp.set_a(1, 0, 1); qp.set_b(0, 7);  //  x + y  <= 7
  qp.set_a(0, 1, -1); qp.set_a(1, 1, 2); qp.set_b(1, 4);  // -x + 2y <= 4
  qp.set_u(1, true, 4);                                   //       y <= 4*/
	/*qp.set_a(0, 45, -2); qp.set_a(1, 45, -3); qp.set_b(45, -4);
	qp.set_l(0, true, 0); 
	qp.set_l(1, true, 0); 
  qp.set_d(0, 0, 6); qp.set_d (1, 1, 2);                  // 3x^2 + y^2*/
	// qp.set_d(1, 0, -268);
 /* qp.set_c(0, 1);    // 1x
	qp.set_c(1, 6); //6y
  qp.set_c0(2);  //+2*/
	/*
    Adding condition, upper lower bounds etc
  */
  int counter = 0;
  for (int j = 0; j < 7; j++) {
		for (int i = 0; i < 7; i++) {
			/*
				Lower Upper Bounds for s parameters, (not neccesary search in whole space)
			*/
			qp.set_l(7*j + i, true, -5); // s +- 1
		 // qp.set_u(7*j + i, true, 5);
			qp.set_l(49 + 7*j + i, true, 0.1); // t +- 0
		 // qp.set_u(49 + 7*j + i, true, 10);
			/* 
				Conditions for counter, 0 <= (B_i + t_i) + D_i*s_i
				also -t_i - s_i * D_i <= B_i
			*/
			qp.set_a(7*7 + 7 * j + i, counter, -1); //ti
			qp.set_a(7 * j + i, counter, -(double)(detailChannels[0].at<float>(j, i)/255.0));
			qp.set_b(counter, (double)(baseChannels[0].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(7*7 + 7 * j + i, counter, -1); //ti
			qp.set_a(7 * j + i, counter, -(double)(detailChannels[1].at<float>(j, i)/255.0));
			qp.set_b(counter, (double)(baseChannels[1].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(7*7 + 7 * j + i, counter, -1); //ti
			qp.set_a(7 * j + i, counter, -(double)(detailChannels[2].at<float>(j, i)/255.0));
			qp.set_b(counter, (double)(baseChannels[2].at<float>(j, i)/255.0));
			counter++;

			/*
				Conditions for counter, (B_i + t_i) + D_i*s_i <= 1
				also t_i + D_i*s_i <= 1 - B_i
			*/
		  /*qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, (double)(detailChannels[0].at<float>(j, i)/255.0));
			qp.set_b(counter, 1 - (double)(baseChannels[0].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, (double)(detailChannels[1].at<float>(j, i)/255.0));
			qp.set_b(counter, 1 - (double)(baseChannels[1].at<float>(j, i)/255.0));
			counter++;
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, (double)(detailChannels[2].at<float>(j, i)/255.0));
			qp.set_b(counter, 1 - (double)(baseChannels[2].at<float>(j, i)/255.0));*/
		//	counter++;*/
		}
	}

	/*
		Setting s
	*/
	for (int j = 0; j < 7; j++) {
		for (int i = 0; i < 7; i++) {
			// std::cout << detailImage.at<float>(j, i) << std::endl;
			double sPowerTerm1 = - 2*(pow((double)(detailImage.at<float>(j, i)/765.0), 2))*1000;
		//	std::cout << sPowerTerm1 << std::endl;
			double sPowerTerm2 = 0.0;
			double sPowerTerm3 = 0.0;
			double sPowerTerm4 = 0.0;
			double sPowerTerm5 = 0.0;
			double sPowerTerm6 = 0.0;
			double sPowerTerm7 = 0.0;
			double sPowerTerm8 = 0.0;
			double sPowerTerm9 = 0.0;
			double sumOfPowerTerms = 0.0;

			if (i > 1) {
				sPowerTerm2 = 2*200*(double)weight1.at<float>(j, i - 1)*1000; //check
				qp.set_d(7 * j + i, 7 * j + (i - 2), -sPowerTerm2); // check
			}
			if (j < (7 - 2)) {
				sPowerTerm3 = 2*200*(double)weight1.at<float>(j + 1, i)*1000;
				qp.set_d(7 * (j + 2) + i, 7 * j + i, -sPowerTerm3);
			}

			if (i < (7 - 2)) {
				sPowerTerm4 = 2*200*(double)weight1.at<float>(j, i + 1)*1000;
				qp.set_d(7 * j + (i + 2), 7 * j + i, -sPowerTerm4);
			}

			if (j > 1) {
				sPowerTerm5 = 2*200*(double)weight1.at<float>(j - 1, i)*1000;
				qp.set_d(7 * j + i, 7 * (j - 2) + i, -sPowerTerm5);
			}
			
		/*	if (i == 0) {
				sPowerTerm6 = (double)(2*200*(double)weight1.at<float>(j, i));
				qp.set_d(7 * j + i + 1, 7 * j + i , -sPowerTerm6);
			}

			if (i == width - 1) {
				sPowerTerm7 = (double)(2*200*weight1.at<float>(j, i));
				 qp.set_d(7 * j + i, 7 * j + i - 1, -sPowerTerm7);
			}

			if (j == 0) {
				sPowerTerm8 = (double)(2*200*weight1.at<float>(j, i));
				qp.set_d((7 + 1) * j + i, 7 * j + i, -sPowerTerm8);
			}

			if (j == height - 1) {
				sPowerTerm9 = (double)(2*200*weight1.at<float>(j, i));
				qp.set_d(7 * j + i, (7 - 1) * j + i, -sPowerTerm9);
			}
*/

			sumOfPowerTerms = sPowerTerm1 +
							sPowerTerm2 +
							sPowerTerm3 +
							sPowerTerm4 +
							sPowerTerm5 +
							sPowerTerm6 + 
							sPowerTerm7 + 
							sPowerTerm8 + 
							sPowerTerm9; 
			
			qp.set_d(7 * j + i, 7 * j + i, sumOfPowerTerms);
		//	qp.set_c(width * j + i, sumOfPowerTerms/4);*/
		}
	}
	

	/*
		Setting t
	*/

	for (int j = 0; j < 7; j++) {
		for (int i = 0; i < 7; i++) {
			double sPowerTerm21 = 0.0;
			double sPowerTerm31 = 0.0;
			double sPowerTerm41 = 0.0;
			double sPowerTerm51 = 0.0;
			double sumOfPowerTerms1 = 0.0;

			if (i > 1) {
				sPowerTerm21 = 2*500*(double)(weight2.at<float>(j, i - 1))*1000; //check
				qp.set_d(7*7 + 7 * j + i, 7 * 7 + 7 * j + (i - 2), -sPowerTerm21); // check
			}

			if (j < (7 - 2)) {
				sPowerTerm31 = 2*500*(double)(weight1.at<float>(j + 1, i))*1000;
				qp.set_d(7*7 + 7 * (j + 2) + i,7*7 +  7 * j + i, -sPowerTerm31);
			}

			if (i < (7 - 2)) {
				sPowerTerm41 = 2*500*(double)(weight2.at<float>(j, i + 1))*1000;
				qp.set_d(7*7 + 7 * j + (i + 2), 7*7 +  7 * j + i, -sPowerTerm41);
			}

			if (j > 1) {
				sPowerTerm51 = 2*500*(double)(weight2.at<float>(j - 1, i))*1000;
				qp.set_d(7*7 + 7 * j + i, 7*7 +  7 * (j - 2) + i, -sPowerTerm51);
			}


		/*if (i == 0) {
				sPowerTerm6 = (double)(2*200*weight1.at<float>(j, i));
				qp.set_d(7*7 + 7 * j + i + 1, 7 * j + i , -sPowerTerm6);
			}

			if (i == width - 1) {
				sPowerTerm7 = (double)(2*200*weight1.at<float>(j, i));
				 qp.set_d(7*7 + 7 * j + i, 7 * j + i - 1, -sPowerTerm7);
			}

			if (j == 0) {
				sPowerTerm8 = (double)(2*200*weight1.at<float>(j, i));
				qp.set_d(7*7 + (7 + 1) * j + i, 7 * j + i, -sPowerTerm8);
			}

			if (j == height - 1) {
				sPowerTerm9 = (double)(2*200*weight1.at<float>(j, i));
				qp.set_d(7*7 + 7 * j + i, (7 - 1) * j + i, -sPowerTerm9);
			}*/
			sumOfPowerTerms1 = sPowerTerm21 +
							sPowerTerm31 +
							sPowerTerm41 +
							sPowerTerm51;
			qp.set_d(7*7 + 7 * j + i, 7*7 + 7 * j + i, sumOfPowerTerms1);
			// std::cout << qp.get_d();
		}
	}
	// qp.set_c0(456456465465);
  CGAL::Quadratic_program_options options;
  // options.set_verbosity(1);                         // verbose mode 
  // options.set_pricing_strategy( CGAL::QP_PARTIAL_FILTERED_DANTZIG); 
	// options.set_auto_validation(true); 
	// solve the program, using ET as the exact type
  Solution solution = CGAL::solve_quadratic_program(qp, ET(), options);
  Solution::Variable_value_iterator it = solution.variable_values_begin();
	Solution::Variable_value_iterator end = solution.variable_values_end();
	int lineCounter = 0;
	int parametersCounterS = 0;
	int parametersCounterT = 0;
	for (; it != end; ++it) {
		std::cout << CGAL::to_double(*it) << std::endl;			
			if (lineCounter >= 7*7) {
				t.at<float>(parametersCounterT / 7, parametersCounterT % width) = CGAL::to_double(*it);
				parametersCounterT++;      
			} else {
				s.at<float>(parametersCounterS / 7, parametersCounterS % width) = CGAL::to_double(*it);
				parametersCounterS++;
			}
			lineCounter++;
	}
  sAndT.push_back(s);
  sAndT.push_back(t);
	return sAndT;
}