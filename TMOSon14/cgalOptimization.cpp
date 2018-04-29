#include "cgalOptimization.h"

void tryIt() {
    // by default, we have a nonnegative QP with Ax <= b
  srand ( time(NULL) );
  /*
    Creating Matrices
  */

  cv::Mat detailLayer;
	detailLayer = cv::Mat::zeros(15, 30, CV_32F);
  cv::Mat w1;
	w1 = cv::Mat::zeros(15, 30, CV_32F);
  cv::Mat w2;
	w2 = cv::Mat::zeros(15, 30, CV_32F);
  Program qp (CGAL::SMALLER, false, 0, false, 0); 
  

  for (int j = 0; j < 15; j++) {
      for (int i = 0; i < 30; i++) {
        detailLayer.at<float>(j, i) = rand() % 255/3;
        w1.at<float>(j, i) = (rand() % 255);
        w2.at<float>(j, i) = (rand() % 255);
      }
  }
 // std::cout << w2 << std::endl;
  // now set the non-default entries

  /*
    Setting s
  */
  /*
  for (int j = 0; j < 15; j++) {
      for (int i = 0; i < 30; i++) {
        double sPowerTerm1 = -2 * pow(detailLayer.at<float>(j, i), 2);
        double sPowerTerm2 = 0.0;
        double sPowerTerm3 = 0.0;
        double sPowerTerm4 = 0.0;
        double sPowerTerm5 = 0.0;
        double sumOfPowerTerms = 0.0;

        if (i > 1) {
          //std::cout << "111" << std::endl;
          sPowerTerm2 = 2*500*w1.at<float>(j, i - 1);
          qp.set_d(30 * j + i, 30 * j + (i - 2), sPowerTerm2);
          
        }

        if (j < (15 - 2)) {
          //std::cout << "222" << std::endl;
          sPowerTerm3 = 2*500*w1.at<float>(j + 1, i);
          qp.set_d(30 * j + i + 60, 30 * j + i, sPowerTerm3);
        }

        if (i < (30 - 2)) {
          //std::cout << "333" << std::endl;
          sPowerTerm4 = 2*500*w1.at<float>(j, i + 1);
          qp.set_d(30 * j + (i + 2), 30 * j + i, sPowerTerm4);
        }

        if (j > 1) {
          sPowerTerm5 = 2*500*w1.at<float>(j - 1, i);
          qp.set_d(30 * j + i, 30 * j + i - 60, sPowerTerm5);
        }

        sumOfPowerTerms = sPowerTerm1 +
                          sPowerTerm2 +
                          sPowerTerm3 +
                          sPowerTerm4 +
                          sPowerTerm5; 
        // setting powers
        qp.set_d(30 * j + i, 30 * j + i, sumOfPowerTerms);
        
      }
  }
  
  for (int j = 0; j < 15; j++) {
      for (int i = 0; i < 30; i++) {
       // double sPowerTerm11 = -2 * pow(detailLayer.at<float>(j, i), 2);
        double sPowerTerm21 = 0.0;
        double sPowerTerm31 = 0.0;
        double sPowerTerm41 = 0.0;
        double sPowerTerm51 = 0.0;
        double sumOfPowerTerms1 = 0.0;

        if (i > 1) {
          //std::cout << "111" << std::endl;
          sPowerTerm21 = 2*500*w2.at<float>(j, i - 1);
          qp.set_d(30*15 + 30 * j + i,30*15 +  30 * j + (i - 2), sPowerTerm21);
          
        }

        if (j < (15 - 2)) {
          //std::cout << "222" << std::endl;
          sPowerTerm31 = 2*500*w2.at<float>(j + 1, i);
          qp.set_d(30*15 + 30 * j + i + 60,30*15 +  30 * j + i, sPowerTerm31);
        }

        if (i < (30 - 2)) {
          //std::cout << "333" << std::endl;
          sPowerTerm41 = 2*500*w2.at<float>(j, i + 1);
          qp.set_d(30*15 + 30 * j + (i + 2),30*15 +  30 * j + i, sPowerTerm41);
        }

        if (j > 1) {
          sPowerTerm51 = 2*500*w2.at<float>(j - 1, i);
          qp.set_d(30*15 + 30 * j + i,30*15 +  30 * j + i - 60, sPowerTerm51);
        }

        sumOfPowerTerms1 = sPowerTerm21 +
                          sPowerTerm31 +
                          sPowerTerm41 +
                          sPowerTerm51; 
        
      }
  }
  */
  qp.set_d(0, 0, 6);
  qp.set_d(1, 1, 2);
  qp.set_d(1, 0, 2);
  qp.set_c(0, 1);
  qp.set_c(1, 6);
  qp.set_c0(2);

 /*  const int X = 0; 
  const int Y = 1;
  qp.set_a(X, 0,  1); qp.set_a(Y, 0, 1); qp.set_b(0, 7);  //  x + y  <= 7
  qp.set_a(X, 1, -1); qp.set_a(Y, 1, 2); qp.set_b(1, 4);  // -x + 2y <= 4
  qp.set_u(Y, true, 4);                                   //       y <= 4
  qp.set_d(X, X, 2); qp.set_d (Y, Y, 8); // !!specify 2D!!    x^2 + 4 y^2
  qp.set_c(Y, -32);                                       // -32y
  qp.set_c0(64); */

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
	
	/*
		Setting s
	*/
  
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			double sPowerTerm1 = -2 * pow(detailImage.at<float>(j, i)/255.0, 2);
			double sPowerTerm2 = 0.0;
			double sPowerTerm3 = 0.0;
			double sPowerTerm4 = 0.0;
			double sPowerTerm5 = 0.0;
			double sumOfPowerTerms = 0.0;

			if (i > 1) {
			//std::cout << "111" << std::endl;
				sPowerTerm2 = 2*200*weight1.at<float>(j, i - 1);
				qp.set_d(width * j + i, width * j + (i - 2), sPowerTerm2);
			
			}

			if (j < (height - 2)) {
			//std::cout << "222" << std::endl;
				sPowerTerm3 = 2*200*weight1.at<float>(j + 1, i);
				qp.set_d(width * j + i + width*2, width * j + i, sPowerTerm3);
			}

			if (i < (width - 2)) {
			//std::cout << "333" << std::endl;
				sPowerTerm4 = 2*200*weight1.at<float>(j, i + 1);
				qp.set_d(width * j + (i + 2), width * j + i, sPowerTerm4);
			}

			if (j > 1) {
				sPowerTerm5 = 2*200*weight1.at<float>(j - 1, i);
				qp.set_d(width * j + i, width * j + i - width*2, sPowerTerm5);
			}

			sumOfPowerTerms = sPowerTerm1 +
							sPowerTerm2 +
							sPowerTerm3 +
							sPowerTerm4 +
							sPowerTerm5; 
			// setting powers
			qp.set_d(width * j + i, width * j + i, sumOfPowerTerms);
			
		}
	}
	
	/*
		Setting t
	*/
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			double sPowerTerm21 = 0.0;
			double sPowerTerm31 = 0.0;
			double sPowerTerm41 = 0.0;
			double sPowerTerm51 = 0.0;
			double sumOfPowerTerms1 = 0.0;

			if (i > 1) {
				//std::cout << "111" << std::endl;
				sPowerTerm21 = 2*500*weight2.at<float>(j, i - 1);
				qp.set_d(width*height + width * j + i, width*height +  width * j + (i - 2), sPowerTerm21);
				
			}

			if (j < (15 - 2)) {
				//std::cout << "222" << std::endl;
				sPowerTerm31 = 2*500*weight2.at<float>(j + 1, i);
				qp.set_d(width*height + width * j + i + width*2, width*height +  30 * j + i, sPowerTerm31);
			}

			if (i < (width - 2)) {
				//std::cout << "333" << std::endl;
				sPowerTerm41 = 2*500*weight2.at<float>(j, i + 1);
				qp.set_d(width*height + width * j + (i + 2), width*height +  width * j + i, sPowerTerm41);
			}

			if (j > 1) {
				sPowerTerm51 = 2*500*weight2.at<float>(j - 1, i);
				qp.set_d(width*height + width * j + i, width*height +  width * j + i - width*2, sPowerTerm51);
			}

			sumOfPowerTerms1 = sPowerTerm21 +
							sPowerTerm31 +
							sPowerTerm41 +
							sPowerTerm51; 
			qp.set_d(width*height + width * j + i, width*height + width * j + i, sumOfPowerTerms1);
		}
	}

  /*
    Adding condition, upper lower bounds etc
  */
  int counter = 0;
  for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			/*
				Lower Upper Bounds for s parameters, (not neccesary search in whole space)
			*/
			qp.set_l(width * j + i, true, 0.001);
			qp.set_u(width * j + i, true, 25);
			/*
				Conditions for counter, 0 <= (B_i + t_i) + D_i*s_i
			*/
			qp.set_a(width*height + width * j + i, counter, -1); //ti
			qp.set_a(width * j + i, counter, -detailChannels[0].at<float>(j, i)/255.0);
			qp.set_b(counter, -0.004 + baseChannels[0].at<float>(j, i)/255.0);
			counter++;
			qp.set_a(width*height + width * j + i, counter, -1); //ti
			qp.set_a(width * j + i, counter, -detailChannels[1].at<float>(j, i)/255.0);
			qp.set_b(counter, -0.004 + baseChannels[1].at<float>(j, i)/255.0);
			counter++;
			qp.set_a(width*height + width * j + i, counter, -1); //ti
			qp.set_a(width * j + i, counter, -detailChannels[2].at<float>(j, i)/255.0);
			qp.set_b(counter, -0.004 + baseChannels[2].at<float>(j, i)/255.0);
			counter++;

			/*
				Conditions for counter, (B_i + t_i) + D_i*s_i <= 1
			*/
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, detailChannels[0].at<float>(j, i)/255.0);
			qp.set_b(counter, 1 - baseChannels[0].at<float>(j, i)/255.0);
			counter++;
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, detailChannels[1].at<float>(j, i)/255.0);
			qp.set_b(counter, 1 - baseChannels[1].at<float>(j, i)/255.0);
			counter++;
			qp.set_a(width*height + width * j + i, counter, 1); //ti
			qp.set_a(width * j + i, counter, detailChannels[2].at<float>(j, i)/255.0);
			qp.set_b(counter, 1 - baseChannels[2].at<float>(j, i)/255.0);
			counter++;
		}
	}
	Solution solution = CGAL::solve_quadratic_program(qp, ET());
  	assert (solution.solves_quadratic_program(qp));
	/*
		Creating tmpDetailFile
	*/
	ofstream fileOutputDetailMaximalization;
	fileOutputDetailMaximalization.open ("outputDetailMaximalization.txt");
	fileOutputDetailMaximalization << solution;
	fileOutputDetailMaximalization.close();

	std::ifstream infile("../outputDetailMaximalization.txt");
	std::string line;
	
	/*
		Getting variable values from tmp file
	*/
	std::vector<double> parameters;
	int lineCounter = 0;
	int parametersCounterS = 0;
	int parametersCounterT = 0;

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
			if (lineCounter > 3 + width*height) {
				t.at<float>(parametersCounterT / width, parametersCounterT % width) = number / number1;
				parametersCounterT++;      
			} else {
				s.at<float>(parametersCounterS / width, parametersCounterS % width) = number / number1;
				parametersCounterS++;
			}
		}
	}
  sAndT.push_back(s);
  sAndT.push_back(t);
	/*
		Removing file
	*/
	if(remove("outputDetailMaximalization.txt" ) != 0 ) {
		std::cout <<  "Error deleting file" << std::endl;
	} else {
		std::cout << "File successfully deleted" << std::endl;
	}
	return sAndT;
}