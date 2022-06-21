/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio                                    *
*                                                                              *
*                       Brno 2022                                              *
*                                                                              *
*                       Implementation of the TMOShan10 class                  *
*                                                                              *
*******************************************************************************/

#include "TMOShan10.h"

#define DIMS 2

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOShan10::TMOShan10()
{
   SetName(L"Shan10");
   SetDescription(L"Globally Optimized Linear Windowed Tone-Mapping");

   levelNum.SetName(L"levelNum");
   levelNum.SetDescription(L"Maltigrid recursion number.");
   levelNum.SetDefault(3);
   levelNum = 3;
   levelNum.SetRange(2, 5);
   this->Register(levelNum);

   winSize.SetName(L"winSize");
   winSize.SetDescription(L"The local window size about each pixel.");
   winSize.SetDefault(3);
   winSize = 3;
   winSize.SetRange(3, 15);
   this->Register(winSize);

   sSat.SetName(L"sSat");
   sSat.SetDescription(L"The saturation factor in Eq. (5).");
   sSat.SetDefault(0.6);
   sSat = 0.6;
   sSat.SetRange(0.4, 0.6);
   this->Register(sSat);

   beta3.SetName(L"beta3");
   beta3.SetDescription(L"Beta3 parameter. Article Eq. (6).");
   beta3.SetDefault(0.0);
   beta3 = 0.0;
   beta3.SetRange(0.0, 0.1);
   this->Register(beta3);

   beta2.SetName(L"beta2");
   beta2.SetDescription(L"Beta2 parameter. Article Eq. (6).");   
   beta2.SetDefault(0.0);
   beta2 = 0.0;
   beta2.SetRange(0.0, 0.4);
   this->Register(beta2);

   beta1.SetName(L"beta1");
   beta1.SetDescription(L"Beta1 parameter. Article Eq. (6).");
   beta1.SetDefault(1.2);
   beta1 = 1.2;
   beta1.SetRange(0.4, 1.4);
   this->Register(beta1);
}

TMOShan10::~TMOShan10()
{
}

/*
 * Convert TMOImage to cv::Mat
 */
Mat TMOImage2Mat(TMOImage *pSrc)
{
   double* pSourceData;
   int rowsCnt, colsCnt;

   pSourceData = pSrc->GetData();
   rowsCnt = pSrc->GetHeight();
   colsCnt = pSrc->GetWidth();

   Mat srcConvMat(rowsCnt, colsCnt, CV_64FC3);

   if (srcConvMat.isContinuous())
   {
      colsCnt *= rowsCnt;
      rowsCnt = 1;
   }

   for (int y = 0; y < rowsCnt; y++)
   {
      for (int x = 0; x < colsCnt; x++)
      {
         for (int c = 0; c < CHANNELSCNT; c++)
         {
            srcConvMat.at<double>(y, x*CHANNELSCNT+c) = pSourceData[CHANNELSCNT-c-1];
         }

         /* Add count of channels (RGB) to pointer */
         pSourceData += CHANNELSCNT;
      }
   }

   return srcConvMat;
}

/*
 * Create guidance map.
 */
Mat TMOShan10::generateGuidanceMap(Mat lumo,
                                   int winSize,
                                   double beta1,
                                   double beta2,
                                   const double beta3,
                                   const double kappa)
{
   int height, width;
   int halfWinSize;

   height = lumo.rows;
   width = lumo.cols;

   halfWinSize = (winSize - 1) / 2;

   Mat meanFilter = Mat::ones(winSize, winSize, lumo.type());

   /*
    * Source: https://stackoverflow.com/questions/21874774/sum-of-elements-in-a-matrix-in-opencv
    * Author: Shai (https://stackoverflow.com/users/1714410/shai)
    */
   meanFilter = meanFilter / cv::sum(meanFilter)[0];

   Mat meanLum(height, width, lumo.type());

   filter2D(lumo, meanLum, -1, meanFilter, Point(-1, -1), 0, BORDER_REPLICATE);


   Mat tLum = matUtil.padArray(lumo, halfWinSize, halfWinSize, PADARR_REPLICATE);

   Mat tmpStdLum = matUtil.colFiltSlidingStd(tLum, winSize, winSize);

   Mat stdLum(height, width, lumo.type());

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         stdLum.at<double>(y, x) = tmpStdLum.at<double>(y + halfWinSize, x + halfWinSize);
      }
   }

   Mat ci(height, width, lumo.type());

   pow(meanLum, beta1, meanLum);
   pow(stdLum, beta2, stdLum);
   pow(lumo, beta3, lumo);

   ci = meanLum.mul(stdLum);
   ci = ci.mul(lumo);   

   ci = ci + kappa;

   ci = ci * 0.2;  
   ci = 1.0 / ci;   

   meanFilter.deallocate();
   meanLum.deallocate();
   tLum.deallocate();
   tmpStdLum.deallocate();
   stdLum.deallocate();

   return ci;
}

/*
 * Returns linear coefficients defined in Eq. (17) in the appendix
 * of the article.
 */
Mat TMOShan10::getLinearCoefficients(Mat lum,
                                     Mat lumo,
                                     Mat map,                                     
                                     uint8_t multiGridFilt,
                                     int winSize,
                                     double epsilon)
{
   int height, width;
   int halfWinSize;
   int xStart, xEnd;
   int yStart, yEnd;
   int m;   

   height = lumo.rows;
   width = lumo.cols;

   halfWinSize = (winSize-1) / 2;

   Mat linearCoefficients = Mat::zeros(height, width, CV_64FC2);
   Mat winLumo;
   Mat winLum;
   Mat invHi;   

   for (int iWidth = 0; iWidth < width; iWidth++)
   {
      xStart = max(0, iWidth-halfWinSize);
      xEnd = min(width-1, iWidth+halfWinSize);

      for (int iHeight = 0; iHeight < height; iHeight++)
      {
         yStart = max(0, iHeight-halfWinSize);
         yEnd = min(height-1, iHeight+halfWinSize);

         m = (yEnd - yStart + 1)*(xEnd - xStart + 1);
         
         winLumo = lumo(Range(yStart, yEnd+1), Range(xStart, xEnd+1));
         winLum = lum(Range(yStart, yEnd+1), Range(xStart, xEnd+1));

         double mean = 0.0;

         // Standard deviation
         for (int y = 0; y < winLumo.rows; y++)
         {
            // Sum of all elements
            for (int x = 0; x < winLumo.cols; x++)
            {
               mean += winLumo.at<double>(y, x);
            }
         }

         // Calculate mean
         mean /= (double)(winLumo.rows * winLumo.cols);

         double sum = 0.0;

         for (int y = 0; y < winLumo.rows; y++)
         {
            for (int x = 0; x < winLumo.cols; x++)
            {
               sum += pow(winLumo.at<double>(y,x) - mean, 2.0);
            }
         }

         // The standard deviation is normalized by the number of observations.
         double tVar = sum / (winLumo.rows * winLumo.cols);
         double tMap = map.at<double>(iHeight, iWidth);

         double mul = pow(tMap, (double)(-multiGridFilt));
         double tDelta = tVar + (epsilon/(double)m) * mul;         


         // Article appendix equasions 13 and 14
         double divider = m * tDelta;         

         double matrixFst = 1.0/divider;
         double matrixSecThd =  (-mean)/divider;
         double matrixFourth = (tDelta + mean * mean)/divider;                       

         invHi = (Mat_<double>(2, 2) << matrixFst, matrixSecThd, matrixSecThd, matrixFourth);       

         // Source: https://stackoverflow.com/questions/21874774/sum-of-elements-in-a-matrix-in-opencv
         // Author: Shai

         // Source: https://stackoverflow.com/questions/43180053/initialize-the-values-into-mat-object-in-opencv
         // Author: Micka
         
         Mat etai = (Mat_<double>(2, 1) << (epsilon / tMap + cv::sum(winLumo.mul(winLum))[0]), cv::sum(winLum)[0]);         
         
         // Mathematical multiplication
         Mat coeffs = (invHi * etai) / 0.8;
         

         linearCoefficients.at<Vec2d>(iHeight, iWidth) = Vec2d(coeffs.at<double>(0, 0), coeffs.at<double>(1, 0));         

         winLumo.deallocate();
         winLum.deallocate();
         invHi.deallocate();
         etai.deallocate();
         coeffs.deallocate();
      }
   }

   return linearCoefficients;
}

/*
 * Upsample by linear coefficients defined in Eq. (17) in the appendix
 * of the article.
 */
Mat TMOShan10::upsampleByLinearCoefficients(Mat lum,
                                            Mat lumo,
                                            Mat bLumo,
                                            Mat map,                                            
                                            uint8_t multiGridFilt,
                                            int winSize,
                                            double epsilon)
{
   
   Mat coeffs = getLinearCoefficients(lum, lumo, map, multiGridFilt, winSize, epsilon);   
   Mat bCoeff(bLumo.rows, bLumo.cols, CV_64FC2);

   vector<Mat> coeffsChannels(2);
   split(coeffs, coeffsChannels);

   vector<Mat> bCoeffChannels(2);
   split(bCoeff, bCoeffChannels);
   
   matUtil.matlabResize(coeffsChannels[0], bCoeffChannels[0], bLumo.rows, bLumo.cols, BILINEAR);
   matUtil.matlabResize(coeffsChannels[1], bCoeffChannels[1], bLumo.rows, bLumo.cols, BILINEAR);
   
   return bCoeffChannels[0].mul(bLumo) + bCoeffChannels[1];   
}

/*
 * Creates matrices S and B defined in the article.
 */
std::pair<Mat, SpMat> TMOShan10::getLaplacian(Mat lumo,                                              
                                              Mat map,
                                              uint8_t multiGridFilt,
                                              int winSize,
                                              double epsilon)
{
   int height, width;
   int xStart, xEnd;
   int yStart, yEnd;
   int m;   

   height = lumo.rows;
   width = lumo.cols;

   int halfWinSize = (winSize-1) / 2;

   Mat indXM(height, width, CV_64FC1);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         indXM.at<double>(y, x) = x*height+y;
      }
   }   

   Mat B = Mat::zeros(height*width, 1, CV_64FC1);

   int cvCount = 0;   
   int vCount = pow((double)winSize, 4.0)*height*width;   

   Mat vy = Mat::zeros(vCount, 1, CV_64FC1);
   Mat vx = Mat::zeros(vCount, 1, CV_64FC1);
   Mat vv = Mat::zeros(vCount, 1, CV_64FC1);
      
   Mat sub;
   Mat addMat;

   Mat powMat(1, 1, CV_64FC1);

   for (int iWidth = 0; iWidth < width; iWidth++)
   {
      xStart = max(0, iWidth-halfWinSize);
      xEnd = min(width-1, iWidth+halfWinSize);
      for (int iHeight = 0; iHeight < height; iHeight++)
      {     
         yStart = max(0, iHeight-halfWinSize);
         yEnd = min(height-1, iHeight+halfWinSize);

         m = (yEnd - yStart + 1)*(xEnd - xStart + 1);        

         Mat winLumo = lumo(Range(yStart, yEnd+1), Range(xStart, xEnd+1));
         Mat indX = indXM(Range(yStart, yEnd+1), Range(xStart, xEnd+1));

         int eCount = (xEnd - xStart + 1)*(yEnd - yStart + 1);

         // WinLum
         Mat tWinLum(winLumo.rows * winLumo.cols, 1, CV_64FC1);

         for (int y = 0; y < winLumo.rows; y++)
         {
            for (int x = 0; x < winLumo.cols; x++)
            {
               tWinLum.at<double>(x * winLumo.rows + y, 0) = winLumo.at<double>(y, x);               
            }
         }         

         Mat tWinLum1;
         Mat tWinLum2;
         
         repeat(tWinLum, 1, eCount, tWinLum1);
         repeat(tWinLum.t(), eCount, 1, tWinLum2);

         // IndX
         Mat tIndX(indX.rows * indX.cols, 1, CV_64FC1);

         for (int x = 0; x < indX.cols; x++)
         {
            for (int y = 0; y < indX.rows; y++)
            {
               tIndX.at<double>(x * indX.rows + y, 0) = indX.at<double>(y, x);               
            }
         }

         Mat tIndX1;
         Mat tIndX2;

         repeat(tIndX, 1, eCount, tIndX1);
         repeat(tIndX.t(), eCount, 1, tIndX2);

         double tMap = map.at<double>(iHeight, iWidth);

         double mean = 0.0;

         // Standard deviation
         for (int i = 0; i < tWinLum.rows; i++)
         {
            // Sum of all elements          
            mean += tWinLum.at<double>(i, 0);       
         }

         // Calculate mean
         mean /= (double)(tWinLum.rows);

         double sum = 0.0;

         for (int i = 0; i < tWinLum.rows; i++)
         {            
            sum += pow(tWinLum.at<double>(i,0) - mean, 2.0);
         }

         // The standard deviation is normalized by the number of observations.         
         double tVar = sum / (tWinLum.rows);

         double mul = pow(tMap, (double)(-multiGridFilt));

         double tDelta = tVar + (epsilon/(double)m) * mul;

         Mat tVal = Mat::zeros(tIndX1.rows, tIndX2.rows, CV_64FC1);

         for (int y = 0; y < tVal.rows; y++)
         {
            for (int x = 0; x < tVal.cols; x++)
            {
               tVal.at<double>(y, x) = (tIndX1.at<double>(y, x) == tIndX2.at<double>(y, x));         
            }
         }         

         sub = (tWinLum1 - mean).mul(tWinLum2 - mean);       
         sub /= tDelta;
         sub += 1.0;
         sub /= (double)m;

         tVal -= sub;         

         double vIncrease = tVal.rows * tVal.cols;

         int position;
         int xCoordinate;
         int yCoordinate;


         // vy
         for (int y = 0; y < tIndX1.rows; y++)
         {
            for (int x = 0; x < tIndX1.cols; x++)
            {
               vy.at<double>(cvCount + x * tIndX1.rows + y, 0) = tIndX1.at<double>(y, x);
            }
         }

         // vx
         for (int y = 0; y < tIndX2.rows; y++)
         {
            for (int x = 0; x < tIndX2.cols; x++)
            {
               vx.at<double>(cvCount + x * tIndX2.rows + y, 0) = tIndX2.at<double>(y, x);
            }
         }

         // vv
         for (int y = 0; y < tVal.rows; y++)
         {
            for (int x = 0; x < tVal.cols; x++)
            {
               vv.at<double>(cvCount + x * tVal.rows + y, 0) = tVal.at<double>(y, x);       
            }
         }
      
         cvCount += vIncrease;      

         addMat = tWinLum - mean;      
         addMat *= (epsilon / ((double)m * tDelta * pow(tMap, multiGridFilt - 1.0)));         

         for (int i = 0; i < tIndX.rows; i++)
         {
            position = tIndX.at<double>(i, 0);
            B.at<double>(position, 0) +=  addMat.at<double>(i, 0);
         }

         indX.deallocate();       
         winLumo.deallocate();
         tWinLum.deallocate();
         tWinLum1.deallocate();
         tWinLum2.deallocate();
         tIndX.deallocate();
         tIndX1.deallocate();
         tIndX2.deallocate();
         tVal.deallocate();         
         sub.deallocate();                  
         addMat.deallocate();       
      }
   }

   Mat tmpvy(cvCount, 1, CV_64FC1);
   Mat tmpvx(cvCount, 1, CV_64FC1);
   Mat tmpvv(cvCount, 1, CV_64FC1);

   for (int i = 0; i < cvCount; i++)
   {
      tmpvy.at<double>(i, 0) = vy.at<double>(i, 0);
      tmpvx.at<double>(i, 0) = vx.at<double>(i, 0);
      tmpvv.at<double>(i, 0) = vv.at<double>(i, 0);
   }

   powMat.deallocate();
   indXM.deallocate();
   vy.deallocate();
   vx.deallocate();
   vv.deallocate();

   int matSize = height * width;

   SpMat S(matSize, matSize);

   for (int i = 0; i < cvCount; i++)
   {
      S.coeffRef(tmpvy.at<double>(i, 0), tmpvx.at<double>(i, 0)) += tmpvv.at<double>(i, 0); 
   }      

  return std::make_pair(B, S);
}

/*
 * Method has two important parts:
 * - creating the matrices S and B
 * - solving the linear system
 */
Mat TMOShan10::solveLinearSystemFire(Mat lumo,
                                     Mat constsMap,
                                     Mat constsValue,
                                     Mat map,
                                     uint8_t multiGridFilt,
                                     int winSize,
                                     double epsilon,
                                     double lambda)
{
   int height, width;

   height = lumo.rows;
   width = lumo.cols;

   // Create linear system - Article appendix equasion 17
   std::pair<Mat, SpMat> laplacian = getLaplacian(lumo, map, multiGridFilt, winSize, epsilon);

   int rLength = height * width;

   Mat cConstsMap(constsMap.rows * constsMap.cols, 1, CV_64FC1);
   
   for (int y = 0; y < constsMap.rows; y++)
   {
      for (int x = 0; x < constsMap.cols; x++)
      {
         cConstsMap.at<double>(x*constsMap.rows + y, 0) = constsMap.at<double>(y, x);
      }
   }
   
   Mat cConstsValue(constsValue.rows * constsValue.cols, 1, CV_64FC1);

   for (int y = 0; y < constsValue.rows; y++)
   {
      for (int x = 0; x < constsValue.cols; x++)
      {
         cConstsValue.at<double>(x*constsValue.rows + y, 0) = constsValue.at<double>(y, x);
      }
   }

   SpMat diag(rLength, rLength);

   for (int i = 0; i < rLength; i++)
   {
      diag.insert(i, i) = cConstsMap.at<double>(i, 0) * lambda;
   }

   SpMat sp = laplacian.second + diag;

   Mat bp = cConstsValue.mul(cConstsMap);   
   bp *= lambda;
   bp += laplacian.first;

   laplacian.first.deallocate();
   cConstsMap.deallocate();
   cConstsValue.deallocate();

   // Solve linear system
   VectorXd bpVec(bp.rows);   

   for (int i = 0; i < bp.rows; i++)
   {
      bpVec[i] = bp.at<double>(i, 0);      
   }

   SimplicialLDLT<SpMat> solver(sp);   

   if(solver.info() != Success) {
      // Decomposition failed.
      cout << "\nError: Decomposition failed ...\n";
   }
   
   VectorXd solverRes = solver.solve(bpVec);

   if(solver.info() != Success) {
      // Decomposition failed.
      cout << "\nError: Decomposition failed ...\n";
   }

   Mat result(height, width, CV_64FC1);

   for (int y = 0; y < height; y++)   
   {
      for (int x = 0; x < width; x++)
      {
         result.at<double>(y, x) = solverRes[x * height + y];
      }
   }
   
   bp.deallocate();   

   return result;
}

/*
 * A multigrid architecture that solves the large linear system.
 */
Mat TMOShan10::solveLinearSystem(Mat lumo,
                                 Mat constsMap,
                                 Mat constsValue,
                                 Mat map,                                 
                                 uint8_t levelNum,
                                 int winSize,
                                 double epsilon,
                                 double lambda)
{
   const uint8_t multiGridFilter = 2;

   assert(levelNum >= 1);

   if (levelNum > 1)
   {
      Mat sLumo((lumo.rows + 1) / 2, (lumo.cols + 1) / 2, CV_64FC1);
      Mat sConstsMap((constsMap.rows + 1) / 2, (constsMap.cols + 1) / 2, CV_64FC1);
      Mat sConstsValue((constsValue.rows + 1) / 2, (constsValue.cols + 1)  / 2, CV_64FC1);
      Mat sMap((map.rows + 1) / 2, (map.cols + 1) / 2, CV_64FC1);      

      matUtil.matlabResize(lumo, sLumo, sLumo.rows, sLumo.cols, BICUBIC);
      matUtil.matlabResize(constsMap, sConstsMap, sConstsMap.rows, sConstsMap.cols, BICUBIC);
      matUtil.matlabResize(constsValue, sConstsValue, sConstsValue.rows, sConstsValue.cols, BICUBIC);
      matUtil.matlabResize(map, sMap, sMap.rows, sMap.cols, BICUBIC);      

      Mat sLum = solveLinearSystem(sLumo,
                                   sConstsMap,
                                   sConstsValue,
                                   sMap,                                   
                                   levelNum - 1,
                                   winSize,
                                   epsilon,
                                   lambda);

      return upsampleByLinearCoefficients(sLum,
                                          sLumo,
                                          lumo,
                                          sMap,                                          
                                          multiGridFilter,
                                          winSize,
                                          epsilon);
   }
   else
   {
      // levelNum is 1 
      for (int y = 0; y < constsMap.rows; y++)
      {
         for (int x = 0; x < constsMap.cols; x++)
         {
            if (constsMap.at<double>(y, x) > 0.87)
            {
               constsMap.at<double>(y, x) = 1.0;
            }
            else
            {
               constsMap.at<double>(y, x) = 0.0;
            }            
         }         
      }

      return solveLinearSystemFire(lumo,
                                   constsMap,
                                   constsValue,
                                   map,                                   
                                   multiGridFilter,
                                   winSize,
                                   epsilon,
                                   lambda);      
   }
}

/*
 * Main method for TMOShan10 operator.
 */
int TMOShan10::Transform()
{
   int height, width;
   int pWinSize;
   int pLevelNum;
   double pBeta1, pBeta2, pBeta3;
   double pSSat;
   
   const double kappa = 0.5;   
   const double epsilon = 0.1;
   const double lambda = 0.0;

   /***********************/
   /* Get user parameters */
   /***********************/

   pBeta1 = beta1.GetDouble();
   pBeta2 = beta2.GetDouble();
   pBeta3 = beta3.GetDouble();
   pSSat = sSat.GetDouble();         
   pWinSize = winSize.GetInt();
   pLevelNum = levelNum.GetInt();

   Mat srcMat = TMOImage2Mat(pSrc);

   height = srcMat.rows;
   width = srcMat.cols;   

   Mat tLum(height, width, CV_64FC1);
   Mat lumo(height, width, CV_64FC1);

   matUtil.setValueHSV(srcMat, tLum);

   Scalar tLumMean =cv::mean(tLum)[0];  

   srcMat *= (200.0 / tLumMean.val[0]);

   matUtil.setValueHSV(srcMat, lumo);
   matUtil.checkLumo(lumo);

   Scalar tmpmean =cv::mean(lumo)[0];    

   Mat lumoMap = lumo.clone();
   Mat map = generateGuidanceMap(lumoMap, pWinSize, pBeta1, pBeta2, pBeta3, kappa); 

   Mat constsMap = Mat::zeros(height, width, CV_64FC1);   
   Mat constsValue = Mat::zeros(height, width, CV_64FC1);   

   Mat lum = solveLinearSystem(lumo,
                               constsMap,
                               constsValue,
                               map,                               
                               pLevelNum,
                               pWinSize,
                               epsilon,
                               lambda);

   matUtil.normImage(lum);
   

   // RGB channels reconstruction   
   Mat retRgb(height, width, CV_64FC3);
   Mat lumoThreeChannels(height, width, CV_64FC3);
   Mat lumThreeChannels(height, width, CV_64FC3);   

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         lumoThreeChannels.at<double>(y, x*CHANNELSCNT+2) = lumo.at<double>(y, x) + 1e-9;
         lumoThreeChannels.at<double>(y, x*CHANNELSCNT+1) = lumo.at<double>(y, x) + 1e-9;
         lumoThreeChannels.at<double>(y, x*CHANNELSCNT)   = lumo.at<double>(y, x) + 1e-9;
      }      
   }
   
   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++)
      {
         lumThreeChannels.at<double>(y, x*CHANNELSCNT+2) = lum.at<double>(y, x) + 1e-9;
         lumThreeChannels.at<double>(y, x*CHANNELSCNT+1) = lum.at<double>(y, x) + 1e-9;
         lumThreeChannels.at<double>(y, x*CHANNELSCNT)   = lum.at<double>(y, x) + 1e-9;
      }      
   }

   Mat rgbRatio = srcMat / (lumoThreeChannels);
   pow(rgbRatio, pSSat, retRgb);

   retRgb = retRgb.mul(lumThreeChannels);

   Mat result = matUtil.cutPeakValley(retRgb, 0.5, 0.5);

   /***********************/
   /*     Save result     */
   /***********************/	

	double *pSourceData = pSrc->GetData();
	double *pDestinationData = pDst->GetData();

	int y = 0;
	for (y = 0; y < pSrc->GetHeight(); y++)
	{
		pSrc->ProgressBar(y, pSrc->GetHeight());
		for (int x = 0; x < pSrc->GetWidth(); x++)
		{
         *pDestinationData++ = result.at<double>(y, x*CHANNELSCNT+2);
         *pDestinationData++ = result.at<double>(y, x*CHANNELSCNT+1);
         *pDestinationData++ = result.at<double>(y, x*CHANNELSCNT);         
		}
	}
	
	pSrc->ProgressBar(y, pSrc->GetHeight());   

	return 0;
}
