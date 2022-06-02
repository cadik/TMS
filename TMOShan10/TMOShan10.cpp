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

   beta2.SetName(L"beta2");
   beta2.SetDescription(L"Beta2 parameter. Article Eq. (6).");   
   beta2.SetDefault(0.2);
   beta2 = 0.2;
   beta2.SetRange(0.1, 0.4);
   this->Register(beta2);

   beta1.SetName(L"beta1");
   beta1.SetDescription(L"Beta1 parameter. Article Eq. (6).");
   beta1.SetDefault(0.4);
   beta1 = 0.4;
   beta1.SetRange(0.4, 0.9);
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
Mat TMOShan10::generateGuidanceMap(Mat &lumo,
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
   Mat tLum(height, width, lumo.type());

   filter2D(lumo, meanLum, -1, meanFilter, Point(-1, -1), 0, BORDER_REPLICATE);

   tLum = matUtil.padArray(lumo, halfWinSize, halfWinSize, PADARR_REPLICATE);

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
Mat TMOShan10::getLinearCoefficients(Mat &lum,
                                     Mat &lumo,
                                     Mat &map,
                                     Mat &epsilonMap,
                                     uint8_t multiGridFilt,
                                     int winSize)
{
   int height, width;
   int halfWinSize;
   int xStart, xEnd;
   int yStart, yEnd;
   int m;
   double epsilon;

   height = lumo.rows;
   width = lumo.cols;

   halfWinSize = (winSize-1) / 2;

   Mat linearCoefficients = Mat::zeros(height, width, CV_64FC2);
   Mat winLumo;
   Mat winLum;
   Mat invHi;
   Mat etai;
   Mat coeffs;

   for (int iWidth = 0; iWidth < width; iWidth++)
   {
      xStart = max(0, iWidth-halfWinSize);
      xEnd = min(width-1, iWidth+halfWinSize);

      for (int iHeight = 0; iHeight < height; iHeight++)
      {
         yStart = max(0, iHeight-halfWinSize);
         yEnd = min(height-1, iHeight+halfWinSize);

         m = (yEnd - yStart)*(xEnd - xStart);

         if (iHeight <= epsilonMap.rows && iWidth <= epsilonMap.cols)
         {
            epsilon = epsilonMap.at<double>(0, 0);
         }
         else
         {
            epsilon = 0.1;
         }         
         
         winLumo = lumo(Range(yStart, yEnd), Range(xStart, xEnd));
         winLum = lum(Range(yStart, yEnd), Range(xStart, xEnd));

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
         double divider = tDelta * m;

         double matrixFst = 1.0/divider;
         double matrixSecThd =  (-mean)/divider;
         double matrixFourth = (tDelta + mean * mean)/divider;

         invHi = (Mat_<double>(2, 2) << matrixFst, matrixSecThd, matrixSecThd, matrixFourth);         

         // Source: https://stackoverflow.com/questions/21874774/sum-of-elements-in-a-matrix-in-opencv
         // Author: Shai

         // Source: https://stackoverflow.com/questions/43180053/initialize-the-values-into-mat-object-in-opencv
         // Author: Micka
         
         etai = (Mat_<double>(2, 1) << (epsilon / tMap + cv::sum(winLumo.mul(winLum))[0]), cv::sum(winLum)[0]);         

         // Mathematical multiplication
         coeffs = (invHi*etai)/0.8;

         linearCoefficients.at<Vec2d>(iHeight, iWidth) = Vec2d(coeffs.at<double>(0, 0),coeffs.at<double>(0, 1));         
         
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
Mat TMOShan10::upsampleByLinearCoefficients(Mat &lum,
                                            Mat &lumo,
                                            Mat &bLumo,
                                            Mat &map,
                                            Mat &epsilonMap,
                                            uint8_t multiGridFilt,
                                            int winSize)
{   
   Mat coeffs = getLinearCoefficients(lum, lumo, map, epsilonMap, multiGridFilt, winSize);   
   Mat bCoeff;   
      
   resize(coeffs, bCoeff, Size(bLumo.cols, bLumo.rows), 0, 0, INTER_LINEAR);
   
   vector<Mat> bCoeffChannels(2);
   split(bCoeff, bCoeffChannels);   

   coeffs.deallocate();
   bCoeff.deallocate();
   
   return bCoeffChannels[0].mul(bLumo) + bCoeffChannels[1];  
}

/*
 * Creates matrices S and B defined in the article.
 */
std::pair<Mat, Mat> TMOShan10::getLaplacian(Mat &lumo,
                                            Mat &epsilonMap,
                                            Mat &map,
                                            uint8_t multiGridFilt,
                                            int winSize)
{
   int height, width;
   int xStart, xEnd;
   int yStart, yEnd;
   int m;
   double epsilon;

   height = lumo.rows;
   width = lumo.cols;

   int halfWinSize = (winSize-1) / 2;

   Mat indXM(height, width, CV_64FC1);

   for (int y = 0; y < height; y++)
   {
      for (int x = 0; x < width; x++) {
         indXM.at<double>(y, x) = x*height+y;
      }
   }

   Mat B = Mat::zeros(height*width, 1, CV_64FC1);


   int cvCount = 0;   
   int vCount = pow((double)winSize, 4.0)*height*width;   

   Mat vy = Mat::zeros(vCount, 1, CV_64FC1);
   Mat vx = Mat::zeros(vCount, 1, CV_64FC1);
   Mat vv = Mat::zeros(vCount, 1, CV_64FC1);
   
   Mat winLumo;
   Mat indX;
   Mat tVal;
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

         m = (yEnd - yStart)*(xEnd - xStart);

         if (iHeight <= epsilonMap.rows && iWidth <= epsilonMap.cols)
         {
            epsilon = epsilonMap.at<double>(0, 0);
         }
         else
         {
            epsilon = 0.1;
         }

         winLumo = lumo(Range(yStart, yEnd), Range(xStart, xEnd));
         indX = indXM(Range(yStart, yEnd), Range(xStart, xEnd));

         int eCount = (xEnd - xStart)*(yEnd-yStart);                  

         // WinLum
         Mat tWinLum(winLumo.rows * winLumo.cols, 1, CV_64FC1);

         for (int y = 0; y < winLumo.rows; y++)
         {
            for (int x = 0; x < winLumo.cols; x++)
            {
               tWinLum.at<double>(x * winLumo.rows + y, 0) = winLumo.at<double>(y, x);
            }
         }

         Mat tWinLum1(tWinLum.rows, tWinLum.cols * eCount, CV_64FC1);
         Mat tWinLum2(tWinLum.cols * eCount, tWinLum.rows, CV_64FC1);
         
         repeat(tWinLum, 1, eCount, tWinLum1);
         repeat(tWinLum.t(), eCount, 1, tWinLum2);

         // IndX
         Mat tIndX(indX.rows * indX.cols, 1, CV_64FC1);

         for (int y = 0; y < indX.rows; y++)
         {
            for (int x = 0; x < indX.cols; x++)
            {
               tIndX.at<double>(x * indX.rows + y, 0) = indX.at<double>(y, x);
            }
         }

         Mat tIndX1(tIndX.rows, tIndX.cols * eCount, CV_64FC1);
         Mat tIndX2(tIndX.cols * eCount, tIndX.rows, CV_64FC1);
         
         repeat(tIndX, 1, eCount, tIndX1);
         repeat(tIndX.t(), eCount, 1, tIndX2);

         double tMap = map.at<double>(iHeight, iWidth);

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

         double mul = pow(tMap, (double)(-multiGridFilt));
         double tDelta = tVar + (epsilon/(double)m) * mul;

         sum = 0.0;
         
         for (int i = 0; i < tWinLum.rows; i++)
         {
            sum += tWinLum.at<double>(i, 0);            
         }
         
         double tMean = sum / tWinLum.rows;

         tVal = Mat::zeros(tIndX1.rows, tIndX2.rows, CV_64FC1);

         for (int y = 0; y < tVal.rows; y++)
         {
            for (int x = 0; x < tVal.cols; x++)
            {
               tVal.at<double>(y, x) = (tIndX1.at<double>(y, x) == tIndX2.at<double>(y, x));
            }
         }
         
         Mat mul1 = tWinLum1 - tMean;
         Mat mul2 = tWinLum2 - tMean;
         
         sub = (tWinLum1 - tMean).mul(tWinLum2 - tMean);
         sub /= tDelta;
         sub += 1.0;
         sub /= m;

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
               position = x * tIndX1.rows + y;               
               xCoordinate = (cvCount + position) / vy.rows;
               yCoordinate = (cvCount + position) % vy.rows;
               vy.at<double>(xCoordinate, yCoordinate) = tIndX1.at<double>(y, x);
            }
         }
         
         // vx
         for (int y = 0; y < tIndX2.rows; y++)
         {
            for (int x = 0; x < tIndX2.cols; x++)
            {
               position = x * tIndX2.rows + y;               
               xCoordinate = (cvCount + position) / vx.rows;
               yCoordinate = (cvCount + position) % vx.rows;
               vx.at<double>(xCoordinate, yCoordinate) = tIndX2.at<double>(y, x);
            }
         }

         // vv
         for (int y = 0; y < tVal.rows; y++)
         {
            for (int x = 0; x < tVal.cols; x++)
            {
               position = x * tVal.rows + y;               
               xCoordinate = (cvCount + position) / vv.rows;
               yCoordinate = (cvCount + position) % vv.rows;
               vv.at<double>(xCoordinate, yCoordinate) = tVal.at<double>(y, x);
            }
         }
      
         cvCount += vIncrease;         
         
         pow(tMap, multiGridFilt - 1.0, powMat);         
         addMat = tWinLum - tMean;
         addMat *= (epsilon / m * tDelta * powMat);

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

   Mat tmpvy = vy(Range(0, cvCount), Range(0, 1));
   Mat tmpvx = vx(Range(0, cvCount), Range(0, 1));
   Mat tmpvv = vv(Range(0, cvCount), Range(0, 1));

   powMat.deallocate();
   indXM.deallocate();          
   vy.deallocate();
   vx.deallocate();
   vv.deallocate();

   int matSize = height * width;

   Mat S = Mat::zeros(matSize, matSize, CV_64FC1);   

   for (int i = 0; i < cvCount; i++)
   {
      S.at<double>(tmpvy.at<double>(i, 0), tmpvx.at<double>(i, 0)) = tmpvv.at<double>(i, 0);      
   }
   
   return std::make_pair(B, S);  
}

/*
 * Method has two important parts:
 * - creating the matrices S and B
 * - solving the linear system
 */
Mat TMOShan10::solveLinearSystemFire(Mat &lumo,
                                     Mat &constsMap,
                                     Mat &constsValue,
                                     Mat &map,
                                     Mat &epsilonMap,
                                     uint8_t multiGridFilt,
                                     int winSize,
                                     double lambda)
{
   int height, width;

   height = lumo.rows;
   width = lumo.cols;

   // Create linear system - Article appendix equasion 17
   std::pair<Mat, Mat> laplacian = getLaplacian(lumo, epsilonMap, map, multiGridFilt, winSize);

   int rLength = height * width;

   Mat cConstsMap = Mat::zeros(constsMap.rows * constsMap.cols, 1, CV_64FC1);
   
   for (int y = 0; y < constsMap.rows; y++)
   {
      for (int x = 0; x < constsMap.cols; x++)
      {
         cConstsMap.at<double>(x*constsMap.rows + y, 0) = constsMap.at<double>(y, x);
      }   
   }   
   
   Mat cConstsValue = Mat::zeros(constsValue.rows * constsValue.cols, 1, CV_64FC1);

   for (int y = 0; y < constsValue.rows; y++)
   {
      for (int x = 0; x < constsValue.cols; x++)
      {
         cConstsValue.at<double>(x*constsValue.rows + y, 0) = constsValue.at<double>(y, x);
      }
   }
   
   Mat diag = Mat::zeros(rLength, rLength, CV_64FC1);   

   for (int i = 0; i < rLength; i++)
   {
      diag.at<double>(i, i) = cConstsMap.at<double>(i, 0);
   }
   
   diag *= lambda;

   Mat sp = laplacian.second + diag;
   laplacian.second.deallocate();
   diag.deallocate();

   Mat bp = cConstsValue.mul(cConstsMap);   
   bp *= lambda;   
   bp += laplacian.first;   

   laplacian.first.deallocate();
   cConstsMap.deallocate();
   cConstsValue.deallocate();

   // Solve linear system
   Mat resLum(sp.rows, sp.cols, CV_64FC1);
   solve(sp, bp, resLum, DECOMP_CHOLESKY);   
         
   sp.deallocate();
   bp.deallocate();   

   return resLum.reshape(0, height);     
}

/*
 * A multigrid architecture that solves the large linear system.
 */
Mat TMOShan10::solveLinearSystem(Mat &lumo,
                                 Mat &constsMap,
                                 Mat &constsValue,
                                 Mat &map,
                                 Mat &epsilonMap,
                                 uint8_t levelNum,
                                 int winSize,
                                 double lambda)
{
   const uint8_t multiGridFilter = 2;

   assert(levelNum >= 1);

   if (levelNum > 1)
   {
      // Divide by 2
      Mat sLumo(lumo.rows * 0.5, lumo.cols * 0.5, CV_64FC1);
      Mat sConstsMap(constsMap.rows * 0.5, constsMap.cols * 0.5, CV_64FC1);
      Mat sConstsValue(constsValue.rows * 0.5, constsValue.cols  * 0.5, CV_64FC1);
      Mat sMap(map.rows * 0.5, map.cols * 0.5, CV_64FC1);
      Mat sEpsilonMap;

      resize(lumo, sLumo, sLumo.size(), 0.5, 0.5);
      resize(constsMap, sConstsMap, sConstsMap.size(), 0.5, 0.5);
      resize(constsValue, sConstsValue, sConstsValue.size(), 0.5, 0.5);
      resize(map, sMap, sMap.size(), 0.5, 0.5);

      if (epsilonMap.rows > 1 && epsilonMap.cols > 1)
      {
         sEpsilonMap = Mat(epsilonMap.rows * 0.5, epsilonMap.cols * 0.5, CV_64FC1);
         resize(epsilonMap, sEpsilonMap, sEpsilonMap.size(), 0.5, 0.5);
      }
      else
      {
         sEpsilonMap = epsilonMap;
      }            

      Mat sLum = solveLinearSystem(sLumo,
                                   sConstsMap,
                                   sConstsValue,
                                   sMap,
                                   sEpsilonMap,
                                   levelNum - 1,
                                   winSize,
                                   lambda);      

      return upsampleByLinearCoefficients(sLum,
                                          sLumo,
                                          lumo,
                                          sMap,
                                          epsilonMap,
                                          multiGridFilter,
                                          winSize);      
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
                                   epsilonMap,
                                   multiGridFilter,
                                   winSize,
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
   double pBeta1, pBeta2;
   double pSSat;

   const uint8_t levelNum = 3;
   const double kappa = 0.5;
   const double beta3 = 0.1;
   const double epsilon = 0.1;
   const double lambda = 0.0;

   /***********************/
   /* Get user parameters */
   /***********************/

   pBeta1 = beta1.GetDouble();
   pBeta2 = beta2.GetDouble();
   pSSat = sSat.GetDouble();         
   pWinSize = winSize.GetInt();   

   Mat srcMat = TMOImage2Mat(pSrc);

   height = srcMat.rows;
   width = srcMat.cols;

   Mat tLum(height, width, CV_64FC1);
   Mat lumo(height, width, CV_64FC1);

   matUtil.setValueHSV(srcMat, tLum);

   Scalar tLumMean = mean(mean(mean(tLum)));

   // Because tLum has only one channel, the tLumMean scalar
   // has only value at first element position -> val[0]
   srcMat = srcMat * (200.0 / tLumMean.val[0]);

   matUtil.setValueHSV(srcMat, lumo);
   matUtil.checkLumo(lumo);

   Mat map = generateGuidanceMap(lumo, pWinSize, pBeta1, pBeta2, beta3, kappa);

   Mat constsMap = Mat::zeros(height, width, CV_64FC1);
   Mat constsValue = Mat::zeros(height, width, CV_64FC1);
   Mat epsilonMap(1, 1, CV_64FC1, epsilon);
  
   Mat lum = solveLinearSystem(lumo,
                               constsMap,
                               constsValue,
                               map,
                               epsilonMap,
                               levelNum,
                               pWinSize,
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

	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight());
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{			
         *pDestinationData++ = result.at<double>(j, i*CHANNELSCNT+2);
         *pDestinationData++ = result.at<double>(j, i*CHANNELSCNT+1);
         *pDestinationData++ = result.at<double>(j, i*CHANNELSCNT);
		}
	}
	
	pSrc->ProgressBar(j, pSrc->GetHeight());

	return 0;
}
