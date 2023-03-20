/* --------------------------------------------------------------------------- *
 * TMOChen05.cpp: implementation of the TMOChen05 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOChen05.h"
#include "wasserstein.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOChen05::TMOChen05()
{
	SetName(L"Chen05");					  // TODO - Insert operator name
	SetDescription(L"Tone reproduction: a perspective from luminance-driven perceptual grouping"); // TODO - Insert description

	Theta.SetName(L"Theta");				// TODO - Insert parameters names
	Theta.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	Theta.SetDefault(1.5);							// TODO - Add default values
	Theta = 1.5;
	Theta.SetRange(1.5, 2.0); // TODO - Add acceptable range if needed
	this->Register(Theta);

   Delta.SetName(L"Delta");				// TODO - Insert parameters names
	Delta.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	Delta.SetDefault(0.5);							// TODO - Add default values
	Delta = 0.5;
	Delta.SetRange(0.5, 1.0); // TODO - Add acceptable range if needed
	this->Register(Delta);
}

#define rgb2luminance(R,G,B) (R*0.2126 + G*0.7152 + B*0.0722)
unsigned int IMAGE_HEIGHT;
unsigned int IMAGE_WIDTH;

typedef vector< vector<float> > PixelMatrix;

typedef vector< vector<double> > PixelDoubleMatrix;

typedef std::vector< vector<int> > PixelIntMatrix;

typedef vector< vector<long> > PixelLongMatrix;

struct Point{
   unsigned short x,y;
};
struct Signature{
   double s;
   double w;
};
typedef std::vector<Point> PointVector;
typedef std::vector<double> HistogramVector;
typedef std::vector<Signature> SignatureVector;
typedef std::vector<int> UnvisitedVector;
typedef std::vector<int> NeighbourContainer;
struct Block_Record{
   PointVector Memebers;
   NeighbourContainer Neighbours;
   float Sum;
   unsigned int Count;
   HistogramVector logHistogram;
   SignatureVector blockSignature;
};
struct Region_Record{
   PointVector Members;
   float Sum;
   unsigned int Count;
   HistogramVector logHistogram;
   SignatureVector regionSignature;
   vector<double> BoundryValues;
   PointVector GridMembers;
   map<double, int> BoundryMap;
   double alpha_r;
   double beta_r;
};

typedef std::vector<Block_Record > Blocks;
typedef std::vector<Region_Record > Regions;

float emdFunctionBB(int firstBlock, int secondBlock, Blocks& pixelBlocks)
{
   //cv::Mat sign1(cv::Size(2,3),CV_32FC1);
   //cv::Mat sign2(cv::Size(2,3),CV_32FC1);
   std::vector<double> av;
   std::vector<double> aw;
   std::vector<double> bv;
   std::vector<double> bw;
   for(int i=0;i < 3; i++)
   {
      //sign1.at<float>(i,0) = pixelBlocks[firstBlock].blockSignature[i].s;
      //sign1.at<float>(i,1) = pixelBlocks[firstBlock].blockSignature[i].w;
      av.push_back(pixelBlocks[firstBlock].blockSignature[i].s);
      aw.push_back(pixelBlocks[firstBlock].blockSignature[i].w);

   }
   for(int k=0; k < 3; k++)
   {
      //sign2.at<float>(k,0) = pixelBlocks[secondBlock].blockSignature[k].s;
      //sign2.at<float>(k,1) = pixelBlocks[secondBlock].blockSignature[k].w;
      bv.push_back(pixelBlocks[secondBlock].blockSignature[k].s);
      bw.push_back(pixelBlocks[secondBlock].blockSignature[k].w);
   } 
   //return cv::EMD(sign2,sign1,cv::DIST_L1);
   return wasserstein(av,aw,bv,bw);
   
}

float emdFunctionBR(int block, int region, Blocks& pixelBlocks, Regions& pixelRegions)
{
   //cv::Mat sign1(cv::Size(2,3),CV_32FC1);
   //cv::Mat sign2(cv::Size(2,3),CV_32FC1);
   std::vector<double> av;
   std::vector<double> aw;
   std::vector<double> bv;
   std::vector<double> bw;
   for(int i=0;i < 3; i++)
   {
      //sign1.at<float>(i,0) = pixelBlocks[block].blockSignature[i].s;
      //sign1.at<float>(i,1) = pixelBlocks[block].blockSignature[i].w;
      av.push_back(pixelBlocks[block].blockSignature[i].s);
      aw.push_back(pixelBlocks[block].blockSignature[i].w);
   }
   for(int k=0; k < 3; k++)
   {
      //sign2.at<float>(k,0) = pixelRegions[region].regionSignature[k].s;
      //sign2.at<float>(k,1) = pixelRegions[region].regionSignature[k].w;
      bv.push_back(pixelRegions[region].regionSignature[k].s);
      bw.push_back(pixelRegions[region].regionSignature[k].w);
   } 
   //return cv::EMD(sign1,sign2,cv::DIST_L1);
   return wasserstein(av,aw,bv,bw);
   
}

bool isValid(int x, int y, int category, PixelIntMatrix& pixels, PixelIntMatrix& pixelCategories)
{
   if(x<0 || x>= IMAGE_HEIGHT || y<0 || y>= IMAGE_WIDTH || pixelCategories[x][y] != category || pixels[x][y] == 1)
   {
      return false;
   }
   return true;
}


void GroupNeighbours(int x, int y, int group, PixelIntMatrix& pixels, PixelIntMatrix& pixelCategories, Blocks& pixelBlocks)
{
   vector<pair<int, int>> queue;
   pair<int, int> p(x,y);
   queue.push_back(p);

   pixels[x][y] = 1;


   while(queue.size() > 0)
   {
      pair<int,int> currPixel = queue[queue.size() - 1];
      queue.pop_back();

      int posX = currPixel.first;
      int posY = currPixel.second;
      if(isValid(posX+1, posY, group, pixels, pixelCategories))
      {
         pixels[posX+1][posY] = 1;
         p.first = posX+1;
         p.second = posY;
         queue.push_back(p);
      }
      if(!(isValid(posX+1, posY, group, pixels, pixelCategories)))
      {
         if(posX+1>=0 && posX+1< IMAGE_HEIGHT && posY>=0 && posY<IMAGE_WIDTH && pixels[posX+1][posY]!=1) 
         {
            if(!(std::find(pixelBlocks[pixelCategories[posX+1][posY]].Neighbours.begin(), pixelBlocks[pixelCategories[posX+1][posY]].Neighbours.end(), pixelCategories[posX][posY]) != pixelBlocks[pixelCategories[posX+1][posY]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX+1][posY]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(pixelBlocks[pixelCategories[posX][posY]].Neighbours.begin(), pixelBlocks[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX+1][posY]) != pixelBlocks[pixelCategories[posX][posY]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX+1][posY]);
            }
         }
      }

      if(isValid(posX-1, posY, group, pixels, pixelCategories))
      {
         pixels[posX-1][posY] = 1;
         p.first = posX-1;
         p.second = posY;
         queue.push_back(p);
      }
      if(!(isValid(posX-1, posY, group, pixels, pixelCategories)))
      {
         if(posX-1>=0 && posX-1< IMAGE_HEIGHT && posY>=0 && posY<IMAGE_WIDTH && pixels[posX-1][posY]!=1)
         {
            if(!(std::find(pixelBlocks[pixelCategories[posX-1][posY]].Neighbours.begin(), pixelBlocks[pixelCategories[posX-1][posY]].Neighbours.end(), pixelCategories[posX][posY]) != pixelBlocks[pixelCategories[posX-1][posY]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX-1][posY]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(pixelBlocks[pixelCategories[posX][posY]].Neighbours.begin(), pixelBlocks[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX-1][posY]) != pixelBlocks[pixelCategories[posX][posY]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX-1][posY]);
            }
         }
         
      }
      if(isValid(posX, posY+1, group, pixels, pixelCategories))
      {
         pixels[posX][posY+1] = 1;
         p.first = posX;
         p.second = posY+1;
         queue.push_back(p);
      }
      if(!(isValid(posX, posY+1, group, pixels, pixelCategories)))
      {
         if(posX>=0 && posX< IMAGE_HEIGHT && posY+1>=0 && posY+1<IMAGE_WIDTH && pixels[posX][posY+1]!=1)
         {
            if(!(std::find(pixelBlocks[pixelCategories[posX][posY+1]].Neighbours.begin(), pixelBlocks[pixelCategories[posX][posY+1]].Neighbours.end(), pixelCategories[posX][posY]) != pixelBlocks[pixelCategories[posX][posY+1]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX][posY+1]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(pixelBlocks[pixelCategories[posX][posY]].Neighbours.begin(), pixelBlocks[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX][posY+1]) != pixelBlocks[pixelCategories[posX][posY]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX][posY+1]);
            }
         }
         
      }
      if(isValid(posX, posY-1, group, pixels, pixelCategories))
      {
         pixels[posX][posY-1] = 1;
         p.first = posX;
         p.second = posY-1;
         queue.push_back(p);
      }
      if(!(isValid(posX, posY-1, group, pixels, pixelCategories)))
      {
         if(posX>=0 && posX< IMAGE_HEIGHT && posY-1>=0 && posY-1<IMAGE_WIDTH && pixels[posX][posY-1]!=1)
         {
            if(!(std::find(pixelBlocks[pixelCategories[posX][posY-1]].Neighbours.begin(), pixelBlocks[pixelCategories[posX][posY-1]].Neighbours.end(), pixelCategories[posX][posY]) != pixelBlocks[pixelCategories[posX][posY-1]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX][posY-1]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(pixelBlocks[pixelCategories[posX][posY]].Neighbours.begin(), pixelBlocks[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX][posY-1]) != pixelBlocks[pixelCategories[posX][posY]].Neighbours.end())){
               pixelBlocks[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX][posY-1]);
            }
         }
         
      }

   }
}

void updateRegionSignature(int regionID, Regions& region, PixelDoubleMatrix& LogLuminancePixels)
{
   double maxlum = LogLuminancePixels[region[regionID].Members[0].y][region[regionID].Members[0].x];
   double minlum = LogLuminancePixels[region[regionID].Members[0].y][region[regionID].Members[0].x];
   region[regionID].logHistogram.clear();
   region[regionID].regionSignature.clear();
   for(int j=0; j < region[regionID].Members.size();j++)
   {
      int x = region[regionID].Members[j].x;
      int y = region[regionID].Members[j].y;
      double tmp = LogLuminancePixels[y][x];
      if(tmp > maxlum)
      {
         maxlum = tmp;
      }
      if(tmp < minlum)
      {
         minlum = tmp;
      }
      region[regionID].logHistogram.push_back(tmp);
   }
   double step = abs((maxlum-minlum)/3.0);
   double first_threshold = minlum + step;
   double second_threshold = first_threshold + step;
   double firstSectionLum =0.0, secondSectionLum=0.0, thirdSectionLum=0.0;
   double firstSectionCount=0.0, secondSectionCount=0.0, thirdSectionCount=0.0;
   for(int k = 0; k < region[regionID].logHistogram.size(); k++)
   {
      if(region[regionID].logHistogram[k] <= first_threshold)
      {
         firstSectionLum += region[regionID].logHistogram[k];
         firstSectionCount += 1.0;
      }
      else if(region[regionID].logHistogram[k] > first_threshold && region[regionID].logHistogram[k] <= second_threshold)
      {
         secondSectionLum += region[regionID].logHistogram[k];
         secondSectionCount += 1.0;
      }
      else{
         thirdSectionLum += region[regionID].logHistogram[k];
         thirdSectionCount += 1.0;
      }
   }
   double sumOfPixels = firstSectionCount + secondSectionCount + thirdSectionCount;
   double s1, w1, s2, w2, s3, w3;
   if(thirdSectionCount == 0)
   {
      s1 = 0.0;
      w1 = 0.0;
   }
   else{
      s1 = thirdSectionLum/thirdSectionCount;
      w1 = thirdSectionCount/sumOfPixels;
   }
   
   Signature signatureOfRegion;
   signatureOfRegion.s = s1;
   signatureOfRegion.w = w1;
   region[regionID].regionSignature.push_back(signatureOfRegion);
   if(secondSectionCount == 0.0){
      s2 = 0.0;
      w2 = 0.0;
   }
   else{
      s2 = secondSectionLum/secondSectionCount;
      w2 = secondSectionCount/sumOfPixels;
   }
   
   signatureOfRegion.s = s2;
   signatureOfRegion.w = w2;
   region[regionID].regionSignature.push_back(signatureOfRegion);
   if(firstSectionCount == 0.0)
   {
      s3 = 0.0;
      w3 = 0.0;
   }
   else{
      s3 = firstSectionLum/firstSectionCount;
      w3 = firstSectionCount/sumOfPixels;
   }
   
   signatureOfRegion.s = s3;
   signatureOfRegion.w = w3;
   region[regionID].regionSignature.push_back(signatureOfRegion);
   //fprintf(stderr,"s1 %f sign1 %f\n",s1,region[regionID].regionSignature[0].s);
}



TMOChen05::~TMOChen05()
{
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOChen05::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information
   //pSrc->Convert(TMO_RGB);
			// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components

   double *imageData = pSrc->GetData();
   int imageHeight = pSrc->GetHeight();
   int imageWidth = pSrc->GetWidth();
   IMAGE_HEIGHT = imageHeight;
   IMAGE_WIDTH = imageWidth;
   double MinimalImageLuminance = 0., MaximalImageLuminance = 0., AverageImageLuminance = 0.;
   pSrc->GetMinMaxAvg(&MinimalImageLuminance, &MaximalImageLuminance, &AverageImageLuminance);
   fprintf(stderr,"Min: %g Max: %g\n",MinimalImageLuminance,MaximalImageLuminance);
   

   //cv::Mat image = cv::imread("/home/matthewlele/images/hdr_images/84y7.tif");
   double stonits = pSrc->GetStonits();
   PixelDoubleMatrix LogLuminancePixels(imageHeight, vector<double>(imageWidth, 0.0));
   double pixelR, pixelG, pixelB;
   cv::Mat image;
   cv::Mat finalImage;
   image = cv::Mat::zeros(imageHeight, imageWidth, CV_32F);
   //Getting log luminance of each pixel
   for(int j=0; j < imageHeight; j++)
   {
      for(int i=0; i<imageWidth; i++)
      {
         pixelR = *imageData++;
         pixelG = *imageData++;
         pixelB = *imageData++;
         LogLuminancePixels[j][i] = log(rgb2luminance(pixelR, pixelG, pixelB));
         image.at<float>(j, i) = log(rgb2luminance(pixelR, pixelG, pixelB));
         //image.at<cv::Vec3f>(j, i)[1] = pixelG;
         //image.at<cv::Vec3f>(j, i)[2] = pixelR;
      }
   }
   cv::normalize(image,image,0.0,1.0,cv::NORM_MINMAX, CV_32F);
   double imgMin, imgMax;
   cv::minMaxLoc(image, &imgMin, &imgMax);
   image.convertTo(finalImage, CV_8U, 255*(imgMax - imgMin));
   //cv::cvtColor(finalImage, finalImage, cv::COLOR_BGR2GRAY);
   cv::Mat contours;
   cv::Mat gray_img;
   cv::Mat testImage = cv::imread("/home/matthewlele/images/cadik/cadik-desk01.hdr");

   //cvtColor(image, gray_img, CV_RGB2GRAY);
   cv::Canny(finalImage, contours, 80, 240);
   cv::namedWindow("Canny");
   cv::imshow("Canny",contours);
   cv::waitKey(0);
   
   
   fprintf(stderr,"Image height: %d Image width: %d\n",imageHeight, imageWidth);
   Blocks pixelBlocks;
   int edgeDetected = 0;
   //creating Blocks of pixels from image
   for(int j=0; j < imageHeight-7; j+=8)
   {
      for(int i=0; i < imageWidth-7; i+=8)
      {
         Block_Record block;
         Point p;
         block.Sum = 0.0;
         block.Count = 0;
         edgeDetected = 0;
         for(int k =j; k < j+8; k++)
         {
            for(int l=i; l < i+8; l++)
            {
               if(contours.at<double>(k, l) != 0.0)
               {
                  edgeDetected = 1;
                  break;
               }
               p.x = l;
               p.y = k;
               block.Count += 1;
               block.Sum += LogLuminancePixels[k][l];
               block.Memebers.push_back(p);
            }
         }
         if(edgeDetected == 1)
         {
            int x_value = i;
            int y_value = j;
            for(int first_row = 0; first_row < 4; first_row++)
            {
               Block_Record small_block;
               Point point;
               small_block.Count = 0;
               small_block.Sum = 0.0;
               for(int k =y_value; k < y_value+2; k++)
               {
                  for(int l=x_value; l < x_value+2; l++)
                  {
                     point.x = l;
                     point.y = k;
                     small_block.Count += 1;
                     small_block.Sum += LogLuminancePixels[k][l];
                     small_block.Memebers.push_back(point);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
            x_value = i;
            y_value += 2;
            for(int second_row = 0; second_row < 4; second_row++)
            {
               Block_Record small_block;
               Point point;
               small_block.Count = 0;
               small_block.Sum = 0.0;
               for(int k =y_value; k < y_value+2; k++)
               {
                  for(int l=x_value; l < x_value+2; l++)
                  {
                     point.x = l;
                     point.y = k;
                     small_block.Count += 1;
                     small_block.Sum += LogLuminancePixels[k][l];
                     small_block.Memebers.push_back(point);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
            x_value = i;
            y_value += 2;
            for(int third_row = 0; third_row < 4; third_row++)
            {
               Block_Record small_block;
               Point point;
               small_block.Count = 0;
               small_block.Sum = 0.0;
               for(int k =y_value; k < y_value+2; k++)
               {
                  for(int l=x_value; l < x_value+2; l++)
                  {
                     point.x = l;
                     point.y = k;
                     small_block.Count += 1;
                     small_block.Sum += LogLuminancePixels[k][l];
                     small_block.Memebers.push_back(point);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
            x_value = i;
            y_value += 2;
            for(int fourth_row = 0; fourth_row < 4; fourth_row++)
            {
               Block_Record small_block;
               Point point;
               small_block.Count = 0;
               small_block.Sum = 0.0;
               for(int k =y_value; k < y_value+2; k++)
               {
                  for(int l=x_value; l < x_value+2; l++)
                  {
                     point.x = l;
                     point.y = k;
                     small_block.Count += 1;
                     small_block.Sum += LogLuminancePixels[k][l];
                     small_block.Memebers.push_back(point);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
         }
         else{
            pixelBlocks.push_back(block);
         }
         
      }
   }
   fprintf(stderr,"Amount of blocks %d\n",pixelBlocks.size());

   //Calculating signature for each pixel
   for(int i=0; i < pixelBlocks.size();i++)
   {
      double maxlum = LogLuminancePixels[pixelBlocks[i].Memebers[0].y][pixelBlocks[i].Memebers[0].x];
      double minlum = LogLuminancePixels[pixelBlocks[i].Memebers[0].y][pixelBlocks[i].Memebers[0].x];
      for(int j=0; j < pixelBlocks[i].Memebers.size();j++)
      {
         int x = pixelBlocks[i].Memebers[j].x;
         int y = pixelBlocks[i].Memebers[j].y;
         double tmp = LogLuminancePixels[y][x];
         if(tmp > maxlum)
         {
            maxlum = tmp;
         }
         if(tmp < minlum)
         {
            minlum = tmp;
         }
         pixelBlocks[i].logHistogram.push_back(tmp);
      }
      double step = abs((maxlum-minlum)/3.0);
      double first_threshold = minlum + step;
      double second_threshold = first_threshold + step;
      double firstSectionLum =0.0, secondSectionLum=0.0, thirdSectionLum=0.0;
      double firstSectionCount=0.0, secondSectionCount=0.0, thirdSectionCount=0.0;
      //fprintf(stderr, "min: %f first_t %f second_t %f max %f\n",minlum,first_threshold,second_threshold,maxlum);
      for(int k = 0; k < pixelBlocks[i].logHistogram.size(); k++)
      {
         if(pixelBlocks[i].logHistogram[k] <= first_threshold)
         {
            firstSectionLum += pixelBlocks[i].logHistogram[k];
            firstSectionCount += 1.0;
         }
         else if(pixelBlocks[i].logHistogram[k] > first_threshold && pixelBlocks[i].logHistogram[k] <= second_threshold)
         {
            secondSectionLum += pixelBlocks[i].logHistogram[k];
            secondSectionCount += 1.0;
         }
         else{
            thirdSectionLum += pixelBlocks[i].logHistogram[k];
            thirdSectionCount += 1.0;
         }
      }
      double sumOfPixels = firstSectionCount + secondSectionCount + thirdSectionCount;
      double s1, w1, s2, w2, s3, w3;
      if(thirdSectionCount == 0)
      {
         s1 = 0.0;
         w1 = 0.0;
      }
      else{
         s1 = thirdSectionLum/thirdSectionCount;
         w1 = thirdSectionCount/sumOfPixels;
      }
      
      Signature signatureOfBlock;
      signatureOfBlock.s = s1;
      signatureOfBlock.w = w1;
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);
      if(secondSectionCount == 0)
      {
         s2 = 0.0;
         w2 = 0.0;
      }
      else{
         s2 = secondSectionLum/secondSectionCount;
         w2 = secondSectionCount/sumOfPixels;
      }
      
      signatureOfBlock.s = s2;
      signatureOfBlock.w = w2;
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);
      if(firstSectionCount == 0.0)
      {
         s3 = 0.0;
         w3 = 0.0;
      }
      else{
         s3 = firstSectionLum/firstSectionCount;
         w3 = firstSectionCount/sumOfPixels;
      }
      signatureOfBlock.s = s3;
      signatureOfBlock.w = w3;
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);
      
      //fprintf(stderr,"s1 %f w1 %f  s2 %f w2 %f  s3 %f w3 %f\n",pixelBlocks[i].blockSignature[0].s,pixelBlocks[i].blockSignature[0].w,pixelBlocks[i].blockSignature[1].s,pixelBlocks[i].blockSignature[1].w,pixelBlocks[i].blockSignature[2].s,pixelBlocks[i].blockSignature[2].w);
      
      
   }

   PixelIntMatrix pixelsGrp(imageHeight, vector<int>(imageWidth,0));
   PixelIntMatrix visitedPixels(imageHeight, vector<int>(imageWidth,0));
   //calculating each pixel's group
   for(int i=0; i < pixelBlocks.size();i++)
   {
      for(int k=0; k < pixelBlocks[i].Memebers.size();k++)
      {
         int x = pixelBlocks[i].Memebers[k].x;
         int y = pixelBlocks[i].Memebers[k].y;
         pixelsGrp[y][x] = i;
      }
   }
   fprintf(stderr,"height: %d width: %d\n",imageHeight,imageWidth);
   //finding neighbouring blocks of each block
   for(int i=0; i < pixelBlocks.size();i++)
   {
      int x = pixelBlocks[i].Memebers[0].x;
      int y = pixelBlocks[i].Memebers[0].y;
      GroupNeighbours(y,x, pixelsGrp[y][x], visitedPixels, pixelsGrp, pixelBlocks);
      if(pixelBlocks[i].Neighbours.size()==0)
      {
         fprintf(stderr,"block %d with %d neighbours\n",i,pixelBlocks[i].Neighbours.size());
      }
   }
   
   
   UnvisitedVector UnvistitedBlocks(pixelBlocks.size(),0);
   PixelIntMatrix pixelsReg(imageHeight, vector<int>(imageWidth,0));
   Regions blocksRegions;
   
   double theta = 1.5; 
   double delta = 1.5;
   int unvisited = pixelBlocks.size();
   int regionID = 0;
   int counterTMP = 0;
   int chosedBlockID = 0;
   vector<int> queue;
   
   while(unvisited > 0)
   {
      queue.clear();
      double biggestS1 = -10.0;
      int brightestBlockID = 0;
      for(int i=0; i < pixelBlocks.size();i++)
      {
         if(UnvistitedBlocks[i]==0 && pixelBlocks[i].blockSignature[0].s > biggestS1)
         {
            biggestS1 = pixelBlocks[i].blockSignature[0].s;
            brightestBlockID = i;
         }
      }
      Region_Record region;
      blocksRegions.push_back(region);
      Signature region_Signature;
      region_Signature.s = pixelBlocks[brightestBlockID].blockSignature[0].s;
      region_Signature.w = pixelBlocks[brightestBlockID].blockSignature[0].w;
      blocksRegions[regionID].regionSignature.push_back(region_Signature);
      region_Signature.s = pixelBlocks[brightestBlockID].blockSignature[1].s;
      region_Signature.w = pixelBlocks[brightestBlockID].blockSignature[1].w;
      blocksRegions[regionID].regionSignature.push_back(region_Signature);
      region_Signature.s = pixelBlocks[brightestBlockID].blockSignature[2].s;
      region_Signature.w = pixelBlocks[brightestBlockID].blockSignature[2].w;
      blocksRegions[regionID].regionSignature.push_back(region_Signature);
      
      queue.push_back(brightestBlockID);
      
      while(queue.size() > 0)
      {
         double smallest  = 5.0;
         int tmpID=0;
         for(int l=0;l < queue.size();l++)
         {
            float tmp = emdFunctionBR(queue[l],regionID,pixelBlocks,blocksRegions);
            if(tmp < smallest )
            {
               smallest = tmp;
               tmpID = l;
            }
         }
         chosedBlockID = queue[tmpID];
         if(smallest < theta)
         {
            queue.erase(queue.begin() + tmpID);
            for(int mem=0;mem < pixelBlocks[chosedBlockID].Memebers.size();mem++)
            {
               blocksRegions[regionID].Members.push_back(pixelBlocks[chosedBlockID].Memebers[mem]);
               int x = pixelBlocks[chosedBlockID].Memebers[mem].x;
               int y = pixelBlocks[chosedBlockID].Memebers[mem].y;
               pixelsReg[y][x] = regionID;
            }
            counterTMP += 1;
            UnvistitedBlocks[chosedBlockID] = 1;
            unvisited -= 1;
            //TODO
            
            for(int n=0; n< pixelBlocks[chosedBlockID].Neighbours.size();n++)
            {
               if(UnvistitedBlocks[pixelBlocks[chosedBlockID].Neighbours[n]] == 0)
               {
                  
                  float tmpEMD = emdFunctionBB(chosedBlockID,pixelBlocks[chosedBlockID].Neighbours[n], pixelBlocks);
                  if(tmpEMD < delta)
                  {
                     if(!(std::find(queue.begin(),queue.end(),pixelBlocks[chosedBlockID].Neighbours[n])!= queue.end()))
                     {
                        queue.push_back(pixelBlocks[chosedBlockID].Neighbours[n]);
                     }
                     
                  }
               }
            }
            //fprintf(stderr,"Chosed block %d region %d iterations %d distance %f queue %d regionID : %d\n",chosedBlockID,blocksRegions[regionID].Members.size(),counterTMP, smallest, queue.size(),regionID);
            fprintf(stderr,"\rBlocks assigned to region %d/%d , regions created %d",counterTMP,pixelBlocks.size(),regionID);
            fflush(stdout);
            updateRegionSignature(regionID, blocksRegions, LogLuminancePixels);
         }
         else{
            
            break;
         }
      }
      regionID += 1;
   }
   fprintf(stderr,"\nCreating regions completed\n");
   /*PixelIntMatrix visitedRegionPixels(imageHeight, vector<int>(imageWidth,0)); 
   int pixelRegionCount =0;
   for(int m=0; m < blocksRegions.size();m++)
   {
      fprintf(stderr,"Region %d members %d \n",m,blocksRegions[m].Members.size());
      pixelRegionCount+= blocksRegions[m].Members.size();
      for(int i =0; i < blocksRegions[m].Members.size();i++)
      {
         int x = blocksRegions[m].Members[i].x;
         int y = blocksRegions[m].Members[i].y;
         visitedRegionPixels[y][x] = 1;
      }
      
      
   }*/
   
   PixelDoubleMatrix localAdaptationPixels(imageHeight, vector<double>(imageWidth, 0.0));
   double sigma_r = 0.4;
   double sigma_rr = 0.5*0.4;
   double sigma_s = (imageHeight*imageWidth)*0.04;
   //double mask = sqrt((imageHeight*imageWidth)*0.04);
   double mask = 70;
   int bilateralIteration = 0;
   for(int i=0; i < blocksRegions.size();i++)
   {
      for(int k=0; k < blocksRegions[i].Members.size();k++)
      {
         fprintf(stderr,"\rBillateral filter iterations %d/%d ",bilateralIteration,imageHeight*imageWidth);
         fflush(stdout);
         bilateralIteration++;
         double inRegion = 0.0;
         double inRegionZ = 0.0;
         double otherRegions = 0.0;
         double otherRegionsZ = 0.0;
         int x = blocksRegions[i].Members[k].x;
         int y = blocksRegions[i].Members[k].y;
         double logLum = LogLuminancePixels[y][x];
         for(int iter_y=y-mask/2.0; iter_y < y+mask/2.0;iter_y++)
         {
            for(int iter_x = x-mask/2.0; iter_x < x+mask/2.0; iter_x++)
            {
               if(iter_x >= 0 && iter_y >=0 && iter_x < imageWidth && iter_y < imageHeight)
               {
                  if(pixelsReg[y][x] == pixelsReg[iter_y][iter_x])
                  {
                     double logLumTmp = LogLuminancePixels[iter_y][iter_x];
                     double functionG = exp(-((iter_x - x)*(iter_x - x) + (iter_y - y)*(iter_y - y))/2*(pow(sigma_s,2.0)));
                     double functionK = exp(-((logLumTmp - logLum)*(logLumTmp - logLum))/2*(pow(sigma_r,2.0)));
                     inRegion += logLumTmp*functionG*functionK;
                     inRegionZ += functionG*functionK;
                  }
                  else{
                     double logLumTmp = LogLuminancePixels[iter_y][iter_x];
                     double functionG = exp(-((iter_x - x)*(iter_x - x) + (iter_y - y)*(iter_y - y))/2*(pow(sigma_s,2.0)));
                     double functionK = exp(-((logLumTmp - logLum)*(logLumTmp - logLum))/2*(pow(sigma_rr,2.0)));
                     otherRegions += logLumTmp*functionG*functionK;
                     otherRegionsZ += functionG*functionK;
                  }
               }
               
            }
         }
         double functionZ = inRegionZ + otherRegionsZ;
         localAdaptationPixels[y][x] = (1.0/functionZ)*(inRegion + otherRegions);
      }
   }
   fprintf(stderr,"\nBillateral filter completed\n");


   //Local tone mapping
   // psi(L,V;p,gamma) = (L/V)^p. 
   MinimalImageLuminance = (MinimalImageLuminance/(MinimalImageLuminance + 1));
   MinimalImageLuminance = pow(MinimalImageLuminance, 0.3);
   MaximalImageLuminance = (MaximalImageLuminance/(MaximalImageLuminance + 1));
   MaximalImageLuminance = pow(MaximalImageLuminance, 0.3);
   double inv_tmp = 1.0/(MaximalImageLuminance - MinimalImageLuminance);
   double alpha = 1.0*inv_tmp;
   double beta = -1.0*MinimalImageLuminance*inv_tmp;
   fprintf(stderr,"Alpha %g  Beta %g\n",alpha,beta);
   double max_scale = alpha*MaximalImageLuminance + beta;
   double min_scale = alpha*MinimalImageLuminance + beta;
   fprintf(stderr,"Max scale %g , Min scale %g\n",max_scale, min_scale);

   PixelIntMatrix regionBoundryPixels(imageHeight, vector<int>(imageWidth,0));
   PixelDoubleMatrix finalValuesPixels(imageHeight, vector<double>(imageWidth, 0.0));
   for(int i = 0; i < blocksRegions.size();i++)
   {
      for(int k = 0; k < blocksRegions[i].Members.size();k++)
      {
         int tmpboundry = 0;
         int x = blocksRegions[i].Members[k].x;
         int y = blocksRegions[i].Members[k].y;
         if(x-1 > 0)
         {
            if(pixelsReg[y][x] != pixelsReg[y][x-1])
            {
               regionBoundryPixels[y][x] = 1;
               tmpboundry = 1;
            }
         }
         if(x+1 < imageWidth)
         {
            if(pixelsReg[y][x] != pixelsReg[y][x+1])
            {
               regionBoundryPixels[y][x] = 1;
               tmpboundry = 1;
            }
         }
         if(y-1 > 0)
         {
            if(pixelsReg[y][x] != pixelsReg[y-1][x])
            {
               regionBoundryPixels[y][x] = 1;
               tmpboundry = 1;
            }
         }
         if(y+1 < imageHeight)
         {
            if(pixelsReg[y][x] != pixelsReg[y+1][x])
            {
               regionBoundryPixels[y][x] = 1;
               tmpboundry = 1;
            }
         }
         if(tmpboundry == 1)
         {
            
            int tmpx = blocksRegions[i].Members[k].x;
            int tmpy = blocksRegions[i].Members[k].y;
            double tmp_val = abs(log(LogLuminancePixels[tmpy][tmpx]/localAdaptationPixels[tmpy][tmpx]));
            blocksRegions[i].BoundryValues.push_back(tmp_val);
            blocksRegions[i].BoundryMap[tmp_val] = k;
         }
      }
   }
   for(int i=0; i < blocksRegions.size();i++)
   {
      if(blocksRegions[i].BoundryValues.size() > 60)
      {
         //vector<double> matrix_a;
         //vector<double> matrix_b;
         int samplesize = (blocksRegions[i].BoundryValues.size() * 5)/100;
         Eigen::MatrixXf a(samplesize,2);
         Eigen::VectorXf b(samplesize);
         for(int n=0; n < samplesize; n++)
         {
            sort(blocksRegions[i].BoundryValues.begin(),blocksRegions[i].BoundryValues.end());
            double tmp = blocksRegions[i].BoundryValues[n];
            auto chosed = blocksRegions[i].BoundryMap.find(tmp);
            int tmpx = blocksRegions[i].Members[chosed->second].x;
            int tmpy = blocksRegions[i].Members[chosed->second].y;
            
            double tmpB = (exp(localAdaptationPixels[tmpy][tmpx])/(exp(localAdaptationPixels[tmpy][tmpx]) + 1));
            tmpB = pow(tmpB, 0.3);
            double finalB = alpha*tmpB + beta;
            //matrix_b.push_back(finalB);
            b(n) = finalB;
            double tmpA = exp(LogLuminancePixels[tmpy][tmpx])/exp(localAdaptationPixels[tmpy][tmpx]);
            
            tmpA = pow(tmpA, 0.3);
            double finalA = tmpA * tmpB;
            //matrix_a.push_back(finalA);
            a(n,0) = finalA;
            a(n,1) = 1.0;

         }
         Eigen::VectorXf result(2);
         //result = a.fullPivHouseholderQr().solve(b);
         result = (a.transpose() * a).ldlt().solve(a.transpose() * b);
         if(result(0)==0.0)
         {
            result(0) = alpha;
         }
         fprintf(stderr,"Alpha-k %g beta-k %g\n",result(0),result(1));
         for(int k = 0; k < blocksRegions[i].Members.size();k++)
         {
            int x = blocksRegions[i].Members[k].x;
            int y = blocksRegions[i].Members[k].y;
            double tmpV = localAdaptationPixels[y][x];
            double tmpL = LogLuminancePixels[y][x];
            double p;
            if(log(exp(tmpL)/exp(tmpV)) <= -1.0)
            {
               p = 0.3;
            }
            else if(log(exp(tmpL)/exp(tmpV)) >= 1.0)
            {
               p = (0.3 + 1.8)/2.0;
            }
            else{
               p = 1.8;
            }
            if(regionBoundryPixels[y][x] == 1)
            {
               p = 0.3;
            }
            tmpL = tmpL/tmpV;
            tmpL = pow(abs(tmpL), p);
            tmpV = (tmpV/(tmpV+1));
            tmpV = pow(abs(tmpV), 0.3);
            //tmpV = result(0)*tmpV + result(1);
            finalValuesPixels[y][x] = result(0)*(tmpL * tmpV) + result(1);
            finalValuesPixels[y][x] = exp(finalValuesPixels[y][x]);
            
         }
      }
      else{
         for(int k = 0; k < blocksRegions[i].Members.size();k++)
         {
            int x = blocksRegions[i].Members[k].x;
            int y = blocksRegions[i].Members[k].y;
            double tmpV = localAdaptationPixels[y][x];
            double tmpL = LogLuminancePixels[y][x];
            double p;
            if(log(exp(tmpL)/exp(tmpV)) <= -1.0)
            {
               p = 0.3;
            }
            else if(log(exp(tmpL)/exp(tmpV)) >= 1.0)
            {
               p = (0.3 + 1.8)/2.0;
            }
            else{
               p = 1.8;
            }
            if(regionBoundryPixels[y][x] == 1)
            {
               p = 0.3;
            }
            tmpL = tmpL/tmpV;
            tmpL = pow(abs(tmpL), p);
            tmpV = (tmpV/(tmpV+1));
            tmpV = pow(abs(tmpV), 0.3);
            tmpV = tmpV;
            finalValuesPixels[y][x] =(tmpL * tmpV);
            finalValuesPixels[y][x] = alpha*exp(finalValuesPixels[y][x]) + beta; 
         }
      }
   }




	double pY, px, py;
   pSrc->Convert(TMO_Yxy);
   pDst->Convert(TMO_Yxy);
   double *pSourceData = pSrc->GetData();
	int j = 0;
   double prev_non_zero = 0.0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

			// Here you can use your transform
			// expressions and techniques...
			//pY *= dParameter; // Parameters can be used like
							  // simple variables
         
         pY *= finalValuesPixels[j][i];
			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
   
	pDst->Convert(TMO_RGB);
   pDst->CorrectGamma(2.2);
	return 0;
}
