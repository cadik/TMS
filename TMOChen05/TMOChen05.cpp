/* --------------------------------------------------------------------------- *
 * TMOChen05.cpp: implementation of the TMOChen05 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOChen05.h"

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

typedef std::vector< vector<int> > PixelIntMatrix;

typedef vector< vector<long> > PixelLongMatrix;

struct Point{
   unsigned short x,y;
};
struct Signature{
   float s;
   float w;
};
typedef std::vector<Point> PointVector;
typedef std::vector<float> HistogramVector;
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
};

typedef std::vector<Block_Record > Blocks;
typedef std::vector<Region_Record> Regions;

float emdFunctionBB(int firstBlock, int secondBlock, Blocks& pixelBlocks)
{
   cv::Mat sign1(cv::Size(2,3),CV_32FC1);
   cv::Mat sign2(cv::Size(2,3),CV_32FC1);
   for(int i=0;i < 3; i++)
   {
      sign1.at<float>(i,0) = pixelBlocks[firstBlock].blockSignature[i].s;
      sign1.at<float>(i,1) = pixelBlocks[firstBlock].blockSignature[i].w;
   }
   for(int k=0; k < 3; k++)
   {
      sign2.at<float>(k,0) = pixelBlocks[secondBlock].blockSignature[k].s;
      sign2.at<float>(k,1) = pixelBlocks[secondBlock].blockSignature[k].w;
   } 
   return cv::EMD(sign1,sign2,cv::DIST_L2);
   
}

float emdFunctionBR(int block, int region, Blocks& pixelBlocks, Regions& pixelRegions)
{
   cv::Mat sign1(cv::Size(2,3),CV_32FC1);
   cv::Mat sign2(cv::Size(2,3),CV_32FC1);
   for(int i=0;i < 3; i++)
   {
      sign1.at<float>(i,0) = pixelBlocks[block].blockSignature[i].s;
      sign1.at<float>(i,1) = pixelBlocks[block].blockSignature[i].w;
   }
   for(int k=0; k < 3; k++)
   {
      sign2.at<float>(k,0) = pixelRegions[region].regionSignature[k].s;
      sign2.at<float>(k,1) = pixelRegions[region].regionSignature[k].w;
   } 
   return cv::EMD(sign1,sign2,cv::DIST_L2);
   
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

void updateRegionSignature(int regionID, Regions& region, PixelMatrix& LogLuminancePixels)
{
   float maxlum = LogLuminancePixels[region[regionID].Members[0].y][region[regionID].Members[0].x];
   float minlum = LogLuminancePixels[region[regionID].Members[0].y][region[regionID].Members[0].x];
   region[regionID].logHistogram.clear();
   region[regionID].regionSignature.clear();
   for(int j=0; j < region[regionID].Members.size();j++)
   {
      int x = region[regionID].Members[j].x;
      int y = region[regionID].Members[j].y;
      float tmp = LogLuminancePixels[y][x];
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
   float step = (maxlum-minlum)/3.0;
   float first_threshold = minlum + step;
   float second_threshold = first_threshold + step;
   float firstSectionLum =0.0, secondSectionLum=0.0, thirdSectionLum=0.0;
   float firstSectionCount=0.0, secondSectionCount=0.0, thirdSectionCount=0.0;
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
   float sumOfPixels = firstSectionCount + secondSectionCount + thirdSectionCount;
   float s1 = thirdSectionLum/thirdSectionCount;
   float w1 = thirdSectionCount/sumOfPixels;
   Signature signatureOfRegion;
   signatureOfRegion.s = s1;
   signatureOfRegion.w = w1;
   region[regionID].regionSignature.push_back(signatureOfRegion);
   float s2 = secondSectionLum/secondSectionCount;
   float w2 = secondSectionCount/sumOfPixels;
   signatureOfRegion.s = s2;
   signatureOfRegion.w = w2;
   region[regionID].regionSignature.push_back(signatureOfRegion);
   float s3 = firstSectionLum/firstSectionCount;
   float w3 = firstSectionCount/sumOfPixels;
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

	double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components

   double *imageData = pSrc->GetData();
   int imageHeight = pSrc->GetHeight();
   int imageWidth = pSrc->GetWidth();
   IMAGE_HEIGHT = imageHeight;
   IMAGE_WIDTH = imageWidth;

   /*
   cv::Mat redM, greenM, blueM;
   redM = cv::Mat::zeros(imageHeight, imageWidth, CV_32F);
   blueM = cv::Mat::zeros(imageHeight, imageWidth, CV_32F);
   greenM = cv::Mat::zeros(imageHeight, imageWidth, CV_32F);
   for(int j=0; j < pSrc->GetHeight(); j++)
   {
      for(int i = 0; i < pSrc->GetWidth(); i++)
      {
         redM.at<float>(j, i) = *imageData++;
         greenM.at<float>(j, i) = *imageData++;
         blueM.at<float>(j, i) = *imageData++;
      }
   }
   cv::normalize(redM,redM, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);
   cv::normalize(blueM,blueM, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);
   cv::normalize(greenM,greenM, 0.0, 1.0, cv::NORM_MINMAX, CV_32F);
   cv::Mat tmpR, tmpB, tmpG;
   double rmin, rmax, bmin, bmax, gmin, gmax;
   cv::minMaxLoc(redM, &rmin, &rmax);
   redM.convertTo(tmpR, CV_8U, 255 * (rmax - rmin));
   cv::minMaxLoc(greenM, &gmin, &gmax);
   greenM.convertTo(tmpG, CV_8U, 255 * (gmax - gmin));
   cv::minMaxLoc(blueM, &bmin, &bmax);
   blueM.convertTo(tmpB, CV_8U, 255 * (bmax - bmin));
   */

   //cv::Mat image = tmpR + tmpG + tmpB;
   cv::Mat image = cv::imread("/home/matthewlele/images/hdr_images/84y7.tif");
   cv::Mat contours;
   cv::Mat gray_img;

   //cvtColor(image, gray_img, CV_RGB2GRAY);
   cv::Canny(image, contours, 80, 240);
   cv::namedWindow("Canny");
   //cv::imshow("Canny",contours);
   //cv::waitKey(0);
   
   double stonits = pSrc->GetStonits();
   PixelMatrix LogLuminancePixels(imageHeight, vector<float>(imageWidth, 0.0));
   double pixelR, pixelG, pixelB;
   for(int j=0; j < imageHeight; j++)
   {
      for(int i=0; i<imageWidth; i++)
      {
         pixelR = *imageData++;
         pixelG = *imageData++;
         pixelB = *imageData++;
         LogLuminancePixels[j][i] = log10(rgb2luminance(pixelR, pixelG, pixelB)+stonits);
      }
   }
   fprintf(stderr,"Image height: %d Image width: %d\n",imageHeight, imageWidth);
   Blocks pixelBlocks;
   int edgeDetected = 0;
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


   for(int i=0; i < pixelBlocks.size();i++)
   {
      float maxlum = LogLuminancePixels[pixelBlocks[i].Memebers[0].y][pixelBlocks[i].Memebers[0].x];
      float minlum = LogLuminancePixels[pixelBlocks[i].Memebers[0].y][pixelBlocks[i].Memebers[0].x];
      for(int j=0; j < pixelBlocks[i].Memebers.size();j++)
      {
         int x = pixelBlocks[i].Memebers[j].x;
         int y = pixelBlocks[i].Memebers[j].y;
         float tmp = LogLuminancePixels[y][x];
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
      float step = (maxlum-minlum)/3.0;
      float first_threshold = minlum + step;
      float second_threshold = first_threshold + step;
      float firstSectionLum =0.0, secondSectionLum=0.0, thirdSectionLum=0.0;
      float firstSectionCount=0.0, secondSectionCount=0.0, thirdSectionCount=0.0;
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
      float sumOfPixels = firstSectionCount + secondSectionCount + thirdSectionCount;
      float s1 = thirdSectionLum/thirdSectionCount;
      float w1 = thirdSectionCount/sumOfPixels;
      Signature signatureOfBlock;
      signatureOfBlock.s = s1;
      signatureOfBlock.w = w1;
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);
      float s2 = secondSectionLum/secondSectionCount;
      float w2 = secondSectionCount/sumOfPixels;
      signatureOfBlock.s = s2;
      signatureOfBlock.w = w2;
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);
      float s3 = firstSectionLum/firstSectionCount;
      float w3 = firstSectionCount/sumOfPixels;
      signatureOfBlock.s = s3;
      signatureOfBlock.w = w3;
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);
 
   }
   
   PixelIntMatrix pixelsGrp(imageHeight, vector<int>(imageWidth,0));
   PixelIntMatrix visitedPixels(imageHeight, vector<int>(imageWidth,0));

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
   
   for(int i=0; i < pixelBlocks.size();i++)
   {
      int x = pixelBlocks[i].Memebers[0].x;
      int y = pixelBlocks[i].Memebers[0].y;
      GroupNeighbours(y,x, pixelsGrp[y][x], visitedPixels, pixelsGrp, pixelBlocks);
   }
   
   
   UnvisitedVector UnvistitedBlocks(pixelBlocks.size(),0);
   Regions blocksRegions;
   vector<int> queue;
   float theta = 0.3; 
   float delta = 0.5;
   int unvisited = 5200;
   int regionID = 0;
   int counterTMP = 0;
   fprintf(stderr,"%g %g\n",theta, delta);
   
   while(unvisited > 0)
   {
      float biggestS1 = 0.0;
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
      region.regionSignature.push_back(pixelBlocks[brightestBlockID].blockSignature[0]);
      region.regionSignature.push_back(pixelBlocks[brightestBlockID].blockSignature[1]);
      region.regionSignature.push_back(pixelBlocks[brightestBlockID].blockSignature[2]);
      blocksRegions.push_back(region);
      queue.push_back(brightestBlockID);
      while(queue.size()>0)
      {
         float smallest  = 5.0;
         int tmpID=0;
         for(int l=0;l < queue.size();l++)
         {
            float tmp = emdFunctionBR(queue[l],regionID,pixelBlocks,blocksRegions);
            if(tmp < smallest)
            {
               smallest = tmp;
               tmpID = l;
            }
         }
         //fprintf(stderr,"Smallest EMD %g\n",smallest);
         int chosedBlockID = queue[tmpID];
         queue.erase(queue.begin() + tmpID);
         if(smallest < theta)
         {
            
            for(int mem=0;mem < pixelBlocks[chosedBlockID].Memebers.size();mem++)
            {
               blocksRegions[regionID].Members.push_back(pixelBlocks[chosedBlockID].Memebers[mem]);
            }
            counterTMP += 1;
            fprintf(stderr,"Chosed block %d region %d iterations %d distance %f\n",chosedBlockID,blocksRegions[regionID].Members.size(),counterTMP, smallest);
            //blocksRegions[regionID].Members = pixelBlocks[chosedBlockID].Memebers;
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
                     if(!(std::find(queue.begin(),queue.end(),pixelBlocks[chosedBlockID].Neighbours[n])!=queue.end()))
                     {
                        //fprintf(stderr,"Neighbour of block %d num %d\n",chosedBlockID,n);
                        queue.push_back(pixelBlocks[chosedBlockID].Neighbours[n]);
                     }
                     
                  }
               }
            }
            updateRegionSignature(regionID, blocksRegions, LogLuminancePixels);
            //fprintf(stderr,"Region %d s1 %f\n",regionID,blocksRegions[regionID].regionSignature[0].s);
            
         }
         else{
            regionID += 1;
            break;
         }
      }
   }
   for(int m=0; m < blocksRegions.size();m++)
   {
      fprintf(stderr,"Region 0 members %d , queue size %d\n",blocksRegions[m].Members.size(),queue.size());
   }
   

   //V(x,y) = 1/Z(x,y).(sum[each pixel in region i,j](LogL(i,j).Gxy(i,j).Kxy(i,j) + sum[everypixel not in region i,j](LogL(ij).Gxy(ij).K'xy(ij))
   //Gxy(i, j) = exp(-((i-x)^2+(j-y)^2)/2sigma_s^2)
   //Kxy(i, j) = exp(-(LogL(i,j) - LogL(x,y))^2/2sigma_r^2)
   //K'xy(i, j) = exp(-(LogL(i, j) - LogL(x,y))^2/2sigma_r'^2)
   //Z(x,y) = sum[pixel in region](Gxy(i, j).Kxy(i, j) + sum[pixel not in region](Gxy(i, j).K'xy(i, j)





   //Local tone mapping
   // psi(L,V;p,gamma) = (L/V)^p. 





	double pY, px, py;

	int j = 0;
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

			// and store results to the destination image
			*pDestinationData++ = pY;
			*pDestinationData++ = px;
			*pDestinationData++ = py;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	pDst->Convert(TMO_RGB);
	return 0;
}
