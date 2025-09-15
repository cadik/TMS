/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                       Bachelor thesis                                             *
*                       Author: Matus Bicanovsky [xbican03 AT stud.fit.vutbr.cz]        *
*                       Brno 2023                                                   *
*                                                                                   *
*                       Implementation of the TMOChen05 class                    *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOChen05.cpp
 * @brief Implementation of the TMOChen05 class
 * @author Matus Bicanovsky
 * @class TMOChen05.cpp
 */

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
	SetName(L"Chen05");					  
	SetDescription(L"Tone reproduction: a perspective from luminance-driven perceptual grouping"); 

	Theta.SetName(L"Theta");				
	Theta.SetDescription(L"ParameterDescription"); 
	Theta.SetDefault(1.5);							
	Theta = 1.5;
	Theta.SetRange(1.5, 2.0); 
	this->Register(Theta);

   Delta.SetName(L"Delta");				
	Delta.SetDescription(L"ParameterDescription"); 
	Delta.SetDefault(0.5);							
	Delta = 0.5;
	Delta.SetRange(0.5, 1.0); 
	this->Register(Delta);
}

#define rgb2luminance(R,G,B) (R*0.2126 + G*0.7152 + B*0.0722)                           // macro for computing luminance from R,G,B values 
unsigned int IMAGE_HEIGHT;
unsigned int IMAGE_WIDTH;

typedef vector< vector<float> > PixelMatrix;

typedef vector< vector<double> > PixelDoubleMatrix;

typedef std::vector< vector<int> > PixelIntMatrix;

typedef vector< vector<long> > PixelLongMatrix;
// structure for pixel , keeping value of x,y of given pixel
struct Point{
   unsigned short x,y;
};
// structre for signatures used in calculating EMD , keeping values s and w of each bin of signature
struct Signature{
   double s;
   double w;
};
typedef std::vector<Point> PointVector;
typedef std::vector<double> HistogramVector;
typedef std::vector<Signature> SignatureVector;
typedef std::vector<int> UnvisitedVector;
typedef std::vector<int> NeighbourContainer;
// structure for Blocks of regions , keeping information about members, neighbours and signature of block
struct Block_Record{
   PointVector Memebers;
   NeighbourContainer Neighbours;
   float Sum;
   unsigned int Count;
   HistogramVector logHistogram;
   SignatureVector blockSignature;
};
// structure for Regions , keeping informations about members, boundry of region , region signature values 
struct Region_Record{
   PointVector Members;
   float Sum;
   unsigned int Count;
   int firstsize;
   double max;
   double min;
   double first_t;
   double second_t;
   double firstLum;
   double secondLum;
   double thirdLum;
   double firstCnt;
   double secondCnt;
   double thirdCnt;
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
// function calculating Earth Movers Distance between 2 given blocks of pixels, used in algorithm for creating regions
float emdFunctionBB(int firstBlock, int secondBlock, Blocks& pixelBlocks)
{
   std::vector<double> av;
   std::vector<double> aw;
   std::vector<double> bv;
   std::vector<double> bw;
   for(int i=0;i < 3; i++)
   {
      av.push_back(pixelBlocks[firstBlock].blockSignature[i].s);
      aw.push_back(pixelBlocks[firstBlock].blockSignature[i].w);

   }
   for(int k=0; k < 3; k++)
   {
      bv.push_back(pixelBlocks[secondBlock].blockSignature[k].s);
      bw.push_back(pixelBlocks[secondBlock].blockSignature[k].w);
   } 
   return wasserstein(av,aw,bv,bw);          // function returns distance between blocks calculated with EMD implementation in wasserstein.h, implemented by Guilherme Nardari , https://github.com/gnardari/wasserstein
   
}
//similar to emdFunctionBB but it calculates the distance between given block and region , used in algorithm for creating regions
float emdFunctionBR(int block, int region, Blocks& pixelBlocks, Regions& pixelRegions)
{
   std::vector<double> av;
   std::vector<double> aw;
   std::vector<double> bv;
   std::vector<double> bw;
   for(int i=0;i < 3; i++)
   {
      av.push_back(pixelBlocks[block].blockSignature[i].s);
      aw.push_back(pixelBlocks[block].blockSignature[i].w);
   }
   for(int k=0; k < 3; k++)
   {
      bv.push_back(pixelRegions[region].regionSignature[k].s);
      bw.push_back(pixelRegions[region].regionSignature[k].w);
   } 
   return wasserstein(av,aw,bv,bw);
   
}
// support function for GroupNeighbours , checks if the neighbouring pixels belong to the same Block
bool isValid(int x, int y, int category, PixelIntMatrix& pixels, PixelIntMatrix& pixelCategories)
{
   if(x<0 || x>= IMAGE_HEIGHT || y<0 || y>= IMAGE_WIDTH || pixelCategories[x][y] != category || pixels[x][y] == 1)
   {
      return false;
   }
   return true;
}

// function for finding neighbouring blocks, similar to one implemented in Yee03, its adjusted floodfill algorithm , that detects if 2 neighbouring pixels are from different blocks
void GroupNeighbours(int x, int y, int group, PixelIntMatrix& pixels, PixelIntMatrix& pixelCategories, Blocks& pixelBlocks)
{
   vector<pair<int, int>> queue;
   pair<int, int> p(x,y);
   queue.push_back(p);

   pixels[x][y] = 1;


   while(queue.size() > 0)
   {
      pair<int,int> currPixel = queue[queue.size() - 1];      //dequeue the front node , and check each of its neighbouring pixels afterwards if they belong to same or different block
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
         // if we find neighbouring pixels which belong to different blocks, we store this information into Neighbours lists of both blocks
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
// function for updating the signature of region after new block was assigned into the region
void updateRegionSignature(int regionID, Regions& region, PixelDoubleMatrix& LogLuminancePixels, int amount)
{
   bool fullupdate = true;
   // checking if the log luminance values of newly added block contains new max or min log luminance value for current region
   // if yes we have to calculate new threshold for signature bins and afterwards calculate new signature from all the log luminance values currently in region
   // if log luminance values from newly added block wont change the thresholds we dont have to calculate new signature, we can just update the current signature with
   // assigning newly added pixels into correct bins 
   if(region[regionID].Members.size()!= 0 && region[regionID].Members.size() > region[regionID].firstsize)
   {
      for(int i= region[regionID].Members.size() - 1 - amount; i < region[regionID].Members.size();i++)
      {
         int x = region[regionID].Members[i].x;
         int y = region[regionID].Members[i].y;
         if(LogLuminancePixels[y][x] < region[regionID].min || LogLuminancePixels[y][x] > region[regionID].max)
         {
            fullupdate = true;
            break;
         }
         else{
            fullupdate = false;
         }
      }
   }
   // if creating new signature is needed based on newly added pixels
   if(fullupdate){
      double maxlum = LogLuminancePixels[region[regionID].Members[0].y][region[regionID].Members[0].x];
      double minlum = LogLuminancePixels[region[regionID].Members[0].y][region[regionID].Members[0].x];
      region[regionID].logHistogram.clear();
      region[regionID].regionSignature.clear();
      // itterate through members of region and find minimal and maximal log luminance value in region
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
      // calculate the step and thresholds for new bins in new signature
      double step = abs((maxlum-minlum)/3.0);
      double first_threshold = minlum + step;
      double second_threshold = first_threshold + step;
      double firstSectionLum =0.0, secondSectionLum=0.0, thirdSectionLum=0.0;                // variables for counting sum of log luminance values in each bin
      double firstSectionCount=0.0, secondSectionCount=0.0, thirdSectionCount=0.0;           // variables for counting amount of pixels in each bin
      region[regionID].max = maxlum;
      region[regionID].min = minlum;
      region[regionID].first_t = first_threshold;
      region[regionID].second_t = second_threshold;
      // itterate through log luminance values of region and assignem them based on their value into bins
      for(int k = 0; k < region[regionID].logHistogram.size(); k++)
      {
         if(region[regionID].logHistogram[k] <= first_threshold)
         {
            firstSectionLum += region[regionID].logHistogram[k];               // update sum of the bin
            firstSectionCount += 1.0;                                          // update number of pixels in the bin
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
      region[regionID].firstLum = firstSectionLum;
      region[regionID].firstCnt = firstSectionCount;
      region[regionID].secondLum = secondSectionLum;
      region[regionID].secondCnt = secondSectionCount;
      region[regionID].thirdLum = thirdSectionLum;
      region[regionID].thirdCnt = thirdSectionCount;
      double sumOfPixels = firstSectionCount + secondSectionCount + thirdSectionCount;
      double s1, w1, s2, w2, s3, w3;
      //calculate the signature values based on given bins, signature = ((s1,w1),(s2,w2),(s3,w3))
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
      region[regionID].regionSignature.push_back(signatureOfRegion);          // store the calculated signature into the struct of current region
   }
   // if its not needed to compute whole new signature we can just update current signature
   else{
      // add log luminance values of newly added pixel into regions histogram
      for(int j=region[regionID].Members.size()-1-amount; j < region[regionID].Members.size();j++)
      {
         int x = region[regionID].Members[j].x;
         int y = region[regionID].Members[j].y;
         double tmp = LogLuminancePixels[y][x];
         region[regionID].logHistogram.push_back(tmp);
      }
      // assign newly added log luminance values into bins based on their value
      for(int k = region[regionID].logHistogram.size()-1-amount; k < region[regionID].logHistogram.size(); k++)
      {
         if(region[regionID].logHistogram[k] <= region[regionID].first_t)
         {
            region[regionID].firstLum += region[regionID].logHistogram[k];
            region[regionID].firstCnt += 1.0;
         }
         else if(region[regionID].logHistogram[k] > region[regionID].first_t && region[regionID].logHistogram[k] <= region[regionID].second_t)
         {
            region[regionID].secondLum += region[regionID].logHistogram[k];
            region[regionID].secondCnt += 1.0;
         }
         else{
            region[regionID].thirdLum += region[regionID].logHistogram[k];
            region[regionID].thirdCnt += 1.0;
         }
      }
      double sumOfPixels = region[regionID].firstCnt + region[regionID].secondCnt + region[regionID].thirdCnt;
      double s1, w1, s2, w2, s3, w3;
      // calculate new signature based on updated values in each bin
      if(region[regionID].thirdCnt == 0)
      {
         s1 = 0.0;
         w1 = 0.0;
      }
      else{
         s1 = region[regionID].thirdLum/region[regionID].thirdCnt;
         w1 = region[regionID].thirdCnt/sumOfPixels;
      }
      
      Signature signatureOfRegion;
      signatureOfRegion.s = s1;
      signatureOfRegion.w = w1;
      region[regionID].regionSignature.push_back(signatureOfRegion);
      if(region[regionID].secondCnt == 0.0){
         s2 = 0.0;
         w2 = 0.0;
      }
      else{
         s2 = region[regionID].secondLum/region[regionID].secondCnt;
         w2 = region[regionID].secondCnt/sumOfPixels;
      }
      
      signatureOfRegion.s = s2;
      signatureOfRegion.w = w2;
      region[regionID].regionSignature.push_back(signatureOfRegion);
      if(region[regionID].firstCnt == 0.0)
      {
         s3 = 0.0;
         w3 = 0.0;
      }
      else{
         s3 = region[regionID].firstLum/region[regionID].firstCnt;
         w3 = region[regionID].firstCnt/sumOfPixels;
      }
      
      signatureOfRegion.s = s3;
      signatureOfRegion.w = w3;
      region[regionID].regionSignature.push_back(signatureOfRegion);              
   }
}



TMOChen05::~TMOChen05()
{
}

int TMOChen05::Transform()
{
	
   double *imageData = pSrc->GetData();
   int imageHeight = pSrc->GetHeight();
   int imageWidth = pSrc->GetWidth();
   IMAGE_HEIGHT = imageHeight;
   IMAGE_WIDTH = imageWidth;
   double MinimalImageLuminance = 0., MaximalImageLuminance = 0., AverageImageLuminance = 0.;
   pSrc->GetMinMaxAvg(&MinimalImageLuminance, &MaximalImageLuminance, &AverageImageLuminance);
   fprintf(stderr,"Min: %g Max: %g\n",MinimalImageLuminance,MaximalImageLuminance);
   

   
   double stonits = pSrc->GetStonits();
   PixelDoubleMatrix LogLuminancePixels(imageHeight, vector<double>(imageWidth, 0.0));
   double pixelR, pixelG, pixelB;
   cv::Mat image;
   cv::Mat finalImage;
   image = cv::Mat::zeros(imageHeight, imageWidth, CV_32F);
   //calculating log luminance value of each pixel in the input images
   for(int j=0; j < imageHeight; j++)
   {
      for(int i=0; i<imageWidth; i++)
      {
         pixelR = *imageData++;
         pixelG = *imageData++;
         pixelB = *imageData++;
         LogLuminancePixels[j][i] = log(rgb2luminance(pixelR, pixelG, pixelB));
         image.at<float>(j, i) = log(rgb2luminance(pixelR, pixelG, pixelB));         // store log luminance values of input images for detecting edges
      }
   }
   // part of the algorithm is usage of Canny edge detector on log luminance values from input image
   cv::normalize(image,image,0.0,1.0,cv::NORM_MINMAX, CV_32F);
   double imgMin, imgMax;
   cv::minMaxLoc(image, &imgMin, &imgMax);
   image.convertTo(finalImage, CV_8U, 255*(imgMax - imgMin));

   cv::Mat contours;
   cv::Mat gray_img;
   cv::Canny(finalImage, contours, 80, 240);                  // finding edges and storing them into the contours matrix

   fprintf(stderr,"Image height: %d Image width: %d\n",imageHeight, imageWidth);
   Blocks pixelBlocks;
   int edgeDetected = 0;
   //creating Blocks of pixels from input image, in the algorithm we are creating blocks of size 8x8 pixels,
   // however if there is edge in the block we split that block into 16 blocks of size 2x2
   // itterating through the image and creating the blocks of size 8x8
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
               if(contours.at<double>(k, l) != 0.0)          // if there is edge in any pixel of current block we dont create it and proceed to create 16 small blocks of 2x2 size
               {
                  edgeDetected = 1;
                  break;
               }
               p.x = l;
               p.y = k;
               block.Count += 1;
               block.Sum += LogLuminancePixels[k][l];
               block.Memebers.push_back(p);                // if the block doesnt contain the edges , we create him and store the pixels into members list
            }
         }
         // if the edge was detected we have to create 16 smaller blocks of 2x2 size
         if(edgeDetected == 1)
         {
            int x_value = i;
            int y_value = j;
            // creating first 4 small blocks in the first row of 8x8 block
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
               pixelBlocks.push_back(small_block);              // storing created block into vector of blocks
               x_value+=2;
            }
            x_value = i;
            y_value += 2;
            // creating second 4 small blocks in the second row of 8x8 block
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
            // creating third 4 small blocks in the third row of 8x8 block
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
            // creating fourth 4 small blocks in the fourth row of 8x8 block
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

   //Calculating signature for each created block
   for(int i=0; i < pixelBlocks.size();i++)
   {
      double maxlum = LogLuminancePixels[pixelBlocks[i].Memebers[0].y][pixelBlocks[i].Memebers[0].x];
      double minlum = LogLuminancePixels[pixelBlocks[i].Memebers[0].y][pixelBlocks[i].Memebers[0].x];
      //iterrate through members of block and find min and max log luminance value in given block
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
      // calculate the steps and thresholds for the 3 bins 
      double step = abs((maxlum-minlum)/3.0);
      double first_threshold = minlum + step;
      double second_threshold = first_threshold + step;
      double firstSectionLum =0.0, secondSectionLum=0.0, thirdSectionLum=0.0;
      double firstSectionCount=0.0, secondSectionCount=0.0, thirdSectionCount=0.0;
      // itterate through the log luminance values of given block and assignem them into the bins based on their value
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
      //calculate the signature values based on given bins, signature = ((s1,w1),(s2,w2),(s3,w3))
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
      pixelBlocks[i].blockSignature.push_back(signatureOfBlock);                // store the calculated signature into block's structure
   }

   PixelIntMatrix pixelsGrp(imageHeight, vector<int>(imageWidth,0));
   PixelIntMatrix visitedPixels(imageHeight, vector<int>(imageWidth,0));
   //calculating each pixel's group for the later checks
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
   //finding neighbouring blocks of each block with help of GroupNeighbours function
   for(int i=0; i < pixelBlocks.size();i++)
   {
      int x = pixelBlocks[i].Memebers[0].x;
      int y = pixelBlocks[i].Memebers[0].y;
      GroupNeighbours(y,x, pixelsGrp[y][x], visitedPixels, pixelsGrp, pixelBlocks);
   }
   
   
   UnvisitedVector UnvistitedBlocks(pixelBlocks.size(),0);                         // marking each block as unvisited
   PixelIntMatrix pixelsReg(imageHeight, vector<int>(imageWidth,0));
   Regions blocksRegions;
   
   double theta = Theta;                                     // variable used to allow block being assigned into region
   double delta = Delta;                                     // variable used to allow neighbouring blocks being put into queue to be assigned into region
   int unvisited = pixelBlocks.size();
   int regionID = 0;
   int counterTMP = 0;
   int chosedBlockID = 0;
   vector<int> queue;
   fprintf(stderr,"Theta: %g , Delta: %g\n",theta,delta);
   // itterate through blocks until all of them had been assigned into region
   while(unvisited > 0)
   {
      queue.clear();                                        // clearing queue for newly created region
      double biggestS1 = -10.0;
      int brightestBlockID = 0;
      for(int i=0; i < pixelBlocks.size();i++)             // itterate through unvisited blocks and find the brightest one , one with the biggest s1 value
      {
         if(UnvistitedBlocks[i]==0 && pixelBlocks[i].blockSignature[0].s > biggestS1)
         {
            biggestS1 = pixelBlocks[i].blockSignature[0].s;
            brightestBlockID = i;
         }
      }
      Region_Record region;                                // create new region
      blocksRegions.push_back(region);
      Signature region_Signature;                          // let the Earth Movers distance between brightest block chosed and newly created region equal 0 , so we copy the signature of block into newly created region 
      region_Signature.s = pixelBlocks[brightestBlockID].blockSignature[0].s;
      region_Signature.w = pixelBlocks[brightestBlockID].blockSignature[0].w;
      blocksRegions[regionID].regionSignature.push_back(region_Signature);
      region_Signature.s = pixelBlocks[brightestBlockID].blockSignature[1].s;
      region_Signature.w = pixelBlocks[brightestBlockID].blockSignature[1].w;
      blocksRegions[regionID].regionSignature.push_back(region_Signature);
      region_Signature.s = pixelBlocks[brightestBlockID].blockSignature[2].s;
      region_Signature.w = pixelBlocks[brightestBlockID].blockSignature[2].w;
      blocksRegions[regionID].regionSignature.push_back(region_Signature);
      blocksRegions[regionID].firstsize = pixelBlocks[brightestBlockID].Memebers.size();
      queue.push_back(brightestBlockID);                 // put the chosed block into the queue
      // while queue of blocks is not empty 
      while(queue.size() > 0)
      {
         double smallest  = 5.0;
         int tmpID=0;
         for(int l=0;l < queue.size();l++)                    // find the block in queue with smallest EMD between the region we are currently expanding and the block
         {
            float tmp = emdFunctionBR(queue[l],regionID,pixelBlocks,blocksRegions);
            if(tmp < smallest )
            {
               smallest = tmp;
               tmpID = l;
            }
         }
         chosedBlockID = queue[tmpID];
         if(smallest < theta)                                // if the distance between chosen block and region is smaller then parameter theta we are allowed to put it into region
         {
            queue.erase(queue.begin() + tmpID);              // remove chosen block from queue
            int amount = 0;
            for(int mem=0;mem < pixelBlocks[chosedBlockID].Memebers.size();mem++)            // insert members of the block into the region
            {
               blocksRegions[regionID].Members.push_back(pixelBlocks[chosedBlockID].Memebers[mem]);
               int x = pixelBlocks[chosedBlockID].Memebers[mem].x;
               int y = pixelBlocks[chosedBlockID].Memebers[mem].y;
               pixelsReg[y][x] = regionID;
               amount++;
            }
            counterTMP += 1;
            UnvistitedBlocks[chosedBlockID] = 1;                // mark the chosed block as visited
            unvisited -= 1;                                    // decrement the value of unvisited blocks
            
            for(int n=0; n< pixelBlocks[chosedBlockID].Neighbours.size();n++)                 // itterate through neighbouring blocks of the chosen block
            {
               if(UnvistitedBlocks[pixelBlocks[chosedBlockID].Neighbours[n]] == 0)
               {
                  
                  float tmpEMD = emdFunctionBB(chosedBlockID,pixelBlocks[chosedBlockID].Neighbours[n], pixelBlocks);       // calculate EMD between chosen block and its current neighbour
                  if(tmpEMD < delta)      // if the distance is smaller than parameter delta, we can put the neighbouring block into the queue
                  {
                     if(!(std::find(queue.begin(),queue.end(),pixelBlocks[chosedBlockID].Neighbours[n])!= queue.end()))        // check if the block isnt already in queue
                     {
                        queue.push_back(pixelBlocks[chosedBlockID].Neighbours[n]);
                     }
                     
                  }
               }
            }
            
            fprintf(stderr,"\rBlocks assigned to region %d/%d , regions created %d",counterTMP,pixelBlocks.size(),regionID);
            fflush(stdout);
            updateRegionSignature(regionID, blocksRegions, LogLuminancePixels, amount);            // update signature of region since new pixels were added into region
         }
         else{        // if the chosen block with smallest EMD between itself and current region doesnt meet the condition that the distance is smaller than theta, then none of other blocks in queue would, so we move on to create new region
            
            break;
         }
      }
      regionID += 1;
   }
   fprintf(stderr,"\nCreating regions completed\n");
   
   // application of Billateral filtering based on equations presented in article https://ieeexplore.ieee.org/document/1467466
   PixelDoubleMatrix localAdaptationPixels(imageHeight, vector<double>(imageWidth, 0.0));
   double sigma_r = 0.4;                                              // weight parameter for pixels in same region
   double sigma_rr = 0.5*0.4;                                         // weight parameter for pixels in different regions
   double sigma_s = (imageHeight*imageWidth)*0.04;
   double mask = 6;
   int bilateralIteration = 0;
   for(int i=0; i < blocksRegions.size();i++)
   {
      for(int k=0; k < blocksRegions[i].Members.size();k++)
      {
         bilateralIteration++;
         fprintf(stderr,"\rBillateral filter iterations %d/%d ",bilateralIteration,imageHeight*imageWidth);
         fflush(stdout);
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
                     double functionG = exp(-((iter_x - x)*(iter_x - x) + (iter_y - y)*(iter_y - y))/2*(pow(sigma_s,2.0)));   // G_xy(i,j) = exp{-((i-x)^2 + (j-y)^2)/2*sigma_s^2}
                     double functionK = exp(-((logLumTmp - logLum)*(logLumTmp - logLum))/2*(pow(sigma_r,2.0)));               // K_xy(i,j) = exp{-(LogLum(i,j) - LogLum(x,y))^2/2*sigma_r^2}
                     inRegion += logLumTmp*functionG*functionK;
                     inRegionZ += functionG*functionK;
                  }
                  else{
                     double logLumTmp = LogLuminancePixels[iter_y][iter_x];
                     double functionG = exp(-((iter_x - x)*(iter_x - x) + (iter_y - y)*(iter_y - y))/2*(pow(sigma_s,2.0)));  // G_xy(i,j) = exp{-((i-x)^2 + (j-y)^2)/2*sigma_s^2}
                     double functionK = exp(-((logLumTmp - logLum)*(logLumTmp - logLum))/2*(pow(sigma_rr,2.0)));             // K_xy(i,j) = exp{-(LogLum(i,j) - LogLum(x,y))^2/2*sigma_rr^2}
                     otherRegions += logLumTmp*functionG*functionK;
                     otherRegionsZ += functionG*functionK;
                  }
               }
               
            }
         }
         double functionZ = inRegionZ + otherRegionsZ;                      // Z_xy = sum(G_xy(i,j) * K_xy(i,j)) + sum(G_xy(i,j) * K_xy(i,j))
         localAdaptationPixels[y][x] = (1.0/functionZ)*(inRegion + otherRegions);     
      }
   }
   fprintf(stderr,"\nBillateral filter completed\n");
   // first part of tone mapping is to find global shaping and scaling parameters alpha and beta, based on solving linear system proposed in https://ieeexplore.ieee.org/document/1467466
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

   // for region wise tone mapping we need to have informations about boundry of each region
   // following block of code is finding pixels at boundry of each region and store them in the region structure
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
            double tmp_val = abs(log(exp(LogLuminancePixels[tmpy][tmpx])/exp(localAdaptationPixels[tmpy][tmpx])));    // we store log(L/V) of each boundry pixel for later calculations, L = log luminance value, V = local adaptation value
            blocksRegions[i].BoundryValues.push_back(tmp_val);
            blocksRegions[i].BoundryMap[tmp_val] = k;    // we create map in which based on its log(L/V) value we can find particular pixel on the boundry
         }
      }
   }
   // this loop itterates through all regions and calculate new luminance values for each pixel based on its region informations
   for(int i=0; i < blocksRegions.size();i++)
   {
      if(blocksRegions[i].BoundryValues.size() > 3)
      {
         int samplesize = 0;
         if(blocksRegions[i].BoundryValues.size() > 59)
         {
            samplesize = (blocksRegions[i].BoundryValues.size() * 5)/100;             // calculating samplesize of boundry pixels for next steps, its 5% of region boundry pixels of given region , as proposed in the article
         }
         else{
            samplesize = 3;                                                           // if the region is smaller than anticipated the samplesize is set to 3
         }
         
         Eigen::MatrixXf a(samplesize,2);
         Eigen::VectorXf b(samplesize);
         for(int n=0; n < samplesize; n++)
         {
            sort(blocksRegions[i].BoundryValues.begin(),blocksRegions[i].BoundryValues.end());    // sorting the values of boundry pixels in ascending order
            double tmp = blocksRegions[i].BoundryValues[n];                                      // storing the value log(L/V) of currently chosen boundry pixel
            auto chosed = blocksRegions[i].BoundryMap.find(tmp);
            int tmpx = blocksRegions[i].Members[chosed->second].x;
            int tmpy = blocksRegions[i].Members[chosed->second].y;
            // calculating the values for least square solution presented in the article https://ieeexplore.ieee.org/document/1467466
            double tmpB = (exp(LogLuminancePixels[tmpy][tmpx])/(exp(LogLuminancePixels[tmpy][tmpx]) + 1));
            double tmpC = (exp(localAdaptationPixels[tmpy][tmpx])/(exp(localAdaptationPixels[tmpy][tmpx]) + 1));
            tmpB = pow(tmpB, 0.3);
            double finalB = alpha*tmpB + beta;
            b(n) = finalB;
            double tmpA = exp(LogLuminancePixels[tmpy][tmpx])/exp(localAdaptationPixels[tmpy][tmpx]);
         
            tmpA = pow(tmpA, 0.3);
            double finalA = tmpA * tmpC;
            a(n,0) = finalA;
            a(n,1) = 1;

         }
         Eigen::VectorXf result(2);
         result = (a.transpose() * a).ldlt().solve(a.transpose() * b);              // solving the least square solution with the Eigen library
         //calculating new pixel luminance values based on the values given by solving the least square solution 
         // equation for new luminance value is newL(x,y) = alpha_k * ( (L(x,y)/V(x,y))^p * (V(x,y)/1 + V(x,y))^gamma ) + beta_k
         // gamma = 0.3, alpha_k is stored in result(0) and beta_k is stored in result(1)
         for(int k = 0; k < blocksRegions[i].Members.size();k++)
         {
            int x = blocksRegions[i].Members[k].x;
            int y = blocksRegions[i].Members[k].y;
            double tmpV = exp(localAdaptationPixels[y][x]);
            double tmpL = exp(LogLuminancePixels[y][x]);
            double p;
            if(log(tmpL/tmpV) <= -1.0)
            {
               p = 0.3;
            }
            else if(log(tmpL/tmpV) >= 1.0)
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
            tmpV = (tmpV/(tmpV+1.0));
            tmpV = pow(abs(tmpV), 0.3);
            finalValuesPixels[y][x] = result(0)*(tmpL * tmpV) + result(1);
            finalValuesPixels[y][x] = finalValuesPixels[y][x];
            if(finalValuesPixels[y][x] < 0 )
            {
               double orig = finalValuesPixels[y][x];
               double tmpX = exp(LogLuminancePixels[y][x]);
               tmpX = tmpX/(tmpX + 1);
               tmpX = pow(tmpX, 0.3);
               finalValuesPixels[y][x] = alpha*tmpX+beta;
            }
            if(isnan(finalValuesPixels[y][x]) || isnan(-1*finalValuesPixels[y][x]))
            {
               double tmpX = exp(LogLuminancePixels[y][x]);
               tmpX = tmpX/(tmpX + 1);
               tmpX = pow(tmpX, 0.3);
               finalValuesPixels[y][x] = alpha*tmpX+beta;
            }
            
         }
      }
      //workaround for really small regions is to just apply the global mapping part of algorithm on them
      else{     
         for(int k = 0; k < blocksRegions[i].Members.size();k++)
         {
            int x = blocksRegions[i].Members[k].x;
            int y = blocksRegions[i].Members[k].y;
            double tmpL = exp(LogLuminancePixels[y][x]);
            tmpL = tmpL/(tmpL + 1);
            tmpL = pow(abs(tmpL), 0.3);
            finalValuesPixels[y][x] = alpha*tmpL + beta;
         }
      }
   }




	double pY, px, py;
   pSrc->Convert(TMO_Yxy);
   pDst->Convert(TMO_Yxy);
   double *pSourceData = pSrc->GetData();
   double *pDestinationData = pDst->GetData();
	int j = 0;
   double prev_non_zero = 0.0;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); 
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			pY = *pSourceData++;
			px = *pSourceData++;
			py = *pSourceData++;

         
         pY *= finalValuesPixels[j][i];          // based on the calculated luminance values from the algorithm we create new output images

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
