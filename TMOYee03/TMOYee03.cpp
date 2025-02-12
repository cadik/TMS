/************************************************************************************
*                                                                                   *
*                       Brno University of Technology                               *
*                       CPhoto@FIT                                                  *
*                                                                                   *
*                       Tone Mapping Studio                                         *
*                                                                                   *
*                       Bachelor thesis                                             *
*                       Author: Matus Bicanovsky [xbican03 AT stud.fit.vutbr.cz]        *
*                       Brno 2022                                                   *
*                                                                                   *
*                       Implementation of the TMOYee03 class                    *
*                                                                                   *
************************************************************************************/
/**
 * @file TMOYee03.cpp
 * @brief Implementation of the TMOYee03 class
 * @author Matus Bicanovsky
 * @class TMOYee03.cpp
 */

/* --------------------------------------------------------------------------- *
 * TMOYee03.cpp: implementation of the TMOYee03 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOYee03.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOYee03::TMOYee03()
{
	SetName(L"Yee03");					 
	SetDescription(L"Segmentation and adaptive assimilation for detail-preserving"); 

	bin_size1.SetName(L"bin_size1");				
	bin_size1.SetDescription(L"bin size determines contrast captured"); 
	bin_size1.SetDefault(0.5);							
	bin_size1 = 0.5;
	bin_size1.SetRange(0.5, 1.0); 
	this->Register(bin_size1);

   bin_size2.SetName(L"bin_size2");				
	bin_size2.SetDescription(L"bin size determines contrast captured"); 
	bin_size2.SetDefault(1.);							
	bin_size2 = 1.;
	bin_size2.SetRange(1., 2.0); 
	this->Register(bin_size2);

   small_threshold.SetName(L"small_threshold");				
	small_threshold.SetDescription(L"threshold determines when larger groups assimilate small groups"); 
	small_threshold.SetDefault(0.00001);						
	small_threshold = 0.00001;
	small_threshold.SetRange(0.00001, 0.03); 
	this->Register(small_threshold);

   big_threshold.SetName(L"big_threshold");				
	big_threshold.SetDescription(L"threshold determines when larger groups assimilate small groups"); 
	big_threshold.SetDefault(0.01);							
	big_threshold = 0.01;
	big_threshold.SetRange(0.01, 0.1); 
	this->Register(big_threshold);   

   max_layers.SetName(L"max_layers");				
	max_layers.SetDescription(L"max_layers controls how many layers to generate, its used to smooth out the boundaries"); 
	max_layers.SetDefault(16.);							
	max_layers = 16.;
	max_layers.SetRange(16., 96.0); 
	this->Register(max_layers);

}
//declaration of macros
#define cdm2ToLambert(C) (C*0.001/3.18)                                 // conversion from cd/m^2 to Lambert
#define lambertToCmd2(L) (L*3.18*1000)                                  // conversion from Lambert to cd/m^2
#define rgb2luminance(R,G,B) (R*0.3 + G*0.6 + B*0.1)                    // conversion of RGB to luminance
#define MAX_DISPLAY_LUMINANCE 125.0                                     // max display luminance value
#define DISPLAY_ADAPTATION_LUMINANCE 25.0
//declaration of global variables
unsigned int IMAGE_HEIGHT;
unsigned int IMAGE_WIDTH;
double bin_size = 0.0;
double MinimumImageLuminance = 0.0;
double stonits = 0.0;
int GroupNumber = 0;
//structers used for algorithm
// structure representing pixels with their coordinates
struct Point{
   unsigned short x,y;
};

typedef std::vector<Point> PointVector;

typedef std::list<unsigned short> NeighbourContainer;
// structure for storing informations about groups of pixels
struct Group_Record{
   PointVector Memebers;
   NeighbourContainer Neighbours;
   double Sum;
   unsigned int Count;
};

typedef std::vector<Group_Record *> Groups;

typedef std::vector<Groups> Layers;

typedef vector< vector<int> > PixelMatrix;

typedef vector< vector<double> > AdaptationMatrix;



Groups CategoryGroups;

TMOYee03::~TMOYee03()
{
}
//function for calculating category of given pixel based on equations presented in article http://www.cs.ucf.edu/~sumant/publications/VisualComp2003.pdf
// category(x,y) = (log10 liminance(x,y) - log10 minImageLuminance) / bin_size
double pixelCategory(double *image, int x, int y)
{
   double pR,pG,pB;
   pR = *(image+((y*IMAGE_WIDTH*3)+(x*3)));
   pG = *(image+((y*IMAGE_WIDTH*3)+(x*3))+1);
   pB = *(image+((y*IMAGE_WIDTH*3)+(x*3))+2);
   double category = ((log10(rgb2luminance(pR,pG,pB)+stonits)-log10(MinimumImageLuminance+stonits))/bin_size);
   category = std::round(category * 10.0) / 10.0;
   return category;
}
// valid function for Flood Fill algorithm to match pixels belonging into same category
bool isValid(double *image, int x, int y, double category, PixelMatrix& pixels)
{
   if(x<0 || x>= IMAGE_WIDTH || y<0 || y>= IMAGE_HEIGHT || pixelCategory(image,x,y) != category || pixels[x][y] == 1)
   {
      return false;
   }
   return true;
}
// implementation of Flood Fill algorithm used for grouping pixels of same category, instead of matching pixels based on their color we match the pixels based on their category
// the Flood Fill is performed in breadth first manner so contiguos pixels of the same category are assigned into same groups
void floodFill(double *image, int x, int y, double category, PixelMatrix& pixels, AdaptationMatrix& pixelCategories)
{
   vector<pair<int, int>> queue;
   pair<int, int> p(x,y);               // appending position of starting pixels in the group , parameter category describes the category we are currently creating
   queue.push_back(p);                  // creating queue for pixels 

   pixels[x][y] = 1;
   pixelCategories[x][y] = GroupNumber;

   Group_Record *group = new Group_Record;                  // creating new group and assigning starting parameters for them
   group->Count = 0;
   group->Sum = 0.0;

   Point point; 
   point.x = x;
   point.y = y;

   group->Memebers.push_back(point);
   double pR,pG,pB;
   pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
   pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
   pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
   double lum = rgb2luminance(pR,pG,pB);
   if(lum == 0.0)
   {
      lum = MinimumImageLuminance;
   }
   group->Sum += log10(lum+stonits);
   group->Count += 1;
   // loop while the queue is not empty in realization of Flood Fill algorithm
   while(queue.size() > 0)
   {
      pair<int,int> currPixel = queue[queue.size() - 1];                 // dequeue the front node
      queue.pop_back();

      int posX = currPixel.first;
      int posY = currPixel.second;
      // checking if the neighbouring pixels are valid , so if they belong into given category
      if(isValid(image, posX+1, posY, category, pixels))
      {
         // if neighbouring pixel is valid in terms of his category we add him into the group and adjust the parameters of group
         pixels[posX+1][posY] = 1;
         pixelCategories[posX+1][posY] = GroupNumber;
         point.x = posX+1;
         point.y = posY;

         group->Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0.0)
         {
            lum = MinimumImageLuminance;
         }
         group->Sum += log10(lum+stonits);
         group->Count += 1;

         p.first = posX+1;
         p.second = posY;
         queue.push_back(p);             // place the pixel into the queue
      }
      // following checks of neighbouring pixels are the same as the previous one just for neighbouring pixels in different directions
      if(isValid(image, posX-1, posY, category, pixels))
      {
         pixels[posX-1][posY] = 1;
         pixelCategories[posX-1][posY] = GroupNumber;
         point.x = posX-1;
         point.y = posY;

         group->Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0.0)
         {
            lum = MinimumImageLuminance;
         }
         group->Sum += log10(lum+stonits);
         group->Count += 1;
         
         p.first = posX-1;
         p.second = posY;
         queue.push_back(p);
      }
      if(isValid(image, posX, posY+1, category, pixels))
      {
         pixels[posX][posY+1] = 1;
         pixelCategories[posX][posY+1] = GroupNumber;
         point.x = posX;
         point.y = posY+1;

         group->Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0.0)
         {
            lum = MinimumImageLuminance;
         }
         group->Sum += log10(lum+stonits);
         group->Count += 1;
         
         p.first = posX;
         p.second = posY+1;
         queue.push_back(p);
      }
      if(isValid(image, posX, posY-1, category, pixels))
      {
         pixels[posX][posY-1] = 1;
         pixelCategories[posX][posY-1] = GroupNumber;
         point.x = posX;
         point.y = posY-1;

         group->Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0.0)
         {
            lum = MinimumImageLuminance;
         }
         group->Sum += log10(lum+stonits);
         group->Count += 1;
         
         p.first = posX;
         p.second = posY-1;
         queue.push_back(p);
      }

   }
   CategoryGroups.push_back(group);
   GroupNumber += 1;

}

// reworked FloodFill algorithm for finding neighbouring groups , the function is an opossite to usuall Flood Fill finding neighbouring pixels of different categories
// if such thing occures, the list of neighbouring groups for each group is afterwards updated
void GroupNeighbours(double *image, int x, int y, double category, PixelMatrix& pixels, AdaptationMatrix& pixelCategories)
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
      if(isValid(image, posX+1, posY, category, pixels))
      {
         pixels[posX+1][posY] = 1;
         
         p.first = posX+1;
         p.second = posY;
         queue.push_back(p);
      }
      if(!(isValid(image, posX+1, posY, category, pixels)))       // if two neighbouring pixels are from different groups then list Neighbours of both groups is updated with this information
      {
         if(posX+1>=0 && posX+1< IMAGE_WIDTH && posY>=0 && posY<IMAGE_HEIGHT && pixels[posX+1][posY]!=1) 
         {
            pixels[posX+1][posY] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX+1][posY]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX+1][posY]]->Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX+1][posY]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX+1][posY]]->Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end(), pixelCategories[posX+1][posY]) != CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]]->Neighbours.push_back(pixelCategories[posX+1][posY]);
            }
         }
      }

      if(isValid(image, posX-1, posY, category, pixels))
      {
         pixels[posX-1][posY] = 1;
         
         p.first = posX-1;
         p.second = posY;
         queue.push_back(p);
      }
      if(!(isValid(image, posX-1, posY, category, pixels)))
      {
         if(posX-1>=0 && posX-1< IMAGE_WIDTH && posY>=0 && posY<IMAGE_HEIGHT && pixels[posX-1][posY]!=1)
         {
            pixels[posX-1][posY] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX-1][posY]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX-1][posY]]->Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX-1][posY]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX-1][posY]]->Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end(), pixelCategories[posX-1][posY]) != CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]]->Neighbours.push_back(pixelCategories[posX-1][posY]);
            }
         }
         
      }
      if(isValid(image, posX, posY+1, category, pixels))
      {
         pixels[posX][posY+1] = 1;
         
         p.first = posX;
         p.second = posY+1;
         queue.push_back(p);
      }
      if(!(isValid(image, posX, posY+1, category, pixels)))
      {
         if(posX>=0 && posX< IMAGE_WIDTH && posY+1>=0 && posY+1<IMAGE_HEIGHT && pixels[posX][posY+1]!=1)
         {
            pixels[posX][posY+1] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY+1]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY+1]]->Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX][posY+1]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY+1]]->Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end(), pixelCategories[posX][posY+1]) != CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]]->Neighbours.push_back(pixelCategories[posX][posY+1]);
            }
         }
         
      }
      if(isValid(image, posX, posY-1, category, pixels))
      {
         pixels[posX][posY-1] = 1;
         
         p.first = posX;
         p.second = posY-1;
         queue.push_back(p);
      }
      if(!(isValid(image, posX, posY-1, category, pixels)))
      {
         if(posX>=0 && posX< IMAGE_WIDTH && posY-1>=0 && posY-1<IMAGE_HEIGHT && pixels[posX][posY-1]!=1)
         {
            pixels[posX][posY-1] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY-1]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY-1]]->Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX][posY-1]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY-1]]->Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]]->Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end(), pixelCategories[posX][posY-1]) != CategoryGroups[pixelCategories[posX][posY]]->Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]]->Neighbours.push_back(pixelCategories[posX][posY-1]);
            }
         }
         
      }

   }
}

int TMOYee03::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst
   pSrc->Convert(TMO_RGB);
   // declaration of variables based on user input or from default values
   int layer = 1;
   double Bin_size1 = bin_size1;
   double Bin_size2 = bin_size2;
   double Max_layers = max_layers;
   double *pSourceData = pSrc->GetData();		
	double *pDestinationData = pDst->GetData(); 
												
												
   
	double pR, pG, pB;
   IMAGE_HEIGHT = pSrc->GetHeight();
   IMAGE_WIDTH = pSrc->GetWidth();
   double Small_threshold = small_threshold;
   double Big_threshold = big_threshold;
   Small_threshold = (IMAGE_HEIGHT*IMAGE_WIDTH)*Small_threshold;
   Big_threshold = (IMAGE_HEIGHT*IMAGE_WIDTH)*Big_threshold;
   stonits = pSrc->GetStonits();
   double MinimalImageLuminance = 0., MaximalImageLuminance = 0., AverageImageLuminance = 0.;
   pSrc->GetMinMaxAvg(&MinimalImageLuminance, &MaximalImageLuminance, &AverageImageLuminance);          // calculating min, max and average luminance in the input image
   MinimumImageLuminance = MinimalImageLuminance;
   MaximalImageLuminance = MaximalImageLuminance*stonits;
   AdaptationMatrix adaptationPixels(IMAGE_WIDTH, vector<double>(IMAGE_HEIGHT,0));

   // loop over all layers of given image, based on parameter max_layers
	for(double currLayer=0.0; currLayer<Max_layers; currLayer++){
      fprintf(stderr,"Layer : %g\n",currLayer+1.0);
      bin_size = log10(Bin_size1)+(Bin_size2-log10(Bin_size1))*(currLayer/(Max_layers-1));
      fprintf(stderr,"Bin_size %g\n",bin_size);
      
      //creating groups of pixels belonging to same category based on Flood Fill algorithm
      AdaptationMatrix pixelCategories(IMAGE_WIDTH, vector<double>(IMAGE_HEIGHT,0));
      PixelMatrix pixels(IMAGE_WIDTH, vector<int>(IMAGE_HEIGHT,0));
      int groups = 0;
      GroupNumber = 0;
      for(int j = 0; j < IMAGE_HEIGHT; j++)
      {
         for(int i = 0; i < IMAGE_WIDTH; i++)
         {
            if(pixels[i][j] != 1)
            {
               double category = pixelCategory(pSourceData, i, j);
               
               floodFill(pSourceData, i, j, category, pixels, pixelCategories);
               groups++;
            }
         }
      }
      fprintf(stderr, "Small threshold: %g Big threshold: %g\n",Small_threshold,Big_threshold);
      fprintf(stderr,"Groups after grouping %d\n",CategoryGroups.size());
      
      // finding neighbouring groups with GroupNeighbours function implemented as adjusted Flood Fill algorithm
      PixelMatrix MembersPixels(IMAGE_WIDTH, vector<int>(IMAGE_HEIGHT,0));
      for(int j = 0; j < IMAGE_HEIGHT; j++)
      {
         for(int i = 0; i < IMAGE_WIDTH; i++)
         {
            
            if(MembersPixels[i][j] != 1)
            {
               double category = pixelCategory(pSourceData, i, j);
            
               GroupNeighbours(pSourceData, i, j, category, MembersPixels, pixelCategories);
               groups++;
            }
         }
      }
      //first part of assimilation process where we get rid of singletons, singleton is group with only one neighbouring group
      // if we want to perform assimilation the amount of pixels in smaller group has to be less than parameter small_threshold
      // and the bigger of the two groups has to have amount of pixels bigger then the parameter big_threshold
      for(int i =0; i < CategoryGroups.size();i++)
      {  
         if(CategoryGroups[i]->Neighbours.size()==1 && CategoryGroups[i]->Count > 0)
         {
            int position = CategoryGroups[i]->Neighbours.front();
            if(CategoryGroups[i]->Count > CategoryGroups[position]->Count)
            {
               if((CategoryGroups[i]->Memebers.size() > Big_threshold) &&(CategoryGroups[position]->Memebers.size()<Small_threshold))
               {
                  // if the conditions for assimilation are met , then we perform assimilation as described in the article http://www.cs.ucf.edu/~sumant/publications/VisualComp2003.pdf
                  CategoryGroups[i]->Memebers.insert(CategoryGroups[i]->Memebers.end(), CategoryGroups[position]->Memebers.begin(), CategoryGroups[position]->Memebers.end());
                  CategoryGroups[i]->Sum = (CategoryGroups[i]->Count + CategoryGroups[position]->Count) * CategoryGroups[i]->Sum/CategoryGroups[i]->Count;
                  CategoryGroups[i]->Count += CategoryGroups[position]->Count;
                  CategoryGroups[position]->Count = 0;
                  
               }
               
            }
            else{
               if((CategoryGroups[position]->Count > Big_threshold) &&(CategoryGroups[i]->Count <Small_threshold))
               {
                  CategoryGroups[position]->Memebers.insert(CategoryGroups[position]->Memebers.end(), CategoryGroups[i]->Memebers.begin(), CategoryGroups[i]->Memebers.end());
                  CategoryGroups[position]->Sum = (CategoryGroups[i]->Count + CategoryGroups[position]->Count) * CategoryGroups[position]->Sum/CategoryGroups[position]->Count;
                  CategoryGroups[position]->Count += CategoryGroups[i]->Count;
                  CategoryGroups[i]->Count = 0;
                  
               }
            }
         }
      }
      
      //second part of assimilation process where we assimilate groups that are smaller than given threshold even though they are not singletons
      // when we find such group , we iterate through their neighbouring groups and chose the biggest one that also meets the condition
      // that number of pixels in that group is larger than parameter big_threshold
      for(int i=0; i < CategoryGroups.size();i++)
      {
         if(CategoryGroups[i]->Count < Small_threshold && CategoryGroups[i]->Count > 0)
         {
            int biggestNeighbour = CategoryGroups[i]->Neighbours.front();
            for(const int & num : CategoryGroups[i]->Neighbours)
            {
               // iterating through neighbouring groups of the located group which is smaller than small_threshold
               if(CategoryGroups[num]->Count > CategoryGroups[biggestNeighbour]->Count)
               {
                  biggestNeighbour = num;
               }
            }
            if(CategoryGroups[biggestNeighbour]->Count > Big_threshold)
            {
               // if conditions are met we perform assimilation in same manner as with assimilating singletons
               CategoryGroups[biggestNeighbour]->Memebers.insert(CategoryGroups[biggestNeighbour]->Memebers.end(),CategoryGroups[i]->Memebers.begin(),CategoryGroups[i]->Memebers.end());
               CategoryGroups[biggestNeighbour]->Sum = (CategoryGroups[biggestNeighbour]->Count + CategoryGroups[i]->Count) * CategoryGroups[biggestNeighbour]->Sum/CategoryGroups[biggestNeighbour]->Count;
               CategoryGroups[biggestNeighbour]->Count += CategoryGroups[i]->Count;
               CategoryGroups[i]->Count = 0;
            }
         }
      }
      
      int groupsAfterAsimilation = 0;
      int pixelsEnd = 0;
      // storing the final values for every pixel from current calculated layer, and clearing the variables and vectors which will be used again in next iteration
      for(int i=0; i < CategoryGroups.size();i++)
      {
         
         if(CategoryGroups[i]->Count > 0)
         {  
            groupsAfterAsimilation++;
            for(int j=0; j < CategoryGroups[i]->Memebers.size();j++)
            {
               pixelsEnd++;
               adaptationPixels[CategoryGroups[i]->Memebers[j].x][CategoryGroups[i]->Memebers[j].y] += CategoryGroups[i]->Sum/CategoryGroups[i]->Count;
            }
            CategoryGroups[i]->Count = 0;
            CategoryGroups[i]->Memebers.clear();
            CategoryGroups[i]->Neighbours.clear();
            CategoryGroups[i]->Sum = 0;
         }
         
      }
      fprintf(stderr,"Groups after asimilation: %d ,pixels after asimilation: %d\n",groupsAfterAsimilation,pixelsEnd);
      fprintf(stderr,"\n");
      CategoryGroups.clear();
      pixelCategories.clear();
      MembersPixels.clear();
      pixels.clear();
   }
   fprintf(stderr,"Stonits : %g\n",stonits);
   
   //usage of TR algorithm to calculate final luminance values for each pixel in image, algorithm is in appendix in article http://www.cs.ucf.edu/~sumant/publications/VisualComp2003.pdf
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++) 
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); 
		for (int i = 0; i < pSrc->GetWidth(); i++) 
		{
			
         pR = *pSourceData++;
			pG = *pSourceData++;
			pB = *pSourceData++;
         
         double L_wa = cdm2ToLambert((adaptationPixels[i][j]/Max_layers));       // final adaptation luminance value of pixel as average from all computed layers
         double L_w = cdm2ToLambert(rgb2luminance(pR, pG, pB)*stonits);
         double f_r = cdm2ToLambert(pR)/L_w, f_g = cdm2ToLambert(pG)/L_w, f_b = cdm2ToLambert(pB)/L_w;

         double S_w = 100.0 + 10.0*log10(L_wa);
         double R_w = 10.0 * log10(L_wa/L_w);

         double L_dMax = cdm2ToLambert(MAX_DISPLAY_LUMINANCE);
         double L_da = cdm2ToLambert(DISPLAY_ADAPTATION_LUMINANCE);
         double S_d = 100.0 + 10.0*log10(L_da);
         double R_d = 8.4 - (S_w-27)*(8.4-R_w)/(S_d-27);

         double L_d = lambertToCmd2(L_da*pow(10,-0.1*R_d))/MAX_DISPLAY_LUMINANCE;
         // calculating new R,G,B values based on computed local adaptation luminance and results of TR algoritmh
         *pDestinationData++ = MIN(1.0,(L_d*f_r));
			*pDestinationData++ = MIN(1.0,(L_d*f_g));
			*pDestinationData++ = MIN(1.0,(L_d*f_b));
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
   pDst->CorrectGamma(2.2);
   pDst->Convert(TMO_RGB);
	return 0;
}
