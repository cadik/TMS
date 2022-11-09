/* --------------------------------------------------------------------------- *
 * TMOYee03.cpp: implementation of the TMOYee03 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOYee03.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOYee03::TMOYee03()
{
	SetName(L"Yee03");					  // TODO - Insert operator name
	SetDescription(L"Segmentation and adaptive assimilation for detail-preserving"); // TODO - Insert description

	bin_size1.SetName(L"bin_size1");				// TODO - Insert parameters names
	bin_size1.SetDescription(L"bin size determines contrast captured"); // TODO - Insert parameter descriptions
	bin_size1.SetDefault(0.5);							// TODO - Add default values
	bin_size1 = 0.5;
	bin_size1.SetRange(0.5, 1.0); // TODO - Add acceptable range if needed
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

#define cdm2ToLambert(C) (C*0.001/3.18)
#define lambertToCmd2(L) (L*3.18*1000)
#define rgb2luminance(R,G,B) (R*0.3 + G*0.6 + B*0.1)
#define MAX_DISPLAY_LUMINANCE 125.0
#define DISPLAY_ADAPTATION_LUMINANCE 25.0

unsigned int IMAGE_HEIGHT;
unsigned int IMAGE_WIDTH;
double bin_size = 0.0;
double MinimumImageLuminance = 0.0;
double stonits = 0.0;
int GroupNumber = 0;

struct Point{
   unsigned short x,y;
};

typedef std::vector<Point> PointVector;

typedef std::list<unsigned short> NeighbourContainer;
struct Group_Record{
   PointVector Memebers;
   NeighbourContainer Neighbours;
   double Sum;
   unsigned int Count;
};

typedef std::vector<Group_Record > Groups;

typedef std::vector<Groups> Layers;

typedef vector< vector<int> > PixelMatrix;

typedef vector< vector<double> > AdaptationMatrix;



Groups CategoryGroups;

TMOYee03::~TMOYee03()
{
}

double pixelCategory(double *image, int x, int y)
{
   double pR,pG,pB;
   pR = *(image+((y*IMAGE_WIDTH*3)+(x*3)));
   pG = *(image+((y*IMAGE_WIDTH*3)+(x*3))+1);
   pB = *(image+((y*IMAGE_WIDTH*3)+(x*3))+2);
   double category = ((log10(rgb2luminance(pR,pG,pB)+stonits)-log10(MinimumImageLuminance+stonits))/bin_size);
   category = std::round(category * 1.0) / 1.0;
   return category;
}

bool isValid(double *image, int x, int y, double category, PixelMatrix& pixels)
{
   if(x<0 || x>= IMAGE_WIDTH || y<0 || y>= IMAGE_HEIGHT || pixelCategory(image,x,y) != category || pixels[x][y] == 1)
   {
      return false;
   }
   return true;
}

void floodFill(double *image, int x, int y, double category, PixelMatrix& pixels, AdaptationMatrix& pixelCategories)
{
   vector<pair<int, int>> queue;
   pair<int, int> p(x,y);
   queue.push_back(p);

   pixels[x][y] = 1;
   pixelCategories[x][y] = GroupNumber;

   Group_Record group;
   group.Count = 0;
   group.Sum = 0.0;

   Point point; 
   point.x = x;
   point.y = y;

   group.Memebers.push_back(point);
   double pR,pG,pB;
   pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
   pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
   pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
   double lum = rgb2luminance(pR,pG,pB);
   if(lum == 0)
   {
      lum = MinimumImageLuminance/stonits;
   }
   group.Sum += log10(lum+stonits);
   group.Count += 1;

   while(queue.size() > 0)
   {
      pair<int,int> currPixel = queue[queue.size() - 1];
      queue.pop_back();

      int posX = currPixel.first;
      int posY = currPixel.second;

      if(isValid(image, posX+1, posY, category, pixels))
      {
         pixels[posX+1][posY] = 1;
         pixelCategories[posX+1][posY] = GroupNumber;
         point.x = posX+1;
         point.y = posY;

         group.Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0)
         {
            lum = MinimumImageLuminance/stonits;
         }
         group.Sum += log10(lum+stonits);
         group.Count += 1;

         p.first = posX+1;
         p.second = posY;
         queue.push_back(p);
      }

      if(isValid(image, posX-1, posY, category, pixels))
      {
         pixels[posX-1][posY] = 1;
         pixelCategories[posX-1][posY] = GroupNumber;
         point.x = posX-1;
         point.y = posY;

         group.Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0)
         {
            lum = MinimumImageLuminance/stonits;
         }
         group.Sum += log10(lum+stonits);
         group.Count += 1;
         
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

         group.Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0)
         {
            lum = MinimumImageLuminance/stonits;
         }
         group.Sum += log10(lum+stonits);
         group.Count += 1;
         
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

         group.Memebers.push_back(point);
         double pR,pG,pB;
         pR = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3)));
         pG = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+1);
         pB = *(image+((point.y*IMAGE_WIDTH*3)+(point.x*3))+2);
         double lum = rgb2luminance(pR,pG,pB);
         if(lum == 0)
         {
            lum = MinimumImageLuminance/stonits;
         }
         group.Sum += log10(lum+stonits);
         group.Count += 1;
         
         p.first = posX;
         p.second = posY-1;
         queue.push_back(p);
      }

   }
   CategoryGroups.push_back(group);
   GroupNumber += 1;

}


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
      if(!(isValid(image, posX+1, posY, category, pixels)))
      {
         if(posX+1>=0 && posX+1< IMAGE_WIDTH && posY>=0 && posY<IMAGE_HEIGHT)
         {
            pixels[posX+1][posY] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX+1][posY]].Neighbours.begin(), CategoryGroups[pixelCategories[posX+1][posY]].Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX+1][posY]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX+1][posY]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]].Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX+1][posY]) != CategoryGroups[pixelCategories[posX][posY]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX+1][posY]);
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
         if(posX-1>=0 && posX-1< IMAGE_WIDTH && posY>=0 && posY<IMAGE_HEIGHT)
         {
            pixels[posX-1][posY] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX-1][posY]].Neighbours.begin(), CategoryGroups[pixelCategories[posX-1][posY]].Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX-1][posY]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX-1][posY]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]].Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX-1][posY]) != CategoryGroups[pixelCategories[posX][posY]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX-1][posY]);
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
         if(posX>=0 && posX< IMAGE_WIDTH && posY+1>=0 && posY+1<IMAGE_HEIGHT)
         {
            pixels[posX][posY+1] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY+1]].Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY+1]].Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX][posY+1]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY+1]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]].Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX][posY+1]) != CategoryGroups[pixelCategories[posX][posY]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX][posY+1]);
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
         if(posX>=0 && posX< IMAGE_WIDTH && posY-1>=0 && posY-1<IMAGE_HEIGHT)
         {
            pixels[posX][posY-1] = 1;
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY-1]].Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY-1]].Neighbours.end(), pixelCategories[posX][posY]) != CategoryGroups[pixelCategories[posX][posY-1]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY-1]].Neighbours.push_back(pixelCategories[posX][posY]);
            }
            if(!(std::find(CategoryGroups[pixelCategories[posX][posY]].Neighbours.begin(), CategoryGroups[pixelCategories[posX][posY]].Neighbours.end(), pixelCategories[posX][posY-1]) != CategoryGroups[pixelCategories[posX][posY]].Neighbours.end())){
               CategoryGroups[pixelCategories[posX][posY]].Neighbours.push_back(pixelCategories[posX][posY-1]);
            }
         }
         
      }

   }
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOYee03::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst
   
   int layer = 1;
   double Bin_size1 = bin_size1;
   double Bin_size2 = bin_size2;
   double Max_layers = max_layers;
   double *pSourceData = pSrc->GetData();		// You can work at low level data
	double *pDestinationData = pDst->GetData(); // Data are stored in form of array
												// of three doubles representing
												// three colour components
	double pR, pG, pB;
   IMAGE_HEIGHT = pSrc->GetHeight();
   IMAGE_WIDTH = pSrc->GetWidth();
   double Small_threshold = small_threshold;
   double Big_threshold = big_threshold;
   Small_threshold = (IMAGE_HEIGHT*IMAGE_WIDTH)*Small_threshold;
   Big_threshold = (IMAGE_HEIGHT*IMAGE_WIDTH)*Big_threshold;
   stonits = pSrc->GetStonits();
   double MinimalImageLuminance = 0., MaximalImageLuminance = 0., AverageImageLuminance = 0.;
   pSrc->GetMinMaxAvg(&MinimalImageLuminance, &MaximalImageLuminance, &AverageImageLuminance);
   MinimumImageLuminance = MinimalImageLuminance;
   MaximalImageLuminance = MaximalImageLuminance*stonits;
   AdaptationMatrix adaptationPixels(IMAGE_WIDTH, vector<double>(IMAGE_HEIGHT,0));
	for(double currLayer=0.0; currLayer<Max_layers; currLayer++){
      fprintf(stderr,"Layer : %g\n",currLayer+1.0);
      bin_size = log10(Bin_size1)+(Bin_size2-log10(Bin_size1))*(currLayer/(Max_layers-1));
      fprintf(stderr,"Bin_size %g\n",bin_size);
      
      
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

      for(int i =0; i < CategoryGroups.size();i++)
      {
         if(CategoryGroups[i].Neighbours.size()==1 && CategoryGroups[i].Count > 0)
         {
            int position = CategoryGroups[i].Neighbours.front();
            if(CategoryGroups[i].Memebers.size() > CategoryGroups[position].Memebers.size())
            {
               if((CategoryGroups[i].Memebers.size() > Big_threshold) &&(CategoryGroups[position].Memebers.size()<Small_threshold))
               {
                  CategoryGroups[i].Memebers.insert(CategoryGroups[i].Memebers.end(), CategoryGroups[position].Memebers.begin(), CategoryGroups[position].Memebers.end());
                  CategoryGroups[i].Sum = (CategoryGroups[i].Count + CategoryGroups[position].Count) * CategoryGroups[i].Sum/CategoryGroups[i].Count;
                  CategoryGroups[i].Count += CategoryGroups[position].Count;
                  CategoryGroups[position].Count = 0;
                  
               }
               
            }
            else{
               if((CategoryGroups[position].Memebers.size() > Big_threshold) &&(CategoryGroups[i].Memebers.size()<Small_threshold))
               {
                  CategoryGroups[position].Memebers.insert(CategoryGroups[position].Memebers.end(), CategoryGroups[i].Memebers.begin(), CategoryGroups[i].Memebers.end());
                  CategoryGroups[position].Sum = (CategoryGroups[i].Count + CategoryGroups[position].Count) * CategoryGroups[position].Sum/CategoryGroups[position].Count;
                  CategoryGroups[position].Count += CategoryGroups[i].Count;
                  CategoryGroups[i].Count = 0;
                  
               }
            }
         }
      }


      for(int i=0; i < CategoryGroups.size();i++)
      {
         if(CategoryGroups[i].Count < Small_threshold && CategoryGroups[i].Count > 0)
         {
            int biggestNeighbour = CategoryGroups[i].Neighbours.front();
            for(const int & num : CategoryGroups[i].Neighbours)
            {
               if(CategoryGroups[num].Count > CategoryGroups[biggestNeighbour].Count)
               {
                  biggestNeighbour = num;
               }
            }
            if(CategoryGroups[biggestNeighbour].Count > Big_threshold)
            {
               CategoryGroups[biggestNeighbour].Memebers.insert(CategoryGroups[biggestNeighbour].Memebers.end(),CategoryGroups[i].Memebers.begin(),CategoryGroups[i].Memebers.end());
               CategoryGroups[biggestNeighbour].Sum = (CategoryGroups[biggestNeighbour].Count + CategoryGroups[i].Count) * CategoryGroups[biggestNeighbour].Sum/CategoryGroups[biggestNeighbour].Count;
               CategoryGroups[biggestNeighbour].Count += CategoryGroups[i].Count;
               CategoryGroups[i].Count = 0;
            }
         }
      }
      int groupsAfterAsimilation = 0;
      int pixelsEnd = 0;
      for(int i=0; i < CategoryGroups.size();i++)
      {
         
         if(CategoryGroups[i].Count > 0)
         {  
            groupsAfterAsimilation++;
            for(int j=0; j < CategoryGroups[i].Memebers.size();j++)
            {
               pixelsEnd++;
               adaptationPixels[CategoryGroups[i].Memebers[j].x][CategoryGroups[i].Memebers[j].y] += CategoryGroups[i].Sum/CategoryGroups[i].Count;
            }
         }
         
      }
      fprintf(stderr,"Groups after asimilation: %d ,pixels after asimilation: %d\n",groupsAfterAsimilation,pixelsEnd);
      fprintf(stderr,"\n");
      CategoryGroups.clear();
   }
   fprintf(stderr,"Stonits : %g\n",stonits);
   
   // Initialy images are in RGB format, but you can
	// convert it into other format
	//pSrc->Convert(TMO_Yxy); // This is format of Y as luminance
	//pDst->Convert(TMO_Yxy); // x, y as color information
   //double aAdaptationLuminance = 800.0;
	int j = 0;
	for (j = 0; j < pSrc->GetHeight(); j++) //
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++) //
		{
			
         pR = *pSourceData++;
			pG = *pSourceData++;
			pB = *pSourceData++;
         
         double L_wa = cdm2ToLambert((adaptationPixels[i][j]/Max_layers));
         double L_w = cdm2ToLambert(rgb2luminance(pR, pG, pB)*stonits);
         double f_r = pR/L_w, f_g = pG/L_w, f_b = pB/L_w;

         double S_w = 100.0 + 10.0*log10(L_wa);
         double R_w = 10.0 * log10(L_wa/L_w);

         double L_dMax = cdm2ToLambert(MAX_DISPLAY_LUMINANCE);
         double L_da = cdm2ToLambert(DISPLAY_ADAPTATION_LUMINANCE);
         double S_d = 100.0 + 10.0*log10(L_da);
         double R_d = 8.4 - (S_w-27)*(8.4-R_w)/(S_d-27);

         double L_d = lambertToCmd2(L_da*pow(10,-0.1*R_d))/MAX_DISPLAY_LUMINANCE;
         
         *pDestinationData++ = MIN(1.0,(L_d*f_r)/1000);
			*pDestinationData++ = MIN(1.0,(L_d*f_g)/1000);
			*pDestinationData++ = MIN(1.0,(L_d*f_b)/1000);
		}
	}
   fprintf(stderr, "\nMinimal image luminance: %g\n", MinimumImageLuminance*stonits);
   fprintf(stderr, "Maximal image luminance: %g\n",MaximalImageLuminance);
	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0;
}
