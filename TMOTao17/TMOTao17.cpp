/* --------------------------------------------------------------------------- *
 * TMOHu14.cpp: implementation of the TMOHu14 class.   *
 * Diploma thesis
 * Author : Vladimir Vlkovic, Brno 2017
 * --------------------------------------------------------------------------- */

#include "TMOTao17.h"
#include <boost/concept_check.hpp>
#include <complex>



/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOTao17::TMOTao17()
{
	SetName(L"Tao17");						
	SetDescription(L"Image and video decolorization in real-time using dominant colors");

	
}

TMOTao17::~TMOTao17()
{
}

void TMOTao17::xyz2lab(double *data)
{
double refX =  95.047; 
double refY =  100.0;
double refZ =  108.883;

double x = *data++ ;
double y = *data++ ;
double z = *data++ ;

x /= refX;
y /= refY;
z /= refZ;

if ( x > 0.008856 ) x = std::pow(x,0.33);
else  x = ( 7.787 * x ) + ( 16 / 116 );

if ( y > 0.008856 ) y = std::pow(y,0.33);
else y = ( 7.787 * y ) + ( 16 / 116 );

if ( z > 0.008856 ) z = std::pow(z,0.33);
else z = ( 7.787 * z ) + ( 16 / 116 );



*--data = 200 * ( y - z ); //b
*--data = 500 * ( x - y ); //a
*--data = ( 116 * y ) - 16; //L


}

void TMOTao17::rgb2xyz(double *data)
{
  double r,g,b,x,y,z;
  r = *data++;
  g = *data++;
  b = *data++;
  
  if (r > 0.04045) r = std::pow((r + 0.055)/1.055, 2.4);
  else r = r /12.92;
  
  if (g > 0.04045) g = std::pow((g + 0.055)/1.055, 2.4);
  else g = g /12.92;
  
  if (b > 0.04045) b = std::pow((b + 0.055)/1.055, 2.4);
  else b = b /12.92;
  
  r = r *100;
  g= g*100;
  b=b*100;
  
  x = r * 0.4124 + g * 0.3576 + b * 0.1805;
  y = r * 0.2126 + g * 0.7152 + b * 0.0722; 
  z = r * 0.0193 + g * 0.1192 + b * 0.9505;

  *--data=z;
  *--data=y;
  *--data=x;
  
  
}

std::vector<double> TMOTao17::getPixelDifferences(std::vector<cv::Vec3d> labVector, int pixelCount)
{
  std::vector<double> pixelDiffs;
  for(int i =0; i<pixelCount; i++)
	{
	  for(int j=i+1;j<pixelCount;j++)
	  {
	    double L1,a1,b1,L2,a2,b2,diff;
	    L1=labVector[i][0];
	    a1=labVector[i][1];
	    b1=labVector[i][2];
	    
	    L2=labVector[j][0];
	    a2=labVector[j][1];
	    b2=labVector[j][2];
	    
	    diff = std::sqrt(std::pow(L1-L2,2)+std::pow(a1-a2,2)+std::pow(b1-b2,2));
	    pixelDiffs.push_back(diff);
	    
	  }
	  
	}
	return pixelDiffs;
}

double TMOTao17::rosen (const column_vector& m)

{
    const double x = m(0); 
    const double y = m(1);
    const double z = m(2);
   

  
      return  pow((x-y)+170,2) + pow((x-z)+114,2) + pow((y-z)+149,2);
}


int TMOTao17::Transform()
{

	
	double* pSourceData = pSrc->GetData();			
	double* pDestinationData = pDst->GetData();			 
													
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	int pixelCount = height * width;
	
	double* pSourceDataLAB=pSrc->GetData();//init
	
	std::vector<double> pixelDiffs;
	std::vector<cv::Vec3d> labVector;
	cv::Vec3d labColor;
	
	column_vector starting_point;
	starting_point.set_size(3);
	starting_point = 0; 
	double a = starting_point(0);
	double b = starting_point(1);
	double c = starting_point(2);
   dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),  
                             dlib::objective_delta_stop_strategy(1e-9),  
                             TMOTao17::rosen, dlib::derivative(TMOTao17::rosen), starting_point, 0, 255);
   
	 a = starting_point(0);
	 b = starting_point(1);
	 c = starting_point(2);

	for (int j = 0; j < pSrc->GetHeight(); j++)
	{
	    pSrc->ProgressBar(j, pSrc->GetHeight());
		
	    for (int i = 0; i < pSrc->GetWidth(); i++) ///result to output, taking only the image correction is discarded
	    {
	     // 
	      
	      //double w=*pSourceData++;
	    //  double v=*pSourceData++;
	   //   double x=*pSourceData++;
	      int d=0;
		/* rgb2xyz(pSourceData);
		 
		 xyz2lab(pSourceData);
		  labColor[0]=*pSourceData++;
		  labColor[1]=*pSourceData++;
		  labColor[2]=*pSourceData++;
		  labVector.push_back(labColor);*/
		 
		 
		 //pSourceData+=3;
	     
	 *pDestinationData++ =*pSourceData++;
		 *pDestinationData++ =*pSourceData++;
		 *pDestinationData++ =*pSourceData++;
	    }
	}
	
	//getPixelDifferences(labVector,pixelCount);
	
	  
	
	return 0;
}


