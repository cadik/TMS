/* --------------------------------------------------------------------------- *
 * TMOHu14.cpp: implementation of the TMOTao17 class.   *
 * Diploma thesis
 * Author : Vladimir Vlkovic, Brno 2017
 * --------------------------------------------------------------------------- */

#include "TMOTao17.h"
#include <boost/concept_check.hpp>
#include <complex>
#include "glog/logging.h"
// Each solver is defined in its own header file.
// include the solver you wish you use:
//#include "pallas/simulated_annealing.h"
//#include "pallas/differential_evolution.h"
//#include "pallas/basinhopping.h"
//#include <dlib/optimization.h>
//#include <include/GridCut/GridGraph_2D_4C.h>
 std::vector<float> pixelDiffs;
 std::vector<cv::Vec3f> labVector2;
 std::vector<cv::Vec3f> labVector;
 std::vector<cv::Vec3f> rgbVec;
 std::vector<float> grayVec;
 std::vector<int> corVec;
cv::Mat flow,cf,pf,pg;
int x,y;
/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOTao17::TMOTao17()
{
	SetName(L"Tao17");						
	SetDescription(L"Video decolorization using visual proximity coherence optimization");

	beta.SetName(L"beta");
	beta.SetDescription(L"Parameter to control the relative weights of spatial and temporal terms");
	beta.SetDefault(0.5);
	beta.SetRange(0.0,1.0);
	this->Register(beta);
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
else  x = ( 7.787 * x ) + ( 16.0 / 116.0 );

if ( y > 0.008856 ) y = std::pow(y,0.33);
else y = ( 7.787 * y ) + ( 16.0 / 116.0 );

if ( z > 0.008856 ) z = std::pow(z,0.33);
else z = ( 7.787 * z ) + ( 16.0 / 116.0 );

double L=( 116.0 * y ) - 16.0;
double a=500.0 * ( x - y );
double b=200.0 * ( y - z );

*--data = 200.0 * ( y - z ); //b
*--data = 500.0 * ( x - y ); //a
*--data = ( 116.0 * y ) - 16.0; //L


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




/*
 * find minimal image energy works quite fine on images , not so much on videos
 * used algorithm from https://www.sciencedirect.com/science/article/pii/S0097849317301024
 * it did kind of work...on images
 * */
int TMOTao17::getWeights(float &wr,float &wg, float &wb)
{
  int weight_r,weight_b,weight_g;	
  double sqrtPi = std::sqrt(2*M_PI) ;
  double grayDiff =0.0;
  float diff,max =0.0;
  for(weight_r = 0; weight_r <= 10; weight_r++)
  {
    for(weight_g = 0; weight_g <= 10; weight_g++)
    {
      weight_b = 10 - (weight_r + weight_g);
      if(weight_b >= 0)
      {
	float weight_r_f = weight_r / 10.0;
	float weight_g_f = weight_g / 10.0;
	float weight_b_f = weight_b / 10.0; 
	
	double maxed=0.0;
	
	 for (int i = 0; i < rgbVec.size(); i++)
	{
	   
		
	    for (int j = i+1; j < rgbVec.size(); j++) 
	    {
	     
	      
		float L1,a1,b1,L2,a2,b2;
		
		  L1=labVector[i][0];
		  a1=labVector[i][1];
		  b1=labVector[i][2];
		  
		  L2=labVector[j][0];
		  a2=labVector[j][1];
		  b2=labVector[j][2];
		  diff = std::sqrt(std::pow(L1-L2,2)+std::pow(a1-a2,2)+std::pow(b1-b2,2)); //lab differences
		  grayDiff = (rgbVec[i][0]*weight_r_f + rgbVec[i][1]*weight_g_f + rgbVec[i][2]*weight_b_f) - 
			      (rgbVec[j][0]*weight_r_f + rgbVec[j][1]*weight_g_f + rgbVec[j][2]*weight_b_f);
		
		
		if (grayDiff < 0.0 )
		{
		  if (diff > 0.0) diff= -diff;  //diff has always diffeent sign
		}
		else 
		{
		  if(diff < 0.0) diff = std::abs(diff);
		}
		  
		  maxed += 1.0 / 0.35*sqrtPi * std::exp(-std::pow(((grayDiff-diff)-0.35),2)/(2*std::pow(0.35,2)));
		 
	    }
	}
	if(maxed> max)
	{
	  max = maxed;
	  wr= weight_r_f;
	  wg=weight_g_f;
	  wb=weight_b_f;
	}
      }
      
    }
  }
  return 0;
}

/*get minal temporal energy, didnt work so i tired same thing as in Transform but the energie should have been minimized
*together now it doesnt work properly
*/

int TMOTao17::getWeightsTemp(float &wr,float &wg, float &wb)
{
  int weight_r,weight_b,weight_g;	
  double sqrtPi = std::sqrt(2*M_PI) ;
  double grayDiff =0.0;
  float diff,min =0.0;
  for(weight_r = 0; weight_r <= 10; weight_r++)
  {
    for(weight_g = 0; weight_g <= 10; weight_g++)  ///66 variants of weights
    {
      weight_b = 10 - (weight_r + weight_g);
      if(weight_b >= 0)
      {
	float weight_r_f = weight_r / 10.0;
	float weight_g_f = weight_g / 10.0;
	float weight_b_f = weight_b / 10.0;
	
	double maxed=0.0;
	
	 for (int i = 0; i < rgbVec.size(); i++)
	{
	
		float L1,a1,b1,L2,a2,b2;
		
		  L1=labVector[i][0];
		  a1=labVector[i][1];
		  b1=labVector[i][2];
		  
		  L2=labVector2[i][0];
		  a2=labVector2[i][1];
		  b2=labVector2[i][2];
		  diff = std::sqrt(std::pow(L1-L2,2)+std::pow(a1-a2,2)+std::pow(b1-b2,2)); //lab color diffs
		  grayDiff =grayVec[i] - (rgbVec[i][0]*weight_r_f + rgbVec[i][1]*weight_g_f + rgbVec[i][2]*weight_b_f);
		
		
		if (grayDiff < 0.0 )
		{
		  if (diff > 0.0) diff= -diff;
		}
		else 
		{
		  if(diff < 0.0) diff = std::abs(diff);
		}
		  
		  maxed += std::pow(grayDiff-diff,2);///1.0 / 0.35*sqrtPi * std::exp(-std::pow(((grayDiff-diff)-0.35),2)/(2*std::pow(0.35,2)));
		 
	    }  ///trid something probably didnt work
	    if (min==0.0) min=maxed;
	
	if(maxed< min)
	{
	  min = maxed;
	  wr= weight_r_f;
	  wg=weight_g_f;
	  wb=weight_b_f;
	}
      }
      
    }
  } 
  return 0;
}



int TMOTao17::Transform()
{

	
  double* pSourceData = pSrc->GetData();			
  double* pDestinationData = pDst->GetData();			 
												  
  int height = pSrc->GetHeight();
  int width = pSrc->GetWidth();
  int pixelCount = height * width;

  
  float wr,wg,wb;

  cv::Mat tmp,rgb, rgb_small; 

rgb=cv::Mat::zeros(height, width, CV_32FC3);

  
  for (int j = 0; j < pSrc->GetHeight(); j++)
  {
      pSrc->ProgressBar(j, pSrc->GetHeight());
	  
      for (int i = 0; i < pSrc->GetWidth(); i++)
      {
      
	  
	    rgb.at<cv::Vec3f>(j,i)[0]=*pSourceData++;
	    rgb.at<cv::Vec3f>(j,i)[1]=*pSourceData++;
	    rgb.at<cv::Vec3f>(j,i)[2]=*pSourceData++;

      }
  }
  
  cv::Size size(SUBSAMPLE,SUBSAMPLE);
  cv::resize(rgb,rgb_small,size);
  cv::cvtColor(rgb_small, tmp, cv::COLOR_RGB2Lab);

  
  
  for (int j = 0; j < SUBSAMPLE; j++)
  {
  
      for (int i = 0; i < SUBSAMPLE; i++) 
      {
	labVector.push_back(tmp.at<cv::Vec3f>(j,i));   ///get samples SUBSAMPLE * SUBSAMPLE of lab image and rgb image
	rgbVec.push_back(rgb.at<cv::Vec3f>(j,i));
      }
  }
  getWeights(wr,wg,wb); 

  double res;
    
    for (int j = 0; j < pSrc->GetHeight(); j++)
  {
      pSrc->ProgressBar(j, pSrc->GetHeight());
	  
      for (int i = 0; i < pSrc->GetWidth(); i++) 
      {
	  
	res = rgb.at<cv::Vec3f>(j,i)[0] * wr + rgb.at<cv::Vec3f>(j,i)[1] * wg + rgb.at<cv::Vec3f>(j,i)[2] * wb;
	  *pDestinationData++ =res;
	  *pDestinationData++ =res;
	  *pDestinationData++ =res;
      }
  }
	

	return 0;
}
/* 
 * Returns the optical flow between frame*/
int TMOTao17::getOpticalFlow(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat &flow)
{
  cv::Mat CFgray,PFgray;
  cv::cvtColor(currentFrame,CFgray,CV_BGR2GRAY);
  cv::cvtColor(previousFrame,PFgray,CV_BGR2GRAY);
  cv::calcOpticalFlowFarneback(CFgray,PFgray, flow,0.5, 3, 15, 3, 5, 1.2, 0);
  return 0;
}
/**
 * Should return coherence refinemt frame but it doesnt, it returns sometihing not worth while
 * @currentFrame current color frame
 * @previousFrame revious color frame
 * @previousGray previos gray frame
 * @CFrame result coherence refinement frame
 */

int TMOTao17::getCoherenceRefinementFrame(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat previousGray, cv::Mat &CFrame)
{
  cv::Mat flow,cf,pf,pg;
  float  new_j,new_i,M_i,aDiff,bDiff,LDiff,tmp;
  double pc_i,firstPart,secondPart;
  new_i=0.0;
  new_j=0.0;
  pc_i=1.0;
  getOpticalFlow(currentFrame,previousFrame,flow);
  cv::normalize(currentFrame,cf,0.0,1.0,cv::NORM_MINMAX,CV_32F);
  cv::normalize(previousFrame,pf,0.0,1.0,cv::NORM_MINMAX,CV_32F);
  cv::normalize(previousGray,pg,0.0,1.0,cv::NORM_MINMAX,CV_32F);
  cv::normalize(flow,flow,0.0,1.0,cv::NORM_MINMAX,CV_32F);
  //typedef GridGraph_2D_4C<short,short,float> Grid;
  //Grid* grid = new Grid(currentFrame.cols,currentFrame.rows);  ///tried grid graph didint work
  for(int j=0; j < CFrame.rows;j++)
  {
    for(int i=0;i<CFrame.cols;i++)
    {
      float f = flow.at<cv::Vec2f>(j,i)[1];
      new_j = std::abs(j - flow.at<cv::Vec2f>(j,i)[1]);
      new_i = std::abs(i - flow.at<cv::Vec2f>(j,i)[0]);
      LDiff = cf.at<cv::Vec3f>(j,i)[0] - pf.at<cv::Vec3f>(new_j,new_i)[0];
      aDiff = cf.at<cv::Vec3f>(j,i)[1] - pf.at<cv::Vec3f>(new_j,new_i)[1];
      bDiff = cf.at<cv::Vec3f>(j,i)[2] - pf.at<cv::Vec3f>(new_j,new_i)[2];
      
      tmp=cf.at<cv::Vec3f>(j,i)[0];
      float ee=pg.at<float>(j,i);
      float eee=pg.at<float>(new_j,new_i); //tried sometihing from the article but didnt work, needs work
      M_i = std::sqrt(std::pow(LDiff,2) + std::pow(aDiff,2) + std::pow(bDiff,2));
      firstPart=std::exp(-(tmp - pg.at<float>(j,i))/std::pow(0.35,2));
      secondPart=std::exp(-1.0/M_i * (tmp - pg.at<float>(new_j,new_i))/std::pow(0.35,2));
      pc_i = firstPart * secondPart;
     CFrame.at<float>(j,i) = pc_i;
   
    /* if (i<CFrame.cols-1)
      {
        

        grid->set_neighbor_cap(grid->node_id(i  ,j),+1,0,pc_i);
       	grid->set_neighbor_cap(grid->node_id(i+1,j),-1,0,pc_i);
      }

      if (y<CFrame.rows-1)
      {
        grid->set_neighbor_cap(grid->node_id(i,j  ),0,+1,pc_i);
       grid->set_neighbor_cap(grid->node_id(i,j+1),0,-1,pc_i);
      }*/
      
    }
  }
 /* grid->compute_maxflow();
    for(int j=0; j < CFrame.rows;j++)
  {
    for(int i=0;i<CFrame.cols;i++)
    {
      std::cerr<< grid->get_segment(grid->node_id(i,j)) <<std::endl;
    }
  }*/
  

  return 0;
  
}

/** gets differential refinement frame
 * @currentFrame previous color frame
 * @previousFrame previous color frame
 * @DFrame result for difference refinemet frame 
 */
int TMOTao17::getDifferentialRefinementFrame(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat &DFrame)
{
  cv::Mat flow;
  getOpticalFlow(currentFrame,previousFrame,flow);
  float sum_a, sum_b, sum_L, new_j,new_i;
  new_i=0.0;
  new_j=0.0;
  int c,r;
  c=flow.cols;
  r=flow.rows;
  int type = flow.type();
  for(int j=0; j < DFrame.rows;j++)
  {
    for(int i=0;i<DFrame.cols;i++)
    {
      new_j = j + flow.at<cv::Vec2f>(j,i)[1];
      new_i = i + flow.at<cv::Vec2f>(j,i)[0];
      sum_a = currentFrame.at<cv::Vec3f>(new_j,new_i)[1] - previousFrame.at<cv::Vec3f>(j,i)[1];
      sum_b = currentFrame.at<cv::Vec3f>(new_j,new_i)[2] - previousFrame.at<cv::Vec3f>(j,i)[2];
      sum_L = currentFrame.at<cv::Vec3f>(j,i)[0] - previousFrame.at<cv::Vec3f>(j,i)[0];
      DFrame.at<float>(j,i) = sum_a + sum_L +sum_b + sum_L;
      
    }
  }
  
  return 0;
}
/**
 * Should get the LPD frame but i dintnt get the minimalization right, my mistake
 */
int TMOTao17::getLowProximityFrame(cv::Mat currentFrame, cv::Mat previousFrame, cv::Mat previousGray, cv::Mat &currentGray)
{
  cv::Mat smallCurrent, smallPrevious,curRGB,smallGray;
  cv::Size size(SUBSAMPLE,SUBSAMPLE);
	  float wrt,wgt,wbt,wr,wg,wb;
	cv::resize(currentFrame,smallCurrent,size);
	cv::resize(previousFrame,smallPrevious,size);
	cv::resize(previousGray,smallGray,size);
	
	cv::cvtColor(currentFrame, curRGB, cv::COLOR_Lab2RGB);
	
	
	 labVector.clear();
	 labVector2.clear();
	 rgbVec.clear();
	 grayVec.clear();
  
	for (int j = 0; j < SUBSAMPLE; j++)   
	{
	
	    for (int i = 0; i < SUBSAMPLE; i++) 
	    {
	      labVector.push_back(smallCurrent.at<cv::Vec3f>(j,i));
	      labVector2.push_back(smallPrevious.at<cv::Vec3f>(j,i));
	      rgbVec.push_back(curRGB.at<cv::Vec3f>(j,i));
	      grayVec.push_back(smallGray.at<float>(j,i)); 
	      
	    }
	}
	getWeights(wr,wg,wb);
	getWeightsTemp(wrt,wgt,wbt);
	wr = (beta * wrt) + ((1.0-beta) * wr);
	wg = (beta * wgt) + ((1.0-beta) * wg);
	wb = (beta * wbt) + ((1.0-beta) * wb);
	
	for (int j = 0; j < currentFrame.rows; j++)
	{
	  for (int i = 0; i < currentFrame.cols; i++) 
	  {
	  
	    currentGray.at<float>(j,i) = curRGB.at<cv::Vec3f>(j,i)[0] * wr + curRGB.at<cv::Vec3f>(j,i)[1] * wg + curRGB.at<cv::Vec3f>(j,i)[2] * wb;
	  }
	}
	
	return 0;
}
/*
 * Should have worked like this get all the proximitz vals, compute count, mean, len complete
 * the algorithm from the artice and done, but i dont the value os lamda and i dont know how L shoulde converge
 * I have tried to make it work but it didnt.
 */
int TMOTao17::classifier(std::vector<cv::Point2f> &proximityVals)
{
      cv::Mat labels;
    cv::Mat1f colors;
    cv::kmeans(proximityVals, 3, labels, cv::TermCriteria(), 3, cv::KMEANS_PP_CENTERS, colors);
    cv::Point2f mean0,mean1,mean2,var0,var1,var2;
    int count0=0;
    int count1=0;
    int count2 =0;
    float len0,len1,len2;


    for (int i = 0; i < proximityVals.size(); ++i)
    {
      if(labels.at<int>(i) == 0)
      {
	mean0.x +=proximityVals[i].x;
	mean0.y +=proximityVals[i].y;
	count0++;
	
      }
      if(labels.at<int>(i) == 1)
      {
	mean1.x +=proximityVals[i].x;
	mean1.y +=proximityVals[i].y;
	count1++;
	
      }
      if(labels.at<int>(i) == 2)
      {
	mean2.x +=proximityVals[i].x;
	mean2.y +=proximityVals[i].y;
	count2++;
	
      }
       std::cerr<<labels.at<int>(i)<< " "<< proximityVals[i]<<std::endl;
    }
    mean0=mean0/(float)count0;
    mean1=mean1/(float)count1;
    mean2=mean2/(float)count2;
    
    len0 = count0/(float)proximityVals.size();
    len1 = count1/(float)proximityVals.size();
    len2 = count2/(float)proximityVals.size();
    for (int i = 0; i < proximityVals.size(); ++i)
    {
      if(labels.at<int>(i) == 0)
      {
	var0.x += std::pow(proximityVals[i].x - mean0.x, 2);
	
      }
      if(labels.at<int>(i) == 1)
      {
	var1.x += std::pow(proximityVals[i].x - mean1.x, 2);
	var1.x += std::pow(proximityVals[i].x - mean1.x, 2);
	
      }
      if(labels.at<int>(i) == 2)
      {
	var2.x += std::pow(proximityVals[i].x - mean2.x, 2);
	var2.y += std::pow(proximityVals[i].y - mean2.y, 2);
	
      }
       
    }
}

/**
 * Get the proximity values
 * @diffMat matrix containig the result of curret frame - previous frame
 * @proximityVlas firts number is L proximity, second is chorimonance proximity
 * */
int TMOTao17::getProximityValues(cv::Mat diffMat,std::vector<cv::Point2f> &proximityVals)
{
  float Lsigma=0.0;
  float asigma=0.0;
  float bsigma =0.0;
  float csigma =0.0;
  cv::Mat Lab[3];
  cv::normalize(diffMat,diffMat,0.0,255.0,cv::NORM_MINMAX,CV_32F);
  cv::split(diffMat,Lab);
  int bins =256;
  float range[]={0,256};
  cv::Point point;
  const float* ranges[]={range};
  int histSize[]={bins};
  cv::Mat Lhist,ahist,bhist;
  int channelss[]={0,1,2};
  
  cv::calcHist( &Lab[0], 1, 0, cv::Mat(), Lhist, 1, histSize, ranges, true, false );
  cv::calcHist( &Lab[1], 1, 0, cv::Mat(), ahist, 1, histSize, ranges, true, false );
  cv::calcHist( &Lab[2], 1, 0, cv::Mat(), bhist, 1, histSize, ranges, true, false );
  
  
  
  for(int i=1; i<Lhist.rows;i++)
  {
    float val =Lhist.at<float>(i);
    float val2 =ahist.at<float>(i);
    float val3 =bhist.at<float>(i);
    if(val != 0)
    {
      Lsigma = Lsigma - val * std::log(val);
    }
    if(val2 != 0)
    {
       asigma = asigma - val2 * std::log(val2);
    }
    if(val3 != 0)
    {
       bsigma = bsigma - val3 * std::log(val3);
    }
  }
  csigma = std::sqrt((std::pow(asigma,2)+std::pow(bsigma,2))/2.0);
  
  point.x=std::abs(Lsigma);
  point.y=csigma;
  proximityVals.push_back(point);
  
  return 0;

}

int TMOTao17::TransformVideo()
{
    int height = vSrc->GetHeight();
    int width = vSrc->GetWidth();
    cv::VideoCapture cap=vSrc->getVideoCaptureObject();
    std::vector<cv::Mat> channels;
    double min, max;
    cv::Mat currentFrame,previousFrame, diff, DFrame, CFrame, previousGray,currentGray, result,small,smallRGB;
    DFrame = cv::Mat::zeros(height, width, CV_32F);
    CFrame = cv::Mat::zeros(height, width, CV_32F);
    previousGray = cv::Mat::zeros(height, width, CV_32F);
    currentGray = cv::Mat::zeros(height, width, CV_32F);
 
    int chosenMethod=0;
    
    float wr,wg,wb;
    
    std::vector<float> cohVector;
    
    std::vector<cv::Point2f> proximityVals; /// fisrt proximity of L. second a, thrid b, fourth c which is an speciality, see method doc
    
    for(int frameCounter = 0; frameCounter<vSrc->GetTotalNumberOfFrames(); frameCounter++)
    {
      vSrc->GetMatVideoFrame(cap,frameCounter,currentFrame); //read frame
      std::cerr<<frameCounter+1 <<" / "<<vSrc->GetTotalNumberOfFrames()<<std::endl;
      
      if(frameCounter >0)
      {
	
	cv::cvtColor(currentFrame,currentFrame,CV_RGB2Lab);  //lab conversion
	cv::absdiff(currentFrame,previousFrame,diff);
	
	getProximityValues(diff,proximityVals);
	
	/*chosenMethod=classifier(proximityVals);  prepared to work
	
	if(chosenMethod ==0)
	{
	  getLowProximityFrame(currentFrame,previousFrame,previousGray,currentGray);
	}
	else if(chosenMethod ==1)
	{
	  getCoherenceRefinementFrame(currentFrame,previousFrame,previousGray,CFrame);
	  currentGray = previousGray + CFrame;
	}
	else if(chosenMethod ==2)
	{
	  getDifferentialRefinementFrame(currentFrame,previousFrame,DFrame);
	  currentGray = previousGray + DFrame;
	}*/
	
	
	//now it is set tu use only LPD
	getLowProximityFrame(currentFrame,previousFrame,previousGray,currentGray);
	currentGray.convertTo(currentGray, CV_32F, 1.0);
	

	cv::cvtColor(currentFrame,currentFrame,CV_Lab2LRGB); 
	for (int j = 0; j < height; j++)
	{
	
	  
	  for (int i = 0; i < width; i++) 
	  {
	  
	    currentGray.at<float>(j,i) = currentFrame.at<cv::Vec3f>(j,i)[0] * wr + currentFrame.at<cv::Vec3f>(j,i)[1] * wg + currentFrame.at<cv::Vec3f>(j,i)[2] * wb;
	    ///gets curent gray by weights
	  }
	}
	cv::cvtColor(currentFrame,currentFrame,CV_RGB2Lab);
	
      }
      else  ///whet it is the first frame
      {
	////same aso Transform go look thre
	cv::Size size(SUBSAMPLE,SUBSAMPLE);
	
	cv::resize(currentFrame,smallRGB,size);
	
	cv::cvtColor(smallRGB, small, cv::COLOR_RGB2Lab);
	
	  
  
	for (int j = 0; j < SUBSAMPLE; j++)
	{
	
	    for (int i = 0; i < SUBSAMPLE; i++) 
	    {
	      labVector.push_back(small.at<cv::Vec3f>(j,i));
	      rgbVec.push_back(currentFrame.at<cv::Vec3f>(j,i));
	    }
	}
	getWeights(wr,wg,wb);
	for (int j = 0; j < height; j++)
	{
	
	  
	  for (int i = 0; i < width; i++) 
	  {
	  
	    currentGray.at<float>(j,i) = currentFrame.at<cv::Vec3f>(j,i)[0] * wr + currentFrame.at<cv::Vec3f>(j,i)[1] * wg + currentFrame.at<cv::Vec3f>(j,i)[2] * wb;
	
	  }
	}
	
	cv::normalize(currentGray,currentGray,0.0,1.0,cv::NORM_MINMAX,CV_32F);
	
	cv::cvtColor(currentFrame,currentFrame,CV_RGB2Lab);
      }
      
      
      previousFrame=currentFrame;  //set frames curret bececomes previous
      previousGray=currentGray;
      channels.clear();
       channels.push_back(currentGray);
      channels.push_back(currentGray);
      channels.push_back(currentGray);
      cv::merge(channels,result);
      vDst->setMatFrame(vDst->getVideoWriterObject(),result);
      
    }
    


}




