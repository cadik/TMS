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
	SetDescription(L"Add your TMO description here"); // TODO - Insert description

	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed
	this->Register(dParameter);
}

#define rgb2luminance(R,G,B) (R*0.2126 + G*0.7152 + B*0.0722)
typedef vector< vector<double> > PixelMatrix;
struct Point{
   unsigned short x,y;
};
typedef std::vector<Point> PointVector;
struct Block_Record{
   PointVector Memebers;
   double Sum;
   unsigned int Count;
};

typedef std::vector<Block_Record > Blocks;

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
   cv::Mat image = cv::imread("/home/matthewlele/images/cadik/cadik-desk01.hdr");
   cv::Mat contours;
   cv::Mat gray_img;

   //cvtColor(image, gray_img, CV_RGB2GRAY);
   cv::Canny(image, contours, 80, 240);
   cv::namedWindow("Canny");
   //cv::imshow("Canny",contours);
   //cv::waitKey(0);
   
   double stonits = pSrc->GetStonits();
   PixelMatrix LogLuminancePixels(imageHeight, vector<double>(imageWidth, 0.0));
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
                     small_block.Memebers.push_back(p);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
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
                     small_block.Memebers.push_back(p);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
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
                     small_block.Memebers.push_back(p);
                  }
               }
               pixelBlocks.push_back(small_block);
               x_value+=2;
            }
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
                     small_block.Memebers.push_back(p);
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
   int pixelCount = 0;
   for(int i =0; i < pixelBlocks.size(); i++)
   {
      fprintf(stderr,"Block num: %d with pixels: %d with sum luminance: %g\n",i,pixelBlocks[i].Memebers.size(),pixelBlocks[i].Sum);
      pixelCount += pixelBlocks[i].Memebers.size();
   }
   fprintf(stderr,"Amount of blocks %d and pixels: %d\n",pixelBlocks.size(), pixelCount);

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
			pY *= dParameter; // Parameters can be used like
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
