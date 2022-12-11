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
typedef vector< vector<float> > PixelMatrix;
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
struct Block_Record{
   PointVector Memebers;
   float Sum;
   unsigned int Count;
   HistogramVector logHistogram;
   SignatureVector blockSignature;
   
};

typedef std::vector<Block_Record > Blocks;

float emdFunction(int firstBlock, int secondBlock, Blocks pixelBlocks)
{
   cv::Mat sign1(cv::Size(2,3),CV_32FC1);
   cv::Mat sign2(cv::Size(2,3),CV_32FC1);
   for(int i=0;i < 3; i++)
   {
      sign1.at<float>(i,0) = pixelBlocks[0].blockSignature[i].s;
      sign1.at<float>(i,1) = pixelBlocks[0].blockSignature[i].w;
   }
   for(int k=0; k < 3; k++)
   {
      sign2.at<float>(k,0) = pixelBlocks[1].blockSignature[k].s;
      sign2.at<float>(k,1) = pixelBlocks[1].blockSignature[k].w;
   } 
   return cv::EMD(sign1,sign2,cv::DIST_L1);
   
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

   float distance = emdFunction(0,1,pixelBlocks);
   fprintf(stderr,"Distance %g\n",distance);
   





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
