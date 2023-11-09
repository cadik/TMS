/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                  *
*                                                                              *
*                       VYF -- term project                                    *
*                       Author: Marek Hlavacka                                 *
*                       Brno 2023                                              *
*                                                                              *
*                       Color to Gray and Back: Color Embedding                *
*                       Into Textured Gray Images                              *
*                       https://ieeexplore.ieee.org/document/1632200           *                                       
********************************************************************************/
/**
 * @file TMOQueiroz06.cpp
 * @brief Color to Gray and Back: Color Embedding Into Textured Gray images
 * @author Marek Hlavacka
 * @class TMOQueiroz06
 */


#include "TMOQueiroz06.h"



/**
  *  @brief Constructor
  */
TMOQueiroz06::TMOQueiroz06()
{
	SetName(L"Queiroz06");					  // TODO - Insert operator name
	SetDescription(L"Add your TMO description here"); // TODO - Insert description
/*
	dParameter.SetName(L"ParameterName");				// TODO - Insert parameters names
	dParameter.SetDescription(L"ParameterDescription"); // TODO - Insert parameter descriptions
	dParameter.SetDefault(1);							// TODO - Add default values
	dParameter = 1.;
	dParameter.SetRange(-1000.0, 1000.0); // TODO - Add acceptable range if needed*/



   dVariant.SetName(L"Variant");
   dVariant.SetDescription(L"Variant of transform. 0 - variant with negativ and positive Cr, Cb");
   dVariant.SetDefault(0);
   dVariant = 0;
   dVariant.SetRange(0,1);

   dSavePartialResults.SetName(L"Save partial results (debug)");
   dSavePartialResults.SetDescription(L"Save partial results of transformations");
   dSavePartialResults.SetDefault(0);
   dSavePartialResults = 0;
   dSavePartialResults.SetRange(0,1);

	//this->Register(dParameter);
   this->Register(dVariant);
   this->Register(dSavePartialResults);
}

/*
@brief destructor
*/
TMOQueiroz06::~TMOQueiroz06()
{
}

/**
 * @brief transformation function
 * @return exit code
 */
int TMOQueiroz06::Transform()
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
	double r, g, b, Y, Cr, Cb;
   int cols = pSrc->GetWidth();
	int rows = pSrc->GetHeight();

   //convert src low level data to cv::Mat with 3 channels
   cv::Mat img = srcToMat(pSrc);  
   std::vector<cv::Mat> channels;
   cv::split(img,channels);
   std::vector<cv::Mat> channelsResized;

   //Resize channels for square 
   for (int c = 0; c < channels.size();c++){
      channelsResized.push_back(resizeToSquare(channels[c]));
      }

   //print resized Y, Cb, Cr channels
   if(this->dSavePartialResults){
         PrintImage(channelsResized[0],"Y");
         PrintImage(channelsResized[1],"Cb");
         PrintImage(channelsResized[2],"Cr");
      } 
   
   //First variant of Color to Gray and Back: Color Embedding Into Textured Gray images with negativ and positiv Cr, Cb channels
   if (this->dVariant == 0){
      cv::Mat CrMinus, CrPlus, CbMinus, CbPlus, ResizedCrMinus, ResizedCrPlus, ResizedCbMinus, ResizedCbPlus, wavelethY, ResizedCb, ResizedCr, ImageForPrint;
      
      //Downsample Cr and Cb channels to 1/2 of origin image
      cv::pyrDown(channelsResized[1], ResizedCb, cv::Size(channelsResized[1].size().height * 0.5,channelsResized[1].size().width * 0.5));
      cv::pyrDown(channelsResized[2], ResizedCr, cv::Size(channelsResized[2].size().height * 0.5,channelsResized[1].size().width * 0.5)); 

      //Create Cr+, Cr-, Cb+ and Cb-
      CrPlus  = ResizedCr.clone().setTo(0,(ResizedCr < 0.0));
      CrMinus =  ResizedCr.clone().setTo(0,(ResizedCr > 0.0));
      CbPlus  =  ResizedCb.clone().setTo(0,(ResizedCb < 0.0));
      CbMinus = ResizedCb.clone().setTo(0,(ResizedCb > 0.0));

      //Downsample Cb- channel to 1/4 of origin image
      cv::pyrDown(CbMinus, ResizedCbMinus, cv::Size(CbMinus.size().height * 0.5,CbMinus.size().width * 0.5)); 
      
      wavelethY = channelsResized[0].clone();
      HaarWavelet(channelsResized[0],wavelethY,2); 
      
      //Copying Cb-, Cb+, Cr+, Cr- to DWT result
      CrPlus.copyTo(wavelethY(cv::Rect(0,wavelethY.cols/2,CrPlus.rows,CrPlus.cols)));
      CbPlus.copyTo(wavelethY(cv::Rect(wavelethY.rows/2,0,CbPlus.rows,CbPlus.cols)));
      CrMinus.copyTo(wavelethY(cv::Rect(wavelethY.rows/2,wavelethY.cols/2,CbPlus.rows,CbPlus.cols)));
      ResizedCbMinus.copyTo(wavelethY(cv::Rect(wavelethY.rows/4,wavelethY.cols/4,ResizedCbMinus.rows,ResizedCbMinus.cols)));
      ImageForPrint = wavelethY.clone();
        
      if (this->dSavePartialResults){
         PrintImage(wavelethY,"DWT_with_Cr_Cb");
      }

      //Gray image to print
      InvHaarWavelet(wavelethY,ImageForPrint,2);
      
      if (this->dSavePartialResults){
         PrintImage(ImageForPrint,"ImageForPrint");
      }

      //Recovery Image colors
      cv::Mat ScannedImage = ImageForPrint.clone();
      
      HaarWavelet(ScannedImage,wavelethY,2);
      if (this->dSavePartialResults){
         PrintImage(wavelethY,"RecoveryImage_DWT");
      }

      //Copying parts of DWT       
      ResizedCbMinus = wavelethY(cv::Rect(wavelethY.rows/4, wavelethY.cols/4, wavelethY.rows/4, wavelethY.cols/4)).clone();
      cv::resize(ResizedCbMinus,CbMinus,cv::Size(wavelethY.rows/2,wavelethY.cols/2));
      CrPlus = wavelethY(cv::Rect(0,wavelethY.cols/2, wavelethY.cols/2,wavelethY.cols/2)).clone();
      CbPlus = wavelethY(cv::Rect(wavelethY.cols/2,0, wavelethY.cols/2,wavelethY.cols/2)).clone();
      CrMinus = wavelethY(cv::Rect(wavelethY.cols/2,wavelethY.cols/2, wavelethY.cols/2,wavelethY.cols/2)).clone();
      
      //Create back Cb and Cr 
      ResizedCb = CbPlus + CbMinus;
      ResizedCr = CrPlus - CrMinus;

      //zeroing parts by Cb and Cr components
      cv::Mat zero1_4 = cv::Mat::zeros(CrPlus.rows, CrPlus.cols, CV_32F);
      cv::Mat zero1_8 = cv::Mat::zeros(ResizedCbMinus.rows, ResizedCbMinus.cols, CV_32F);
      zero1_4.copyTo(wavelethY(cv::Rect(0,wavelethY.cols/2, zero1_4.rows, zero1_4.cols)));
      zero1_4.copyTo(wavelethY(cv::Rect(wavelethY.rows/2,wavelethY.cols/2, zero1_4.rows, zero1_4.cols)));
      zero1_4.copyTo(wavelethY(cv::Rect(wavelethY.rows/2, 0, zero1_4.rows, zero1_4.cols)));
      zero1_8.copyTo(wavelethY(cv::Rect(wavelethY.rows/4, wavelethY.rows/4, zero1_8.rows, zero1_8.cols)));

      cv::Mat Cr_, Cb_, Y_, InvWavelethY;
      InvWavelethY = wavelethY.clone();
      InvHaarWavelet(wavelethY,InvWavelethY,2);
         
      if (this->dSavePartialResults){
         PrintImage(InvWavelethY,"RecoveryImage_IDWT");
      }
      //resize Image to original size and crop if is necessary
      cv::resize(ResizedCr,Cr_,cv::Size(channelsResized[0].rows,channelsResized[0].cols));
      cv::resize(ResizedCb,Cb_,cv::Size(channelsResized[0].rows,channelsResized[0].cols));
      Cb_ = Cb_(cv::Rect(0,0, cols, rows)).clone();
      Cr_ = Cr_(cv::Rect(0,0, cols, rows)).clone();
      Y_ = InvWavelethY(cv::Rect(0,0,cols, rows)).clone();
      if (this->dSavePartialResults){
         PrintImage(Cb_,"RecoveryImage_ResizedCb");
         PrintImage(Cr_,"RecoveryImage_ResizedCr");
         PrintImage(Y_ ,"RecoveryImage_Y");
      }

      float delta = 0.5;
      for (int x = 0; x < rows; x++)
         for (int y = 0; y < cols; y++){
            *pDestinationData++ = Y_.at<float>(x,y) + 1.403 * (Cr_.at<float>(x,y) - delta);
            *pDestinationData++ = Y_.at<float>(x,y) - 0.344 * (Cb_.at<float>(x,y) - delta) - 0.714 * (Cr_.at<float>(x,y) - delta);
            *pDestinationData++ = Y_.at<float>(x,y) + 1.773 * (Cb_.at<float>(x,y) - delta);
         }
   
   }
   //less robust variant
   else if (this->dVariant == 1){
      cv::Mat ResizedCr, ResizedCb, Cr_, Cb_, Y_;
      //downsample Cr,Cb to 1/2 of Image
      cv::pyrDown(channelsResized[1], ResizedCb, cv::Size(channelsResized[1].size().height * 0.5,channelsResized[1].size().width * 0.5));
      cv::pyrDown(channelsResized[2], ResizedCr, cv::Size(channelsResized[2].size().height * 0.5,channelsResized[1].size().width * 0.5)); 

      cv::Mat wavelethY = cv::Mat::zeros(channelsResized[0].rows,channelsResized[0].cols,CV_32F);
      //DWT
      HaarWavelet(channelsResized[0],wavelethY,2);

      //Copy Cr and Cb to DWT result 
      ResizedCr.copyTo(wavelethY(cv::Rect(0,wavelethY.cols/2,ResizedCr.rows,ResizedCr.cols)));
      ResizedCb.copyTo(wavelethY(cv::Rect(wavelethY.rows/2,0,ResizedCr.rows,ResizedCr.cols)));
      if (this->dSavePartialResults){
         PrintImage(wavelethY,"DWT_with_Cr_Cb");
      }

      //Gray image to print
      cv::Mat InvWavelethY = cv::Mat::zeros(channelsResized[0].rows,channelsResized[0].cols,CV_32F); 
      InvHaarWavelet(wavelethY,InvWavelethY,2);
      if (this->dSavePartialResults){
         PrintImage(InvWavelethY,"ImageForPrint");
      }
      
      //Recovery Image colors
      cv::Mat revoveryWav = cv::Mat::zeros(channelsResized[0].rows,channelsResized[0].cols,CV_32F);
      HaarWavelet(InvWavelethY,revoveryWav,2);
      if (this->dSavePartialResults){
         PrintImage(revoveryWav,"RecoveryImage_DWT");
      }

      //Get Cr and Cb from DWT results 
      ResizedCr = revoveryWav(cv::Rect(0,wavelethY.cols/2, wavelethY.cols/2,wavelethY.cols/2)).clone();
      ResizedCb = revoveryWav(cv::Rect(wavelethY.cols/2,0, wavelethY.cols/2,wavelethY.cols/2)).clone();
      cv::Mat zero = cv::Mat::zeros(wavelethY.cols/2, wavelethY.cols/2, CV_32F);
      
      //zeroing parts by Cb and Cr components
      zero.copyTo(revoveryWav(cv::Rect(0,wavelethY.cols/2,ResizedCr.rows,ResizedCr.cols)));
      zero.copyTo(revoveryWav(cv::Rect(wavelethY.rows/2,0,ResizedCr.rows,ResizedCr.cols)));
      
      //IDWT
      InvHaarWavelet(revoveryWav,InvWavelethY,2);

      //resize Image to original size and crop if is necessary
      cv::resize(ResizedCr,Cr_,cv::Size(channelsResized[0].rows,channelsResized[0].cols));
      cv::resize(ResizedCb,Cb_,cv::Size(channelsResized[0].rows,channelsResized[0].cols));
      Cb_ = Cb_(cv::Rect(0,0, cols, rows)).clone();
      Cr_ = Cr_(cv::Rect(0,0, cols, rows)).clone();
      Y_ = InvWavelethY(cv::Rect(0,0,cols, rows)).clone();

      if (this->dSavePartialResults){
         PrintImage(Cb_,"RecoveryImage_ResizedCb");
         PrintImage(Cr_,"RecoveryImage_ResizedCr");
         PrintImage(Y_,"RecoveryImage_IDWT");
      }


      //r = Y + 1.403 * (Cr - 0.5)
      //g = Y - 0.344 * (Cb - 0.5) - 0.714 * (Cr - 0.5)
      //b = Y + 1.773 * (Cb - 0.5)
      float delta = 0.5;
      for (int x = 0; x < rows; x++)
         for (int y = 0; y < cols; y++){
            *pDestinationData++ = Y_.at<float>(x,y) + 1.403 * (Cr_.at<float>(x,y) - delta);
            *pDestinationData++ = Y_.at<float>(x,y) - 0.344 * (Cb_.at<float>(x,y) - delta) - 0.714 * (Cr_.at<float>(x,y) - delta);
            *pDestinationData++ = Y_.at<float>(x,y) + 1.773 * (Cb_.at<float>(x,y) - delta);
         }
      }
      //pSrc->ProgressBar(j, pSrc->GetHeight());
      //pDst->Convert(TMO_RGB);
	return 0;
}


/* @brief Convert DWT with Haar waveleth to IMG
 * 
 * @param src - input matrix
 * @param dst - output matrix
 * @int level - level of DWT
 **/
void TMOQueiroz06::InvHaarWavelet(cv::Mat &src,cv::Mat &dst,int level){
    float c,dh,dv,dd;
    for (int l=level;l>0;l--){
        for (int y=0;y<(src.rows>>l);y++){
            for (int x=0; x<(src.cols>>l);x++){
                c=src.at<float>(y,x);
                dh=src.at<float>(y,x+(src.cols>>l));
                dv=src.at<float>(y+(src.rows>>l),x);
                dd=src.at<float>(y+(src.rows>>l),x+(src.cols>>l));

                dst.at<float>(y*2,x*2)=0.5*(c+dh+dv+dd);
                dst.at<float>(y*2,x*2+1)=0.5*(c-dh+dv-dd);
                dst.at<float>(y*2+1,x*2)=0.5*(c+dh-dv-dd);
                dst.at<float>(y*2+1,x*2+1)=0.5*(c-dh-dv+dd);            
            }
        }
        cv::Mat C=src(cv::Rect(0,0,src.cols>>(l-1),src.rows>>(l-1)));
        cv::Mat D=dst(cv::Rect(0,0,src.cols>>(l-1),src.rows>>(l-1)));
        D.copyTo(C);
    }   
}


/* @brief Convert img to DWT with Haar waveleth
 * 
 * @param src - input matrix
 * @param dst - output matrix
 * @int level - level of DWT
 **/
void TMOQueiroz06::HaarWavelet(cv::Mat &src,cv::Mat &dst,int level){
    float c,dh,dv,dd;
    for (int l=0;l<level;l++){
        for (int y=0;y<(src.rows/pow(2,l+1));y++){
            for (int x=0; x<(src.cols/pow(2,l+1));x++){
               c=src.at<float>(2*y,2*x)+src.at<float>(2*y,2*x+1)+src.at<float>(2*y+1,2*x)+src.at<float>(2*y+1,2*x+1);
               dh=src.at<float>(2*y,2*x)+src.at<float>(2*y+1,2*x)-src.at<float>(2*y,2*x+1)-src.at<float>(2*y+1,2*x+1);
               dv=src.at<float>(2*y,2*x)+src.at<float>(2*y,2*x+1)-src.at<float>(2*y+1,2*x)-src.at<float>(2*y+1,2*x+1);
               dd=src.at<float>(2*y,2*x)-src.at<float>(2*y,2*x+1)-src.at<float>(2*y+1,2*x)+src.at<float>(2*y+1,2*x+1);

               dst.at<float>(y,x)=c*0.5;
               dst.at<float>(y,x+(src.cols/pow(2,l+1)))=dh*0.5;
               dst.at<float>(y+(src.rows/pow(2,l+1)),x)=dv*0.5;               
               dst.at<float>(y+(src.rows/pow(2,l+1)),x+(src.cols/pow(2,l+1)))=dd*0.5;
            }
        }
        dst.copyTo(src);
    }   
};


/* @brief Resize image from NxM to square NxN
 * 
 * @param input - input matrix
 * @param name  - name of img
 **/
void TMOQueiroz06::PrintImage(const cv::Mat input, std::string name){
   cv::Mat out = input.clone();
   out.convertTo(out, CV_8U, 255.0);
   auto ImgName = "./" + name + ".png";
   cv::imwrite(ImgName,out);

}

/* @brief Resize image from NxM to square NxN
 * 
 * @param input - input matrix
 * @return output matrix with NxN size
 **/
cv::Mat TMOQueiroz06::srcToMat(TMOImage *image)
{
	double *data = image->GetData();
	int cols = image->GetWidth();
	int rows = image->GetHeight();
   double r,g,b;
	cv::Mat Y = cv::Mat::zeros(rows,cols, CV_32F);
   cv::Mat Cr = cv::Mat::zeros(rows,cols, CV_32F);   
   cv::Mat Cb = cv::Mat::zeros(rows,cols, CV_32F);   
   double delta = 0.5;
	for (int x = 0; x < rows; x++)
	{
		for (int y = 0; y < cols; y++)
		{
			r = *data++;
			g = *data++;
			b = *data++;
			Y.at<float>(x,y) = 0.299 * r + 0.587 * g + 0.114 * b;
			Cr.at<float>(x,y) = 0.500 * r - 0.419 * g - 0.081 * b + delta; 
			Cb.at<float>(x,y) = -0.169 * r - 0.331 * g + 0.500 * b + delta;
         
		}
	}
   std::vector<cv::Mat> channels;
   channels.push_back(Y);
   channels.push_back(Cb);
   channels.push_back(Cr);
   cv::Mat fin;
   cv::merge(channels,fin);
   return fin;
}

/* @brief Resize image from NxM to square NxN
 * 
 * @param input - input matrix
 * @return output matrix with NxN size
 **/
cv::Mat TMOQueiroz06::resizeToSquare(const cv::Mat &input)
{
   int square_size = cv::max(input.cols, input.rows);
   cv::Mat output = cv::Mat::zeros(cv::Size(square_size, square_size), input.type());
   cv::Rect roi(0, 0, input.cols, input.rows);
   input.copyTo(output(roi));
   return output;
}

