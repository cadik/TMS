/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                                   *
*                                                                              *
*                       Diploma thesis                                         *
*                       Author: Matus Bicanovsky                               *
*                       Brno 2024                                              *
*                                                                              *
*                       Implementation of the TMOSong14 class                  *
*                Real-time video decolorization using bilateral filtering      *
*******************************************************************************/
/* --------------------------------------------------------------------------- *
 * TMOSong14.cpp: implementation of the TMOSong14 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOSong14.h"

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOSong14::TMOSong14()
{
	SetName(L"Song14");					  
	SetDescription(L"Color to grayscale operator for images and video, method from paper Real-time video decolorization using bilateral filtering"); 
}

TMOSong14::~TMOSong14()
{
}

//main function for processing image and video
int TMOSong14::Transform()
{
	int height = pSrc->GetHeight();
	int width = pSrc->GetWidth();
	double *pSourceData = pSrc->GetData();		
	double *pDestinationData = pDst->GetData(); 
										
				
	//create the copy of input image (I) with just luminance information (L)
	cv::Mat I(height, width, CV_64FC3, pSourceData);
	cv::Mat L(height, width, CV_64F);
	for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            cv::Vec3d pixel = I.at<cv::Vec3d>(j, i);
            L.at<double>(j, i) = 0.299 * pixel[2] + 0.587 * pixel[1] + 0.114 * pixel[0];
        }
    }
	//convert the (I) and (L) to CV_32F format for bilateral filtering function by OpenCV
	I.convertTo(I, CV_32F);
	L.convertTo(L, CV_32F);
	//apply bilateral filter to input (I) store the result in I_bilateral
	cv::Mat I_bilateral;
	double sigma_s = 2.0;                           //spatial parameter according to the paper
	double sigma_r = 0.02;                          //range parameter according to the paper                      
	cv::bilateralFilter(I, I_bilateral, -1, sigma_r, sigma_s);

	// apply joint bilateral filter to input with luminance image as the guidance image
	cv::Mat I_joint_bilateral;
	cv::ximgproc::jointBilateralFilter(L, I, I_joint_bilateral, -1, sigma_r, sigma_s);

	//convert the original,luminance and filtered images back to CV_64F (doubles) for better precision
	I_bilateral.convertTo(I_bilateral, CV_64F);
	I_joint_bilateral.convertTo(I_joint_bilateral, CV_64F);
	I.convertTo(I, CV_64F);
	L.convertTo(L, CV_64F);
	
	//estimate the detail loss (D), D = I_bilateral - I_joint_bilateral
	cv::Mat D = I_bilateral - I_joint_bilateral;
	fprintf(stderr, "Detail loss estimated\n");
	
	//compute gradients of the detail loss (D) and the input image (I) using Sobel operator
    std::vector<cv::Mat> grad_Dx_channels(3), grad_Dy_channels(3);               //vector of cv::Mat for gradient of D in x and y direction for each channel
    std::vector<cv::Mat> grad_Ix_channels(3), grad_Iy_channels(3);               //vector of cv::Mat for gradient of I in x and y direction for each channel
    for (int c = 0; c < 3; c++) {
        cv::Mat D_channel, I_channel;
        cv::extractChannel(D, D_channel, c);
        cv::extractChannel(I, I_channel, c);
        cv::Sobel(D_channel, grad_Dx_channels[c], CV_64F, 1, 0, 3);         //compute gradient of D in x direction for given channel
        cv::Sobel(D_channel, grad_Dy_channels[c], CV_64F, 0, 1, 3);         //compute gradient of D in y 
        cv::Sobel(I_channel, grad_Ix_channels[c], CV_64F, 1, 0, 3);         //compute gradient of I in x 
        cv::Sobel(I_channel, grad_Iy_channels[c], CV_64F, 0, 1, 3);         //compute gradient of I in y 
    }
	//construct matrices A and B using the calculated gradients
    cv::Mat A(height * width, 3, CV_64F); //each row is a 3D gradient vector
    cv::Mat B(height * width, 1, CV_64F); //each row is a scalar value
    cv::Mat x(3, 1, CV_64F); // 3x1 vector

    int idx = 0;
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
			//gradient of D at pixel (j, i) for each channel
            cv::Vec3d grad_D(
                grad_Dx_channels[0].at<double>(j, i) + grad_Dy_channels[0].at<double>(j, i),
                grad_Dx_channels[1].at<double>(j, i) + grad_Dy_channels[1].at<double>(j, i),
                grad_Dx_channels[2].at<double>(j, i) + grad_Dy_channels[2].at<double>(j, i)
            );
            A.at<double>(idx, 0) = grad_D[0];               //gradient of D for red channel
            A.at<double>(idx, 1) = grad_D[1];               //gradient of D for green channel
            A.at<double>(idx, 2) = grad_D[2];               //gradient of D for blue channel

            //find the channel with the largest amount of contrast loss
            int tmp = 0;
            double max_contrast_loss = std::abs(D.at<cv::Vec3d>(j, i)[0]);
            for (int c = 1; c < 3; c++) {
                double contrast_loss = std::abs(D.at<cv::Vec3d>(j, i)[c]);
                if (contrast_loss > max_contrast_loss) {
                    max_contrast_loss = contrast_loss;
                    tmp = c;
                }
            }

            //set B to the gradient of I at the channel with the largest contrast loss
            B.at<double>(idx, 0) = grad_Ix_channels[tmp].at<double>(j, i) + grad_Iy_channels[tmp].at<double>(j, i);
            idx++;
        }
    }
	fprintf(stderr, "Linear system constructed\n");

	//compute the gradient of the luminance image (L) and subtract it from matrix B
    cv::Mat grad_Lx, grad_Ly;
    cv::Sobel(L, grad_Lx, CV_64F, 1, 0, 3);
    cv::Sobel(L, grad_Ly, CV_64F, 0, 1, 3);
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            B.at<double>(j * width + i, 0) -= grad_Lx.at<double>(j, i) + grad_Ly.at<double>(j, i);   //computing (B - ∇L) from equation (8)
        }
    }
	//solve the mapping function
    cv::Mat AtA = A.t() * A;    //calculate A^TA
    cv::Mat AtB = A.t() * B;    //calculate A^T(B - ∇L), matrix B is already in the form (B - ∇L)
	cv::Mat AtA_inv;            
	cv::invert(AtA, AtA_inv, cv::DECOMP_SVD);   //inverse of A^TA
    
	x = AtA_inv * AtB;        //equation to solve x = (A^TA)^-1 * A^T(B - ∇L), equation (8) in the paper

    fprintf(stderr, "Mapping function computed\n");

    //calculate residual image (R) using the mapping function (x) and the detail loss (D), R(x,y) = D(x,y) * x
    cv::Mat R(height, width, CV_64F);
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            cv::Vec3d d_pixel = D.at<cv::Vec3d>(j, i);
            R.at<double>(j, i) = d_pixel[0] * x.at<double>(0, 0) + d_pixel[1] * x.at<double>(1, 0) + d_pixel[2] * x.at<double>(2, 0);                   //dot product to map detail loss to residual image
        }
    }

    fprintf(stderr, "Residual image computed\n");
	//reconstruct the final image as sum of residual image(R) and the luminance image(L)
	cv::Mat G = R + L;
    int j;
	for (j = 0; j < pSrc->GetHeight(); j++)
	{
		pSrc->ProgressBar(j, pSrc->GetHeight()); // You can provide progress bar
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			double pixel = G.at<double>(j, i);
            *pDestinationData++ = pixel;
            *pDestinationData++ = pixel;
            *pDestinationData++ = pixel;
		}
	}
	pSrc->ProgressBar(j, pSrc->GetHeight());
	//pDst->Convert(TMO_RGB);
	return 0;
}
