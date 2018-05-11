/* --------------------------------------------------------------------------- *
 * TMOKuhn08.cpp: implementation of the TMOKuhn08 class.   *
 * --------------------------------------------------------------------------- */

#include "TMOKuhn08.h"
#include "particles.h"
#include <time.h>    

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOKuhn08::TMOKuhn08()
{
	SetName(L"Kuhn08");						
	SetDescription(L"Add your TMO description here");	// TODO - Insert description

	iK.SetName(L"iK");				
	iK.SetDescription(L"Cluster count");
	iK.SetDefault(15);							
	iK=15;
	iK.SetRange(1,255);				
	this->Register(iK);

	iAttempts.SetName(L"iAttempts");				
	iAttempts.SetDescription(L"Cluster attempts");
	iAttempts.SetDefault(5);							
	iAttempts=5;
	iAttempts.SetRange(1,20);				
	this->Register(iAttempts);

	bColor.SetName(L"bColor");				
	bColor.SetDescription(L"Enable color");
	bColor.SetDefault(false);							
	bColor=false;
	this->Register(bColor);

	bChrominance.SetName(L"bChrominance");				
	bChrominance.SetDescription(L"Only chrominance");
	bChrominance.SetDefault(false);							
	bChrominance=false;
	this->Register(bChrominance);

	iStep.SetName(L"iStep");				
	iStep.SetDescription(L"Maximum steps for optimizing gray");
	iStep.SetDefault(10);							
	iStep=10;
	iStep.SetRange(0,2000);				
	this->Register(iStep);

	bOriginal.SetName(L"bOriginal");				
	bOriginal.SetDescription(L"Nothink done only display input");
	bOriginal.SetDefault(false);							
	bOriginal=false;
	this->Register(bOriginal);


	bInterpolate.SetName(L"bInterpolate");				
	bInterpolate.SetDescription(L"Interpolate final gray image.");
	bInterpolate.SetDefault(true);							
	bInterpolate=true;
	this->Register(bInterpolate);
		
}

TMOKuhn08::~TMOKuhn08()
{
}

/* ------------------------------- *
 * Convert to OpenCV matrix format *
 * ------------------------------- */
cv::Mat TMOKuhn08::convertToMat(TMOImage* image)
{
	double* data = image->GetData();
	int rows = image->GetWidth();
	int cols = image->GetHeight();
	double pL, pa, pb;
	cv::Mat mat(rows * cols, 3, CV_32F);
  	for( int y = 0; y < rows; y++ )
	{
    	for( int x = 0; x < cols; x++ )
		{
			pL = *data++;
			pa = *data++;
			pb = *data++;

			mat.at<float>(y + x*rows, 0) = pL;
			mat.at<float>(y + x*rows, 1) = pa;
			mat.at<float>(y + x*rows, 2) = pb;
		}
	}
	return mat;
}
void TMOKuhn08::copy(TMOImage* in,TMOImage* out)
{
	double* pSourceData = in->GetData();				// You can work at low level data
	double* pDestinationData = out->GetData();			// Data are stored in form of array 
														// of three doubles representing
														// three colour components
	int rows = pSrc->GetWidth();
	int cols = pSrc->GetHeight();

	double pL, pa, pb;
		for (int j = 0; j < cols; j++)
		{
			for (int i = 0; i < rows; i++)
			{
				pL = *pSourceData++;
				pa = *pSourceData++;
				pb = *pSourceData++;
				*pDestinationData++ = pL;
				*pDestinationData++ = pa;
				*pDestinationData++ = pb;
			}
		}
}
/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
int TMOKuhn08::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can 
	// convert it into other format
	pSrc->Convert(TMO_LAB);								// This is format of Y as luminance
	pDst->Convert(TMO_LAB);								// x, y as color information


	int rows = pSrc->GetWidth();
	int cols = pSrc->GetHeight();

	double pL, pa, pb;
	if(bOriginal)
	{
		copy(pSrc, pDst);
	}
	else
	{
		  /* initialize random seed: */
  		srand (time(NULL));
		cv::Mat samples = convertToMat(pSrc);
		int clusterCount = iK;
		cv::Mat labels;
		int attempts = iAttempts;
		cv::Mat centers;
		cv::kmeans(samples, clusterCount, labels, cv::TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, cv::KMEANS_PP_CENTERS, centers );


		ParticlesPool particles(rows,cols,bColor,bChrominance);
		particles.initialize(labels,centers);
		double* pDestinationData = pDst->GetData();	

		if(!bColor)
		{
			particles.computeMaxDistance();
			particles.createSprings();
			particles.compute(iStep);
		

			if(bInterpolate)
			{
				double* pSourceData = pSrc->GetData();
				
				for( int y = 0; y < rows; y++ )
				{
					for( int x = 0; x < cols; x++ )
					{ 
						if(x == 0 && y == 0)
						{
							//std::cerr << "*" <<std::endl;
						}
						Particle source;

						source.L = *pSourceData++;
						source.a = *pSourceData++;
						source.b = *pSourceData++;

						double out = 0;
						Particle* p = particles.pool[y][x];
						if(source.L >= p->L)
						{
							
							out = p->gray + particles.getSkFactor(p) *  particles.getDistance(p,&source);
							if(x == 0 && y == 0)
						{

							/*std::cerr << "# " <<particles.getDistance(p,&source) << " s:" << source.L <<" " 
							<< source.a <<  "" << source.b << " p:" << p->L  <<" " << p->a << " " << p->b   <<std::endl;*/
						}
						}
						else
						{
							if(x == 0 && y == 0)
						{
							//std::cerr << "-" <<std::endl;
						}
							out = p->gray - particles.getSkFactor(p) * particles.getDistance(p,&source);
						}

						*pDestinationData++ = out;
						*pDestinationData++ = 0;
						*pDestinationData++ = 0;
					}
				}
			}
			else
			{
				pDestinationData = pDst->GetData();			// Data are stored in form of array 

				particles.toDestination(pDestinationData);
			}
		}
		else
		{
			pDestinationData = pDst->GetData();			// Data are stored in form of array 

			particles.toDestination(pDestinationData);
		}
	}
		
	pDst->Convert(TMO_RGB);

	return 0;
}

