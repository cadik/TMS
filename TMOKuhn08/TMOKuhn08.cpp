/*******************************************************************************
*                                                                              *
*                       Brno University of Technology                          *
*                       CPhoto@FIT                                             *
*                                                                              *
*                       Tone Mapping Studio	                               *
*                                                                              *
*                       Brno 2018                                              *
*                                                                              *
*                       Implementation of the TMOKuhn08 class                  *
*                                                                              *
*******************************************************************************/
/**
 * @file TMOKuhn08.cpp
 * @brief Implementation of the TMOKuhn08 class
 * @class TMOKuhn08.cpp
 * 
 * @todo Insert description
 */

#include "TMOKuhn08.h"
#include "particles.h"
#include <time.h>

/**
  *  @brief Constructor
  */
TMOKuhn08::TMOKuhn08()
{
	SetName(L"Kuhn08");
	SetDescription(L"Add your TMO description here"); /** TODO - Insert description */

	iK.SetName(L"iK");
	iK.SetDescription(L"Cluster count");
	iK.SetDefault(15);
	iK = 15;
	iK.SetRange(1, 255);
	this->Register(iK);

	iAttempts.SetName(L"iAttempts");
	iAttempts.SetDescription(L"Cluster attempts");
	iAttempts.SetDefault(1);
	iAttempts = 1;
	iAttempts.SetRange(1, 20);
	this->Register(iAttempts);

	bColor.SetName(L"bColor");
	bColor.SetDescription(L"Enable color");
	bColor.SetDefault(false);
	bColor = false;
	this->Register(bColor);

	bChrominance.SetName(L"bChrominance");
	bChrominance.SetDescription(L"Only chrominance");
	bChrominance.SetDefault(false);
	bChrominance = false;
	this->Register(bChrominance);

	dTime.SetName(L"dTime");
	dTime.SetDescription(L"Maximum time for optimizing gray");
	dTime.SetDefault(1.0);
	dTime = 1.0;
	dTime.SetRange(0, 10000);
	this->Register(dTime);

	bOriginal.SetName(L"bOriginal");
	bOriginal.SetDescription(L"Nothink done only display input");
	bOriginal.SetDefault(false);
	bOriginal = false;
	this->Register(bOriginal);

	bInterpolate.SetName(L"bInterpolate");
	bInterpolate.SetDescription(L"Interpolate final gray image.");
	bInterpolate.SetDefault(true);
	bInterpolate = true;
	this->Register(bInterpolate);
}

/**
  *  @brief Destructor
  * /
TMOKuhn08::~TMOKuhn08()
{
}

/**
  *  @brief Convert to OpenCV matrix format
  * 
  *  @param image
  */
cv::Mat TMOKuhn08::convertToMat(TMOImage *image)
{
	double *data = image->GetData();
	int rows = image->GetWidth();
	int cols = image->GetHeight();
	double pL, pa, pb;
	cv::Mat mat(rows * cols, 3, CV_32F);
	for (int y = 0; y < rows; y++)
	{
		for (int x = 0; x < cols; x++)
		{
			pL = *data++;
			pa = *data++;
			pb = *data++;

			mat.at<float>(y + x * rows, 0) = pL;
			mat.at<float>(y + x * rows, 1) = pa;
			mat.at<float>(y + x * rows, 2) = pb;
		}
	}
	return mat;
}

/**
  *  @brief Copy image
  * 
  *  @param in Input
  *  @param out Output
  */
void TMOKuhn08::copy(TMOImage *in, TMOImage *out)
{
	double *pSourceData = in->GetData();	   /** You can work at low level data */
	double *pDestinationData = out->GetData(); /** Data are stored in form of array 
														  of three doubles representing
														   three colour components */
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

/**
  *  @brief Converts image
  */
int TMOKuhn08::Transform()
{
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	/** Covert to LAB color space */
	pSrc->Convert(TMO_LAB);
	pDst->Convert(TMO_LAB);

	int rows = pSrc->GetWidth();
	int cols = pSrc->GetHeight();

	double pL, pa, pb;

	/** display original image */
	if (bOriginal)
	{
		copy(pSrc, pDst);
	}
	else
	{
		/** initialize random seed: */
		srand(time(NULL));
		/** Convert to openCV mat format */
		cv::Mat samples = convertToMat(pSrc);
		/** Number of clustering colors */
		int clusterCount = iK;
		cv::Mat labels;
		int attempts = iAttempts;
		cv::Mat centers;
		/** quantization of image by kmeans in opencv */
		cv::kmeans(samples, clusterCount, labels, cv::TermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 10000, 0.0001), attempts, cv::KMEANS_PP_CENTERS, centers);

		/** Initialize particles */
		ParticlesManager particles(rows, cols, bColor, bChrominance);
		particles.initialize(labels, centers);

		double *pDestinationData = pDst->GetData();

		/** Show gray image */
		if (!bColor)
		{
			particles.computeMaxDistance();
			particles.createSprings();
			particles.compute(dTime);

			/** Final interpolation of gray color */
			/** Compare with src color and modify gray level */
			if (bInterpolate)
			{
				double *pSourceData = pSrc->GetData();

				/** Interpolate for all pixels */
				for (int y = 0; y < rows; y++)
				{
					for (int x = 0; x < cols; x++)
					{
						Particle source;
						source.L = *pSourceData++;
						source.a = *pSourceData++;
						source.b = *pSourceData++;

						double out = 0;
						Particle *p = particles.getParticle(x, y);
						/** Source L is smaller then quantized */
						if (source.L >= p->L)
						{
							out = p->gray + particles.getSkFactor(p) * particles.CalculateDistance(p, &source);
						}
						else
						{
							out = p->gray - particles.getSkFactor(p) * particles.CalculateDistance(p, &source);
						}
						*pDestinationData++ = out;
						*pDestinationData++ = 0;
						*pDestinationData++ = 0;
					}
				}
			}
			else
			{
				pDestinationData = pDst->GetData();
				particles.toDestination(pDestinationData);
			}
		}
		/** Show colored quantized image */
		else
		{
			pDestinationData = pDst->GetData();
			particles.toDestination(pDestinationData);
		}
	}

	pDst->Convert(TMO_RGB);
	return 0;
}
