/*******************************************************************************
*                                                                              *
*                         Brno University of Technology                        *
*                       Faculty of Information Technology                      *
*                                                                              *
*                         Color-to-Grayscale Conversions                       *
*                                                                              *
*                                 diploma thesis                               *
*             Author: Petr Pospisil [xpospi68 AT stud.fit.vutbr.cz]            *
*                                    Brno 2016                                 *
*                                                                              *
*******************************************************************************/

/*---------------------------------------------------------------------------- *
 * TMOKim09.cpp: implementation of the TMOKim09 class.                         *
 *               Robust Color-to-gray via Nonlinear Global Mapping             *
 * Method number: 1                                                            *
 * --------------------------------------------------------------------------- */

#include "TMOKim09.h"
#include "../tmolib/matrix.h"							// matrix library

/**
 * constructor, prepare parameters
 */
TMOKim09::TMOKim09()
{
	SetName(L"Kim09");
	SetDescription(L"Robust Color-to-gray via Nonlinear Global Mapping");

	alpha.SetName(L"alpha");
	alpha.SetDescription(L"Influence of the chromatic contrast on feature discriminability");
	alpha.SetDefault(1);
	alpha=1.;
	alpha.SetRange(-100.0,100.0);
	this->Register(alpha);
	
	////Added by V.Vlkovic in 2018 for video conversion
	beta.SetName(L"beta");
	beta.SetDescription(L"Parameter to control the relative weights of spatial and temporal terms");
	beta.SetDefault(0.5);
	beta.SetRange(0.0,1.0);
	this->Register(beta);
	
	// show/hide debug info
	verbose.SetName(L"v");
	verbose.SetDescription(L"Verbose output");
	verbose.SetDefault(false);
	verbose=false;	
	this->Register(verbose);
}

TMOKim09::~TMOKim09()
{
}

// this needs to be global variable, bec. matrix.h contains function's bodies
// => it can not be included to TMOKim09.h file
mtx::Vector X(9);

/**
 * global mapping function
 * 
 * @param theta angle of input point from Lch color space
 * @return chrominance influance to resulting color
 */
double TMOKim09::FunctionF(double theta){
	double result = 0.0;
	double h = TMOImage::DegreesToRadians(theta);
	
	if (verbose) std::cerr.precision(10);
		
	for (int k = 1; k <= (2*N + 1); k++){
		if (k <= N){							// k = 1, 2, 3, 4
			result += cos(k * h) * X[k-1];			
		}else if ((k > N) && (k < 2*N + 1)){				// k = 5, 6, 7, 8	
			result += sin((k - N) * h) * X[k-1];			
		}else{								// k = 9	
			result += X[k-1];			
		}		
	}	
	
	return result;
}

/**
 * Helmholtz–Kohlrausch Effect
 * 
 * @param input - point in Luv color space
 * @return HK predictor
 */
double TMOKim09::HkEffectPredictor(double * input){
	double u = input[1] / 255 * 354.0 - 134.0;
	double v = input[2] / 255.0 * 256.0 - 140.0;

	double luvU = LUV_WHITE_U;
	double luvV = LUV_WHITE_V;
	
	double theta = atan2(v - luvV, u - luvU);
	
	double q_theta = -0.01585
		- 0.03017 * cos(theta)     - 0.04556 * cos(2 * theta)
		- 0.02667 * cos(3 * theta) - 0.00295 * cos(4 * theta)
		+ 0.14592 * sin(theta)     + 0.05084 * sin(2 * theta)
		- 0.01900 * sin(3 * theta) - 0.00764 * sin(4 * theta);
		
	double K_Br = 0.2717 * ((6.469 + 6.362 * pow(L_A, 0.4495)) / (6.469 + pow(L_A, 0.4495)));
	double s_uv = 13.0 * sqrt(pow(u - luvU, 2) + pow(v - luvV, 2));

	return input[0] + (-0.1340 * q_theta + 0.0872 * K_Br) * s_uv * input[0];
}

/**
 * compute gradient of 2 pixels
 * 
 * @param first - first pixel in Lab 
 * @param second - second pixel in Lab 
 * @return gradient of 2 points
 */
double TMOKim09::Gradient(double* first, double* second){
	double delta_L = (first[0] - second[0]) / LAB_TO_LUV;
	double delta_a = first[1] - second[1];
	double delta_b = first[2] - second[2];
	
	double luv_first[3], luv_second[3];
	
	PixelLabToLuv(first, luv_first);
	PixelLabToLuv(second, luv_second);
	
	/*if (verbose) std::cerr << "L: " << first[0] << ", a: " << first[1] << ", b: " << first[2] <<
	", l: " << luv_first[0] << ", u: " << luv_first[1] << ", v: " << luv_first[2] << std::endl;*/
	
	double delta_lhk = HkEffectPredictor(luv_first) - HkEffectPredictor(luv_second);	
	const double R = 2.54 * sqrt(2.0);
	int sign = 0;
		
	if (delta_lhk != 0.0){
		sign = (delta_lhk > 0) ? 1 : -1;
	} else if (delta_L != 0.0){
		sign = (delta_L > 0) ? 1 : -1;
	} else {
		sign = (pow(delta_L, 3) + pow(delta_a, 3) + pow(delta_b, 3) > 0) ? 1 : -1;
	}	 					
	
	return sign * sqrt(pow(delta_L, 2) + pow(alpha * (sqrt((pow(delta_a, 2) + pow(delta_b, 2))) / R), 2));
}

/**
 * wrapper for converting Lab > XYZ and XYZ > Luv
 * 
 * @param Lab - input parameter: pixel in Lab color space
 * @param Luv - output parameter: pixel in Luv color space
 */
void TMOKim09::PixelLabToLuv(double* Lab, double* Luv){
	double x, y, z, L, u, v;
	
	TMOImage::LabToXyz(Lab[0], Lab[1], Lab[2], &x, &y, &z);		
	TMOImage::XyzToLuv(x, y, z, &L, &u, &v);		
	
	Luv[0] = L;
	Luv[1] = u;
	Luv[2] = v;
}

/**
 * converts one pixel in Lch to Lab color space
 * 
 * @param Lch pixel in lch (input)
 * @param Lab pixel in Lab (output)
 */
void TMOKim09::PixelLchToLab(double* Lch, double* Lab){
	double h = TMOImage::DegreesToRadians(Lch[2]);	
	
	Lab[0] = Lch[0];
	Lab[1] = Lch[1] * cos(h);
	Lab[2] = Lch[1] * sin(h);
}

/**
 * transformation function
 * @return exit code
 */
int TMOKim09::Transform(){
	pSrc->Convert(TMO_LCH);
	pDst->Convert(TMO_LCH);

	double* pSourceData = pSrc->GetData();
	double* pDestinationData = pDst->GetData();

	double L, c, h, g, f, p, q;
	double c_shift_left, h_shift_left, c_shift_right, h_shift_right, 
		c_shift_up, h_shift_up, c_shift_down, h_shift_down;		// variables for u and v
	double Gx, Gy, Lx, Ly;							// variables for p and q
	int lambda = pSrc->GetWidth() * pSrc->GetHeight();
	
	mtx::Matrix Ms(9, 9);
	mtx::Vector bs(9), u(9), v(9);	
	//Vector X(9);

	// get Ms and bs
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			//if (v) std::cerr << "L: " << L << ", c: " << c << ", h: " << h << std::endl;
			
			// prepare shifted c and h
			c_shift_left = (i == 0) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i - 1, j)[1];
			h_shift_left = (i == 0) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i - 1, j)[2];
			c_shift_right = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i + 1, j)[1];
			h_shift_right = (i == pSrc->GetWidth() - 1) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i + 1, j)[2];
			c_shift_down = (j == 0) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i, j - 1)[1];
			h_shift_down = (j == 0) ?  pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i, j - 1)[2];
			c_shift_up = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[1] : pSrc->GetPixel(i, j + 1)[1];
			h_shift_up = (j == pSrc->GetHeight() - 1) ? pSrc->GetPixel(i, j)[2] : pSrc->GetPixel(i, j + 1)[2];	
			
			h_shift_left = TMOImage::DegreesToRadians(h_shift_left);
			h_shift_right = TMOImage::DegreesToRadians(h_shift_right);
			h_shift_down = TMOImage::DegreesToRadians(h_shift_down);
			h_shift_up = TMOImage::DegreesToRadians(h_shift_up);
			
			// 1. compute u and v
			// u = (C * t)_x	v = (C * t)_y
			for (int k = 1; k <= (2*N + 1); k++){				
				// compute u
				if (k <= N){					// k = 1, 2, 3, 4
					u[k-1] = c_shift_right * cos(k * h_shift_right) - c_shift_left * cos(k * h_shift_left);					
				}else if ((k > N) && (k < 2*N + 1)){		// k = 5, 6, 7, 8
					u[k-1] = c_shift_right * sin((k-N) * h_shift_right) - c_shift_left * sin((k-N) * h_shift_left);					
				}else{						// k = 9
					u[k-1] = c_shift_right - c_shift_left;					
				}
				
				// compute v
				if (k <= N){					// k = 1, 2, 3, 4
					v[k-1] = c_shift_up * cos(k * h_shift_up) - c_shift_down * cos(k * h_shift_down);
				}else if ((k > N) && (k < 2*N + 1)){		// k = 5, 6, 7, 8
					v[k-1] = c_shift_up * sin((k-N) * h_shift_up) - c_shift_down * sin((k-N) * h_shift_down);
				}else{						// k = 9
					v[k-1] = c_shift_up - c_shift_down;
				}
			}
			
			// 2. compute p, q
			// p = G^x - L_x	q = G^y - L_y			
			double lab1[3], lab2[3];
			int minus = (i == 0) ? 0 : 1;
			int plus = (i == pSrc->GetWidth() - 1) ? 0 : 1;
			PixelLchToLab(pSrc->GetPixel(i + plus, j), lab1);
			PixelLchToLab(pSrc->GetPixel(i - minus, j), lab2);						
			Gx = Gradient(lab1, lab2);									
			Lx = pSrc->GetPixel(i + plus, j)[0] - pSrc->GetPixel(i - minus, j)[0];
			
			minus = (j == 0) ? 0 : 1;
			plus = (j == pSrc->GetHeight() - 1) ? 0 : 1;			
			PixelLchToLab(pSrc->GetPixel(i, j + plus), lab1);
			PixelLchToLab(pSrc->GetPixel(i, j - minus), lab2);
			Gy = Gradient(lab1, lab2);			
			Ly = pSrc->GetPixel(i, j + plus)[0] - pSrc->GetPixel(i, j - minus)[0];
			
			p = Gx - Lx;
			q = Gy - Ly;
			
			// 3. compute Ms
			// Ms = E(u*u^T + v*v^T)			
			Ms += ((u * transpose(u)) + (v * transpose(v)));
			
			// 4. compute bs
			// bs = E(pu + qv)
			bs += ((p * u) + (q * v));
		}
	}
		
	// 5. compute x	
	// X_image = (M_s + lambda * I)^(-1) * bs			
	mtx::Matrix ident(9,9);	
	ident.identity();	
	ident = ident * lambda;	
	
	if (verbose) std::cout << "bs:" << std::endl;	
	if (verbose) std::cerr << bs << std::endl << std::endl;	
	if (verbose) std::cout << "Ms:" << std::endl;
	if (verbose) std::cerr << Ms << std::endl << std::endl;	
	
	Ms = Ms + ident;	
	X = pseudoinverse(Ms) * bs;	
	
	if (verbose) std::cout << "X:" << std::endl;
	if (verbose) std::cerr << X << std::endl << std::endl;	
	
	if (verbose) std::cerr << "Function F test" << std::endl;
	if (verbose) std::cerr.precision(10);
	if (verbose){
		for (double i = 0.0; i < 360; i++){
			//std::cerr << fixed << "f(" << i << "): " << functionF(i) << std::endl;
			std::cerr << fixed << FunctionF(i) << std::endl;
		}
	}
	
	// restore pointer to source data
	pSourceData = pSrc->GetData();	
	
	// 6. compute resulted grayscale image
	for (int j = 0; j < pSrc->GetHeight(); j++)
	{		
		for (int i = 0; i < pSrc->GetWidth(); i++)
		{
			L = *pSourceData++;
			c = *pSourceData++;
			h = *pSourceData++;			
			
			// 6. resolve f(theta)					
			g = L + FunctionF(h)*c;					// global mapping

			// store results to the destination image			
			*pDestinationData++ = g;
			*pDestinationData++ = 0.0;
			*pDestinationData++ = LCH_H_WHITE;			
		}
	}	
	
	pDst->Convert(TMO_RGB);
	return 0;
}
/**
 * Disclaimer for the sake of code continuity and for espect to my colleague
 * I din't want to change any of his code and I made a decision just to tweak his code as little as possible
 * That is why I addet just bits to the code that I needed for finishing the method
 * In this way I would like to thank my colleague Mr. Pospisil for writing the foundations of this method
 *  This method was written by Vladimir Vlkovic for diploma thesis in 2018
 * */
int TMOKim09::TransformVideo()
{
  int height = vSrc->GetHeight();
  int width = vSrc->GetWidth();
  
  
  
  TMOImage currentFrame, previousFrame, currentGrayFrame, previousGrayFrame;
  previousFrame.New(width,height,TMO_LCH);
  currentGrayFrame.New(width,height,TMO_LCH);
  previousGrayFrame.New(width,height,TMO_RGB);
  
  double L, c, h, g, f, p, q;
  double c_shift_left, h_shift_left, c_shift_right, h_shift_right, 
	  c_shift_up, h_shift_up, c_shift_down, h_shift_down, c_current, h_current,
	  l_current, gray_previous, color_diff;		
  double Gx, Gy, Lx, Ly;						
  int lambda = width * height;
  
  mtx::Matrix Ms(9, 9),Mv(9,9),Mt(9,9);
  mtx::Vector bs(9),bv(9),bt(9), u(9), v(9),t(9);	
  
  cv::VideoCapture cap=vSrc->getVideoCaptureObject();
  for(int frameCounter = 0; frameCounter<vSrc->GetTotalNumberOfFrames(); frameCounter++)
  {
     
      vSrc->getTMOImageVideoFrame(cap,frameCounter,currentFrame);
      currentFrame.Convert(TMO_LCH);
      
      for (int j = 0; j < height; j++){
	
		for (int i = 0; i < width; i++){
			//if (v) std::cerr << "L: " << L << ", c: " << c << ", h: " << h << std::endl;
			
			// prepare shifted c and h
			c_shift_left = (i == 0) ? currentFrame.GetPixel(i, j)[1] : currentFrame.GetPixel(i - 1, j)[1];
			h_shift_left = (i == 0) ? currentFrame.GetPixel(i, j)[2] : currentFrame.GetPixel(i - 1, j)[2];
			c_shift_right = (i == currentFrame.GetWidth() - 1) ? currentFrame.GetPixel(i, j)[1] : currentFrame.GetPixel(i + 1, j)[1];
			h_shift_right = (i == currentFrame.GetWidth() - 1) ? currentFrame.GetPixel(i, j)[2] : currentFrame.GetPixel(i + 1, j)[2];
			c_shift_down = (j == 0) ? currentFrame.GetPixel(i, j)[1] : currentFrame.GetPixel(i, j - 1)[1];
			h_shift_down = (j == 0) ?  currentFrame.GetPixel(i, j)[2] : currentFrame.GetPixel(i, j - 1)[2];
			c_shift_up = (j == currentFrame.GetHeight() - 1) ? currentFrame.GetPixel(i, j)[1] : currentFrame.GetPixel(i, j + 1)[1];
			h_shift_up = (j == currentFrame.GetHeight() - 1) ? currentFrame.GetPixel(i, j)[2] : currentFrame.GetPixel(i, j + 1)[2];	
			
			h_shift_left = TMOImage::DegreesToRadians(h_shift_left);
			h_shift_right = TMOImage::DegreesToRadians(h_shift_right);
			h_shift_down = TMOImage::DegreesToRadians(h_shift_down);
			h_shift_up = TMOImage::DegreesToRadians(h_shift_up);
			
			l_current = currentFrame.GetPixel(i, j)[0];
			c_current = currentFrame.GetPixel(i, j)[1];
			h_current = currentFrame.GetPixel(i, j)[2];
			
			// 1. compute u and v and t
			
			for (int k = 1; k <= (2*N + 1); k++){				
				// compute u, v ,t
				if (k <= N){					// k = 1, 2, 3, 4
					u[k-1] = c_shift_right * cos(k * h_shift_right) - c_shift_left * cos(k * h_shift_left);	
					v[k-1] = c_shift_up * cos(k * h_shift_up) - c_shift_down * cos(k * h_shift_down);
					t[k-1] = cos(k * h_current);
				}else if ((k > N) && (k < 2*N + 1)){		// k = 5, 6, 7, 8
					u[k-1] = c_shift_right * sin((k-N) * h_shift_right) - c_shift_left * sin((k-N) * h_shift_left);
					v[k-1] = c_shift_up * sin((k-N) * h_shift_up) - c_shift_down * sin((k-N) * h_shift_down);
					t[k-1] = sin((k-N) * h_current);
				}else{						// k = 9
					u[k-1] = c_shift_right - c_shift_left;	
					v[k-1] = c_shift_up - c_shift_down;
					t[k-1] = 1.0;
				}
	
			}
			
			// 2. compute p, q
			// p = G^x - L_x	q = G^y - L_y			
			double lab1[3], lab2[3];
			int minus = (i == 0) ? 0 : 1;
			int plus = (i == width - 1) ? 0 : 1;
			PixelLchToLab(currentFrame.GetPixel(i + plus, j), lab1);
			PixelLchToLab(currentFrame.GetPixel(i - minus, j), lab2);						
			Gx = Gradient(lab1, lab2);									
			Lx = currentFrame.GetPixel(i + plus, j)[0] - currentFrame.GetPixel(i - minus, j)[0];
			
			minus = (j == 0) ? 0 : 1;
			plus = (j == height - 1) ? 0 : 1;			
			PixelLchToLab(currentFrame.GetPixel(i, j + plus), lab1);
			PixelLchToLab(currentFrame.GetPixel(i, j - minus), lab2);
			Gy = Gradient(lab1, lab2);			
			Ly = currentFrame.GetPixel(i, j + plus)[0] - currentFrame.GetPixel(i, j - minus)[0];
			
			p = Gx - Lx;
			q = Gy - Ly;
			
			// 3. compute Ms
			// Ms = E(u*u^T + v*v^T)			
			Ms += ((u * transpose(u)) + (v * transpose(v)));
			
			// 4. compute bs
			// bs = E(pu + qv)
			bs += ((p * u) + (q * v));
			
			if(frameCounter > 0)  /// if int is not the first frame
			{
			  PixelLchToLab(currentFrame.GetPixel(i , j), lab1);
			  PixelLchToLab(previousFrame.GetPixel(i, j), lab2);						
			  color_diff = Gradient(lab1, lab2);
			  gray_previous = previousGrayFrame.GetPixel(i, j)[0];
			  Mt += std::pow(c_current,2) * t * transpose(t);              ////get Mt, see paper
			
			  bt += c_current * (gray_previous + color_diff - l_current) * t;  ///get bt, see paper
			}
			
			
			
		}
	}
	mtx::Matrix ident(9,9);	
	ident.identity();	
	ident = ident * lambda;
	if(frameCounter > 0)  ///if it si not the first frame
	{
	  Mv = (1 - beta) * Ms + beta * Mt;   ///computation of Mv, bv see papaer
	  bv = (1 - beta) * bs + beta * bt;
	  
	  Mv = Mv + ident;
	  X = pseudoinverse(Mv) * bv;
	}
	else  ///if it is first frame
	{

	  Ms = Ms + ident;	
	  X = pseudoinverse(Ms) * bs;
	}
	double* color= new double[3];
	double* imageData = new double[width*height*3];
	int count= 0;
	for (int j = 0; j < height; j++)
	{		
		for (int i = 0; i < width; i++)
		{
		  
		  color = currentFrame.GetPixel(i, j);
		  L = color[0];
		  c = color[1];
		  h = color[2];			
		  
		  // 6. resolve f(theta)					
		  g = L + FunctionF(h)*c;					

		  // store results to the destination image			
		  imageData[count++] = g;	
		  imageData[count++] = 0.0;
		  imageData[count++] = LCH_H_WHITE;
		
		}
	}
	currentGrayFrame.SetData(imageData);   ///set data to current image
	currentGrayFrame.Convert(TMO_RGB);
	vDst->setTMOImageFrame(vDst->getVideoWriterObject(),currentGrayFrame);  //write to video
	
	previousFrame = currentFrame;
	previousGrayFrame = currentGrayFrame;  //set previous frames
	currentGrayFrame.Convert(TMO_LCH);
	
	//std::cerr<<"Frames: "<< frameCounter+1<<" / " <<vSrc->GetTotalNumberOfFrames()<<std::endl;
      
  }
}



