/*--------------------------------------------------------------------------- *
 * TMOrc2gngm.cpp: implementation of the TMOrc2gngm class.   *
 * 	rc2gngm = Robust Color-to-gray via Nonlinear Global Mapping
 * --------------------------------------------------------------------------- */

#include "TMOrc2gngm.h"
#include "../matrix.h"								// matrix library

/* --------------------------------------------------------------------------- *
 * Constructor serves for describing a technique and input parameters          *
 * --------------------------------------------------------------------------- */
TMOrc2gngm::TMOrc2gngm()
{
	SetName(L"rc2gngm");
	SetDescription(L"Robust Color-to-gray via Nonlinear Global Mapping");

	alpha.SetName(L"alpha");
	alpha.SetDescription(L"Influence of the chromatic contrast on feature discriminability");
	alpha.SetDefault(1);
	alpha=1.;
	alpha.SetRange(-100.0,100.0);
	this->Register(alpha);
}

TMOrc2gngm::~TMOrc2gngm()
{
}

// this needs to be global variable, bec. matrix.h contains function's bodies
// => it can not be included to TMOrc2gngm.h file
Vector X(9);

/**
 * global mapping function
 * 
 * @param theta angle of input point from Lch color space
 * @return chrominance influance to resulting color
 */
double TMOrc2gngm::FunctionF(double theta){
	double result = 0.0;
	double h = TMOImage::DegreesToRadians(theta);		
	
	std::cerr.precision(10);
	
	for (int k = 1; k <= 9; k++){
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
 * Helmholtz–Kohlrausch Effect VAC method
 * 
 * @param input - point in Luv color space
 * @return HK predictor
 */
double TMOrc2gngm::Hk_effect_predictor(double * input){
	double u = input[1];
	double v = input[2];		
	
	double theta = atan2(v - LUV_WHITE_V, u - LUV_WHITE_U);
	
	double q_theta = -0.01585
		- 0.03017 * cos(theta)     - 0.04556 * cos(2 * theta)
		- 0.02667 * cos(3 * theta) - 0.00295 * cos(4 * theta)
		+ 0.14592 * sin(theta)     + 0.05084 * sin(2 * theta)
		- 0.01900 * sin(3 * theta) - 0.00764 * sin(4 * theta);
		
	double K_Br = 0.2717 * ((6.469 + 6.362 * pow(L_A, 0.4495)) / (6.469 + pow(L_A, 0.4495)));
	double s_uv = 13.0 * sqrt(pow(u - LUV_WHITE_U, 2) + pow(v - LUV_WHITE_V, 2));

	return 1 + (-0.1340 * q_theta + 0.0872 * K_Br) * s_uv;
}

/**
 * compute gradient of 2 pixels
 * 
 * @param first - first pixel in Lab 
 * @param second - second pixel in Lab 
 * @return gradient of 2 points
 */
double TMOrc2gngm::Gradient(double* first, double* second){
	double delta_L = first[0] - second[0];
	double delta_a = first[1] - second[1];
	double delta_b = first[2] - second[2];
	
	double luv_first[3], luv_second[3];
	
	PixelLabToLuv(first, luv_first);
	PixelLabToLuv(second, luv_second);
	
	/*std::cerr << "L: " << first[0] << ", a: " << first[1] << ", b: " << first[2] <<
	", l: " << luv_first[0] << ", u: " << luv_first[1] << ", v: " << luv_first[2] << std::endl;*/
	
	double delta_lhk = Hk_effect_predictor(luv_first) - Hk_effect_predictor(luv_second);
	int sign;
	const double R = 2.54 * sqrt(2);
		
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
void TMOrc2gngm::PixelLabToLuv(double* Lab, double* Luv){
	double x, y, z, L, u, v;
	
	TMOImage::LabToXyz(Lab[0], Lab[1], Lab[2], &x, &y, &z);		
	TMOImage::XyzToLuv(x, y, z, &L, &u, &v);
	
	//if ((counter > 50000) && (counter < 100000)) std::cerr << "x: " << x << ", y: " << y << ", z: " << z << "   L: " << L << ", u: " << u << ", v: " << v << std::endl;	
	//std::cerr << "l: " << Lab[0] << ", a: " << Lab[1] << ", b: " << Lab[2] << "   L: " << L << ", u: " << u << ", v: " << v << std::endl;	
	//std::cerr << "l: " << Lab[0] << ", a: " << Lab[1] << ", b: " << Lab[2] << "   x: " << x << ", y: " << y << ", z: " << z << std::endl;				
	
	Luv[0] = L;
	Luv[1] = u;
	Luv[2] = v;
}

/**
 * converts one pixel in Lch to Lab color space
 * 
 * @param Lch pixel in lch
 * @return pixel in Lab
 */
// TODO move to TMOImage
void TMOrc2gngm::PixelLchToLab(double* Lch, double* Lab){		
	double h = TMOImage::DegreesToRadians(Lch[2]);			
	
	Lab[0] = Lch[0];
	Lab[1] = Lch[1] * cos(h);
	Lab[2] = Lch[1] * sin(h);				
}

/* --------------------------------------------------------------------------- *
 * This overloaded function is an implementation of your tone mapping operator *
 * --------------------------------------------------------------------------- */
/**
 * transoformation function
 */
int TMOrc2gngm::Transform()
{
	
	// Source image is stored in local parameter pSrc
	// Destination image is in pDst

	// Initialy images are in RGB format, but you can convert it into other format
	pSrc->Convert(TMO_LCH);
	pDst->Convert(TMO_LCH);

	double* pSourceData = pSrc->GetData();					// You can work at low level data
	double* pDestinationData = pDst->GetData();				// Data are stored in form of array 
										// of three doubles representing
										// three colour components
	double L, c, h, g, f, p, q;
	double c_shift_left, h_shift_left, c_shift_right, h_shift_right, 
		c_shift_up, h_shift_up, c_shift_down, h_shift_down;		// variables for u and v
	double Gx, Gy, Lx, Ly;							// variables for p and q
	int lambda = pSrc->GetWidth() * pSrc->GetHeight();
	
	Matrix Ms(9, 9);
	Vector bs(9), u(9), v(9);	
	//Vector X(9);

	// get Ms and bs
	for (int j = 0; j < pSrc->GetHeight(); j++){
		for (int i = 0; i < pSrc->GetWidth(); i++){
			//std::cerr << "L: " << L << ", c: " << c << ", h: " << h << std::endl;
			
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
			for (int k = 1; k <= 9; k++){				
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
			//pSrc->Convert(TMO_LAB);						
			double lab1[3], lab2[3];
			int minus = (i == 0) ? 0 : 1;
			int plus = (i == pSrc->GetWidth() - 1) ? 0 : 1;
			PixelLchToLab(pSrc->GetPixel(i + plus, j), lab1);
			PixelLchToLab(pSrc->GetPixel(i - minus, j), lab2);			
			//std::cerr << "LCH2LAB L: " << pSrc->GetPixel(i + plus, j)[0] << ", c: " << pSrc->GetPixel(i + plus, j)[1] << ", h: " << pSrc->GetPixel(i + plus, j)[2] << ", l: " << lab1[0] << ", a: " << lab1[1] << ", b: " << lab1[2] << std::endl;			
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
	Matrix ident(9,9);	
	ident.identity();	
	ident = ident * lambda;	
	
	std::cout << "bs:" << std::endl;
	std::cerr << bs << std::endl << std::endl;
	
	std::cout << "Ms:" << std::endl;
	std::cerr << Ms << std::endl << std::endl;	
	
	Ms = Ms + ident;	
	X = pseudoinverse(Ms) * bs;	
	
	std::cout << "X:" << std::endl;
	std::cerr << X << std::endl << std::endl;	
	
	std::cerr << "Function F test" << std::endl;
	std::cerr.precision(10);
	for (double i = 0.0; i < 360; i++){
		//std::cerr << fixed << "f(" << i << "): " << functionF(i) << std::endl;
		std::cerr << fixed << FunctionF(i) << std::endl;
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
			*pDestinationData++ = 296.812926236627;
		}
	}	
	
	pDst->Convert(TMO_RGB);
	return 0;
}



