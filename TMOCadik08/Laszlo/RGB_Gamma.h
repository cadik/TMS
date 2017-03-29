
//RGB linear - RGB gamma corrected - XYZ 

//void XYZ_RGB709(double X, double Y, double Z,
//				  double*R, double*G, double*B);

//void RGB709_XYZ(double R, double G, double B,
//   	          double*X, double*Y, double*Z);

//double Rlinear_Rdisplay_01(double Rlinear);
//double Rdisplay_01_Rlinear(double Rdisplay);

//unsigned char Rdisplay_clipping(double Rdisplay_01);
//double Rdisplay_normalization(unsigned char Rdisplay);

//void RGBlinear_RGBdisplay(
//		double Rlinear, 
//		double Glinear, 
//		double Blinear,
//		BYTE * Rdisplay,
//		BYTE * Gdisplay,
//		BYTE * Bdisplay);

void GRAM(long n, double vmat[][4], double w[][4])
{
//global v-bol w vektorokat csinal

long i,j,k;
double s,t;

	for(j=1;j<=n;j++)
		w[1][j] = vmat[1][j];

	for(j=1,s=0;j<=n;j++)
	{
		t = w[1][j];
		s += t*t;
	}

	s = sqrt((double)s);
	for(j=1;j<=n;j++)
		w[1][j] /= s;
		// w egysegvektor

	for(k=2;k<=n;k++)
	{//k

     for(i=1;i<=n;i++)
		w[k][i] = vmat[k][i];

        for(j=1;j<k;j++)
        {//j

        //wv skalarszorzat
        for(i=1,s=0;i<=n;i++)
			s+= w[j][i]*vmat[k][i];

        for(i=1;i<=n;i++)
			w[k][i] -= s * w[j][i];

        }//j

        for(i=1,s=0;i<=n;i++)
        {
			t = w[k][i];
			s+= t*t;
        }
			s = sqrt((double)s);

        for(i=1;i<=n;i++)
			w[k][i] /= s;
       
	}//k

}//end GRAM

double min_distance_step(long n,
double vmat[][4],double w[][4],double cons[4])
{
double s1,s2,t,TT;
long j,k;

// a global x tombon valtozon hajt vegre modositast

	TT=0;

	for(k=1;k<=n;k++)
	{

        for(j=1,s1=s2=0;j<=n;j++)
        {
			t = vmat[k][j];
			s1 += t * x[j];
			s2 += t * w[k][j];
        }

        t = ((double)cons[k] - s1) / s2;

        for(j=1;j<=n;j++)
			x[j] += t * w[k][j];

			TT += t*t;
	}

	return(TT);

}

long INVERZIO(long n)
{
//input  a global  vmat
//output a global  inverz

    double t,s,q1,q2;

    long i,j,k;	

	GRAM(n, vmat, wmat);

	for(k=1;k<=n;k++)
	{

		for(i=0;i<=n;i++)
		b[i] = 0;

		b[k] = 1;


		for(i=0;i<=n;i++)
		x[i] = 0;

		{
		min_distance_step(n,vmat,wmat,b);

		t = 0;
		for(i=1;i<=n;i++)
		{
			 s=0;
			 for(j=1;j<=n;j++)
			 s+=vmat[i][j] * x[j];
 
			s = fabs(s - b[i]);
			if(t<s) t = s;
		}

		}//l

		for(i=1;i<=n;i++)
		inverz[i][k] = x[i];

	}//k


	q1=q2=0;
	for(i=1;i<=n;i++)
	for(j=1;j<=n;j++)
	{
		s=0;
		for(k=1;k<=n;k++)
		s += vmat[i][k]*inverz[k][j]; 

		if(i==j)
		{
		  if(fabs(s-1.) > q1) q1 = fabs(s-1.);
		}
		else
		{
		  if(fabs(s) > q2) q2 = fabs(s);
		}

	}//ij

	if(q1 > EPS_9 || q2 > EPS_9) 
	return(0);
	//ekkor szingularis eset van, nincs inverz, vagy nem eleg pontos

	return(1);
}

long XYZ_rgb_conversion_matrices(double x_R,
								 double y_R,
								 double x_G,
								 double y_G,
								 double x_B,
								 double y_B,
								 double x_white,
								 double y_white)
{
//default D65 white point
//#define D65_x (double) 0.3127269
//#define D65_y (double) 0.3290232
//			x       y      
//R        0.6400  0.3300  
//G        0.3000  0.6000  
//B        0.1500  0.0600  
//The CIE chromaticities for the red, green, and blue
//ITU-R BT.709 reference primaries,
//and for CIE Standard Illuminant D65

	long i, j;
	double s;

	s  = 100. / y_white ;
	white[1] = s * x_white;
	white[2] = 100.;
	white[3] = s * (1. - (x_white + y_white));

	vmat[1][1] = x_R;
	vmat[1][2] = x_G;
	vmat[1][3] = x_B;

	vmat[2][1] = y_R;
	vmat[2][2] = y_G;
	vmat[2][3] = y_B;

	vmat[3][1] = 1-(x_R+y_R);
	vmat[3][2] = 1-(x_G+y_G);
	vmat[3][3] = 1-(x_B+y_B);

	if(INVERZIO(3)==0) return(0);
	//problem with white ballance
	
	for(i=1;i<=3;i++)
	{
		x[i] = 0;
		for(j=1;j<=3;j++)
			x[i] += inverz[i][j] * white[j];
	}

	for(i=1;i<=3;i++)
		for(j=1;j<=3;j++)
		{
			vmat[i][j] *= x[j];
			rgb_XYZ[i][j] = vmat[i][j];
			/////////////
		}
		
	if(INVERZIO(3)==0) return(0);
	//problem with white ballance
	
	for(i=1;i<=3;i++)
		for(j=1;j<=3;j++)
			XYZ_rgb[i][j] = inverz[i][j];
			/////////////

	return(1); 
	//white ballance Ok
}

void _XYZ_rgb(double X,
			 double Y,
			 double Z,
			 double *r,
			 double *g,
			 double *b)
{

	*r = XYZ_rgb[1][1]*X + XYZ_rgb[1][2]*Y + XYZ_rgb[1][3]*Z;

	*g = XYZ_rgb[2][1]*X + XYZ_rgb[2][2]*Y + XYZ_rgb[2][3]*Z;

	*b = XYZ_rgb[3][1]*X + XYZ_rgb[3][2]*Y + XYZ_rgb[3][3]*Z;

}

void _rgb_XYZ(double r,
			  double g,
			  double b,
			  double *X,
			  double *Y,
			  double *Z)
{

	*X = rgb_XYZ[1][1]*r + rgb_XYZ[1][2]*g + rgb_XYZ[1][3]*b;

	*Y = rgb_XYZ[2][1]*r + rgb_XYZ[2][2]*g + rgb_XYZ[2][3]*b;

	*Z = rgb_XYZ[3][1]*r + rgb_XYZ[3][2]*g + rgb_XYZ[3][3]*b;

}

void XYZ_RGB709(double X, double Y, double Z,
				double*R, double*G, double*B)
{
// white Y = 100
// white R,G,B = 1,  D65
// [ R709 ] [ 3.240479 -1.53715  -0.498535 ] [ X ] 
// [ G709 ]=[-0.969256  1.875991  0.041556 ]*[ Y ] 
// [ B709 ] [ 0.055648 -0.204043  1.057311 ] [ Z ] 

	*R	=	0.03240479*X - 0.0153715 *Y - 0.00498535*Z;
    *G	=  -0.00969256*X + 0.01875991*Y + 0.00041556*Z;
    *B  =   0.00055648*X - 0.00204043*Y + 0.01057311*Z;

}

void RGB709_XYZ(double R, double G, double B,
				double*X, double*Y, double*Z)
{
// white Y = 100
// white R,G,B = 1,  D65
// [ X ] [ 0.412453  0.35758   0.180423 ] [ R709 ] 
// [ Y ]=[ 0.212671  0.71516   0.072169 ]*[ G709 ] 
// [ Z ] [ 0.019334  0.119193  0.950227 ] [ B709 ] 

	*X = 41.2453*R + 35.758 *G + 18.0423*B;
	*Y = 21.2671*R + 71.516 *G +  7.2169*B;
    *Z =  1.9334*R + 11.9193*G + 95.0227*B;

}

double flair_value(double Gamma_ambient)
{
//ambient or viewing Gamma dependent flair (fog) value
//the display lightness for zero input would be 0.5 - 5 % 
//the gamma function works with this offset

	double a, b;
	//Gamma_ambient az [1., 1.5] intervallumban lehet

	if(Gamma_ambient < 1.1498987)
	{
		a = 5;
		b = -4./0.1498987;
		return((a + b * (Gamma_ambient - 1.))/100.);
	}
	else
	{
		b = - 0.5 / (0.5 - 0.1498987);
		a = 0.5 * (1. - b);
		return((a + b * (Gamma_ambient - 1.))/100.);
	}

}

void Gamma_precomputation(double Gamma_R,
						  double Gamma_G,
						  double Gamma_B,
						  double Gamma_ambient)
{
//if Gamma = 2.2 and Gamma_ambient = 1.1498987,
//than delta = 0.01 => gamma function with 1% offset
//beginnig part linear: value and derivated are contineous

	double f,F,gr,gg,gb,delta;

	f = flair_value(Gamma_ambient);

	gr = Gamma_R / Gamma_ambient;
	F = pow(f,1./gr);
	_Rdelta = delta = F/(1-F);
	_R0 = delta /(gr * (1.+ delta) - 1.);
	_RC = gr * pow((_R0+delta) / (1.+ delta) , gr - 1.);

	gg = Gamma_G / Gamma_ambient;
	F = pow(f,1./gg);
	_Gdelta = delta = F/(1-F);
	_G0 = delta /(gg * (1.+ delta) - 1.);
	_GC = gg * pow((_G0+delta) / (1.+ delta) , gg - 1.);

	gb = Gamma_B / Gamma_ambient;
	F = pow(f,1./gb);
	_Bdelta = delta = F/(1-F);
	_B0 = delta /(gb * (1.+ delta) - 1.);
	_BC = gb * pow((_B0+delta) / (1.+ delta) , gb - 1.);

//	printf("\n\tgamma_R=%lf delta=%lf  R0=%lf  C=%lf",
//			Gamma_R,delta,_R0,_RC);

}
/*
void rgb_RGB(double r, 
			 double g,
			 double b,
			 double * R,
			 double * G,
			 double * B,
			 double Gamma_R,
			 double Gamma_G,
			 double Gamma_B,
			 double Gamma_ambient)
{

	double r0, g0, b0, rg, gg, bg;

	Gamma_precomputation(Gamma_R,Gamma_G,Gamma_B,Gamma_ambient);

	r0 = _R0 / _RC;
	rg = Gamma_ambient/Gamma_R;
	if(r < 0) r = 0;
	if(r > 1) r = 1;
	//clipping

	if(r < r0)
		*R = 255. * r / _RC;
	else
		*R = 255. *( - _Rdelta + (1.+_Rdelta)*pow(r,rg) );

	g0 = _G0 / _GC;
	gg = Gamma_ambient/Gamma_G;
	if(g < 0) g = 0;
	if(g > 1) g = 1;
	//clipping

	if(g < g0)
		*G = 255. * g / _GC;
	else
		*G = 255. *( - _Gdelta + (1.+_Gdelta)*pow(g,gg) );

	b0 = _B0 / _BC;
	bg = Gamma_ambient/Gamma_B;
	if(b < 0) b = 0;
	if(b > 1) b = 1;
	//clipping

	if(b < b0)
		*B = 255. * b / _BC;
	else
		*B = 255. *( - _Bdelta + (1.+_Bdelta)*pow(b,bg) );

}

void RGB_rgb(double R, 
			 double G,
			 double B,
			 double * r,
			 double * g,
			 double * b,
			 double Gamma_R,
			 double Gamma_G,
			 double Gamma_B,
			 double Gamma_ambient)
{

	double gr, gg, gb;

	Gamma_precomputation(Gamma_R,Gamma_G,Gamma_B,Gamma_ambient);

	gr = Gamma_R / Gamma_ambient;
	R /= 255.;
	if(R < 0) R = 0;
	if(R > 1) R = 1;

	if(R < _R0) 
		*r = _RC * R;
	else
		*r = pow((R+_Rdelta)/(1.+_Rdelta) , gr);

	gg = Gamma_G / Gamma_ambient;
	G /= 255.;
	if(G < 0) G = 0;
	if(G > 1) G = 1;

	if(G < _G0) 
		*g = _GC * G;
	else
		*g = pow((G+_Gdelta)/(1.+_Gdelta) , gg);

	gb = Gamma_B / Gamma_ambient;
	B /= 255.;
	if(B < 0) B = 0;
	if(B > 1) B = 1;

	if(B < _B0) 
		*b = _BC * B;
	else
		*b = pow((B+_Bdelta)/(1.+_Bdelta) , gb);

}
*/

void rgb_RGB(double r, 
			 double g,
			 double b,
			 double * R,
			 double * G,
			 double * B,
			 double Gamma_R,
			 double Gamma_G,
			 double Gamma_B,
			 double Gamma_ambient)
{

	double rg, gg, bg;

	rg = Gamma_ambient/Gamma_R;
	if(r < 0) r = 0;
	if(r > 1) r = 1;
	//clipping
	*R = 255. * pow(r,rg);

	gg = Gamma_ambient/Gamma_G;
	if(g < 0) g = 0;
	if(g > 1) g = 1;
	//clipping
	*G = 255. * pow(g,gg);

	bg = Gamma_ambient/Gamma_B;
	if(b < 0) b = 0;
	if(b > 1) b = 1;
	//clipping
	*B = 255. * pow(b,bg);

}

void RGB_rgb(double R, 
			 double G,
			 double B,
			 double * r,
			 double * g,
			 double * b,
			 double Gamma_R,
			 double Gamma_G,
			 double Gamma_B,
			 double Gamma_ambient)
{

	double gr, gg, gb;

	gr = Gamma_R / Gamma_ambient;
	R /= 255.;
	if(R < 0) R = 0;
	if(R > 1) R = 1;
	*r = pow(R , gr);

	gg = Gamma_G / Gamma_ambient;
	G /= 255.;
	if(G < 0) G = 0;
	if(G > 1) G = 1;
    *g = pow(G , gg);

	gb = Gamma_B / Gamma_ambient;
	B /= 255.;
	if(B < 0) B = 0;
	if(B > 1) B = 1;
	*b = pow(B , gb);

}

double Rlinear_Rdisplay_01(double Rlinear)
{
//both are normalized to [0,1]
//default 709 gamma = 2.2
	double Rdisplay_01;

	Rdisplay_01 =  Rlinear <= 0.018 ?	4.5 * Rlinear : 
									-0.099 + 1.099 * pow(Rlinear, 0.45); 
	return(Rdisplay_01);								

}

double Rdisplay_01_Rlinear(double Rdisplay)
{
//both are normalized to [0,1]
//default 709 gamma = 2.2

	double Rlinear;

	Rlinear =  Rdisplay <= 0.081 ?  Rdisplay / 4.5 : 
									pow((Rdisplay + 0.099) / 1.099, 2.2); 

	return(Rlinear);

}

unsigned char Rdisplay_clipping(double Rdisplay_01)
{
//from Rdisplay in order [0,1] to [0,255]
//clipping examples: -0.01 -> 0 ; 1.12 -> 255

	double R;

	if(Rdisplay_01 < 0) Rdisplay_01 = 0;
	if(Rdisplay_01 > 1) Rdisplay_01 = 1;

	R = 0.5 + Rdisplay_01 * 255.;
	
	return( (unsigned char) R );

}

double Rdisplay_normalization(unsigned char Rdisplay)
{
//from [0,255] to [0,1]

	double R;

	R = (double) Rdisplay;
	
	return(R / 255.);

}


void RGBlinear_RGBdisplay(
	double Rlinear, 
	double Glinear, 
	double Blinear,
	BYTE * Rdisplay,
	BYTE * Gdisplay,
	BYTE * Bdisplay)
{
//===================================================
		
	double Rdisplay_01, Gdisplay_01, Bdisplay_01;

		Rdisplay_01 = Rlinear_Rdisplay_01(Rlinear);
		*Rdisplay   = Rdisplay_clipping(Rdisplay_01);

		Gdisplay_01 = Rlinear_Rdisplay_01(Glinear);
		*Gdisplay   = Rdisplay_clipping(Gdisplay_01);

		Bdisplay_01 = Rlinear_Rdisplay_01(Blinear);
		*Bdisplay   = Rdisplay_clipping(Bdisplay_01);

}


void RGBdisplay_RGBlinear(
	 BYTE Rdisplay,
	 BYTE Gdisplay,
	 BYTE Bdisplay,
	 double *Rlinear, 
	 double *Glinear, 
	 double *Blinear)
{
//===================================================
		
	double Rdisplay_01, Gdisplay_01, Bdisplay_01;

	Rdisplay_01 = Rdisplay_normalization(Rdisplay);
    *Rlinear = Rdisplay_01_Rlinear(Rdisplay_01);

	Gdisplay_01 = Rdisplay_normalization(Gdisplay);
    *Glinear = Rdisplay_01_Rlinear(Gdisplay_01);

	Bdisplay_01 = Rdisplay_normalization(Bdisplay);
    *Blinear = Rdisplay_01_Rlinear(Bdisplay_01);

}
