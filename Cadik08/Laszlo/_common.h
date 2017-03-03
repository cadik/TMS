//#include <conio.h>
#include <stdio.h>

//common cross transformation functions

#include "color_datae.h"
#include "classic_transformations.h"
#include "RGB_Gamma.h"
#include "spectrum_functions.h"
#include "COLOROID.h"
#include "dither.h"

//#include "magassagi_szin.h"


//from and to XYZ

//	the permitted INPUT:

//	1:	XYZ
//	2:	xyY

//	3:	Lab
//	4:	Luv

//	5:	ATV
//	6:	fiTV
//	7:	lambdaTV

//	8:	rgb_linear(709)
//	9:	RGB_monitor(709)
// 10:  new selection

//	OUTPUT:
//	xy purity lambda lambda_complementer XYZ
//  LabhC  LuvhsC  A,fi-T-V  fi_A_harmony_series rgb RGB


void _read_datae(void)
{




}

void XYZ_to_ALL(double X, 
				double Y, 
				double Z, 
				double Gamma_R, 
				double Gamma_G, 
				double Gamma_B, 
				double Gamma_ambient,
				FILE * out)
{

	double x, y, z, Xn, Yn, Zn,
		   L, a, b, h_ab, Chroma_ab, u, v, h_uv, s_uv, Chroma_uv,
		   A, T, V, white, black, color, lambda, fi, purity,
		   HL, Ha, Hb, A34, A130, A180, A230, A326,
		   R, G, B, r, g, bb;

	BYTE   Ok, Rdisplay, Gdisplay, Bdisplay;

	Xn = D65_X;
	Yn = 100.;
	Zn = D65_Z;

	XYZ_xyz(X, Y, Z, &x, &y, &z);

	Ok = xy_fi_p_lambda(x, y, &fi, &purity, &lambda);
	
	if(Ok == 1)
	//existing color
	{
		XYZ_Lab(X, Y, Z, Xn, Yn, Zn,
				&L, &a, &b, &h_ab, &Chroma_ab);

		XYZ_Hunter_Lab(X, Y, Z, &HL, &Ha, &Hb);
		//Hunter_Lab_XYZ(HL, Ha, Hb, &X, &Y, &Z);
		//printf("\n\tHUNTER XYZ   %lf  %lf  %lf",X,Y,Z);

		XYZ_Luv(X, Y, Z, Xn, Yn, Zn,
				&L, &u, &v, &h_uv, &s_uv, &Chroma_uv);

		XYZ_ATV(X, Y, Z,
				&A, &T, &V, &white, &black, &color, &lambda, &fi);

		HARMONY_SERIES(A, &A34, &A130, &A180, &A230, &A326);

		_XYZ_rgb(X, Y, Z, &r, &g, &bb);
		//printf("\n\tXYZ %lf %lf %lf",X,Y,Z);

		//_rgb_XYZ(r, g, bb, &X, &Y, &Z);
		//printf("\n\tXYZ %lf %lf %lf",X,Y,Z);

		rgb_RGB(r, g, bb, &R, &G, &B,
				Gamma_R, Gamma_G, Gamma_B, Gamma_ambient);

		//RGB_rgb(R, G, B, &r, &g, &bb,
		//		Gamma_R, Gamma_G, Gamma_B, Gamma_ambient);

		fprintf(out,"\n\tr = %6.3lf	g = %6.3lf	b = %6.3lf\n",r,g,bb);
		fprintf(out,"\n\tR = %lf  G = %lf  B = %lf",R,G,B);

		//XYZ_RGB709(X, Y, Z, &R, &G, &B);	
		//RGBlinear_RGBdisplay(R, G, B, &Rdisplay, &Gdisplay, &Bdisplay);
	}

	fprintf(out,"\n\n\n\tXYZ - xyz - p   [ p(D65) = 0, p(spectrum) = 1]\n\t-------------------------------------------------------");
	fprintf(out,"\n\tX = %7.3lf	Y = %7.3lf	Z = %7.3lf\n",X,Y,Z);
	fprintf(out,"\n\tx = %7.3lf	y = %7.3lf	p = %9.6lf\n\n",x,y,purity);
	if(Ok == 0)
	{
		fprintf(out,"NOT EXISTING COLOR");
		return;
	}

	fprintf(out,"\n\t\CIE Lab - Luv\n\t-------------------------------------------------------");
	fprintf(out,"\n\tL = %7.3lf	a = %7.3lf	b = %7.3lf\n",L,a,b);
	fprintf(out,"\n\th_ab = %7.3lf	Chroma_ab = %7.3lf\n",h_ab,Chroma_ab);

	fprintf(out,"\n\tHunter L = %7.3lf  a = %7.3lf   b = %7.3lf\n",HL,Ha,Hb);

	fprintf(out,"\n\tL = %7.3lf	u = %7.3lf	v = %7.3lf\n",L,u,v);
	fprintf(out,"\n\th_uv = %7.3lf	s_uv = %7.3lf	Chroma_uv = %7.3lf\n\n",h_uv,s_uv,Chroma_uv);

	fprintf(out,"\n\tColoroid\n\t-------------------------------------------------------");
	fprintf(out,"\n\tA = %7.2lf	T = %7.2lf	V = %7.2lf",A,T,V);
	fprintf(out,"\n\tfi = %7.2lf	lambda = %7.3lf nm",fi,lambda);
	fprintf(out,"\n\twhite = %7.3lf	black = %7.3lf	color = %7.3lf",white,black,color);
	fprintf(out,"\n\tA34=%7.2lf A130=%7.2lf A180=%7.2lf A230=%7.2lf A326=%7.2lf\n\n",
		             A34, A130, A180, A230, A326);

	//fprintf(out,"\n\trgb linear [0,1] - RGB monitor [0,255] with gamma 2.2\n\t-------------------------------------------------------");
	//fprintf(out,"\n\tr = %6.3lf	g = %6.3lf	b = %6.3lf\n",r,g,bb);
	//fprintf(out,"\n\tR = %3u		G = %3u		B = %3u",R,G,B);

}

BYTE ATV_to_XYZLab(long A, 
				   long T, 
				   long V, 
				   FILE * out)
{

	double X, Y, Z, x, y, z, Xn, Yn, Zn,
		   L, a, b, h_ab, Chroma_ab, fi, lambda, purity;

	BYTE   Ok;
	char   _;

	_ = '_';

	Xn = D65_X;
	Yn = 100.;
	Zn = D65_Z;

	if(ATV_XYZ((double)A, (double)T, (double)V, 
							   &X, &Y, &Z) == 0) return(0);

	if(XYZ_xyz(X, Y, Z, &x, &y, &z) == 0) return(0);

	//if(ITERACIO(X,Y,Z, 0.005, 0.98, 30, 600)==0) return(0);

	Ok = xy_fi_p_lambda(x, y, &fi, &purity, &lambda);
	
	if(Ok == 1)
	//existing color
	{

		XYZ_Lab(X, Y, Z, Xn, Yn, Zn,
				&L, &a, &b, &h_ab, &Chroma_ab);

	if(T >= 10 && V >= 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%2ld%c%2ld",_,A,_,T,_,V);

	if(T >= 10 && V < 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%2ld%c%1ld",_,A,_,T,_,V);

	if(T < 10 && V >= 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%1ld%c%2ld",_,A,_,T,_,V);

	if(T < 10 && V < 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%1ld%c%1ld",_,A,_,T,_,V);

	fprintf(out,"\t%.3lf\t%.3lf\t%.3lf",X,Y,Z);

	fprintf(out,"\t%.3lf\t%.3lf\t%.3lf",L,a,b);

	}
	else return(0);
	
return(1);
	
}

void Coloroid_lap(long A, 
				  long T_also, 
				  long d_T, 
				  long T_felso,
				  long V_also, 
				  long d_V,
				  long V_felso,
				  FILE * out)
{

	long i, j;

	fprintf(out,"\nBEGIN_DATA_FORMAT\nSample_Name\tXYZ_X\tXYZ_Y\tXYZ_Z\tLab_L\tLab_a\tLab_b\nEND_DATA_FORMAT\nBEGIN_DATA",A);

	for(i = V_also; i<= V_felso; i+= d_V)
	{
		for(j = T_also; j<= T_felso; j+= d_T)
			if(ATV_to_XYZLab(A, j, i, out)==0) break;

	}

	fprintf(out,"\nEND_DATA");

}

BYTE ATV___to_XYZLab(long A, 
				   long T, 
				   double V, 
				   FILE * out)
{

	double X, Y, Z, x, y, z, Xn, Yn, Zn,
		   L, a, b, h_ab, Chroma_ab, fi, lambda, purity;

	BYTE   Ok;
	char   _;

	_ = '_';

	Xn = D65_X;
	Yn = 100.;
	Zn = D65_Z;

	if(ATV_XYZ((double)A, (double)T, (double)V, 
							   &X, &Y, &Z) == 0) return(0);

	if(XYZ_xyz(X, Y, Z, &x, &y, &z) == 0) return(0);

	Ok = xy_fi_p_lambda(x, y, &fi, &purity, &lambda);
	
	if(Ok == 1)
	//existing color
	{

		XYZ_Lab(X, Y, Z, Xn, Yn, Zn,
				&L, &a, &b, &h_ab, &Chroma_ab);

	if(T >= 10 && V >= 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%2ld",_,A,_,T);

	if(T >= 10 && V < 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%2ld",_,A,_,T);

	if(T < 10 && V >= 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%1ld",_,A,_,T);

	if(T < 10 && V < 10)
	fprintf(out,"\nColoroid_ATV%c%2ld%c%1ld",_,A,_,T);

	fprintf(out,"\t%.3lf\t%.3lf\t%.3lf",X,Y,Z);

	fprintf(out,"\t%.3lf\t%.3lf\t%.3lf",L,a,b);

	}
	else return(0);
	
return(1);
	
}

void Coloroid_szinkor(FILE * out)
{

	short A_, Tt;
	long i, j, k;
	double X,Y,Z,E_;
	FILE * in;

	in = fopen("Kor_v70.txt","rt");

	fprintf(out,"\nBEGIN_DATA_FORMAT\nSample_Name\tXYZ_X\tXYZ_Y\tXYZ_Z\tLab_L\tLab_a\tLab_b\nEND_DATA_FORMAT\nBEGIN_DATA",A);

	E_ = 10.*sqrt(70.);
	for(k=0;k<48;k++)
	{
		fscanf(in,"%d%d",&A_,&Tt);
		printf("\n   %ld  %ld  %ld", A_,Tt,(long)E_);
		//getch();

		while(1)
		{
			if(ATV_XYZ((double)A_, (double)Tt, E_,  &X,&Y,&Z)==0)
			{
				printf("\n  k = %ld  ***",k);
				Tt *= 0.99;
			}
			else break;
		}
		
		ATV___to_XYZLab((double)A_, (double)Tt, E_, out);
	}

	fclose(in);
	
	fprintf(out,"\nEND_DATA");

}

void _7_basic_fi_computation(void)
{
	//7_basic_fi[8]
	int i;
	double _FI;

	for(i=1; i<=7; i++)
	{
		Coloroid_A_fi(10.*i, &_FI);
		_7_basic_fi[i]= _FI;
	//printf("\n\n\tA = %5.2lf     FI = %5.2lf", 10.*i, _7_basic_fi[i]);
	//	getch();
	}

}

void GRAY_EQUI(double X,	double Y,	double Z, 
			   double *_A,	double *_T,	double *_V,
			   double *rel_T, /*only for hue formula */ int *T_lower,  int *T_upper,   double *mu_T_lower,
			   int *fi_lower, int *fi_upper, double *mu_fi_lower,
			   int *V_lower, int *V_upper, double *mu_V_lower) 
{
// for given XYZ after comuting the Coloroid ATV we give the values of A10_prev, A10_next, mu_prev, mu_next, rel_T

	BYTE feasible;
	int A1, AA1, AA2, i, j, j1, k, k1;
	double dA, T1, T2, TT, V1, V2, VV, Vil, YY, fi, fi1, fi2, mu, f, T_max;
	double wh, bl, col, lambda;

	feasible = XYZ_ATV(X, Y, Z, _A, _T, _V, &wh, &bl, &col, &lambda, &fi);
	//ismert A, T, fi, interpolalunk ket szomszedos alap-szinezeti lap kozott
	//j, j1 a szomszedos hue lapok szamai 1...7
	//k, k1 a szomszedos vilagossagertkek szamai, k = 0 a fekete, k vagy k1 = 1,2,3 a harom nevezetes V=45,65,85 ertek es k1 = 4 a feher
	
	if(feasible == 0) 
	{	
		*_A = 10;
		fi = _7_basic_fi[1];
		*_T = 0;
		*_V = 10.0 * sqrt(Y);
		//for security a gray pixel intsead of an invalid color
	}

	A1 = (int) *_A;
	dA = *_A - A1;  //az A2 sulya
	AA1 = Aa[A1];
	AA2 = AA1+1;
	if(AA2 == 49) AA2 = 1;
	//1...48 indexek

	T1 = MAX_TEL[AA1].T;
	T2 = MAX_TEL[AA2].T;
	TT = dA * T2 + (1.0 - dA) * T1;

	V1 = MAX_TEL[AA1].V;
	V2 = MAX_TEL[AA2].V;
	VV = dA * V2 + (1.0 - dA) * V1;

	//printf("\n\t A1= %ld AA1 = %ld   AA2 = %ld  T1=%g   T2=%g  TT=%g --------",A1, AA1, AA2, T1, T2, TT);

	Vil = *_V;
	if(Vil < 45) { k = 1; k1 = 2; goto cim1; }
		if(Vil < 65) { k = 2; k1 = 3; goto cim1; }
			if(Vil < 85) { k = 3; k1 = 4; goto cim1; }
						  k = 4; k1 = 5; // a feherrel zarodo eset

cim1:

	*V_lower = k;
	*V_upper = k1;
	*mu_V_lower = (VIL_GRAY[k1] - Vil)/ (VIL_GRAY[k1] - VIL_GRAY[k]);

/*			VIL_GRAY[5] = 100.;  //Y = 100
			VIL_GRAY[4] =  85.;  //Y = 72.25
			VIL_GRAY[3] =  65.;  //Y = 42.25
			VIL_GRAY[2] =  45.;  //Y = 20.25
			VIL_GRAY[1] =   0.;  //Y = 0
*/

	YY = 0.01 * VV * VV;

	if(Y >= YY) T_max = (100. - Y) / (100. - YY) * TT;
	else        T_max =         Y  /         YY  * TT;
	//itt T_max az adott folytonosan tetszoleges hue es adott az V ill. Y melletti legnagyobb telitettseg
	//_T a tenyleges telitettseg
	f = *rel_T = *_T/T_max * 5.;
	if(f >= 5) 
		{
			f = 5;
			*T_lower = 5; *T_upper = 5;
			*mu_T_lower = 0.;
		}
	else
		{
			i = (int) f;
			*T_lower = i;
			*T_upper = i+1;
			f = f - i;
			*mu_T_lower = 1.0 - f;
		}
//----------------------------------------------------------------

	if(*_A >= 60 && *_A < 70)
	{
		fi1 = _7_basic_fi[6] + 360.;
		fi2 = _7_basic_fi[7];
		if(fi < 0)  fi += 360.;
		j = 6;
		j1 = 7;

		goto cim2;
	}
	else
		for(i = 1 ; i <= 7; i++)
		{
				j = i;
				j1 = i+1;
				if(j1>7) j1=1;
			if(i != 6 && fi <= _7_basic_fi[j]+1E-5 &&  fi > _7_basic_fi[j1])
			{
				fi1 = _7_basic_fi[j];
				fi2 = _7_basic_fi[j1];
				goto cim2;
			}
		}

cim2:
		//j, j1 es f, fi1, fi2 ismertek

		*fi_lower = j;
		*fi_upper = j1;
		*mu_fi_lower = (fi - fi2) / (fi1 - fi2);

		//printf("\n j=%ld j1=%ld fi=%g fi1=%g   fi2= %g   mu_fi = %g",j,j1,fi, fi1, fi2, *mu_fi_lower);

		 //TELITETTSEG_GRAY[8][6][6],  7 szinezetre, 85, 65, 45 vilagossagra,  (1...5) * TELITETTSEG_LEPCSO-re
		 //VIL_GRAY[6]  novekvo 5 vilagossag-ertek
         //LEPCSO_SZAM[8][6] hany T-egyseg fer a max testszin telitettsegig
}


double GRAY_T_EQUI(double rel_T, int T_lower,  int T_upper,  double mu_T_lower,
							     int fi_lower, int fi_upper, double mu_fi_lower,
							     int V_lower,  int V_upper,  double mu_V_low)
{

//10*i,5. + j*20., TELITETTSEG_GRAY[i][k][j]

double relT_lower1, relT_lower2, relT_upper1, relT_upper2, 
	   relT_lower,  relT_upper /* V szerint upper - lower */, T_gray_interpolated;


	//int *T_lower,  int *T_upper,  double *mu_T_lower
	//int *fi_lower, int *fi_upper, double *mu_fi_lower
	//int *V_lower,  int *V_upper,  double *mu_V_low
	//tri-linear interpolation
   
		relT_lower1 = mu_T_lower  * TELITETTSEG_GRAY[fi_lower][T_lower][V_lower] +
		       (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_lower][T_upper][V_lower];
			   //lower V level, lower fi level

		relT_lower2 = mu_T_lower  * TELITETTSEG_GRAY[fi_upper][T_lower][V_lower] +
		       (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_upper][T_upper][V_lower];
			   	//lower V level, upper fi level

		relT_upper1 = mu_T_lower  * TELITETTSEG_GRAY[fi_lower][T_lower][V_upper] +
		       (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_lower][T_upper][V_upper];
			   //upper V level, lower fi level

		relT_upper2 = mu_T_lower  * TELITETTSEG_GRAY[fi_upper][T_lower][V_upper] +
		       (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_upper][T_upper][V_upper];
			   //upper V level, upper fi level

	   //from the 4 linear interpolations with bilinear will be the value obtained


	
	//after this point mu_T_lower will be NOT used more
	//now first at the upper and lower V levels

	relT_lower = mu_fi_lower  * relT_lower1 + 
		  (1.0 - mu_fi_lower) * relT_lower2;

	relT_upper = mu_fi_lower  * relT_upper1 + 
		  (1.0 - mu_fi_lower) * relT_upper2;

	//after this only a linear interpolation according to V levels

	T_gray_interpolated = mu_V_low  * relT_lower +
		           (1.0 - mu_V_low) * relT_upper;
	
	if(fabs(T_gray_interpolated) > 11)
	{
		//printf("\n    T_gray_interpolated = %g", T_gray_interpolated);     
		return(0);
	}

	return(T_gray_interpolated);

}


double quasi_sqrt_rel_T(double rel_T)
{
	//x = T_rel < 1:   0.5 * x,  else   sqrt(x) - 0.5;  ertek es derivaltfolytonos
	double x;

	if(rel_T < 0.) return (-1.0);
	//e r r o r

	x = 2 * rel_T;  //10 is the solid color max

	if(x < 1) return (0.5 * x);

	return(sqrt(x) - 0.5);

}

double GRAY_T_DIFF_EQUI_FOR_2_COLORS(double rel_T1, 
									 int T_lower1,  int T_upper1,  double mu_T_lower1,
									 int fi_lower1, int fi_upper1, double mu_fi_lower1,
							         int V_lower1,  int V_upper1,  double mu_V_lower1, 
									 double rel_T2, 
									 int T_lower2,  int T_upper2,  double mu_T_lower2,
									 int fi_lower2, int fi_upper2, double mu_fi_lower2,
							         int V_lower2,  int V_upper2,  double mu_V_lower2)
{
	double s, t;

	 s = 0.113681;
	
	 t  = GRAY_T_EQUI(rel_T2, T_lower2,  T_upper2,  mu_T_lower2,
							   fi_lower2, fi_upper2, mu_fi_lower2,
							   V_lower2,  V_upper2,  mu_V_lower2) -
							   
		  GRAY_T_EQUI(rel_T1, T_lower1,  T_upper1,  mu_T_lower1,
							   fi_lower1, fi_upper1, mu_fi_lower1,
							   V_lower1,  V_upper1,  mu_V_lower1);

	 return(s * t);
}

double GRAY_HUE_DIFF_EQUI(double rel_T1, double rel_T2,
					      int fi_lower1, int fi_upper1, double mu_fi_lower1,
						  int fi_lower2, int fi_upper2, double mu_fi_lower2)
{

	double h_gray, h11, h12, s;

	s = 0.146095;

	h11 = HUE_PAIRS[fi_lower1][fi_lower2] * mu_fi_lower2 +
		  
		  HUE_PAIRS[fi_lower1][fi_upper2] * (1.0 - mu_fi_lower2); 

	h12 = HUE_PAIRS[fi_upper1][fi_lower2] * mu_fi_lower2 +
		  
		  HUE_PAIRS[fi_upper1][fi_upper2] * (1.0 - mu_fi_lower2); 

	h_gray = h11 * mu_fi_lower1 + h12 * (1.0 - mu_fi_lower1);

	//printf("\n\th11=%g h12=%g h-gray=%g\n", h11, h12, h_gray);

	h_gray *= sqrt ( quasi_sqrt_rel_T(rel_T1) *  quasi_sqrt_rel_T(rel_T2));

	return(s * h_gray);

}

double LUMINANCE_GRAD(double X1,  double Y1,  double Z1,
				      double X2,  double Y2,  double Z2,
		    	      double *dA, double *dT, double *dV)
{

	int T_lower1, T_upper1, fi_lower1, fi_upper1, V_lower1, V_upper1,  
		T_lower2, T_upper2, fi_lower2, fi_upper2, V_lower2, V_upper2;
	double A1, T1, V1, rel_T1, mu_T_lower1, mu_fi_lower1, mu_V_lower1,
	       A2, T2, V2, rel_T2, mu_T_lower2, mu_fi_lower2, mu_V_lower2,
		   s1, s2;

	GRAY_EQUI(X1,Y1,Z1,&A1,&T1,&V1,
			   &rel_T1, &T_lower1,  &T_upper1,  &mu_T_lower1,
			            &fi_lower1, &fi_upper1, &mu_fi_lower1,
			            &V_lower1,  &V_upper1,  &mu_V_lower1);

	GRAY_EQUI(X2,Y2,Z2,&A2,&T2,&V2,
			   &rel_T2, &T_lower2,  &T_upper2,  &mu_T_lower2,
			            &fi_lower2, &fi_upper2, &mu_fi_lower2,
			            &V_lower2,  &V_upper2,  &mu_V_lower2);

	*dT = GRAY_T_DIFF_EQUI_FOR_2_COLORS(  rel_T1, 
										  T_lower1, T_upper1, mu_T_lower1,
										  fi_lower1, fi_upper1, mu_fi_lower1,
										  V_lower1,  V_upper1, mu_V_lower1, 
										  rel_T2, 
										  T_lower2, T_upper2, mu_T_lower2,
										  fi_lower2, fi_upper2, mu_fi_lower2,
										  V_lower2, V_upper2, mu_V_lower2 );

	*dA = GRAY_HUE_DIFF_EQUI( rel_T1, rel_T2,
					          fi_lower1, fi_upper1, mu_fi_lower1,
						      fi_lower2, fi_upper2, mu_fi_lower2);

	*dV = V2 - V1;

	return(*dV + WEIGHT_LIGHTNESS_CHROMINANCE * (*dA + *dT));
}


void regression_hue_saturation_weigths(void)
{
	double a11, a12, a21, a22, b1, b2, D, D1, D2, 
		   A_, T, V , TT, HH, MM, weightT, weightH,
		   X1, Y1, Z1, X2, Y2, Z2, dA, dT, dV;
	int i;

	a11 = a12 = a21 = a22 = b1 = b2 = 0;

		for(i = 1; i <= 16; i++)
		{
			MM = - parok_rel_gray[i];

			ATV_XYZ(parok[i].A1,parok[i].T1,parok[i].V1,
				    &X1,&Y1,&Z1);

			
			ATV_XYZ(parok[i].A2,parok[i].T2,parok[i].V2,
				    &X2,&Y2,&Z2);

			//printf("\n X1=%g Y1=%g Z1=%g\n", X1, Y1, Z1);
			//printf("\n X2=%g Y2=%g Z2=%g\n", X2, Y2, Z2);

		   LUMINANCE_GRAD(X1,  Y1,  Z1,
					      X2,  Y2,  Z2,
		    		      &dA, &dT, &dV);
		   HH = dA;
		   TT = dT;

		   a11 += TT*TT;
		   a12 = a21 += HH*TT;
		   a22 += HH*HH;

		   b1 += MM * TT;
		   b2 += MM * HH;

			D  = a11*a22 - a12*a21;
			D1 = b1 *a22 - a12* b2;
			D2 = a11* b2 -  b1*a21;

		   
		   //printf("\n i = %ld TT = %g  HH = %g  MM = %g\n", i, dT, dA, MM);
		   //getch();

		}
		   weightT = D1 / D;
		   weightH = D2 / D;

		   X1 = a11*weightT + a12 * weightH - b1;
		   Y1 = a21*weightT + a22 * weightH - b2;
		   printf("\n\n  weightT = %g  weightH = %g  hiba1 = %g  hiba2 = %g\n", weightT, weightH, X1, Y1);

		   	for(i = 1; i <= 16; i++)
			{
				MM = - parok_rel_gray[i];

				ATV_XYZ(parok[i].A1,parok[i].T1,parok[i].V1,
				    &X1,&Y1,&Z1);

			
				ATV_XYZ(parok[i].A2,parok[i].T2,parok[i].V2,
				    &X2,&Y2,&Z2);

				LUMINANCE_GRAD(X1,  Y1,  Z1,
					           X2,  Y2,  Z2,
		    		           &dA, &dT, &dV);

				printf("\n formula = %g  observation = %g", dA + dT + dV, MM);
				//getch();
			}

}
/*

for gray - invariance 

it doesnt exist in the human vision system or in the nature, it is an artificial abstarction,
but some partial effects can be perteptually based describe (the best tool e.g. the Coloroid)

Thereby the chrominance has a less weight than his color difference
The result image can not have bigger diffrence than the black-white difference of the media,
the col diff can be even 2 times larger

We have to conserve the white and black pixels (they are fixed)

The problem can not derive directy from some col-diff formulas, especially
not from ds-based ones, but from Coloroid view conditions 

practically all image contains before and after the transformation a set of isoluminant pixel pairs
which are undistinguishible, and had visible difference on the color image
The existance of this set is necessary for global and local methods,
but they have different characters

In a complex spatio-chromatic approach most probably exist a perceptually best approach,
but giving options by some free parameters, it can be obtained a large set of 
overemphesized visaulizations of the original image both returning in the color world or using gray-scale images
			
*/
