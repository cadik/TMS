
#include "BMP_.h"

//COLOROID functions

/*
double _p(long i1, long i2, double dx, double dy, 
double *p, double *q, double *lambda);

void _p_C(long i1, long i2, double dx, double dy, double *q);

void _p_D65(long i1, long i2, double dx, double dy, double *q);

void LINEG_MIX(double q,
			 double x1, double y1, double Y1,
			 double x2, double y2, double Y2,
			 double *x, double *y,
			 double *X, double *Y, double *Z,
			 double *p);
//this function is used to linear interpolation by fi
//or in 1 nm interval for lambda
//the XYZ values have p != q factor

void _xy_XYZ_48_based_on_fi_C(long index, long i1, long i2, double x, double y);
void _fi_lambda_C(long index, double fi);

void fi_Coloroid_A(double fi, double * A_hue);
BYTE Coloroid_A_fi(double A_hue, double * fi);

BYTE keres(short k1,short k2, double fi, double * lambda);

void fi_lambda(double fi, double * lambda);
void fi_lambda_complementer(double fi, double * lambda_complementer);
BYTE lambda_fi(double lambda, double * fi);
BYTE Coloroid_A_lambda(double A_hue, double * lambda);
BYTE lambda_Coloroid_A(double lambda, double * A_hue);

void fi_Coloroid_limes_color(double fi,
							 double * X, double * Y, double * Z,
							 double * x, double * y);
void fi_CIE1931_limes_color(double fi,
							double * X, double * Y, double * Z,
							double * x, double * y);
BYTE xy_fi_p_lambda(double x, double y, 
					double *fi, double *p, double *lambda);

BYTE Coloroid_decomposition(
	double X,
	double Y,
	double X_limes,
	double Y_limes,
	double *white,
	double *black,
	double *color);

BYTE XYZ_ATV(
		double X,
		double Y,
		double Z,
		double *A,
		double *T,
		double *V,
		double *white,
		double *black,
		double *color,
		double *lambda,
		double *fi);

BYTE ATV_XYZ(double  A, double  T, double V,
             double *X, double *Y, double *Z);

BYTE HARMONY_SERIES(double A,
					double *A34,
					double *A130,
					double *A180,
					double *A230,  //-130
					double *A326); // -34

void READ_COLOROID_DATAE(void);

*/

double _p(long i1, long i2, double dx, double dy, 
double *p, double *q, double *lambda)
{
// p szintisztasagot es lambda jellemzo hullamhosszat szamitunk
// p > 1 esetben a szin nem letezik

double a11, a12, a21, a22, b1, b2, D, D1, D2, t;

a11 = dx;
a21 = dy;

a12 = CIE1931_x[i1] - CIE1931_x[i2];
a22 = CIE1931_y[i1] - CIE1931_y[i2];

b1 = CIE1931_x[i1] - D65_x;
b2 = CIE1931_y[i1] - D65_y;

D  = a11*a22 - a12*a21;
D1 = b1 *a22 - a12* b2;
D2 = a11* b2 -  b1*a21;

 t = D1 / D;
*q = D2 / D;

*lambda = (double)i1 + *q;
*p = 1./t;

//q * (i2) + (1-q) * (i1)
 
return(1);
}

void _p_C(long i1, long i2, double dx, double dy, double *q)
//==========================================================
{
double a11, a12, a21, a22, b1, b2, D, D1, D2, t;

a11 = dx;
a21 = dy;

a12 = CIE1931_x[i1] - CIE1931_x[i2];
a22 = CIE1931_y[i1] - CIE1931_y[i2];

b1 = CIE1931_x[i1] - C_x;
b2 = CIE1931_y[i1] - C_y;

D  = a11*a22 - a12*a21;
D1 = b1 *a22 - a12* b2;
D2 = a11* b2 -  b1*a21;

 t = D1 / D;
*q = D2 / D;

return;
}

void _p_D65(long i1, long i2, double dx, double dy, double *q)
//============================================================
{
double a11, a12, a21, a22, b1, b2, D, D1, D2, t;

a11 = dx;
a21 = dy;

a12 = CIE1931_x[i1] - CIE1931_x[i2];
a22 = CIE1931_y[i1] - CIE1931_y[i2];

b1 = CIE1931_x[i1] - D65_x;
b2 = CIE1931_y[i1] - D65_y;

D  = a11*a22 - a12*a21;
D1 = b1 *a22 - a12* b2;
D2 = a11* b2 -  b1*a21;

 t = D1 / D;
*q = D2 / D;

return;
}


void LIN_MIX(double q,
			 double x1, double y1, double Y1,
			 double x2, double y2, double Y2,
			 double *x, double *y,
			 double *X, double *Y, double *Z,
			 double *p)
//===========================================
{
//(x,y) =   (1-q)*(x1,y1)    + q*(x2,y2),   q INPUT is known e.g. from fi
//(X,Y,X) = (1-p)*(X1,Y1,Z1) + p*(X2,Y2,Z2) 

double d, dX, dY, s, X1, Z1, X2, Z2;

*x = (1 - q) * x1 + q * x2;
*y = (1 - q) * y1 + q * y2;

xyY_XZ(x1,y1,Y1,&X1,&Z1);
xyY_XZ(x2,y2,Y2,&X2,&Z2);

s  = X1 + Y1 + Z1;
dX = X2 - X1;
dY = Y2 - Y1;
d  = dX + dY + Z2 - Z1;
    
if(fabs(dX) > fabs(dY))
*p = (*x * s - X1) / (dX - *x * d);   
else
*p = (*y * s - Y1) / (dY - *y * d);

*Y = (1 - *p) * Y1 +  *p * Y2;
xyY_XZ(*x,*y,*Y,X,Z);

}

void _xy_XYZ_48_based_on_fi_C(long index, long i1, long i2, double x, double y)
//=============================================================================
{
double xx,yy,Y,X,Z,p,q;

   _p_C(i1,i2,x,y,&q);

   LIN_MIX(q,CIE1931_x[i1],CIE1931_y[i1],100.*CIE1931_Y[i1],
			 CIE1931_x[i2],CIE1931_y[i2],100.*CIE1931_Y[i2],			 
			 &xx, &yy, &X, &Y, &Z, &p);

	xc[index] = xx;
	yc[index] = yy;

	Xc[index] = X;
	Yc[index] = Y;
	Zc[index] = Z;

	//basic datae of Coloroid system

}

void _fi_lambda_C(long index, double fi)
//======================================
{
long i,i1,i2;
double f,x,y;

x = cos(GRAD_RAD*fi);//  /10.;
y = sin(GRAD_RAD*fi);//  /10.;

f = fi;

if(fabs(f - fi_C_[450]) < EPS_6)
{
	xc[index] = CIE1931_x[450];
	yc[index] = CIE1931_y[450];

	Xc[index] = 100. * CIE1931_X[450];
	Yc[index] = 100. * CIE1931_Y[450];
	Zc[index] = 100. * CIE1931_Z[450];
	return;
}

if(fabs(f - fi_C_[625]) < EPS_6)
{
	xc[index] = CIE1931_x[625];
	yc[index] = CIE1931_y[625];

	Xc[index] = 100. * CIE1931_X[625];
	Yc[index] = 100. * CIE1931_Y[625];
	Zc[index] = 100. * CIE1931_Z[625];;
	return;
}

if(f<fi_C_[490] || f>fi_C_[491])
{
//fi elojelvaltasi helye
	i1 = 490;
	i2 = 491;
	_xy_XYZ_48_based_on_fi_C(index,i1,i2,x,y);
	return;
}

if(f<fi_C_[625] && f>fi_C_[450])
{
//biborszin tartomany
	i1 = 450;
	i2 = 625;
    _xy_XYZ_48_based_on_fi_C(index,i1,i2,x,y);
	return;
}

if(f<=fi_C_[450])
  {
	  for(i=451;i<=490;i++)
	  if(f>fi_C_[i])
	  {
		i1=i-1;
		i2=i;
		_xy_XYZ_48_based_on_fi_C(index,i1,i2,x,y);
		return;
	  }
  }


  for(i=491;i<=625;i++)
  if(f>fi_C_[i])
  {
    i1=i-1;
    i2=i;
    _xy_XYZ_48_based_on_fi_C(index,i1,i2,x,y);
    return;
   }

}

void fi_Coloroid_A(double fi, double * A_hue)
//===========================================
{

short  k, k1, k2;
double f, fi1, fi2, dfi, dfi12, q;

//  fi given in (-360, 360) interval
	f = fi;
	if(f <= -180) f += 360;
	if(f >   180) f -= 360;
	//normalization to (-180,180] interval

	if(f <= Coloroid_fi_D65[35] || f >= Coloroid_fi_D65[36])
	{
			fi1 = 360. + Coloroid_fi_D65[35];
			fi2 =		 Coloroid_fi_D65[36];
			if(f < 0)   f += 360;

			dfi12 =	fi1 - fi2;
				  
			dfi   =	fi1 - f;	

			q = dfi / dfi12; 

			*A_hue = 60. + q;
			
			return;
	}

	k1 = 35+48;
	k2 = 36;
	
	while(1)
	{
		k = (k1 + k2)/2;

		if(f < Coloroid_fi_D65[k>48?(k-48):k])
			k2 = k;
		else
			k1 = k;

		if(k1-k2 == 1) break;
	}
	
		dfi12 =	Coloroid_fi_D65[k1>48?(k1-48):k1]-
				Coloroid_fi_D65[k2>48?(k2-48):k2];
				  

		dfi   =	f - Coloroid_fi_D65[k2>48?(k2-48):k2];	

		q = dfi / dfi12;

		if(fabs(q - 1) < EPS_5)
		{
			k2++;
			q = 0;
		}
		
		*A_hue = (double)A[k2>48?(k2-48):k2] + q;
		
		return;
}

BYTE Coloroid_A_fi(double A_hue, double * fi)
//===========================================
{
	//not permitted A_hue examples: -56, 9, 38, 59, 78, 103
	//permitted A_hue values:        10, 11, 12, 12.46, 75.84, 76, 76.99
short k, k1, k2;
double fi1, fi2, q;

	if(A_hue < 10 || A_hue > 77) return(0);

	if(fabs(floor(A_hue) - A_hue) < EPS_5 || fabs(ceil(A_hue) - A_hue) < EPS_5)
	{
		k  = (short)(A_hue + 0.5);
		if(Aa[k] == -1) return(0);
		//this A number isnt permitted
		k  = Aa[k];
		*fi = Coloroid_fi_D65[k];
		return(1);
	}

	else

	{

		if((short)A_hue == 60)
		{
			fi1 = 360. + Coloroid_fi_D65[35];
			fi2 =		 Coloroid_fi_D65[36];
			q = A_hue - 60;
			*fi = (1-q) * fi1 + q * fi2;
			if(*fi > 180) *fi -= 360;
			//60 - 61 is the place of sign change

		}
		else
		{
			k1 = (short) A_hue;
			q = A_hue - k1;
			if(Aa[k1] == -1) return(0);
			//this A number isnt permitted
			k1 = Aa[k1];
			k2 = k1 + 1;
			if(k2 > 48) k2 = 1;
			fi1 = Coloroid_fi_D65[k1];
			fi2 = Coloroid_fi_D65[k2];
		
			*fi = (1-q) * fi1 + q * fi2;
		}

	}
	
	return(1);

}


BYTE keres(short k1,short k2, double fi, double * lambda)
//=======================================================
{
	//assumption: fi(k1) > fi(k2)

	double p, q, fi1, fi2, dx, dy;
	short k;

	fi1 = fi_D65[k1];
	fi2 = fi_D65[k2];

	if(fi1 < fi || fi2 > fi) return(0);

	while(1)
	{
		k = (k1 + k2)/2;

		if(fi < fi_D65[k])
			k1 = k;
		else
			k2 = k;

		if(k2-k1 == 1) break;
	}

	dx = cos(GRAD_RAD*fi);
	dy = sin(GRAD_RAD*fi);
	_p(k1,k2,dx,dy,&p,&q,lambda);
				
	return(1);
}


void fi_lambda(double fi, double * lambda)
//========================================
{
	//white: D65
	//for purple colors lambda is negativ:
	//with the vawelength of complementer color

	//maximum *lambda vawelength is 699 nm:
	//over 699 nm the fi is constant

	double f, fi0, fi1, fi2, fi699_l, fi699_u, q, p, dx, dy;

	//  fi given in (-360, 360) interval
	f = fi;
	if(f <= -180) f += 360;
	if(f >   180) f -= 360;
	//normalization to (-180,180] interval

	//[699-830]: small variance, not exactly constant CIE datae
	fi699_l = -8.5863733643;
	fi699_u = -8.5863553570;
	//699:    -8.586369
	//698:	  -8.584490      grad

	if(fi699_l < f && f < fi699_u)
	{
		*lambda = 699.;
		return;
	}

	if(f > fi_D65[360] && f < fi699_l)
	{
	//range of purple colors
		fi1 = f + 180;
		keres(493,567,fi1,lambda);
		*lambda = - *lambda;
		return;
	}

	if(f < fi_D65[491] || f > fi_D65[492])
	{
	//change of sign interval
		
		dx = cos(GRAD_RAD*fi);
		dy = sin(GRAD_RAD*fi);
		_p(491,492,dx,dy,&p,&q,lambda);
		//491,492 nm is the place of sign change
		return;
	}

	if(f <= fi_D65[360] && f >= fi_D65[491])
	{
		keres(360,491,f,lambda);
		return;
	}

	if(f <= fi_D65[492] && f > fi699_u)
	{
		keres(492,699,f,lambda);
		return;
	}

		
}

void fi_lambda_complementer(double fi, double * lambda_complementer)
//==================================================================
{
	//white: D65
	//for purple colors lambda is negativ:
	//with the vawelength of complementer color

	//maximum *lambda vawelength is 699 nm:
	//over 699 nm the fi is constant

	double fi_c;

	fi_c = fi + 180;
	if(fi_c > 180) fi_c -= 360;
	fi_lambda(fi_c, lambda_complementer);

	return;

}

BYTE lambda_fi(double lambda, double * fi)
//========================================
{
	//white: D65
	//for purple colors lambda is negativ:
	//with the vawelength of complementer color

	//maximum *lambda vawelength is 699 nm:
	//over 699 nm the fi is constant

	double fi1, fi2, X, Y, Z, x, y, z, q, p;
	short lam1,lam2;

	*fi = 0;

	if(lambda >=0) 
	{

		if(lambda < 360. || lambda > 830.) return(0);

		if(lambda >= 699. && lambda <=830.)
		{
			//constant hue in red
			*fi = fi_D65[699];
			return(1);
		}

		if(lambda >= 360. && lambda < 699.)
		{
			//normal spectrum colors
			lam1 = (short)lambda;
			lam2 = lam1 + 1;
			q = lambda - lam1;
			LIN_MIX(q,CIE1931_x[lam1],CIE1931_y[lam1],100.*CIE1931_Y[lam1],
				  CIE1931_x[lam2],CIE1931_y[lam2],100.*CIE1931_Y[lam2],			 
				  &x, &y, &X, &Y, &Z, &p);
			*fi = RAD_GRAD*atan2(y-D65_y,x-D65_x);

			return(1);
		}
	}

	if(lambda < 0)
	{
		if(lambda < - 566.42886 || lambda > - 493.331499/*493.33773*/)
			return(0);

		//lambda is negativ with possible complementer
		//according to the purpur colors
		lambda = - lambda;

		lam1 = (short)lambda;
		lam2 = lam1 + 1;
		q = lambda - lam1;
		LIN_MIX(q,CIE1931_x[lam1],CIE1931_y[lam1],100.*CIE1931_Y[lam1],
			  CIE1931_x[lam2],CIE1931_y[lam2],100.*CIE1931_Y[lam2],			 
			  &x, &y, &X, &Y, &Z, &p);

        *fi = RAD_GRAD*atan2(y-D65_y,x-D65_x) - 180;

		return(1);
	}

	return(1);
}

BYTE Coloroid_A_lambda(double A_hue, double * lambda)
{
//vawelengths to Coloroid A values
	double fi;

	if(Coloroid_A_fi(A_hue,&fi) == 0) return(0);
	fi_lambda(fi,lambda);

	return(1);
}

BYTE lambda_Coloroid_A(double lambda, double * A_hue)
{
//Coloroid A hue values to given vawelength
	double fi;

	if(lambda_fi(lambda,&fi) == 0) return(0);
	fi_Coloroid_A(fi,A_hue);

	return(1);

}

void fi_Coloroid_limes_color(double fi,
							 double * X, double * Y, double * Z,
							 double * x, double * y)
//=========================================================================
{

	double lambda, xx, yy, q, p, f;
	long i1,i2;

	//  fi given in (-360, 360) interval
	f = fi;
	if(f <= -180) f += 360;
	if(f >   180) f -= 360;
	//normalization to (-180,180] interval

	if(f <= fi_D65[625] && f >= fi_D65[450])
	{
	//purple color
		xx = cos(GRAD_RAD*f);
		yy = sin(GRAD_RAD*f);
		_p_D65(625,450,xx,yy,&q);
		i1 = 625;
		i2 = 450;
	}
	else
	{
	//[450,625] nm
		fi_lambda(f,&lambda);
		i1 = (long) lambda;
		i2 = i1+1;
		q  = lambda - i1;
	}

		LIN_MIX(q,CIE1931_x[i1],CIE1931_y[i1],100.*CIE1931_Y[i1],
				  CIE1931_x[i2],CIE1931_y[i2],100.*CIE1931_Y[i2],			 
				  x, y, X, Y, Z, &p);
}

void fi_CIE1931_limes_color(double fi,
							double * X, double * Y, double * Z,
							double * x, double * y)
//=========================================================================
{

	double xx, yy, lambda, q, t, f;
	long i1,i2;

	//  fi given in (-360, 360) interval
	f = fi;
	if(f <= -180) f += 360;
	if(f >   180) f -= 360;
	//normalization to (-180,180] interval

	if(f <= fi_D65[699] && f >= fi_D65[360])
	{
	//purple color
		xx = cos(GRAD_RAD*f);
		yy = sin(GRAD_RAD*f);
		_p_D65(699,360,xx,yy,&q);
		i1 = 699;
		i2 = 360;
	}
	else
	{
	//[360,699] nm
		fi_lambda(f,&lambda);
		i1 = (long) lambda;
		i2 = i1+1;
		q  = lambda - i1;
	}

		LIN_MIX(q,CIE1931_x[i1],CIE1931_y[i1],100.*CIE1931_Y[i1],
				  CIE1931_x[i2],CIE1931_y[i2],100.*CIE1931_Y[i2],			 
				  x, y, X, Y, Z, &t);

}


BYTE xy_fi_p_lambda(double x, double y, 
double *fi, double *p, double *lambda)
//=====================================
{
//D65-re fi szoget, p (purity) erteket es jellemzo hullamhosszat
//szamolunk. Ha p > 1, nem letezik a szin. Ha lambda negativ -> 
//az adott biborszin komplementerenek hullamhosszat adjuk meg negativ elojellel
//xs,ys a jellemzo spektrumszin koordinatai

double x1,y1,X,Y,Z,xL,yL;

if(x<0 || y<0)
	{
		 *p = -1;
		 return(0);
	}

     x1 = x-D65_x;
     y1 = y-D65_y;

     if(fabs(x1) > EPS_6)
          *fi = RAD_GRAD * atan2(y1,x1);
	 else
          *fi = y1 >= 0 ? 90. : -90.;

	 fi_lambda(*fi,lambda);

	 fi_CIE1931_limes_color(*fi,&X,&Y,&Z,&xL,&yL);

	 xL -= D65_x;
	 yL -= D65_y;

	 
	 if(fabs(xL) > fabs(yL))
		 *p = fabs(x1)/fabs(xL);
	 else
		 *p = fabs(y1)/fabs(yL);

	 if(*p<0 || *p > 1.0001) 
	 {
		 *p = -1;
		 return(0);
	 }
	 
	
	 return(1);

}

BYTE Coloroid_decomposition(
	double X,
	double Y,
	double X_limes,
	double Y_limes,
	double *white,
	double *black,
	double *color)
{
// assumtion: limes color of XY is X_limes,Y_limes
// return = 1 =>   feasible Colorid color
// return = 0 => infeasible Coloroid color (but perhaps existing color)

// T = 100 * c (Coloroid saturation)
// white + black + color =1   ADDITIV mixing
// valid for D65 white point

double D, D1, D2, X_white, Y_white;

X_white = (100./D65_y) * D65_x; //95.0449218;   //X_D65
Y_white = 100.;         //Y_D65
//Ideal black: X = Y = Z = 0

D  = X_white * Y_limes - X_limes * Y_white;
D1 = X       * Y_limes - X_limes * Y;

*white = D1 / D;
if(fabs(*white) < EPS_8) *white = 0;

D2 = X_white * Y       - X       * Y_white;
*color = D2 / D;
if(fabs(*color) < EPS_8) *color = 0;

*black = 1 - (*white + *color);
if(fabs(*black) < EPS_8) *black = 0;

return(1);
}

BYTE XYZ_ATV(
		double X,
		double Y,
		double Z,
		double *A,
		double *T,
		double *V,
		double *white,
		double *black,
		double *color,
		double *lambda,
		double *fi)
{
	double x, y, z, X_limes, Y_limes, Z_limes, x_limes, y_limes;

	if(XYZ_xyz(X,Y,Z,&x,&y,&z)==0) return(0);
	
	*fi = RAD_GRAD * atan2(y-D65_y, x-D65_x);
	fi_lambda(*fi,lambda);
	
	fi_Coloroid_limes_color(*fi,&X_limes,&Y_limes,&Z_limes,
				                &x_limes,&y_limes);

	Coloroid_decomposition(X, Y, X_limes, Y_limes,
							white, black, color);

	fi_Coloroid_A(*fi, A);

	*T = 100. * *color;
	*V = 10. * sqrt(Y);

	return(1);

}


BYTE ATV_XYZ(double  A, double  T, double V,
             double *X, double *Y, double *Z)
{

	double x, y, X_limes, Y_limes, Z_limes, x_limes, y_limes,
		A1,T1,V1,white,black,color,lambda,p;
    double fi, s, t;
 

		if(V < EPS_5  || T < -EPS_5 || V > 100.01 || T > 120) return(0);
		//for high V and T is Coloroid not valid
		//but T and V can the valu 100 exceed

		*Y = V * V * 0.01;

		if(Coloroid_A_fi(A, &fi) == 0) return(0);

		fi_Coloroid_limes_color(fi,&X_limes,&Y_limes,&Z_limes,
				                &x_limes,&y_limes);

		if(T > EPS_8)
        {
		        s = epsw * (V * V - T * Y_limes);
                t = T * (X_limes + Y_limes + Z_limes);

                x = (D65_x * s + x_limes * t) / (s + t);

                y = (V * V + y_limes * t - T * Y_limes) / (s + t);

        }
        else
        {
                x = D65_x;
                y = D65_y;
        }

		if(x < EPS_8 || y < EPS_8) return(0);

        xyY_XZ(x,y,*Y,X,Z);    
				
		xy_fi_p_lambda(x,y,&fi,&p,&lambda);

	
	if(fabs (p + 1.0) < 1E-5 ) return(0);

	return(1);
}

BYTE fiTV_XYZ(double fi, 
			  double T, 
			  double V,
			  double X_limes, 
			  double Y_limes, 
			  double Z_limes,
			  double x_limes,
			  double y_limes,
              double *X, 
			  double *Y, 
			  double *Z)
{
//ismert és megengedett fi esetén

	double x, y, A1,T1,V1;
    double s, t;
 

		if(V < EPS_5  || T < -EPS_5 || V > 100.01 || T > 120) return(0);
		//for high V and T is Coloroid not valid
		//but T and V can the valu 100 exceed

		*Y = V * V * 0.01;

		if(T > EPS_8)
        {
		        s = epsw * (V * V - T * Y_limes);
                t = T * (X_limes + Y_limes + Z_limes);

                x = (D65_x * s + x_limes * t) / (s + t);

                y = (V * V + y_limes * t - T * Y_limes) / (s + t);

        }
        else
        {
                x = D65_x;
                y = D65_y;
        }

		if(x < EPS_8 || y < EPS_8) return(0);

        xyY_XZ(x,y,*Y,X,Z);    
				
	return(1);
}


BYTE HARMONY_SERIES(double A,
					double *A34,
					double *A130,
					double *A180,
					double *A230, //-130
					double *A326) // -34
{

	double fi;

	if(Coloroid_A_fi(A, &fi) == 0) return(0);

	fi_Coloroid_A(fi+ 34, A34);

	fi_Coloroid_A(fi+130,A130);

	fi_Coloroid_A(fi+180,A180);
	
	fi_Coloroid_A(fi-130,A230);

	fi_Coloroid_A(fi- 34,A326);

	return(1);

}

void _hue(double A_48_as, double * A_rendes)
{

	long AA, i;
	double AAtort;

	if(A_48_as < 1) A_48_as += 48;

	AA		=	(int)A_48_as;
	AAtort	=	A_48_as - AA; 

	if(AAtort > 0.99) 
	{
		AA++;
		AAtort = 0;
	}

	i = AA%48;
	if(i==0) i=48;

	*A_rendes = A[i] + AAtort;

	if(*A_rendes > 76.99) *A_rendes = 10;

	return;
	
}

double	hue_diff(double A,
				 double diff)
{

	double result, Atort, AA, alap;
	long Aint;

	alap = A;
	Aint	= Aa[(int)alap];
	Atort	= alap - (int)alap;
	AA		= Aint + Atort; //[1,...,48.999)

	AA	+=	diff;

	_hue(AA, &result);

	return(result);
	
}

double	hue_complement_diff(double A,
							double diff1,
							double diff2)
{

	double result, Atort, AA, A180, alap, fi;
	long Aint;

	alap = A;
	alap = hue_diff(alap, diff1);
	Coloroid_A_fi(alap, &fi);
	fi_Coloroid_A(fi+180,&A180);

	Aint	= Aa[(int)A180];
	Atort	= A180 - (int)A180;
	AA		= Aint + Atort; //[1,...,48.999)

	AA	+=	diff2;

	_hue(AA, &result);

	return(result);

}
			
void READ_COLOROID_DATAE(void)
{
		
		double f, fi, lambda;
		long i,j;
		FILE *inp, *tab;
  
        for(i=0;i<77;i++)
        Aa[i] = -1;
              
        inp =  fopen("Cadik08/Laszlo/Coloroid.txt","rt");
        tab =  fopen("Cadik08/Laszlo/coloroid_datae.txt","wt");
		//coltab=fopen("Coloroid.dat","wb");
        for(i=1; i<=48; i++)
        {

        //fscanf(inp,"%d %f\n",&j,&f);
		A[i]  = j = (long) _read_data(inp);
        Aa[j] = i;
        fi_C[i] = f = _read_data(inp);

			_fi_lambda_C(i,f);

			fi = Coloroid_fi_D65[i] =
				RAD_GRAD * atan2(yc[i]-D65_y,xc[i]-D65_x);
			fi_lambda(fi,&lambda);
        
			fprintf(tab,"\n A=%2ld %8.3lf nm fi_D=%10.5lf x=%7.5lf y=%7.5lf Y=%8.5lf fi_C=%10.5lf",
				A[i],lambda,fi,xc[i],yc[i],Yc[i],fi_C[i]);

				/*		
				k = A[i];
				fwrite(&k,1,2,coltab);

				s = lambda;
				fwrite(&s,1,8,coltab);

				s = xc[i];
				fwrite(&s,1,8,coltab);
				s = yc[i];
				fwrite(&s,1,8,coltab);

				s = Xc[i];
				fwrite(&s,1,8,coltab);
				s = Yc[i];
				fwrite(&s,1,8,coltab);
				s = Zc[i];
				fwrite(&s,1,8,coltab);

				s = fi_C[i];
				fwrite(&s,1,8,coltab);

				s = Coloroid_fi_D65[i];
				fwrite(&s,1,8,coltab);
				/**/		
		
        }

        fclose(inp);
		fclose(tab);
		//fclose(coltab);
}

void READ_BASIC_COLORS(void)
{

	FILE *inp, *tab;
	double i;
	double _A, _T, _V, X, Y, Z;

        inp =  fopen("basic_color.txt","rt");
        tab =  fopen("basic_colors_53.txt","wt");

        for(i=1; i<=53; i++)
        {
			
		_A = _read_data(inp);
		_T = _read_data(inp);
		_V = _read_data(inp);

		ATV_XYZ(_A, _T, _V, &X, &Y, &Z);
        
		fprintf(tab,"\n\t%6.2lf\t%6.2lf\t%6.2lf\t%7.3lf\t%7.3lf\t%7.3lf",
				_A,_T,_V,X,Y,Z);

		}

		fclose(inp);
		fclose(tab);
}


//INNEN KEZDODIK AZ UJ RESZ
//RGB_trangle...  EGYSZER az eletben (csak display fuggo)
//global -> R_fi,R_X,R_Y, ... (ld. color_datae.h)

void RGB_trangle_precomputation(void)
{

		double X, Y, Z, x ,y ,z, xx, yy;

		_rgb_XYZ(1,0,0,&R_X,&R_Y,&R_Z);
		xx = x_r - D65_x;
		yy = y_r - D65_y;
		R_fi = RAD_GRAD*atan2(yy,xx);

		_rgb_XYZ(0,1,0,&G_X,&G_Y,&G_Z);
		xx = x_g - D65_x;
		yy = y_g - D65_y;
		G_fi = RAD_GRAD*atan2(yy,xx);

		_rgb_XYZ(0,0,1,&B_X,&B_Y,&B_Z);
		xx = x_b - D65_x;
		yy = y_b - D65_y;
		B_fi=RAD_GRAD*atan2(yy,xx);

}

//stretch egy COLOROID hue (A) lapra csak egyszer kell 
//return-je Tmax lehetseges max Coloroid-i telitettseg

double stretch(double X,
			   double Y,
			   double Z,
			   double *X_lim,
			   double *Y_lim,
			   double *Z_lim,
			   double *fi_lim)
{

		double  x, y, z, xx, yy, fi, A_, T, V, wh, bl, co, la,
				x1, x2, y1, y2, a11, a21, a12, a22, D, D2, //D1,
				b1, b2, p, q, X1, Y1, Z1, X2, Y2, Z2,
				s, s1, s2, d, dX, dY;

		XYZ_xyz(X,Y,Z,&x,&y,&z);
		xx = x - D65_x;
		yy = y - D65_y;

		if(fabs(xx) < EPS_6)
		{
			if(yy >=0) fi =  90.;
			else       fi = -90.;

		}
		else

			fi = RAD_GRAD*atan2(yy,xx);

		if(fabs(fi - R_fi) < EPS_6) 
		{
			XYZ_ATV(R_X,R_Y,R_Z,&A_,&T,&V,&wh,&bl,&co,&la,&fi);
			return(T);
		}

		if(fabs(fi - G_fi) < EPS_6) 
		{
			XYZ_ATV(G_X,G_Y,G_Z,&A_,&T,&V,&wh,&bl,&co,&la,&fi);
			return(T);
		}

		if(fabs(fi - B_fi) < EPS_6) 
		{
			XYZ_ATV(B_X,B_Y,B_Z,&A_,&T,&V,&wh,&bl,&co,&la,&fi);
			return(T);
		}

		if(fi > R_fi && fi < G_fi)
		{
			x1 = x_r;
			y1 = y_r;
			Y1 = R_Y;

			x2 = x_g;
			y2 = y_g;
			Y2 = G_Y;

			goto eq;
		}

		if(fi < R_fi && fi > B_fi)
		{
			x1 = x_r;
			y1 = y_r;
			Y1 = R_Y;

			x2 = x_b;
			y2 = y_b;
			Y2 = B_Y;

			goto eq;
		}
		else
		{
			x1 = x_b;
			y1 = y_b;
			Y1 = B_Y;

			x2 = x_g;
			y2 = y_g;
			Y2 = G_Y;

			goto eq;
		}

eq:
		a11 = xx;
		a21 = yy;

		a12 = x1-x2;
		a22 = y1-y2;

		b1 = x1 - D65_x;
		b2 = y1 - D65_y;

		D  = a11*a22 - a12*a21;
//		D1 = b1 *a22 - a12* b2;
		D2 = a11* b2 -  b1*a21;

//		p = D1 / D;
		q = D2 / D;

		x = (1-q) * x1 + q * x2;
		y = (1-q) * y1 + q * y2;

		s = Y1/y1;
		X1 = s * x1;
		Z1 = s * (1 - x1 - y1);

		s = Y2/y2;
		X2 = s * x2;
		Z2 = s * (1 - x2 - y2);		

		s  = X1 + Y1 + Z1;
		dX = X2 - X1;
		dY = Y2 - Y1;
		d  = dX + dY + Z2 - Z1;
    
		if(fabs(dX) > fabs(dY))
		p = (x * s - X1) / (dX - x * d);   
		else
		p = (y * s - Y1) / (dY - y * d);

		s1 = 1/p;
		s2 = 1/(1-p);
		s = s1 < s2 ? s1 : s2;

		Y = (1 - p) * Y1 +  p * Y2;
		*Y_lim = s * Y;
		
		xyY_XZ(x,y,*Y_lim,X_lim,Z_lim);

		printf("\n\tX = %lf  Y = %lf  Z=%lf",
			*X_lim, *Y_lim, *Z_lim);


		XYZ_ATV(*X_lim,*Y_lim,*Z_lim,&A_,&T,&V,&wh,&bl,&co,&la,fi_lim);
		printf("\n\tdisplay lim color: A=%lf  T=%lf  V=%lf",
			A_,T,V);

		return(T);

}

void T_V_x_y_koord(long		y_size,	//0,...,y_size  //0,...,100,
				   double	T_max,  //ha fixen 100, akkor nem nyujtunk
				   double   T,
				   double   V,
				   long		*x_koord,    //0,...,x_size     
				   long		*y_koord)
				  

{

	*x_koord = (long) (T * (y_size / T_max) + 0.5);

	*y_koord = (long) (V * (y_size /  100.) + 0.5);

}

void Coloroid_BMP(long		y_size,	//0,...,y_size  //0,...,100,
				  double	T_max,  //ha fixen 100, akkor nem nyujtunk
				  double	Coloroid_A,	//melyik szinezet lapja
				  long		x_koord,    //0,...,x_size     
				  long		y_koord,
				  double	*r,
				  double	*g,
				  double	*b,
				  double    *T,
				  double    *V)
{
//output : a linearis (r,g,b) az adott pixelen, 
//hatterpont eseten D65 feher 

//T_max == 100, akkor nem nyujtunk,
//egyik lap csak 38 - ig, egy masik pedig 114-ig mehet telitettesegben
//ekkor x_size legyen >= 1.15*y_size
//ezt nem ellenorizzuk

//T_max != 100, akkor nyujtunk, negyzet alaku hasznos terulet
//ekkor "x_size" = y_size

double	s1, s2, X, Y, Z,
		fi, X_limes, Y_limes, Z_limes, x_limes, y_limes;

	Coloroid_A_fi(Coloroid_A, &fi);

	fi_Coloroid_limes_color(fi,&X_limes,&Y_limes,&Z_limes,
				                &x_limes,&y_limes);

	s1 = 100. / y_size;
	s2 = T_max / 100;

	*T = s1 * s2 * x_koord;
	*V = s1 *      y_koord;
	if(*V<1E-2) *V=1E-2;

		if(fiTV_XYZ(fi, *T, *V, X_limes, Y_limes, Z_limes,
			x_limes, y_limes, &X, &Y, &Z)==0)
	{
		*r = *g = *b = 1.;
		return;
	}

	_XYZ_rgb(X,Y,Z,r,g,b);

	if(*r > 1 || *g > 1 || *b > 1 || *r < 0  || *g < 0 || *b < 0)

		*r = *g = *b = 1.;

		
// Tmax != 100 esetben
// a BMP szamitasok elott egyszer a kovetkezo fuggvenyeket le
// kell futtatni
//			RGB_trangle_precomputation();
//			T_max = stretch(X,Y,Z,&X_lim,&Y_lim,&Z_lim,&fi_lim);
//ez utobbi egy szinezeti lapra, az elobbi univerzalisan ervenyes

}


void nyujtott_proba_BMP(long y_res,
						double Coloroid_A)
{

FILE * fileki;
double T_max,X_lim,Y_lim,Z_lim,fi_lim,r,g,b,X,Y,Z,T,V,
		A_hue,col,bl,wh,lambda,fi,
		dR, dG, dB;
long x, y, x_res, dummy;
pixel pix;
BYTE cc;

cc=0;
x_res = y_res;

	RGB_trangle_precomputation();
	ATV_XYZ(Coloroid_A,5,50,&X,&Y,&Z);
	T_max = stretch(X,Y,Z,&X_lim,&Y_lim,&Z_lim,&fi_lim);

    fileki = fopen("Coloroid.bmp","wb");

	BMP_WRITE_HEADER_24(x_res+1,y_res+1,&dummy,fileki);

	printf("\n COL_A=%lf  T_max=%lf  dummy=%ld\n",
	Coloroid_A,T_max,dummy);

	for(y=0;y<=y_res;y++)
    {
        for(x=0;x<=x_res;x++)
        {
	
			Coloroid_BMP(	y_res,
							T_max,
							Coloroid_A,	
							x,
							y, 
							&r,
							&g,
							&b,
							&T,
							&V);

			rgb_RGB(r, g, b, &dR, &dG, &dB, 1.5, 1.5, 1.5, 1.);

			pix.r = (BYTE) (dR + 0.5);
			pix.g = (BYTE) (dG + 0.5);
			pix.b = (BYTE) (dB + 0.5);

			fwrite(&pix.b, 1, 1, fileki);
			fwrite(&pix.g, 1, 1, fileki);
			fwrite(&pix.r, 1, 1, fileki);

        }//x

	    fwrite(&cc,1,dummy,fileki);
              
   }//y
 
   fclose(fileki);
}
