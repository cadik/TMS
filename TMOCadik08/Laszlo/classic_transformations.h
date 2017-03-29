
//basic color transformations
//XYZ xyY Lab Luv

//BYTE XYZ_xyz(double  X, double  Y, double  Z,
//             double *x, double *y, double *z);

//BYTE xyY_XZ(double  x, double  y, double Y,
//            double *X, double *Z);

//long CIE_L_from_Y(double Y, double Yn, double * L);
//long CIE_Y_from_L(double L, double Yn, double * Y);

//BYTE XYZ_Lab(double  X,  double Y,  double Z,
//             double  Xn, double Yn, double Zn,
//             double *L,  double *a, double *b,
//             double *h_ab,          double *Chroma_ab);

//long Lab_XYZ(double  L,  double a,  double b,
//             double  Xn, double Yn, double Zn,
//             double  *X, double *Y, double *Z);

//long XYZ_Luv(double X,  double Y,  double Z,
//             double Xn, double Yn, double Zn,
//             double *L, double *u, double *v,
//             double *h_uv, double *s_uv, double *Chroma_uv);

//long Luv_XYZ(double L,  double u,  double v,
//             double Xn, double Yn, double Zn,
//             double *X, double *Y, double *Z);

//BYTE LABNHU(double x, double y,  double Y,
//            double*L, double *A, double *B);

//BYTE Hunter_Lab(double X, double Y, double Z,
//				  double*L, double*a, double*b);

//---------------------------------------------------------------------

BYTE XYZ_xyz(double  X, double  Y, double  Z,
             double *x, double *y, double *z)
//===========================================
{

	double t;

	if(X < EPS_5 || Y < EPS_5 || Z < EPS_5) return(0);

	t = X + Y + Z;
	
	*x = X / t;
	*y = Y / t;
	*z = Z / t;

	return(1);

}

BYTE xyY_XZ(double  x, double  y, double Y,
            double *X, double *Z)
//=========================================
{
	double s, z;

	if(x < EPS_8 || y < EPS_8 || Y < EPS_5) return (0);

	z = 1-(x+y);

	s = Y / y;

	*X = s * x;
	*Z = s * z;

	return(1);

}

long CIE_L_from_Y(double Y, double Yn, double * L)
{
//================================================
	double y_rel;

	//0 <= y_rel = Y/Yn <= 100
	//hier permitted y_rel > 1

	if(Y < 0 || Yn < EPS_8) return(0);

	y_rel = Y / Yn;

	if(y_rel < 0.008856)
			*L = 903.3 * y_rel;
	else
			*L = 116. * pow(y_rel,_third) - 16.;

	return(1);
}

long CIE_Y_from_L(double L, double Yn, double * Y)
{
	double y_rel;

	if(L < 0) return(0);

	if(L < 0.008856*903.3)
	   y_rel = L / 903.3;

	else
	   y_rel = pow((L + 16.)/116. , 3.);

	*Y = y_rel * Yn;

	return(0);
}

BYTE XYZ_Lab(double  X,  double Y,  double Z,
             double  Xn, double Yn, double Zn,
             double *L,  double *a, double *b,
             double *h_ab,          double *Chroma_ab)
{

	double x_rel, y_rel, z_rel, xx, yy, zz;

	if(Xn < EPS_8 || Yn < EPS_8 || Zn < EPS_8) return(0);
	//bad values for nominally white

	if(CIE_L_from_Y(Y,Yn,L) == 0) return(0);
	//Y < 0 || Yn < EPS_8

	x_rel = X/Xn;
	y_rel = Y/Yn;
	z_rel = Z/Zn;

	xx = pow(x_rel,_third);
	yy = pow(y_rel,_third);
	zz = pow(z_rel,_third);

	*a = 500. * (xx - yy);
	*b = 200. * (yy - zz);

	*h_ab = RAD_GRAD * atan2(*b,*a);
	//hue

	*Chroma_ab = sqrt(*a**a + *b**b);
	//chroma

	return(1);
}

long Lab_XYZ(double  L,  double a,  double b,
             double  Xn, double Yn, double Zn,
             double  *X, double *Y, double *Z)
{
//============================================
	
	double x_rel, y_rel, z_rel, xx, yy, zz;

	if(L < 0 || Yn < EPS_8) return(0);

	CIE_Y_from_L(L,Yn,Y);
	
	y_rel = *Y/Yn;
	yy = pow(y_rel,_third);

	xx = a / 500. + yy;
	zz = yy - b / 200.;

	x_rel = pow(xx,3.);
	*X = x_rel * Xn;

	z_rel = pow(zz,3.);
	*Z = z_rel * Zn;

	return(1);
}

long XYZ_Luv(double X,  double Y,  double Z,
             double Xn, double Yn, double Zn,
             double *L, double *u, double *v,
             double *h_uv, double *s_uv, double *Chroma_uv)
{
//=========================================================

	double u_, u_n, v_, v_n, d_, d_n, t;

	if(Xn < EPS_8 || Yn < EPS_8 || Zn < EPS_8) return(0);
	//bad values for nominally white

	if(CIE_L_from_Y(Y,Yn,L) == 0) return(0);
	//Y < 0 || Yn < EPS_8
	
	d_  = X  + 15 * Y  + 3 * Z;
	d_n = Xn + 15 * Yn + 3 * Zn;

	u_  = 4 * X  / d_;
	u_n = 4 * Xn / d_n;

	v_  = 9 * Y  / d_;
	v_n = 9 * Yn / d_n;

	t = 13 **L;

	*u = t * (u_ - u_n);
	*v = t * (v_ - v_n);

	*h_uv = RAD_GRAD * atan2(*v,*u);
	//hue

	*Chroma_uv = sqrt(*u**u + *v**v);
	//chroma

	*s_uv = *Chroma_uv / *L;     
	//psychometric saturation

	return(1);
}

long Luv_XYZ(double L,  double u,  double v,
             double Xn, double Yn, double Zn,
             double *X, double *Y, double *Z)
{
//===========================================

	double d_n,u_n,v_n,u_,v_,t;

	if(L < 0 || Yn < EPS_8) return(0);

	CIE_Y_from_L(L,Yn,Y);
	
	d_n = Xn + 15 * Yn + 3 * Zn;
	u_n = 4 * Xn / d_n;
	v_n = 9 * Yn / d_n;

	t = 13. * L;

	u_ = u / t + u_n;
	v_ = v / t + v_n;

	*X = 9./4. **Y * u_ / v_;
	*Z = (4. **X / u_ - *X - 15 **Y) / 3.;

	return(1);

}

long XYZ_sqrtLuv(double X,  double Y,  double Z,
             double Xn, double Yn, double Zn,
             double *L, double *u, double *v,
             double *h_uv, double *s_uv, double *Chroma_uv)
{
//===========================================================
//gyokos, Coloroidi vilagossag es u,v koordinatak rendrszere

	double u_, u_n, v_, v_n, d_, d_n, t, y_rel;

	if(Xn < EPS_8 || Yn < EPS_8 || Zn < EPS_8) return(0);
	//bad values for nominally white

	if(Y < 0) return(0);

	y_rel = Y / Yn;

	*L = 100. * sqrt(Y);
	
	d_  = X  + 15 * Y  + 3 * Z;
	d_n = Xn + 15 * Yn + 3 * Zn;

	u_  = 4 * X  / d_;
	u_n = 4 * Xn / d_n;

	v_  = 9 * Y  / d_;
	v_n = 9 * Yn / d_n;

	t = 13 **L;

	*u = t * (u_ - u_n);
	*v = t * (v_ - v_n);

	*h_uv = RAD_GRAD * atan2(*v,*u);
	//hue

	*Chroma_uv = sqrt(*u**u + *v**v);
	//chroma

	*s_uv = *Chroma_uv / *L;     
	//psychometric saturation

	return(1);
}


BYTE LABNHU(double x, double y,  double Y,
            double*L, double *A, double *B)
{
//=========================================

		double  Av,Bv,kY,q,t,z;
    
		if(x < 0 || y < EPS_8 || Y < 0) return(0);
        t = 1./6.;
        q = - 1/6;

        z = 1. - x - y;

        Av = 0.25 * pow(x/y + t, _third);
        Bv =   q  * pow(z/y + t, _third);

        kY = 500. * pow(Y,_third);

        *L = 116. * pow( Y/100. , _third) - 16.;
        *A = kY * (Av - ANv);
        *B = kY * (Bv - BNv);

		return(1);

}

BYTE XYZ_Hunter_Lab(double X, double Y, double Z,
    		    double*L, double*a, double*b)
{
	if(Y < EPS_8) return(0);

	*L = 10. * sqrt(Y);

	*a = 175. * (1.02*X - Y )/(*L);
	*b =  70. * (Y - 0.847*Z)/(*L);

	return(1);

}

BYTE Hunter_Lab_XYZ(double L, double a, double b,
					double *X, double *Y, double *Z)
{
	double s;

	s = (L / 10.);
	*Y = s * s;

	*X = (a * s / 17.5 + *Y)/1.02;
	*Z = (*Y - b * s / 7.)/0.847;

	return (1); 
}

