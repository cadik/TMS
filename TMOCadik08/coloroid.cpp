// original authors: Laszlo Neumann,
//                   Martin Cadik
//
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "coloroid.h"

//______________________________________________________________________________
double Coloroid::luminanceGrad(double X1,  double Y1,  double Z1,
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

	GRAY_EQUI(X2, Y2, Z2, &A2, &T2, &V2,
	          &rel_T2, &T_lower2,  &T_upper2,  &mu_T_lower2,
	          &fi_lower2, &fi_upper2, &mu_fi_lower2,
	          &V_lower2,  &V_upper2,  &mu_V_lower2);

	*dT = GRAY_T_DIFF_EQUI_FOR_2_COLORS(rel_T1, 
	                                    T_lower1, T_upper1, mu_T_lower1,
	                                    fi_lower1, fi_upper1, mu_fi_lower1,
	                                    V_lower1,  V_upper1, mu_V_lower1, 
	                                    rel_T2, 
	                                    T_lower2, T_upper2, mu_T_lower2,
	                                    fi_lower2, fi_upper2, mu_fi_lower2,
	                                    V_lower2, V_upper2, mu_V_lower2 );

	*dA = GRAY_HUE_DIFF_EQUI(rel_T1, rel_T2,
	                         fi_lower1, fi_upper1, mu_fi_lower1,
	                         fi_lower2, fi_upper2, mu_fi_lower2);

	*dV = V2 - V1;

	constexpr const double WEIGHT_LIGHTNESS_CHROMINANCE = 3.;
	return(*dV + WEIGHT_LIGHTNESS_CHROMINANCE * (*dA + *dT));
}

//==========================================================
void Coloroid::_p_C(long i1, long i2, double dx, double dy, double* q)
{
	double a11, a12, a21, a22, b1, b2, D, D1, D2, t;

	a11 = dx;
	a21 = dy;

	a12 = CIE1931_x[i1] - CIE1931_x[i2];
	a22 = CIE1931_y[i1] - CIE1931_y[i2];

	b1 = CIE1931_x[i1] - C_x;
	b2 = CIE1931_y[i1] - C_y;

	D  = a11 * a22 - a12 * a21;
	D1 = b1 * a22 - a12 * b2;
	D2 = a11 * b2 -  b1 * a21;

	 t = D1 / D;
	*q = D2 / D;
}

//============================================================
void Coloroid::_p_D65(long i1, long i2, double dx, double dy, double *q)
{
	double a11, a12, a21, a22, b1, b2, D, D1, D2, t;

	a11 = dx;
	a21 = dy;

	a12 = CIE1931_x[i1] - CIE1931_x[i2];
	a22 = CIE1931_y[i1] - CIE1931_y[i2];

	b1 = CIE1931_x[i1] - D65_x;
	b2 = CIE1931_y[i1] - D65_y;

	D  = a11 * a22 - a12 * a21;
	D1 = b1 * a22 - a12 * b2;
	D2 = a11 * b2 -  b1 * a21;

	 t = D1 / D;
	*q = D2 / D;
}

//=========================================================
bool Coloroid::xyY_XZ(double x, double y, double Y,
            double *X, double *Z)
{
	double s, z;

	if(x < EPS_8 || y < EPS_8 || Y < EPS_5)
		return false;

	z = 1 - (x + y);

	s = Y / y;

	*X = s * x;
	*Z = s * z;

	return true;
}

//===========================================
void Coloroid::LIN_MIX(double q, double x1, double y1, double Y1,
             double x2, double y2, double Y2,
             double *x, double *y, double *X, double *Y,
             double *Z, double *p)
{
	//(x,y) =   (1-q)*(x1,y1)    + q*(x2,y2),   q INPUT is known e.g. from fi
	//(X,Y,X) = (1-p)*(X1,Y1,Z1) + p*(X2,Y2,Z2) 

	double d, dX, dY, s, X1, Z1, X2, Z2;

	*x = (1 - q) * x1 + q * x2;
	*y = (1 - q) * y1 + q * y2;

	xyY_XZ(x1, y1, Y1, &X1, &Z1);
	xyY_XZ(x2, y2, Y2, &X2, &Z2);

	s  = X1 + Y1 + Z1;
	dX = X2 - X1;
	dY = Y2 - Y1;
	d  = dX + dY + Z2 - Z1;
	    
	if(fabs(dX) > fabs(dY))
		*p = (*x * s - X1) / (dX - *x * d);   
	else
		*p = (*y * s - Y1) / (dY - *y * d);

	*Y = (1 - *p) * Y1 +  *p * Y2;
	xyY_XZ(*x, *y, *Y, X, Z);
}

//=============================================================================
void Coloroid::_xy_XYZ_48_based_on_fi_C(long index, long i1, long i2, double x, double y)
{
	double xx, yy, Y, X, Z, p, q;

	_p_C(i1, i2, x, y, &q);

	LIN_MIX(q, CIE1931_x[i1], CIE1931_y[i1], 100. * CIE1931_Y[i1],
	        CIE1931_x[i2], CIE1931_y[i2], 100. * CIE1931_Y[i2],			 
	        &xx, &yy, &X, &Y, &Z, &p);

	xc[index] = xx;
	yc[index] = yy;

	Xc[index] = X;
	Yc[index] = Y;
	Zc[index] = Z;
}

//=============================================================================
double Coloroid::_read_data(FILE * file)
{
	double szam;

	long k, volt_szamjegy;

	for (k = 0; k < 800; ++k)
		jelsor[k] = 0;

	k = volt_szamjegy = 0;
		
	while (true) {
		jel = (unsigned char) fgetc(file);
		jelsor[k] = jel;
		k++;
		if ('0' <= jel && jel <= '9')
			++volt_szamjegy;
		if (volt_szamjegy > 0 && !(jel == '.' ||
		    (jel >= '0' && jel <= '9') ||
		    jel == '-' || jel == '+')) 
			break;
	}

	szam = atof(jelsor);

	return(szam);
}

//=============================================================================
void Coloroid::readSpectrumDatae()
{
	FILE * spectrum_in;
	short i;
	double s,t;
	//char c;

	spectrum_in = fopen("TMOCadik08/Xyz31_1.txt",  "rt");

	for (i = 360; i <= 830; ++i) {
		//fscanf(spectrum_in,"%d,%f,%f,%f",&vawelength,&x,&y,&z);

		_read_data(spectrum_in);
		CIE1931_X[i] = _read_data(spectrum_in);
		CIE1931_Y[i] = _read_data(spectrum_in);
		CIE1931_Z[i] = _read_data(spectrum_in);

		//ezek 100.-szorosai a hatarszinek XYZ koordinatai

		t = CIE1931_X[i] + CIE1931_Y[i] + CIE1931_Z[i];
		t = 1./t;

		CIE1931_x[i] = t * CIE1931_X[i];
		CIE1931_y[i] = t * CIE1931_Y[i];
		CIE1931_z[i] = t * CIE1931_Z[i];

			
		t = (CIE1931_y[i]-D65_y) / (CIE1931_x[i]-D65_x);
		fi_D65[i]=RAD_GRAD * atan2( (CIE1931_y[i]-D65_y), (CIE1931_x[i]-D65_x) );

		fi_C_[i]=RAD_GRAD * atan2(CIE1931_y[i]-C_y, CIE1931_x[i]-C_x);

	}

	fclose(spectrum_in);

	spectrum_in = fopen("TMOCadik08/LIGHTD65.TXT","rt");

	for(i=300;i<=830;i++)
	{
		_read_data(spectrum_in);
		light_D65[i] = _read_data(spectrum_in);
	}
	
	fclose(spectrum_in);

	for (i = 360, s = 0; i <= 830; ++i)
		s += light_D65[i] * CIE1931_Y[i];

		t = 100. / s;

	for (i = 360; i <= 830; ++i)
		light_D65_normalized[i] = light_D65[i] * t;

	//for(i=360,s=0;i<=830;i++)
	//	s += light_D65_normalized[i] * CIE1931_y[i];
	//ennek 100 az Y erteke konstans 1 reflektancia mellett

	spectrum_in = fopen("TMOCadik08/light_A.txt","rt");

	for (i = 300; i <= 830; ++i) {
		_read_data(spectrum_in);
		light_A[i] = _read_data(spectrum_in);        
	}

	for (i = 360, s = 0; i <= 830; ++i)
		s += light_A[i] * CIE1931_Y[i];

		t = 100. / s;

	for (i = 360 ; i <= 830; ++i)
		light_A_normalized[i] = light_A[i] * t;
	
	fclose(spectrum_in);
}

//===========================================
bool Coloroid::Coloroid_A_fi(double A_hue, double* fi)
{
	//not permitted A_hue examples: -56, 9, 38, 59, 78, 103
	//permitted A_hue values:        10, 11, 12, 12.46, 75.84, 76, 76.99
	short k, k1, k2;
	double fi1, fi2, q;

	if (A_hue < 10 || A_hue > 77)
		return false;

	if (fabs(floor(A_hue) - A_hue) < EPS_5 ||
	    fabs(ceil(A_hue) - A_hue) < EPS_5) {
		k  = (short)(A_hue + 0.5);
		if (Aa[k] == -1)
			return false;
		//this A number isnt permitted
		k = Aa[k];
		*fi = Coloroid_fi_D65[k];
		return true;
	}
	else {
		if (static_cast<short>(A_hue) == 60) {
			fi1 = 360. + Coloroid_fi_D65[35];
			fi2 = Coloroid_fi_D65[36];
			q = A_hue - 60;
			*fi = (1-q) * fi1 + q * fi2;
			if (*fi > 180)
				*fi -= 360;
			//60 - 61 is the place of sign change

		}
		else {
			k1 = (short) A_hue;
			q = A_hue - k1;
			if (Aa[k1] == -1)
				return false;
			//this A number isnt permitted
			k1 = Aa[k1];
			k2 = k1 + 1;
			if (k2 > 48)
				k2 = 1;
			fi1 = Coloroid_fi_D65[k1];
			fi2 = Coloroid_fi_D65[k2];
		
			*fi = (1 - q) * fi1 + q * fi2;
		}

	}
	
	return true;
}

//=======================================================
void Coloroid::_7_basic_fi_computation()
{
	//7_basic_fi[8]
	int i;
	double _FI;

	for(i = 1; i <= 7; i++) {
		Coloroid_A_fi(10. * i, &_FI);
		_7_basic_fi[i] = _FI;
	}

}

//=======================================================
void Coloroid::_p(long i1, long i2, double dx, double dy, 
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

	D  = a11 * a22 - a12 * a21;
	D1 = b1 * a22 - a12 * b2;
	D2 = a11 * b2 -  b1 * a21;

	 t = D1 / D;
	*q = D2 / D;

	*lambda = (double) i1 + *q;
	*p = 1. / t;

	//q * (i2) + (1-q) * (i1)
}

//=======================================================
bool Coloroid::search(short k1, short k2, double fi,
            double* lambda)
{
	//assumption: fi(k1) > fi(k2)

	double p, q, fi1, fi2, dx, dy;
	short k;

	fi1 = fi_D65[k1];
	fi2 = fi_D65[k2];

	if (fi1 < fi || fi2 > fi)
		return false;

	while (true) {
		k = (k1 + k2) / 2;

		if(fi < fi_D65[k])
			k1 = k;
		else
			k2 = k;

		if (k2 - k1 == 1)
			break;
	}

	dx = cos(GRAD_RAD*fi);
	dy = sin(GRAD_RAD*fi);
	_p(k1, k2, dx, dy, &p, &q, lambda);

	return true;
}

//========================================
void Coloroid::fi_lambda(double fi, double* lambda)
{
	// white: D65
	// for purple colors lambda is negative:
	// with the vawelength of complementer color

	// maximum *lambda vawelength is 699 nm:
	// over 699 nm the fi is constant

	double f, fi0, fi1, fi2, fi699_l, fi699_u,
	       q, p, dx, dy;

	//  fi given in (-360, 360) interval
	f = fi;
	if (f <= -180)
		f += 360;
	if (f > 180)
		f -= 360;
	//normalization to (-180,180] interval

	//[699-830]: small variance, not exactly constant CIE datae
	fi699_l = -8.5863733643;
	fi699_u = -8.5863553570;
	//699:    -8.586369
	//698:	  -8.584490      grad

	if (fi699_l < f && f < fi699_u) {
		*lambda = 699.;
		return;
	}

	if (f > fi_D65[360] && f < fi699_l) {
		//range of purple colors
		fi1 = f + 180;
		search(493, 567, fi1, lambda);
		*lambda = - *lambda;
		return;
	}

	if (f < fi_D65[491] || f > fi_D65[492]) {
		//change of sign interval
		dx = cos(GRAD_RAD * fi);
		dy = sin(GRAD_RAD * fi);
		_p(491, 492, dx, dy, &p, &q, lambda);
		//491,492 nm is the place of sign change
		return;
	}

	if (f <= fi_D65[360] && f >= fi_D65[491]) {
		search(360, 491, f, lambda);
		return;
	}

	if (f <= fi_D65[492] && f > fi699_u)
		search(492, 699, f, lambda);
}


//======================================
void Coloroid::_fi_lambda_C(long index, double fi)
{
	long i,i1,i2;
	double f,x,y;

	x = cos(GRAD_RAD * fi);//  /10.;
	y = sin(GRAD_RAD * fi);//  /10.;

	f = fi;

	if (fabs(f - fi_C_[450]) < EPS_6) {
		xc[index] = CIE1931_x[450];
		yc[index] = CIE1931_y[450];

		Xc[index] = 100. * CIE1931_X[450];
		Yc[index] = 100. * CIE1931_Y[450];
		Zc[index] = 100. * CIE1931_Z[450];
		return;
	}

	if (fabs(f - fi_C_[625]) < EPS_6) {
		xc[index] = CIE1931_x[625];
		yc[index] = CIE1931_y[625];

		Xc[index] = 100. * CIE1931_X[625];
		Yc[index] = 100. * CIE1931_Y[625];
		Zc[index] = 100. * CIE1931_Z[625];;
		return;
	}

	if (f<fi_C_[490] || f>fi_C_[491]) {
		//fi elojelvaltasi helye
		i1 = 490;
		i2 = 491;
		_xy_XYZ_48_based_on_fi_C(index, i1, i2, x, y);
		return;
	}

	if(f < fi_C_[625] && f > fi_C_[450]) {
		//biborszin tartomany
		i1 = 450;
		i2 = 625;
		_xy_XYZ_48_based_on_fi_C(index, i1, i2, x, y);
		return;
	}

	if(f <= fi_C_[450]) {
		for(i = 451; i <= 490; ++i)
			if(f > fi_C_[i]) {
				i1 = i - 1;
				i2 = i;
				_xy_XYZ_48_based_on_fi_C(index, i1, i2, x, y);
				return;
			}
	  }


	for (i = 491; i<=625; ++i)
		if(f > fi_C_[i]) {
			i1 = i - 1;
			i2 = i;
			_xy_XYZ_48_based_on_fi_C(index, i1, i2, x, y);
			return;
		}
}

void Coloroid::readColoroidDatae()
{
	double f, fi, lambda;
	long i,j;
	FILE *inp, *tab;
  
        for (i = 0; i < 77; ++i)
        	Aa[i] = -1;
              
        inp =  fopen("TMOCadik08/Coloroid.txt","rt");
        tab =  fopen("TMOCadik08/coloroid_datae.txt","wt");
	//coltab=fopen("Coloroid.dat","wb");
        for (i = 1; i <= 48; ++i) {
        	//fscanf(inp,"%d %f\n",&j,&f);
		A[i]  = j = (long) _read_data(inp);
		Aa[j] = i;
		fi_C[i] = f = _read_data(inp);

		_fi_lambda_C(i,f);

		fi = Coloroid_fi_D65[i] = RAD_GRAD * atan2(yc[i] - D65_y, xc[i] - D65_x);
		fi_lambda(fi, &lambda);

		fprintf(tab,"\n A=%2ld %8.3lf nm fi_D=%10.5lf x=%7.5lf y=%7.5lf Y=%8.5lf fi_C=%10.5lf",
		        A[i], lambda, fi, xc[i], yc[i], Yc[i], fi_C[i]);

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
		fwrite(&s,1,8,coltab);*/
        }

	fclose(inp);
	fclose(tab);
	//fclose(coltab);
}

//===========================================
void Coloroid::fi_Coloroid_A(double fi, double * A_hue)
{
	short  k, k1, k2;
	double f, fi1, fi2, dfi, dfi12, q;

	//  fi given in (-360, 360) interval
	f = fi;
	if (f <= -180)
		f += 360;
	if (f > 180)
		f -= 360;
	//normalization to (-180,180] interval

	if (f <= Coloroid_fi_D65[35] || f >= Coloroid_fi_D65[36]) {
			fi1 = 360. + Coloroid_fi_D65[35];
			fi2 = Coloroid_fi_D65[36];
			if (f < 0)
				f += 360;

			dfi12 =	fi1 - fi2;
			dfi = fi1 - f;	

			q = dfi / dfi12; 

			*A_hue = 60. + q;

			return;
	}

	k1 = 35 + 48;
	k2 = 36;

	while (true) {
		k = (k1 + k2) / 2;

		if (f < Coloroid_fi_D65[k > 48 ? (k - 48) : k])
			k2 = k;
		else
			k1 = k;

		if (k1 - k2 == 1)
			break;
	}

	dfi12 =	Coloroid_fi_D65[k1 > 48 ? (k1 - 48) : k1] -
	        Coloroid_fi_D65[k2 > 48 ? (k2 - 48) : k2];

	dfi = f - Coloroid_fi_D65[k2 > 48 ? (k2 - 48) : k2];	

	q = dfi / dfi12;

	if (fabs(q - 1) < EPS_5) {
		k2++;
		q = 0;
	}

	*A_hue = (double) A[k2 > 48 ? (k2 - 48) : k2] + q;
	
	return;
}

//______________________________________________________________________________
bool Coloroid::readColor2GrayDatae()
{
	FILE *inp, *tab;
	long i, j, k, i1, j1, tel;
	double g, s, _A, _T, _V, _V45, _V65, _V85, _T45, _T65, _T85;

        inp =  fopen("TMOCadik08/COLOR_GRAY_UNIF.txt","rt");
        tab =  fopen("TMOCadik08/GOLOR2GRAY_datae.txt","wt");


	for(i = 1; i <= 16; i++) {
		pair_rel_gray[i] = _read_data(inp);
		pairs[i].A1 = _read_data(inp);
		pairs[i].T1 = _read_data(inp);
		pairs[i].V1 = _read_data(inp);

		pairs[i].A2 = _read_data(inp);
		pairs[i].T2 = _read_data(inp);
		pairs[i].V2 = _read_data(inp);
		g = _read_data(inp);

		if (fabs(g + pair_rel_gray[i]) > 1E-8)
			return true;

		fprintf(tab,"\n\n\ti=%2ld\trel szurke =  %5.1lf\tA1=%5.1lf\tT1=%5.1lf\tV1=%5.1lf\tA2=%5.1lf\tV2=%5.1lf\tT2=%5.1lf",
		        i, pair_rel_gray[i], pairs[i].A1, pairs[i].T1, pairs[i].V1, pairs[i].A2, pairs[i].V2, pairs[i].T1);
	}

	//--------------------------------------------------------------------------------------------------
	for (i = 1; i <= 7; ++i)
		for (j = 1; j <= 7; ++j) {
			_read_data(inp);
			_read_data(inp);
			HUE_PAIRS[i][j] = s = _read_data(inp);
			if (s < 0)
				sig[i][j] = -1;
			else
				sig[i][j] = 1;
			if (i == j)
				HUE_PAIRS[i][j] = 0.;
		}

	for (i = 1; i <= 7; i++)
		for (j=1; j < i; j++) {
			s = fabs(HUE_PAIRS[i][j]) + fabs(HUE_PAIRS[j][i]);
			s /= 2.;
			HUE_PAIRS[i][j] = s * sig[i][j];
			HUE_PAIRS[i][j] = - HUE_PAIRS[j][i];

			//fuggveny hue parokra a T_rel relativ telitettsegek 0...5 es azon t u l i  ertekeire
			//x = T_rel < 1:   0.5 * x,  else   sqrt(x) - 0.5;  ertek es derivaltfolytonos ************************
		}

	for (i = 1; i <= 7; ++i) {
		fprintf(tab,"\n");
		//for(j=1; j <=7; j++)
		//	fprintf(tab,"\n\t\tA1=%2ld\t\tA2=%2ld\t\trel szurke =  %5.2lf",
		//		10*i,10*j,HUE_PAIRS[i][j]);
		
		fprintf(tab,"\n\t\tA = %2ld\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf\t%5.2lf",
			10 * i, HUE_PAIRS[i][1], HUE_PAIRS[i][2], HUE_PAIRS[i][3], HUE_PAIRS[i][4], HUE_PAIRS[i][5],
			HUE_PAIRS[i][6], HUE_PAIRS[i][7]);
	}

	//--------------------------------------------------------------------------------------------------
	fprintf(tab,"\n");
	for (i = 1; i <= 48; ++i) {
			//int A[49],Aa[77];
			MAX_TEL[i].T = _read_data(inp);
			MAX_TEL[i].V = s = _read_data(inp);
			MAX_TEL[i].Y = 0.01 * s * s;
			fprintf(tab,"\n\t\ MAX  A =%3ld\tT =%5.1lf\tV =%5.1lf\tY =%5.1lf",
			        A[i], MAX_TEL[i].T, MAX_TEL[i].V, MAX_TEL[i].Y);
	}

	fprintf(tab,"\n");
	//TELITETTSEG_GRAY[8][6][6];

	for (i = 1; i <= 7; ++i)
		for (j = 1; j <= 5; ++j)
			for (k = 1; k <= 5; ++k)
				TELITETTSEG_GRAY[i][j][k] = 0;

	for (i = 1; i <= 7; ++i) {
		fprintf(tab,"\n");
		i1 = 10 * i;
		for (j = 4; j >= 2; j--) {
			fprintf(tab,"\n");
			for (k = 1; k <= 5; ++k) {
				_read_data(inp);
				_read_data(inp);

				s = _read_data(inp);
				//if(k==1) TELITETTSEG_LEPCSO[i][j] = s;
				//helyette a max solid color telitettseg 5-ode az egyseg

				TELITETTSEG_GRAY[i][k][j] = _read_data(inp);
				fprintf(tab,"\n\t\tA =%3ld\tV = %5.2lf\t rel szurke = %5.2lf",
				        10 * i, 5. + j * 20., TELITETTSEG_GRAY[i][k][j] /*ATV sorrend*/);
			}//k
		}//j
	}//i

	VIL_GRAY[5] = 100.;  //Y = 100
	VIL_GRAY[4] =  85.;  //Y = 72.25
	VIL_GRAY[3] =  65.;  //Y = 42.25
	VIL_GRAY[2] =  45.;  //Y = 20.25
	VIL_GRAY[1] =   0.;  //Y = 0


	fclose(inp);
	fclose(tab);

	return false;
}

//=========================================================================
void Coloroid::fi_Coloroid_limes_color(double fi, double * X, double * Y, double * Z,
                             double * x, double * y)
{
	double lambda, xx, yy, q, p, f;
	long i1,i2;

	//  fi given in (-360, 360) interval
	f = fi;
	if (f <= -180)
		f += 360;
	if (f > 180)
		f -= 360;
	//normalization to (-180,180] interval

	if (f <= fi_D65[625] && f >= fi_D65[450]) {
		//purple color
		xx = cos(GRAD_RAD * f);
		yy = sin(GRAD_RAD * f);
		_p_D65(625, 450, xx, yy, &q);
		i1 = 625;
		i2 = 450;
	}
	else {
		//[450,625] nm
		fi_lambda(f, &lambda);
		i1 = (long) lambda;
		i2 = i1 + 1;
		q  = lambda - i1;
	}

	LIN_MIX(q, CIE1931_x[i1], CIE1931_y[i1], 100. * CIE1931_Y[i1],
	        CIE1931_x[i2], CIE1931_y[i2], 100. * CIE1931_Y[i2],
	        x, y, X, Y, Z, &p);
}

//===========================================
void Coloroid::Coloroid_decomposition(double X, double Y, double X_limes, double Y_limes, double *white, double *black, double *color)
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
	D1 = X * Y_limes - X_limes * Y;

	*white = D1 / D;
	if (fabs(*white) < EPS_8)
		*white = 0;

	D2 = X_white * Y - X * Y_white;
	*color = D2 / D;
	if (fabs(*color) < EPS_8)
		*color = 0;

	*black = 1 - (*white + *color);
	if (fabs(*black) < EPS_8)
		*black = 0;
}

//===========================================
bool Coloroid::XYZ_xyz(double  X, double  Y, double  Z,
             double *x, double *y, double *z)
{
	double t;

	if (X < EPS_5 || Y < EPS_5 || Z < EPS_5)
		return false;

	t = X + Y + Z;
	
	*x = X / t;
	*y = Y / t;
	*z = Z / t;

	return true;
}

//===========================================
bool Coloroid::XYZ_ATV(double X, double Y, double Z,
             double *A, double *T, double *V,
             double *white, double *black,
             double *color, double *lambda, double *fi)
{
	double x, y, z, X_limes, Y_limes, Z_limes, x_limes, y_limes;

	if (XYZ_xyz(X, Y, Z, &x, &y, &z) == 0)
		return false;
	
	*fi = RAD_GRAD * atan2(y - D65_y, x - D65_x);
	fi_lambda(*fi,lambda);

	fi_Coloroid_limes_color(*fi, &X_limes, &Y_limes, &Z_limes,
	                        &x_limes,&y_limes);

	Coloroid_decomposition(X, Y, X_limes, Y_limes, white, black, color);

	fi_Coloroid_A(*fi, A);

	*T = 100. * *color;
	*V = 10. * sqrt(Y);

	return true;
}

//===========================================
void Coloroid::rgb709Xyz(double R, double G, double B,
                double*X, double*Y, double*Z)
{
	// white Y = 100
	// white R,G,B = 1,  D65
	// [ X ] [ 0.412453  0.35758   0.180423 ] [ R709 ] 
	// [ Y ]=[ 0.212671  0.71516   0.072169 ]*[ G709 ] 
	// [ Z ] [ 0.019334  0.119193  0.950227 ] [ B709 ] 

	*X = 41.2453 * R + 35.758 * G + 18.0423 * B;
	*Y = 21.2671 * R + 71.516 * G +  7.2169 * B;
	*Z = 1.9334 * R + 11.9193 * G + 95.0227 * B;
}

//===========================================
void Coloroid::GRAY_EQUI(double X, double Y, double Z, 
               double *_A, double *_T, double *_V,
               double *rel_T, /*only for hue formula */ int *T_lower,
               int *T_upper,   double *mu_T_lower,
               int *fi_lower, int *fi_upper, double *mu_fi_lower,
               int *V_lower, int *V_upper, double *mu_V_lower) 
{
	// for given XYZ after comuting the Coloroid ATV we give the values of A10_prev, A10_next, mu_prev, mu_next, rel_T
	bool feasible;
	int A1, AA1, AA2, i, j, j1, k, k1;
	double dA, T1, T2, TT, V1, V2, VV, Vil, YY, fi, fi1, fi2, mu, f, T_max;
	double wh, bl, col, lambda;

	feasible = XYZ_ATV(X, Y, Z, _A, _T, _V, &wh, &bl, &col, &lambda, &fi);
	//ismert A, T, fi, interpolalunk ket szomszedos alap-szinezeti lap kozott
	//j, j1 a szomszedos hue lapok szamai 1...7
	//k, k1 a szomszedos vilagossagertkek szamai, k = 0 a fekete, k vagy k1 = 1,2,3 a harom nevezetes V=45,65,85 ertek es k1 = 4 a feher
	
	if (!feasible) {	
		*_A = 10;
		fi = _7_basic_fi[1];
		*_T = 0;
		*_V = 10.0 * sqrt(Y);
		//for security a gray pixel intsead of an invalid color
	}

	A1 = (int) *_A;
	dA = *_A - A1;  //az A2 sulya
	AA1 = Aa[A1];
	AA2 = AA1 + 1;
	if (AA2 == 49)
		AA2 = 1; //1...48 indexek

	T1 = MAX_TEL[AA1].T;
	T2 = MAX_TEL[AA2].T;
	TT = dA * T2 + (1.0 - dA) * T1;

	V1 = MAX_TEL[AA1].V;
	V2 = MAX_TEL[AA2].V;
	VV = dA * V2 + (1.0 - dA) * V1;

	//printf("\n\t A1= %ld AA1 = %ld   AA2 = %ld  T1=%g   T2=%g  TT=%g --------",A1, AA1, AA2, T1, T2, TT);

	Vil = *_V;
	if (Vil < 45) {
		k = 1;
		k1 = 2;
		goto cim1;
	}

	if (Vil < 65) {
		k = 2;
		k1 = 3;
		goto cim1;
	}
		
	if (Vil < 85) {
		k = 3;
		k1 = 4;
		goto cim1;
	}
	k = 4;
	k1 = 5; // a feherrel zarodo eset

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

	if (Y >= YY)
		T_max = (100. - Y) / (100. - YY) * TT;
	else
		T_max = Y  / YY * TT;
	//itt T_max az adott folytonosan tetszoleges hue es adott az V ill. Y melletti legnagyobb telitettseg
	//_T a tenyleges telitettseg
	f = *rel_T = *_T / T_max * 5.;
	if(f >= 5) {
		f = 5;
		*T_lower = 5;
		*T_upper = 5;
		*mu_T_lower = 0.;
	}
	else {
		i = (int) f;
		*T_lower = i;
		*T_upper = i + 1;
		f = f - i;
		*mu_T_lower = 1.0 - f;
	}

	//----------------------------------------------------------------
	if (*_A >= 60 && *_A < 70) {
		fi1 = _7_basic_fi[6] + 360.;
		fi2 = _7_basic_fi[7];
		if (fi < 0)
			fi += 360.;
		j = 6;
		j1 = 7;

		goto cim2;
	}
	else
		for(i = 1 ; i <= 7; ++i) {
			j = i;
			j1 = i + 1;
			if (j1 > 7)
				j1 = 1;
			if (i != 6 && fi <= _7_basic_fi[j] + 1E-5 && fi > _7_basic_fi[j1]) {
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

//______________________________________________________________________________
double Coloroid::GRAY_T_EQUI(double rel_T, int T_lower,  int T_upper,  double mu_T_lower,
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

	// TODO segfault
	// this was added >
	if (fi_lower >= 8)
		fi_lower = 7;
	else if (fi_lower < 0)
		fi_lower = 0;

	if (fi_upper >= 8)
		fi_upper = 7;
	else if (fi_upper < 0)
		fi_upper = 0;

	if (T_lower >= 6)
		T_lower = 5;
	else if (T_lower < 0)
		T_lower = 0;

	if (T_upper >= 6)
		T_upper = 5;
	else if (T_upper < 0)
		T_upper = 0;

	if (V_lower >= 6)
		V_lower = 5;
	else if (V_lower < 0)
		V_lower = 0;

	if (V_upper >= 6)
		V_upper = 5;
	else if (V_upper < 0)
		V_upper = 0;
	// < here it ends
   
	relT_lower1 = mu_T_lower  * TELITETTSEG_GRAY[fi_lower][T_lower][V_lower] +
	              (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_lower][T_upper][V_lower]; //lower V level, lower fi level

	relT_lower2 = mu_T_lower  * TELITETTSEG_GRAY[fi_upper][T_lower][V_lower] +
	              (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_upper][T_upper][V_lower]; //lower V level, upper fi level

	relT_upper1 = mu_T_lower  * TELITETTSEG_GRAY[fi_lower][T_lower][V_upper] +
	              (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_lower][T_upper][V_upper]; //upper V level, lower fi level

	relT_upper2 = mu_T_lower  * TELITETTSEG_GRAY[fi_upper][T_lower][V_upper] +
	              (1.0 - mu_T_lower) * TELITETTSEG_GRAY[fi_upper][T_upper][V_upper]; //upper V level, upper fi level

	//from the 4 linear interpolations with bilinear will be the value obtained
	//after this point mu_T_lower will be NOT used more
	//now first at the upper and lower V levels

	relT_lower = mu_fi_lower * relT_lower1 + (1.0 - mu_fi_lower) * relT_lower2;

	relT_upper = mu_fi_lower * relT_upper1 + (1.0 - mu_fi_lower) * relT_upper2;

	//after this only a linear interpolation according to V levels

	T_gray_interpolated = mu_V_low  * relT_lower + (1.0 - mu_V_low) * relT_upper;
	
	if (fabs(T_gray_interpolated) > 11)
		//printf("\n    T_gray_interpolated = %g", T_gray_interpolated);     
		return 0.;

	return T_gray_interpolated;
}

//______________________________________________________________________________
double Coloroid::GRAY_T_DIFF_EQUI_FOR_2_COLORS(double rel_T1, int T_lower1,  int T_upper1,  double mu_T_lower1,
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

//______________________________________________________________________________
double Coloroid::quasi_sqrt_rel_T(double rel_T)
{
	//x = T_rel < 1:   0.5 * x,  else   sqrt(x) - 0.5;  ertek es derivaltfolytonos
	double x;

	if (rel_T < 0.)
		return -1.;
	//e r r o r

	x = 2 * rel_T;  //10 is the solid color max

	if (x < 1)
		return .5 * x;

	return sqrt(x) - .5;
}

//______________________________________________________________________________
double Coloroid::GRAY_HUE_DIFF_EQUI(double rel_T1, double rel_T2,
                          int fi_lower1, int fi_upper1, double mu_fi_lower1,
                          int fi_lower2, int fi_upper2, double mu_fi_lower2)
{
	double h_gray, h11, h12, s;

	s = .146095;

	h11 = HUE_PAIRS[fi_lower1][fi_lower2] * mu_fi_lower2 +
	      HUE_PAIRS[fi_lower1][fi_upper2] * (1.0 - mu_fi_lower2); 
	h12 = HUE_PAIRS[fi_upper1][fi_lower2] * mu_fi_lower2 +
	      HUE_PAIRS[fi_upper1][fi_upper2] * (1.0 - mu_fi_lower2); 

	h_gray = h11 * mu_fi_lower1 + h12 * (1.0 - mu_fi_lower1);

	//printf("\n\th11=%g h12=%g h-gray=%g\n", h11, h12, h_gray);

	h_gray *= sqrt(quasi_sqrt_rel_T(rel_T1) * quasi_sqrt_rel_T(rel_T2));

	return s * h_gray;
}
