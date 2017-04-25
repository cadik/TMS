// original authors: Laszlo Neumann,
//                   Martin Cadik

#ifndef TMOCADIK08_COLOROID_H
#define TMOCADIK08_COLOROID_H

#define GRAD_RAD  (double) (M_PI / 180.)
#define RAD_GRAD  (double) (180. / M_PI)
#define EPS_5 2E-5
#define EPS_8 1E-8
#define EPS_9 1E-9
#define EPS_6 1E-6

#define C_x (double) 0.310063
#define C_y (double) 0.316158
#define C_z (double) (1. - (C_x + C_y))

#define C_X (double) 98.0721664
#define C_Y (double) 100.
#define C_Z (double) 118.2253810
//           x       y      
// R       0.6400  0.3300  
// G       0.3000  0.6000  
// B       0.1500  0.0600  
//The CIE chromaticities for the red, green, and blue
//ITU-R BT.709 reference primaries,
//and for CIE Standard Illuminant D65
//default D65 white point
#define D65_x (double) 0.3127269
#define D65_y (double) 0.3290232
#define D65_z (double) (1. - (D65_x + D65_y))

class Coloroid {
	public:
	double luminanceGrad(double X1,  double Y1,  double Z1,
                             double X2,  double Y2,  double Z2,
                             double *dA, double *dT, double *dV);
	void readSpectrumDatae();
	void readColoroidDatae();
	bool readColor2GrayDatae();
	void _7_basic_fi_computation();
	void rgb709Xyz(double R, double G, double B,
	                double*X, double*Y, double*Z);
	
	private:
	struct maxtel {
		double T, V, Y;
	};

	struct ATV {
		double A1, T1, V1,
		       A2, T2, V2;
	};

	void _p_C(long i1, long i2, double dx, double dy, double* q);

	void _p_D65(long i1, long i2, double dx, double dy, double *q);

	bool xyY_XZ(double x, double y, double Y, double *X, double *Z);

	void LIN_MIX(double q, double x1, double y1, double Y1,
		     double x2, double y2, double Y2,
		     double *x, double *y, double *X, double *Y,
		     double *Z, double *p);

	void _xy_XYZ_48_based_on_fi_C(long index, long i1, long i2, double x, double y);

	double _read_data(FILE * file);

	bool Coloroid_A_fi(double A_hue, double* fi);

	void _p(long i1, long i2, double dx, double dy, double* p, double* q, double* lambda);

	bool search(short k1, short k2, double fi, double* lambda);

	void fi_lambda(double fi, double* lambda);

	void _fi_lambda_C(long index, double fi);

	void fi_Coloroid_A(double fi, double* A_hue);

	void fi_Coloroid_limes_color(double fi, double* X, double* Y,
	                             double* Z, double* x, double* y);

	void Coloroid_decomposition(double X, double Y, double X_limes,
	                            double Y_limes, double* white, double* black,
	                            double* color);
	bool XYZ_xyz(double X, double Y, double Z, double* x, double* y, double* z);

	bool XYZ_ATV(double X, double Y, double Z, double *A, double *T,
	             double* V, double* white, double* black, double* color,
	             double* lambda, double* fi);

	void GRAY_EQUI(double X, double Y, double Z,
	               double* _A, double* _T, double* _V,
	               double* rel_T, int* T_lower,
	               int* T_upper, double* mu_T_lower,
	               int* fi_lower, int* fi_upper, double* mu_fi_lower,
	               int* V_lower, int* V_upper, double* mu_V_lower);

	double GRAY_T_EQUI(double rel_T, int T_lower, int T_upper, double mu_T_lower,
	                   int fi_lower, int fi_upper, double mu_fi_lower,
	                   int V_lower, int V_upper, double mu_V_low);

	double GRAY_T_DIFF_EQUI_FOR_2_COLORS(double rel_T1, int T_lower1,  int T_upper1,  double mu_T_lower1,
                                     int fi_lower1, int fi_upper1, double mu_fi_lower1,
                                     int V_lower1,  int V_upper1,  double mu_V_lower1, 
                                     double rel_T2, 
                                     int T_lower2,  int T_upper2,  double mu_T_lower2,
                                     int fi_lower2, int fi_upper2, double mu_fi_lower2,
                                     int V_lower2,  int V_upper2,  double mu_V_lower2);

	double quasi_sqrt_rel_T(double rel_T);

	double GRAY_HUE_DIFF_EQUI(double rel_T1, double rel_T2,
                          int fi_lower1, int fi_upper1, double mu_fi_lower1,
                          int fi_lower2, int fi_upper2, double mu_fi_lower2);

	int A[49],Aa[77];
	double _7_basic_fi[8];

	double	spectrum[831], fi_D65[831], fi_C_[831];
	double	fi_C[49];
	double	Coloroid_fi_D65[49],lam[49],
		xc[49],yc[49],Xc[49],Yc[49],Zc[49],epsc[49];

	double	CIE1931_X[831], CIE1931_Y[831], CIE1931_Z[831],
		CIE1931_x[831], CIE1931_y[831], CIE1931_z[831];
	double	light_D65[831], light_D65_normalized[831],
		light_A[831], light_A_normalized[831];

	ATV pairs[23];
	double pair_rel_gray[23];
	maxtel MAX_TEL[49];

	double HUE_PAIRS[8][8], TELITETTSEG_GRAY[8][6][6], VIL_GRAY[6];
	int sig[8][8];
	char jel, jelsor[800];
};

#endif
