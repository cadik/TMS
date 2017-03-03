//basic coefficients and transformations of color science

#define WEIGHT_LIGHTNESS_CHROMINANCE 3.0

#define EPS_5 2E-5
#define EPS_8 1E-8
#define EPS_9 1E-9
#define EPS_6 1E-6

//			x       y      
//R        0.6400  0.3300  
//G        0.3000  0.6000  
//B        0.1500  0.0600  
//The CIE chromaticities for the red, green, and blue
//ITU-R BT.709 reference primaries,
//and for CIE Standard Illuminant D65
//default D65 white point
#define D65_x (double) 0.3127269
#define D65_y (double) 0.3290232
#define D65_z (double) (1. - (D65_x + D65_y))

#define x_r 0.64
#define y_r 0.33
#define x_g 0.30
#define y_g 0.60
#define x_b 0.15
#define y_b 0.06

#define epsw (1./D65_y)

#define D65_X (double) (100./D65_y) * D65_x
#define D65_Y (double)  100.
#define D65_Z (double) (100./D65_y) * D65_z
//approximative: Planck 6504 K

#define A_x (double) 0.4475735				//Nemcsics: Szindinamika
#define A_y (double) 0.4074394
#define A_z (double) (1. - (A_x + A_y))

#define A_X (double) 109.8503147
#define A_Y (double) 100.
#define A_Z (double)  35.5849300
//standard A source
//approximative: Planck 2856 K

/*
#define B_X (double) ?
#define B_Y (double) 100.
#define B_Z (double) ?
//standard B source, not used
//approximative Planck 4874 K
*/

#define BYTE unsigned char

#define UINT1 unsigned short
#define UINT  unsigned short
#define DWORD unsigned long
#define LONG  long
#define WORD  unsigned short

#define C_x (double) 0.310063
#define C_y (double) 0.316158
#define C_z (double) (1. - (C_x + C_y))

#define C_X (double)  98.0721664
#define C_Y (double) 100.
#define C_Z (double) 118.2253810
//standard C source
//approximative: Planck 6774 K

#define ANv  0.25940166
#define BNv -0.089901568
//LABHNU konstansok

#define hP (double) 6.626176E-34      /* Planck allando [Js] */
#define kP (double) 1.380662E-23      /* Boltzmann allando [J/K] */
#define cP (double) 2.99792458E+8     /* fenysebesseg [m/s] */
#define c1 (double) 3.741832E-16      //Wm^2
#define c2 (double) 1.438786E-2       //mK
//Planck formulahoz 

/*
Y = 0.212671 * R + 0.715160 * G + 0.072169 * B;
         x       y       z
R        0.6400  0.3300  0.0300
G        0.3000  0.6000  0.1000
B        0.1500  0.0600  0.7900
The CIE chromaticities for the red, green, and blue
ITU-R BT.709 reference primaries, and for CIE Standard Illuminant D65
/**/

#define PI (double) 3.1415926535
#define GRAD_RAD  (double) (PI/180.)
#define RAD_GRAD  (double) (180./PI)
#define _third (double) (1./3.)

#define EPS1 1E-8
#define EPS2 1E-4

#define ANv  0.25940166
#define BNv -0.089901568
//LABHNU konstansok

double	light_D65[831], light_D65_normalized[831],
		light_A[831],   light_A_normalized[831];
double	spectrum[831],fi_D65[831],fi_C_[831];


double	CIE1931_X[831],CIE1931_Y[831],CIE1931_Z[831],
		CIE1931_x[831],CIE1931_y[831],CIE1931_z[831];
//2 fokos latomezo

double	CIE1964_X[831],CIE1964_Y[831],CIE1964_Z[831],
		CIE1964_x[831],CIE1964_y[831],CIE1964_z[831];
//10 fokos latomezo

double	fi_C[49];
double	Coloroid_fi_D65[49],lam[49],
		xc[49],yc[49],Xc[49],Yc[49],Zc[49],epsc[49];

int A[49],Aa[77];

double _7_basic_fi[8];

double vmat[4][4], wmat[4][4], inverz[4][4],
	   XYZ_rgb[4][4], rgb_XYZ[4][4], white[4], b[4], x[4],
	   XL[10], YL[10], ZL[10], Lab_error[10];

double	_Rdelta, _R0, _RC,
		_Gdelta, _G0, _GC,
		_Bdelta, _B0, _BC;



typedef struct
{
double L, a, b;
}
Lab;

typedef struct
{
	double T, V, Y;
}
maxtel;


typedef struct
{
double A1, T1, V1, A2, T2, V2;
}
ATV;

typedef struct
{
double r, g, b; //linearis
}
Lin_rgb;

typedef struct
{
unsigned char r, g, b;	
}
pixel;

typedef struct
{
long tipus;
double X, Y, Z;	
}
metameria;
metameria metamer[6000];

ATV parok[23];
double parok_rel_gray[23];

maxtel MAX_TEL[49];

double HUE_PAIRS[8][8], TELITETTSEG_GRAY[8][6][6], VIL_GRAY[6], LEPCSO_SZAM[8][6];
int sig[8][8];


Lab		Selected_Colors[256],
		Color_Lab_sequence[2600],
		unit_vector, actual_point, k_vector;
long	Selected_Heights[256];
double	Relative_Delta[255];
//ezek magassagi vagy pseudo-coloring skalakhoz kellenek
//Lin_rgb sequence_rgb[2600];
double  Height_delta[2600];
//meterenkent a magassagi rgb linearis szinek es a relativ delta col diff
double AA[256], TT[256], VV[256], delta[255]; 

double Xave[20], Yave[20], Zave[20];

long _2D_uniform_matrix[256][256];

double	eredeti[30][472], ortonormalt[30][472], 
		light_source_normalized[30][472],
		konstans[30], _spectrum[472], vector[472],
		grad[472], grad_prev[472], conjugated_direction[472];

double	R_fi,	G_fi,	B_fi, 
		R_X,	R_Y,	R_Z,
		G_X,	G_Y,	G_Z,
		B_X,	B_Y,	B_Z;

double dL[256], da[256], db[256];
unsigned char Rps[256], Gps[256], Bps[256];

FILE * out_metamer;

