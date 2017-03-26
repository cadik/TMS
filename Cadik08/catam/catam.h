#ifndef CATAM_H
#define CATAM_H

/*
    The C CATAM software library     (CCATSL)
    Faculty of Mathematics, University of Cambridge

    David Harris
    V2.0 August 2000
    V2.1 March  2001

*/

/* Include two commonly used include files */

#include <stdio.h>		/* Input/output functions e.g. printf    */ 
#include <math.h>		/* Simple mathematical functions in libm */

/* C does not have a logical (boolean) type like Pascal so define one */
 
typedef unsigned char boolean;
typedef unsigned char Boolean;	/* alternate spelling */


/* Maximum length of returned strings */
#define STRINGMAXLENGTH 300
typedef char StringCT[STRINGMAXLENGTH+1];

/* Values in float.h preferred */
#define MAXREAL         1.7e+308
#define MINREAL         5.0e-324
#define EPSILON         1.2e-16

/* Useful constants.  */
#define M_E	2.7182818284590452354   /* exp(1.0)      */
#define M_PI	3.14159265358979323846  /* 4.0*atan(1.0) */

/* 
The parameter names should be made more descriptive and match
the documentation. The definitions should be reformatted to make
them more legible. See Rk4CL. 
*/

/* Ordinary Differential Equations */
typedef void (*ODEFunctionCT)(double x, double y[], double ydot[]);
 
void Rk4CL(
       int n,		/*  size of system                                   */
       double *xi,	/*  POINTER to initial x value, gets increased by    */
			/*  the stepsize dx                                  */
       double dx,	/*  stepsize                                         */
       ODEFunctionCT F, /* compute ydot(x,y)                                 */
       double y[],	/*  array of initial conditions, updated to y(x+dx)  */
       double ydot[]);	/*  array of ydot values                             */

/* 
   Note the definition of the fourth parameter, F, in the above. When calling
   Rk4CL you should supply the name of a function which computes vector ydot
   from x and vector y.    See the CCATSL example program introwin.c for
   further details.
*/

boolean RkfCL(int n, double aber, double reler, double *t, double *dt, double dtmin,
	      ODEFunctionCT F, double y[], double yp[], int *nleft);
void    Rkf45CL(int n, double aber, double reler, double tout, double *t, ODEFunctionCT F, 
		double y[], double yp[], int *rkf);

/* Integration */
double QuadCL(double (*F)(double x), double a, double b, double aber, double reler,
	      double *errest, double *flag);
double RombergCL(double (*F)(double x), double a, double b, double *error, int n);

/* Note on two-dimensional arrays - dnh - 12 December 00

The Windows gcc compiler allows incomplete array declarations of
function parameters such as double a[][] (I think that ANSI C does
allow this). MS Visual C++ 5.0 and Borland C++ 5.01 do not: you have
to say double a[][N]. The former is useful for passing general
two-dimensional array parameters as the size of the second dimension
does not have to be specified, and the function can be written to
process arrays of any size.  The called function does of course have
to be given the size of the second dimension some other way, usually
as another parameter, and it has to calculate array indexing
explicitly. All the CCATSL functions which take two-dimensional array
parameters do this. If your compiler does not allow [][] pass the
parameter as a pointer to double using a (double*) cast. See
examples/2d.array.c */

#ifdef CCATSL2D
typedef double two_d_array[][];
typedef int    two_d_int_array[][];
#else
typedef double *two_d_array;
typedef int    *two_d_int_array;
#endif

/* Matrix routines */
double BandCL(int n, int n0, int lbw, int ubw, two_d_array a, double z[], double x[]);
double DecomposeCL(int n, int n0, two_d_array a);
double EigenMaxCL(int n, int n0, double acc, two_d_array a, double eigvec[]);
double EigenValCL(int n, int n0, double acc, two_d_array a, double eigvec[], double theta);
double GaussElimCL(int n, int n0, two_d_array a, double z[], double x[]);
double InvertCL(int n, int n0, two_d_array a);
void   JacobiCL(int n, int n0, two_d_array a, double eigval[]);
void   SolveCL(int n, int n0, two_d_array a, double z[]);
int    SvdCL(int m, int n, int n0, two_d_array a, two_d_array u, two_d_array v, double w[]);
void   TridiagCL(int m, double a[], double b[], double c[], double r[], double x[]);

/* Sorting */
void   XSortCL(double xp[], int nx);
void   XYSortCL(double xp[], double yp[], int nx);
void   XYZSortCL(double xp[], double yp[], two_d_array zp, int nx, int ny, int na0);

/* Special functions */
double BesselCL(int n, double x);	   /* Jn(x) */
void   CubicRootsCL(double a, double b, double c, int *nroots, double *r1, double *r2, double *r3);
void   FftCL(double f[][2]/*V2.1*/, int m, double sign);
void   FftSinCL(int pwr, double fvs[], int v);
double MinCL(double (*F)(double x), double a, double b, double tol);
double PhiCL(double x);	      /* Cumulative normal density function */
double InvPhiCL(double pa);   /* Inverse       "                    */
void   PoissonCL(two_d_array psp, two_d_array ztp, int nnx, int nny, double dlx, double dly);
void   SplineCL(int n, double x[], double y[], double b[], double c[], double d[]);
double SplineValCL(int n, double xx, double x[], double y[], double b[], double c[], double d[]);
double ZeroCL(double (*F)(double x), double a, double b, double tol);

/* Random number generation */
double RandomCL(void);
int    RandomIntCL(int range);
void   SetRandomCL(int seed);

/* functions which are available only in Windows-mode */
#define EscapeCL	EscapeCL_is_not_available_in_character_mode
#define PauseCL		PauseCL_is_not_available_in_character_mode
#define WaitCL		WaitCL_is_not_available_in_character_mode
/* 2.0g HaltCL simply exits in character_mode */
#define HaltCL()	exit(0)

/* program timing */                          
void   TickCL(void);      /* reset timer                     */
long   TockCL(void);      /* elapsed centiseconds            */

/* GnuplotfCL formats a list of expressions (like printf) and transmits
   the resulting string as a command to the gnuplot program.  */
int    GnuplotfCL(char *format, ...);

/* Numerical function error reporting */

extern boolean ErrorFlagCD;
extern char *ErrorMessageCD;

#ifdef CCATSLWIN
#include <catamwin.h>
#endif

#endif /* CATAM_H */
