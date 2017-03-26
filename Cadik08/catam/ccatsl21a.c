/*CATAM Software Library*/
/*Copyright  DAMTP  1994, 1998, 1999, 2000, 2001*/
/*Department of Applied Mathematics and Theoretical Physics*/
/*University of Cambridge*/
/* Ancestry:
       Scientific Programmers Toolkit - authors
       CATSL Pascal version 2.41 - John Evans
       CCATSL 1.1 
       CCATSL 2.0
*/

#include <stdlib.h>
#include <time.h> 

//#define BUILDING /* fix ErrorFlagCD declaration */
#include "catam.h"
#define NUM_C
//#include "numdef.h"
//#include "numdata.h"

/* these external functions are for use only by other library modules */ 
#define _EXTERN_ extern
//#include "num.int.funs.h"
#include "ccatsl21a.h"

int _Imin(int a, int b)
{
  if (a > b)
    return b;
  else
    return a;
}

int _Imax(int a, int b)
{
  if (a >= b)
    return a;
  else
    return b;
}

void _Set_error(int errcode, char *errstr)
{
  Errorflag = true;
  Errorcode = errcode;
  strcpy(Errorstring, errstr);
}


double *_New_array(long nn)
{
  char *p;
  long m;

  m = nn * sizeof(double);
  p = _NewPtr(m);
  if (p == NULL)
    return NULL;
  else
    return ((double *)p);
}


void _Free_array(double **p, long nn)
{
  long M;

  M = nn * sizeof(double);
  _DisposePtr((char *)(*p), M);
  *p = NULL;
}




double _Float(int i)
{
  return i;
}  /* float */


double _Nraise(double value, int np)
{
  int n;
  double val;
  int FORLIM;

  if (np == 0)
    return 1.0;

  n = abs(np);
  val = value;
  FORLIM = n;
  for (n = 2; n <= FORLIM; n++)
    val *= value;
  if (np < 0)
    val = 1.0 / val;
  return val;
}

double _Raise(double value, double power)
{
  double Result;
  boolean npi;
  double Rs, f, fp, d, t;
  long nfp;

  if (value == 0.0)
    return 0.0;

  if (power == 0.0)
    return 1.0;

  Rs = exp(power * log(fabs(value)));
  Result = Rs;

  if (value >= 0.0)
    return Result;
  Result = 0.0;
  f = 1.0;
  npi = false;

  /*loop based on finding if sin(2n+1)r¹) is zero - imaginary part zero*/
  /*students seem to be finding problems with original Raise*/

  do {
    fp = f * power;
    nfp = (long)floor(fp + 0.5);
    d = fabs(nfp - fp);
    if (d < 0.001) {
      t = 1.0;
      if (nfp != nfp / 2 * 2)
	t = -1.0;
      npi = true;
    }
    f += 2.0;
  } while (npi != true && f <= 30.0);   /*limited search - 15 roots*/
  if (npi)
    return (Rs * t);
  return Result;
}


int _Sgn(double x)
{
  if (x > 0.0)
    return 1;
  else if (x == 0.0)
    return 0;
  else
    return -1;
}

double _Tan(double x)
{
  double c, s;

  c = cos(x);
  s = sin(x);
  if (c == 0.0)
    return (MAXREAL * _Sgn(s));
  else
    return (s / c);
}


/* return arcsin(x) - result in range [-pi/2,pi/2] */

double _Arcsin(double x)
{
  if (fabs(x) > 1.0 - EPSILON)
    return (_Sgn(x) * 0.5 * M_PI);
  else if (x > 0.5)
    return (0.5 * M_PI - atan(sqrt(1 - x * x) / x));
  else if (x > -0.5)
    return atan(x / sqrt(1 - x * x));
  else
    return (-0.5 * M_PI - atan(sqrt(1 - x * x) / x));
}


double _Min(double a, double b)
{
  if (a > b)
    return b;
  else
    return a;
}


double _Max(double a, double b)
{
  if (a > b)
    return a;
  else
    return b;
}


double _log10_(double x)
{
  /*--- D11aug99 - dnh*/
  return (log(x) / log(10.0));
}



void _CTerror(char *s)
{
  puts(s);
  Errorflag = false;
  Errorcode = 0;
  *Errorstring = '\0';
}

#define MAXITER         100
#define MAXITERSTRING   "100"


/* last interval used for spline evaluation */
static int search_interval;


/*dnh - there is a local seterr in ctcntrl too so we'll rename this one
  setmatherror
*/
static void
_setmatherror(int code)
{
  switch (code) {

  case 1:
    strcpy(Errorstring, "Error in function evaluation");
    break;

  case 2:
    strcpy(Errorstring, "f(a) and f(b) not of opposite sign");
    break;

  case 3:
    strcpy(Errorstring, "Bad interval: high end <= low end");
    break;

  case 4:
    strcpy(Errorstring, "Spline needs >=3 points");
    break;

  case 5:
    strcpy(Errorstring, "x[i+1] not > x[i] for some i");
    break;

  case 6:
    strcpy(Errorstring, "Must halve interval at least once");
    break;

  case 7:
    strcpy(Errorstring, "Not enough heap space");
    break;

  case 8:
    strcpy(Errorstring, "Tolerance must not be negative");
    break;

  case 9:
    strcpy(Errorstring, "Order must be >= 1");
    break;

  case 10:
    strcpy(Errorstring, "Order must be <= 30");
    break;

  case 11:
    strcpy(Errorstring, "Time Step must be positive");
    break;

  case 12:
    strcpy(Errorstring, "Minimum time step must be positive");
    break;

  case 13:
    strcpy(Errorstring, "Power of 2 must be >= 1");
    break;

  case 14:
    sprintf(Errorstring, "No convergence after %s iterations.", MAXITERSTRING);
    break;
  }/*of case*/
  _Set_error(code + 20900, Errorstring);
  /*_CTerror(Errorstring);  ---------------MAC********/
}

static void
_swop(double *x, double *y)
{
  /*If x>y then swop x and y*/
  double temp;

  if (*x <= *y)
    return;
  temp = *x;
  *x = *y;
  *y = temp;
}


/*--- C there is a weird gcc error*/
/*--- C  ctmaths.c:109: warning: static declaration for cbrt follows non-static*/
/*--- C fixed by simply changing cbrt to cubrt*/
static double
_cubrt(double a)
{
  if (a == 0.0)
    return 0.0;
  else
    return (_Sgn(a) * _Raise(fabs(a), 1.0 / 3.0));
}


/* Local variables for cubic: */
struct LOC_cubic {
  double *r1, *r2, *r3, e, f, g;
} ;

static void
_notzero(struct LOC_cubic *LINK)
{
  double t;

  LINK->g = LINK->f * _Raise(-LINK->e, -3.0 / 2.0) / -2.0;
  if (fabs(LINK->g) >= 1.0)
    LINK->g = _Sgn(LINK->g) * 0.999999999;
  LINK->g = atan(sqrt(1.0 - LINK->g * LINK->g) / LINK->g);
  t = M_PI;
  if (LINK->g < 0.0)
    LINK->g += t;
  *LINK->r1 = cos(LINK->g / 3.0);
  *LINK->r2 = cos((LINK->g + 2.0 * t) / 3.0);
  *LINK->r3 = cos((LINK->g + 4.0 * t) / 3.0);
}


void CubicRootsCL(double a, double b, double c, int *nroots, double *r1_, double *r2_, double *r3_)
{

  /*This procedure finds the roots of x**3 + a*x**2 + b*x +c = 0.
It returns the number of real roots, nroots, and the roots, r1,r2 and r3,
in order of increasing value.
If there is only one root, it is r1, and the other two roots are set to 0.

     The standard form is y**3 + 3*e*y + f, where y = x +a/3*/
  struct LOC_cubic V;
  double temp, d;
  boolean warn;

  V.r1 = r1_;
  V.r2 = r2_;
  V.r3 = r3_;
  *V.r1 = 0.0;
  *V.r2 = 0.0;
  *V.r3 = 0.0;
  temp = a / 3.0;
  V.e = b / 3.0 - temp * temp;
  V.f = 2.0 * temp * temp * temp - b * a / 3.0 + c;

  /*Evaluate the discriminant*/
  d = V.f * V.f + V.e * V.e * V.e * 4.0;

  warn = false;

  /*If d close to zero, make it zero, but issue a warning*/
  if (d != 0.0) {
    if (V.e <= 0.0 && fabs(d) <= EPSILON) {
      warn = true;
      d = 0.0;
    }
  }

  /*if d>0, there is only one real root (Cardan's solution)*/
  if (d > 0) {
    if (V.e == 0.0)
      *V.r1 = _cubrt(-V.f) - a / 3.0;
    else {
      *V.r1 = _cubrt((sqrt(d) - V.f) / 2.0);
      *V.r1 += -(V.e / *V.r1) - a / 3.0;
    }
    *nroots = 1;
    return;
  }
  /*If d=0, two roots are equal*/
  if (d == 0) {
    *V.r1 = -_Sgn(V.f);
    *V.r2 = *V.r1 / -2.0;
    *V.r3 = *V.r2;
  } else {
    /*If d<0, use the general solution*/
    if (V.f == 0.0) {
      *V.r1 = sqrt(3.0) / 2.0;
      *V.r2 = 0.0;
      *V.r3 = -*V.r1;
    } else
      _notzero(&V);
  }

  /*This part of the program applies to d<=0*/
  d = 2.0 * sqrt(-V.e);
  V.g = a / 3.0;
  *V.r1 = d * *V.r1 - V.g;
  *V.r2 = d * *V.r2 - V.g;
  *V.r3 = d * *V.r3 - V.g;

  /*Order the roots*/
  do {
    _swop(V.r1, V.r2);
    _swop(V.r2, V.r3);
  } while (*V.r1 > *V.r2 || *V.r2 > *V.r3);
  *nroots = 3;
  if (!warn)
    return;
  printf("\nWARNING: It is possible that there is only one root.\n");
  printf("Due to the limitation on the accuracy of the machine,\n");
  printf("the two equal roots may be spurious.\n\n");
}

static double
_Zeroin(double a, double b, double tol, double (*F)(double x))
{

  /* Find a zero of f(x) in interval (a,b) with tolerance tol >=0 */
  double Result, c, d, e, eps, fa, fb, fc, xm, tol1, p, q, r, s;
  boolean converged;
  int iter;

  Errorflag = false;
  /* add this to prevent C compilation warnings*/
  Result = 0.0;
  if (Errorflag)
    return Result;

  /* relative machine accuracy */
  eps = EPSILON;
  tol = fabs(tol / 2);

  /* initialise ... */

  fa = F(a);
  if (Errorflag) {   /*'Error in function evaluation'*/
    _setmatherror(1);
    return Result;
  }

  fb = F(b);
  if (Errorflag) {   /*'Error in function evaluation'*/
    _setmatherror(1);
    return Result;
  }

  if (fa * fb > 0) {  /*'f(a) and f(b) not of opposite sign'*/
    _setmatherror(2);
    return Result;
  }

  /* begin step ... */

  c = a;
  fc = fa;
  d = b - a;
  e = d;

  iter = 0;
  do {

    if (fabs(fc) < fabs(fb)) {   /* ensure fa,fc>=fb */
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    /* convergence test ... */

    tol1 = eps * fabs(b) + tol;
    xm = (c - b) / 2;

    converged = (fabs(xm) <= tol1 || fb == 0);

    if (!converged) {
      if (fabs(e) < tol1 || fabs(fa) < fabs(fb)) {
	/* bisection */
	d = xm;
	e = d;
      } else {
	if (a == c) {   /* linear interpolation */
	  s = fb / fa;
	  p = 2 * xm * s;
	  q = 1 - s;
	} else {  /* inverse quadratic interpoln */
	  q = fa / fb;
	  r = fb / fc;
	  s = fb / fa;
	  p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
	  q = (q - 1) * (r - 1) * (s - 1);
	}

	/* adjust signs ... */

	if (p > 0)
	  q = -q;
	p = fabs(p);

	/* test acceptability ... */

	if (2 * p > 3 * xm * q - fabs(tol1 * q) || 2 * p > fabs(e * q)) {
	  /* not acceptable, .. bisect */
	  d = xm;
	  e = d;
	} else {  /* acceptable */
	  e = d;
	  d = p / q;
	}
      }

      /* complete the step ... */

      a = b;
      fa = fb;
      if (fabs(d) > tol1)
	b += d;
      else
	b += _Sgn(xm) * tol1;
      fb = F(b);
      if (Errorflag) {   /*'Error in function evaluation'*/
	_setmatherror(1);
	return Result;
      }

      /* move left or right? */
      if (fb * fc > 0) {
	c = a;
	fc = fa;
	d = b - a;
	e = d;
      }
    }

    iter++;
  } while (!(converged || iter > MAXITER));

  Result = b;

  /* check any remaining error */
  if (iter > MAXITER)
    _setmatherror(14);

  return Result;
} /* zeroin */

static double
_Fmin(double a, double b, double tol, double (*fp)(double x))
{
  double Result, c, d, e, eps, xm, p, q, r, tol1, tol2, u, v, w, fu, fv, fw, fx, x;
  boolean done;
  int iter;


  /* add this to prevent C compilation warnings*/
  Result = 0.0;
  if (Errorflag)
    return Result;

  /* golden proportion */
  c = (3 - sqrt(5.0)) / 2;

  /* m/c precision */
  eps = sqrt(EPSILON);

  /* initialise ... */

  e = 0.0;
  if (b <= a) {  /*'Value for high end of interval less than low end'*/
    _setmatherror(3);
    return Result;
  }
  v = a + c * (b - a);
  w = v;
  x = v;
  fx = fp(x);
  if (Errorflag) {   /*'Error in function evaluation'*/
    _setmatherror(1);
    return Result;
  }
  fv = fx;
  fw = fx;

  done = false;

  /* main loop ... */

  iter = 0;
  do {

    xm = (a + b) / 2;

    /* check if search done ... */

    tol1 = eps * fabs(x) + tol / 3;
    tol2 = 2 * tol1;
    done = (fabs(x - xm) <= tol2 - (b - a) / 2);

    if (!done) {
      /* choose and use method ... */

      if (fabs(e) > tol1) {   /* try fit parabola */
	r = (x - w) * (fx - fv);
	q = (x - v) * (fx - fw);
	p = (x - v) * q - (x - w) * r;
	q = 2 * (q - r);
	if (q > 0)
	  p = -p;
	q = fabs(q);
	r = e;
	e = d;

	/* check if parabola acceptable ... */

	if (2 * fabs(p) < fabs(q * r) && p > q * (a - x) && p > q * (b - x)) {
	  d = p / q;
	  u = x + d;
	  if (u - a < tol2 || b - u < tol2)
	    d = tol2 * _Sgn(xm - x);
	} else {  /* ... if not, then use golden section */
	  if (x < xm)
	    e = b - x;
	  else
	    e = a - x;
	  d = c * e;
	}
      } else {  /* otherwise use golden section anyway */
	if (x < xm)
	  e = b - x;
	else
	  e = a - x;
	d = c * e;
      }


      /* keep away from x ... */

      if (fabs(d) < tol1)
	u = x + tol1 * _Sgn(d);
      else
	u = x + d;
      fu = fp(u);
      if (Errorflag) {   /*'Error in function evaluation'*/
	_setmatherror(1);
	return Result;
      }

      /* update ... */

      if (fu <= fx) {
	if (u < x)
	  b = x;
	else
	  a = x;
	v = w;
	fv = fw;
	w = x;
	fw = fx;
	x = u;
	fx = fu;
      } else {
	if (u < x)
	  a = u;
	else
	  b = u;
	if (fu <= fw || w == x) {
	  v = w;
	  fv = fw;
	  w = u;
	  fw = fu;
	} else if (fu <= fv || v == x || v == w) {
	  v = u;
	  fv = fu;
	}
      }
    }


    iter++;
  } while (!(done || iter > MAXITER));

  Result = x;
  if (iter > MAXITER)
    _setmatherror(14);

  return Result;


}  /*fmin*/


void SplineCL(int n, double *x, double *y, double *b, double *c, double *d)
{
  int i;
  double t;



  if (Errorflag)
    return;

  if (n < 3) {  /*'Spline needs >=3 points'*/
    _setmatherror(4);
    return;
  }

  /* tri-diagonal coefficients : b[] diagonal, d[] off-d, c[] RHS ... */

  d[0] = x[1] - x[0];
  if (d[0] <= 0.0) {   /*'x[i+1] not > x[i] for some i'*/
    _setmatherror(5);
    return;
  }
  c[1] = (y[1] - y[0]) / d[0];
  for (i = 2; i < n; i++) {
    d[i - 1] = x[i] - x[i - 1];
    if (d[i - 1] <= 0.0) {   /*'x[i+1] not > x[i] for some i'*/
      _setmatherror(5);
      return;
    }
    b[i - 1] = 2 * (d[i - 2] + d[i - 1]);
    c[i] = (y[i] - y[i - 1]) / d[i - 1];
    c[i - 1] = c[i] - c[i - 1];
  }

  /* end condition ... */

  b[0] = -d[0];
  b[n - 1] = -d[n - 2];

  if (n == 3) {
    c[0] = 0.0;
    c[n - 1] = 0.0;
  } else {
    c[0] = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
    c[n - 1] = c[n - 2] / (x[n - 1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4]);
    c[0] = c[0] * d[0] * d[0] / (x[3] - x[0]);
    c[n - 1] = c[n - 1] * d[n - 2] * d[n - 2] / (x[n - 4] - x[n - 1]);
  }

  /* forward elimination ... */

  for (i = 1; i < n; i++) {
    t = d[i - 1] / b[i - 1];
    b[i] -= t * d[i - 1];
    c[i] -= t * c[i - 1];
  }

  /* back substitution ... */

  c[n - 1] /= b[n - 1];
  for (i = n - 2; i >= 0; i--)
    c[i] = (c[i] - d[i] * c[i + 1]) / b[i];

  /* ... now find polynomial coefficients ... */

  b[n - 1] = (y[n - 1] - y[n - 2]) / d[n - 2] +
	     d[n - 2] * (c[n - 2] + 2 * c[n - 1]);
  for (i = 0; i <= n - 2; i++) {
    b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2 * c[i]);
    d[i] = (c[i + 1] - c[i]) / d[i];
    c[i] = 3 * c[i];
  }
  c[n - 1] = 3 * c[n - 1];
  d[n - 1] = d[n - 2];

  search_interval = 1;

  /* no errors */
  Errorflag = false;

}  /*of spline*/


double SplineValCL(int n, double xx, double *x, double *y, double *b, double *c, double *d)
{
  double Result;
  int i, j, k;
  double t;


  /* add this to prevent C compilation warnings*/
  Result = 0.0;
  if (Errorflag)
    return Result;

  i = search_interval;
  if (i < 1 || i > n)
    i = 1;

  /* find interval ... */
  if (xx <= x[i - 1] || xx > x[i]) {
    i = 1;
    j = n;
    do {
      k = (i + j) / 2;
      if (xx < x[k - 1])
	j = k;
      else
	i = k;
    } while (j != i + 1);
  }

  t = xx - x[i - 1];
  search_interval = i;
  Result = y[i - 1] + t * (b[i - 1] + t * (c[i - 1] + t * d[i - 1]));

  /* no errors */
  Errorflag = false;

  return Result;
}  /*of seval*/


/* Local variables for Romb: */
struct LOC_Romb {
  int n;
  double *q;
} ;

static void
_errstop(struct LOC_Romb *LINK)
{
  _Free_array(&LINK->q, LINK->n + 1);
  /*'Error in function evaluation'*/
  _setmatherror(1);
}


static double
_Romb(double a, double b, int n_, double (*fp)(double x), double *error)
{
  struct LOC_Romb V;
  double Result;
  int i, j;
  double c, d, ff, h, h2, sum, x;
  int FORLIM, FORLIM1;


  V.n = n_;
  /* add this to prevent C compilation warnings*/
  Result = 0.0;
  /* Must pre-set in case user routines don't */
  Errorflag = false;

  if (V.n <= 0) {  /*'Must halve interval at least once'*/
    _setmatherror(6);
    return Result;
  }

  *error = 0.0;
  h = b - a;
  /* h current interval size */
  if (h <= 0) {  /*'Bad interval'*/
    _setmatherror(3);
    return Result;
  }

  /* generate local heap array q */
  V.q = _New_array(V.n + 1);
  if (V.q == NULL) {  /*'Not enough heap space'*/
    _setmatherror(7);
    return Result;
  }

  /* Trapezium rule sum ... */

  sum = fp(a);
  if (Errorflag) {
    _errstop(&V);
    return Result;
  }
  sum += fp(b);
  if (Errorflag) {
    _errstop(&V);
    return Result;
  }

  /* zeroth estimate of integral .. */

  V.q[0] = h * sum / 2;

  FORLIM = V.n;
  /* interval halving loop ... */

  for (i = 1; i <= FORLIM; i++) {
    h2 = h;
    h /= 2;

    /* trapezium sum ...*/

    x = a + h;
    do {
      ff = fp(x);
      if (Errorflag) {
	_errstop(&V);
	return Result;
      }

      /* use error to track function maximum ... */
      if (fabs(ff) > *error)
	*error = fabs(ff);

      sum += 2 * ff;
      x += h2;
    } while (x <= b);

    /* store next estimate ... */
    V.q[i] = h * sum / 2;
  }

  FORLIM = V.n;
  /* Romberg error elimination ... */

  for (i = 1; i <= FORLIM; i++) {
    c = 1.0;   /* d=2^(2*i)-1 */
    for (j = 1; j <= i; j++)
      c += c;
    c *= c;
    d = c - 1;

    x = V.q[0];
    FORLIM1 = V.n;
    for (j = 1; j <= FORLIM1; j++)
      V.q[j - 1] = (c * V.q[j] - V.q[j - 1]) / d;
  }

  *error = 2e-10 * *error + fabs(V.q[0] - x);
  Result = V.q[0];
  _Free_array(&V.q, V.n + 1);

  return Result;
}  /*of Romb*/


/* Local variables for Quanc8: */
struct LOC_Quanc8 {
  double a, b, *flag;
  double qf[17];
  double qw[5];
  double st, xfail;
} ;

static void
_errstop_(struct LOC_Quanc8 *LINK)
{
  *LINK->flag = (long)(*LINK->flag) + (LINK->b - LINK->xfail) / (LINK->b - LINK->a);
  /*'Error in function evaluation'*/
  _setmatherror(1);
}

static double _qrule(int i, struct LOC_Quanc8 *LINK)
{
  int j;
  double t;

  t = 0.0;
  for (j = 0; j <= 4; j++)
    t += LINK->qw[j] * (LINK->qf[i + j] + LINK->qf[i - j + 8]);
  return (t * LINK->st);
}

static double
_Quanc8(double a_, double b_, double aber, double reler, double (*fp)(double x), double *errest, double *flag_, int *nf)
{
  struct LOC_Quanc8 V;
  double Result;
  double qr[31];
  double qx[17];
  double qfs[8][30], qxs[8][30];

  double ar, st1, cor, qp, q, qd, ql, est, t, tl, r;
  int lm, lx, lo, nx, fi, l, ni, i, j, k, m;
  boolean finished;

  /*qrule*/
  V.a = a_;
  V.b = b_;
  V.flag = flag_;
  /* add this to prevent C compilation warnings*/
  Result = 0.0;

  if (Errorflag)
    return Result;

  if (V.b <= V.a) {  /*'Bad interval'*/
    _setmatherror(3);
    return Result;
  }

  /* Test tolerances */
  if (aber < 0.0 || reler < 0.0) {
    /*'Tolerance must not be negative'*/
    _setmatherror(8);
    return Result;
  }

  ar = 0.0;
  r = 0.0;
  cor = 0.0;

  lm = 1;
  lx = 30;
  lo = 6;
  nx = 5000;

  m = 1;
  for (k = 0; k <= lo; k++)
    m += m;
  fi = nx - (lx - lo + m) * 8;

  t = 14175.0;
  V.qw[0] = 3956 / t;
  V.qw[1] = 23552 / t;
  V.qw[2] = -3712 / t;
  V.qw[3] = 41984.0 / t;
  V.qw[4] = -9080 / t;

  *V.flag = 0.0;
  *errest = 0.0;
  *nf = 0;
  Result = 0.0;
  if (V.a == V.b)
    return Result;

  /* initialise for 1st interval ...*/

  l = 0;
  ni = 1;
  qx[0] = V.a;
  qx[16] = V.b;
  qp = 0.0;
  st1 = (V.b - V.a) / 16;
  i = 8;
  do {
    j = i;
    do {
      qx[j] = (qx[j - i] + qx[j + i]) / 2;
      j += i + i;
    } while (j <= 16);
    i /= 2;
  } while (i != 1);
  *nf = 9;

  j = 0;
  finished = false;
  do {
    V.qf[j] = fp(qx[j]);
    j += 2;
  } while (!(j > 16 || Errorflag));
  if (Errorflag) {
    V.xfail = qx[j - 2];
    _errstop_(&V);
    return Result;
  }
  *nf = j / 2;

  /* central calculation ... */

  do {

    j = 0;
    do {
      qx[j + 1] = (qx[j] + qx[j + 2]) / 2;
      V.qf[j + 1] = fp(qx[j + 1]);
      j += 2;
      finished = (finished || Errorflag);
    } while (!(j > 14 || finished));
    if (Errorflag) {
      V.xfail = qx[j - 2];
      _errstop_(&V);
      return Result;
    }

    *nf += j / 2;
    V.st = (qx[16] - qx[0]) / 16;
    ql = _qrule(0, &V);
    qr[l] = _qrule(8, &V);
    q = ql + qr[l];
    qd = q - qp;
    ar += qd;

    /* convergence test ... */
    est = fabs(qd) / 1023;
    tl = reler * fabs(ar);
    if (aber > tl)
      tl = aber;
    tl = tl * V.st / st1;

    if (l < lm || l < lx && *nf <= fi && est >= tl) {
      /* next interval - save data */
      ni += ni;
      for (i = 0; i <= 7; i++) {
	qfs[i][l] = V.qf[i + 9];
	qxs[i][l] = qx[i + 9];
      }
      l++;

      /* set up LH data for next step ...*/
      qp = ql;
      for (j = -1; j >= -8; j--) {
	V.qf[j + j + 18] = V.qf[j + 9];
	qx[j + j + 18] = qx[j + 9];
      }
      finished = false;
    } else {
      if (l >= lx)
	(*V.flag)++;
      else if (*nf > fi) {
	fi += fi;
	lx = lo;
	*V.flag += (V.b - qx[0]) / (V.b - V.a);
      }

      /* add contribns into sums ...*/
      r += q;
      *errest += est;
      cor += qd / 1023;

      /* next interval ... */
      while ((ni & 1) == 1) {
	ni /= 2;
	l--;
      }
      ni++;
      if (l > 0) {
	/* set up data for next interval ...*/
	qp = qr[l - 1];
	qx[0] = qx[16];
	V.qf[0] = V.qf[16];
	for (i = 1; i <= 8; i++) {
	  V.qf[i + i] = qfs[i - 1][l - 1];
	  qx[i + i] = qxs[i - 1][l - 1];
	}
	finished = false;
      } else
	finished = true;

    }

  } while (!finished);

  r += cor;
  if (*errest != 0) {
    while (fabs(r) + *errest == *errest)
      *errest = 2 * *errest;
  }
  return r;


  /* no. of fn evals near limit */
}  /*of Quanc8*/


#define MAXN            30


void Rk4CL(int n, double *t, double dt, void (*F)(double t, double *y, double *ydash), double *y, double *yp)
{

  /* Fourth-order Runge-Kutta procedure to effect one step*/
  /*     in the solution of the nth order system*/
  /*      y'j = fj(t, y1, ..... , yn)      j=1,...n, n<=30    */
  /* <<<<<<<<< maximum 30 variables */
  int i, j, u, v;
  double k[MAXN], y0[MAXN];


  if (Errorflag)
    return;

  if (n < 1) {  /*'Order must be >= 1'*/
    _setmatherror(9);
    return;
  }

  if (n > MAXN)  /*'Order must be <= 30'*/
    _setmatherror(10);

  /* copy values */
  for (i = 0; i < n; i++)
    y0[i] = y[i];
  for (i = 0; i < n; i++)
    k[i] = 0.0;

  u = 1;
  for (i = 0; i <= 3; i++) {
    v = 2 - i / 2;
    *t += dt / 2 * (i & 1);
    /*------------------------------------------*/
    /*user-func uses globals t,y[] to set yp[]*/

    /*dnh dummy := F;*/
    F(*t, y, yp);

    /*-------------------------------------------*/
    if (Errorflag) {   /*'Error in function evaluation'*/
      _setmatherror(1);
      return;
    }
    for (j = 0; j < n; j++) {
      k[j] += yp[j] * u;
      y[j] = y0[j] + yp[j] * dt / v;
    }
    u = v;
  }

  for (i = 0; i < n; i++)
    y[i] = y0[i] + k[i] * dt / 6;

}  /*of RK4*/

#undef MAXN


#define MAXN            30


/* <<<<<<<<<<< maximum 30 variables */

typedef double odearray[MAXN];


/* Local variables for Rkf: */
struct LOC_Rkf {
  int n;
  double *t;
  void (*F)(double t, double *y, double *ydash);
  double *y, *yp;
  double t0;
} ;

static void _getyp(double *fstep, double dt, struct LOC_Rkf *LINK)
{
  int j, FORLIM;

  *LINK->t = LINK->t0 + dt;
  /*dnh dummy := F;*/
  LINK->F(*LINK->t, LINK->y, LINK->yp);
  if (Errorflag) {   /*'Error in function evaluation'*/
    _setmatherror(1);
    return;
  }
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    fstep[j] = LINK->yp[j];
}


static void
_Rkf(int n_, double aber, double reler, double *t_, double *dt, double dtmin, void (*F_)(double t, double *y, double *ydash), double *y_, double *yp_, int *nleft, boolean *flag)
{

  /* Variable step Runge-Kutta-Fehlberg method to complete one step*/
  /*      in the solution of the nth order system*/
  /*     y'j = fj(t, y1, ..... , yn)      j=1,...n,  n<=30   */
  struct LOC_Rkf V;
  odearray y0, yy, f1, f2, f3, f4, f5;
  int j;
  double d, e, s, te;
  int FORLIM;


  V.n = n_;
  V.t = t_;
  V.F = F_;
  V.y = y_;
  V.yp = yp_;
  if (Errorflag)
    return;
  *flag = false;

  /* test order*/
  if (V.n < 1) {  /*'Order must be >= 1'*/
    _setmatherror(9);
    return;
  }

  /* test max order */
  if (V.n > MAXN)  /*'Order must be <= 30'*/
    _setmatherror(10);

  /* test dtmin */
  if (dtmin <= 0.0) {  /*'Minimum time step must be positive'*/
    _setmatherror(12);
    return;

  }
  /* Test tolerances */
  if (aber < 0.0 || reler < 0.0) {
    /*'Tolerance must not be negative'*/
    _setmatherror(8);
    return;
  }

  /* save t,y at start of step */
  V.t0 = *V.t;
  FORLIM = V.n;
  for (j = 0; j < FORLIM; j++)
    y0[j] = V.y[j];

  /* store derivatives at start of step..... */

  _getyp(yy, 0.0, &V);
  if (Errorflag)
    return;

  /* loop until error criteria is satisfied...... */

  do {
    d = *dt / 4;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++)
      V.y[j] = y0[j] + d * yy[j];
    _getyp(f1, d, &V);
    if (Errorflag)
      return;

    d = 3 * *dt / 32;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++)
      V.y[j] = y0[j] + d * (yy[j] + 3 * f1[j]);
    _getyp(f2, 3 * *dt / 8, &V);
    if (Errorflag)
      return;

    d = *dt / 2197;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++)
      V.y[j] = y0[j] + d * (1932 * yy[j] + 7296 * f2[j] - 7200 * f1[j]);
    _getyp(f3, 12 * *dt / 13, &V);
    if (Errorflag)
      return;

    d = *dt / 4104;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++)
      V.y[j] = y0[j] +
	  d * (8341 * yy[j] - 845 * f3[j] + 29440 * f2[j] - 32832.0 * f1[j]);
    _getyp(f4, *dt, &V);
    if (Errorflag)
      return;


    d = *dt / 20520;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++)
      V.y[j] = y0[j] + d * (9295 * f3[j] - 6080 * yy[j] - 5643 * f4[j] +
			    41040.0 * f1[j] - 28352 * f2[j]);
    _getyp(f5, *dt / 2, &V);
    if (Errorflag)
      return;


    /* compute largest error relative to tolerance..... */

    s = 100.0;
    FORLIM = V.n;
    for (j = 0; j < FORLIM; j++) {
      te = fabs(yy[j] / 360 - 128 * f2[j] / 4275 - 2197 * f3[j] / 75240.0 +
		f4[j] / 50 + 2 * f5[j] / 55);
      if (te > 0)
	e = (aber + reler * fabs(y0[j])) / te;
      if (e < s)
	s = e;
    }

    /* s>1 means accept this step, else reduce dt ...*/

    if (s > 1) {
      *V.t = V.t0 + *dt;
      (*nleft)--;
      d = *dt / 7618050.0;
      FORLIM = V.n;
      for (j = 0; j < FORLIM; j++)
	V.y[j] = y0[j] + d * (902880.0 * yy[j] + 3855735.0 * f3[j] - 1371249.0 *
				f4[j] + 3953664.0 * f2[j] + 277020.0 * f5[j]);
    } else {
      *dt /= 2;
      *nleft *= 2;
      if (fabs(*dt) < dtmin) {
	*flag = false;
	*V.t = V.t0;
	FORLIM = V.n;
	for (j = 0; j < FORLIM; j++)
	  V.y[j] = y0[j];
	return;
      }
    }

    /* s>2^5 means dt could be doubled if nleft even */

    if (s > 32 && ((*nleft) & 1) == 0) {
      *dt = 2 * *dt;
      *nleft /= 2;
    }

  } while (s <= 1);

  *flag = true;

}  /*of Rkf*/

#undef MAXN


#define MAXN            30


/* <<<<<<<<<< maximum 30 variables */

typedef double odearray_[MAXN];

/* local parameters for rkf45 */


/* Local variables for Rkf45: */
struct LOC_Rkf45 {
  int n;
  double *t;
  void (*F)(double t, double *y, double *ydash);
  double *y, *yp;

  odearray_ y0, yy, f1, f2, f3, f4, f5;
} ;

static void _getyp_(double *fstep, double dt, struct LOC_Rkf45 *LINK)
{
  int j;
  double t0;
  int FORLIM;

  t0 = *LINK->t;
  *LINK->t = t0 + dt;
  /*dnh dummy := F;*/
  LINK->F(*LINK->t, LINK->y, LINK->yp);
  *LINK->t = t0;
  if (Errorflag) {   /*'Error in function evaluation'*/
    _setmatherror(1);
    return;
  }
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    fstep[j] = LINK->yp[j];
}

static void _rkfehl(double dt, struct LOC_Rkf45 *LINK)
{
  double d;
  int j, FORLIM;

  d = dt / 4;
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    LINK->y[j] = LINK->y0[j] + d * LINK->yy[j];
  _getyp_(LINK->f1, d, LINK);
  if (Errorflag)
    return;

  d = 3 * dt / 32;
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    LINK->y[j] = LINK->y0[j] + d * (LINK->yy[j] + 3 * LINK->f1[j]);
  _getyp_(LINK->f2, 3 * dt / 8, LINK);
  if (Errorflag)
    return;

  d = dt / 2197;
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    LINK->y[j] = LINK->y0[j] +
	d * (1932 * LINK->yy[j] + 7296 * LINK->f2[j] - 7200 * LINK->f1[j]);
  _getyp_(LINK->f3, 12 * dt / 13, LINK);
  if (Errorflag)
    return;

  d = dt / 4104;
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    LINK->y[j] = LINK->y0[j] + d * (8341 * LINK->yy[j] - 845 * LINK->f3[j] +
				 29440 * LINK->f2[j] - 32832.0 * LINK->f1[j]);
  _getyp_(LINK->f4, dt, LINK);
  if (Errorflag)
    return;

  d = dt / 20520;
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    LINK->y[j] = LINK->y0[j] +
	d * (9295 * LINK->f3[j] - 6080 * LINK->yy[j] - 5643 * LINK->f4[j] +
	     41040.0 * LINK->f1[j] - 28352 * LINK->f2[j]);
  _getyp_(LINK->f5, dt / 2, LINK);
  if (Errorflag)
    return;

  d = dt / 7618050.0;
  FORLIM = LINK->n;
  for (j = 0; j < FORLIM; j++)
    LINK->f1[j] = LINK->y0[j] + d * (902880.0 * LINK->yy[j] +
		      3855735.0 * LINK->f3[j] - 1371249.0 * LINK->f4[j] +
		      3953664.0 * LINK->f2[j] + 277020.0 * LINK->f5[j]);
}  /*rkfehl*/

/* initialise array yp */
static void _inityp(struct LOC_Rkf45 *LINK)
{
  _getyp_(LINK->yy, 0.0, LINK);
}  /*inityp*/


/*dnh changed fp to F - to be consistent!*/

void Rkf45CL(int n_, double aber, double reler, double tout, double *t_, void (*F_)(double t, double *y, double *ydash), double *y_, double *yp_, int *rkf)
{

  /* Solves ordinary differential equations in non-stiff or mildly*/
  /*stiff cases. Integration stepsize is automatically*/
  /*chosen by the routine to meet error criteria set by the user.*/

  /*The ODE system is of nth order*/
  /*y'j = fj(t, y1, ..... , yn)    j=1,...n, n<=30  */
  struct LOC_Rkf45 V;
  int nf, k0, inn, jf, kf;
  double h, re1, ae1, u26, eps;
  int k, maxnf, mf;
  boolean op, hf;
  double re0, ae, dt, est, hmin, rer, s, scale, tol, w1, w2;
  int FORLIM;

  V.n = n_;
  V.t = t_;
  V.F = F_;
  V.y = y_;
  V.yp = yp_;
  /*------------------------------------initialsation*/
  nf = 0;
  k0 = 0;
  inn = 0;
  jf = 0;
  kf = 0;
  h = 0.0;
  re1 = 0.0;
  ae1 = 0.0;
  u26 = 0.0;
  eps = 0.0;
  /*------------------------------------*/
  if (Errorflag)
    return;

  /* test order*/
  if (V.n < 1) {  /*'Order must be >= 1'*/
    _setmatherror(9);
    return;
  }

  /* test max order */
  if (V.n > MAXN)  /*'Order must be <= 30'*/
    _setmatherror(10);

  if (*V.t > tout) {  /*'Time Step must be positive'*/
    _setmatherror(11);
    return;
  }

  /* Test tolerances */
  if (aber < 0.0 || reler < 0.0) {
    /*'Tolerance must not be negative'*/
    _setmatherror(8);
    return;
  }

  FORLIM = V.n;
  /* save y at start of step */
  for (k = 0; k < FORLIM; k++)
    V.y0[k] = V.y[k];

  re0 = 1e-12;
  maxnf = 3000;
  mf = abs(*rkf);

  if (V.n < 0 || reler < 0 || aber < 0 || mf == 0 || mf > 8 ||
      mf != 1 && *V.t == tout && kf != 3) {
    *rkf = 8;
    return;   /* bad input */
  }

  if (mf == 1)   /* 'machine accuracy'= epsilon/2 */
    u26 = 26 * EPSILON / 2;
  else if (mf == 2) {
    if (kf == 4)   /* zero fn eval count */
      nf = 0;
    if (kf == 3 || inn == 0 || kf == 4 && mf != 2) {
      *rkf = jf;
      if (kf == 3)   /* set flag from earlier call */
	mf = abs(*rkf);
    } else if (kf == 5 && aber == 0 || kf == 6 && reler <= re1 && aber <= ae1)
      return;
  } else {
    if (*rkf == 4)
      nf = 0;
    /* zero fn eval count */
    if (*rkf != 3 && (*rkf != 4 || mf == 2) && (*rkf != 5 || aber <= 0))
      return;
    *rkf = jf;
    if (kf == 3)
      mf = abs(*rkf);   /* set flag from earlier call */
  }

  /* save inputs, set contin flag ... */

  jf = *rkf;
  re1 = reler;
  ae1 = aber;
  kf = 0;

  /* test error tol ... */

  rer = 2 * eps + re0;
  if (reler < rer) {
    reler = rer;
    *rkf = 3;
    kf = 3;
    return;
  }
  dt = tout - *V.t;

  /* initialise ...*/

  if (mf == 1) {
    inn = 0;
    k0 = 0;
    nf = 0;
    _inityp(&V);
    if (Errorflag)
      return;
    if (*V.t == tout) {
      *rkf = 2;
      return;
    }
  }
  if (inn == 0) {
    inn = 1;
    h = fabs(dt);
    FORLIM = V.n;
    for (k = 0; k < FORLIM; k++) {
      tol = reler * fabs(V.y0[k]) + aber;
      if (tol > 0) {
	w1 = tol;
	w2 = fabs(V.yy[k]);
	if (w2 > 0) {
	  if (log(w2) + 5 * log(h) > log(tol))
	    h = exp(0.2 * log(tol / w2));
	}
      }
    }
    w2 = fabs(*V.t);
    if (w2 < fabs(dt))
      w2 = fabs(dt);
    if (w1 <= 0)
      h = 0.0;
    w2 = u26 * w2;
    if (w2 > h)
      h = w2;
    jf = _Sgn(*rkf) * 2;
  }

  h = _Sgn(dt) * h;

  /* test if too stiff ... */

  if (fabs(h) >= 2 * fabs(dt))
    k0++;
  if (k0 >= 100) {
    k0 = 0;
    *rkf = 7;
    return;
  }
  if (fabs(dt) <= u26 * fabs(*V.t)) {
    FORLIM = V.n;
    for (k = 0; k < FORLIM; k++) {
      V.y0[k] += dt * V.yy[k];
      V.y[k] = V.y0[k];
    }
    nf++;
    _inityp(&V);
    if (Errorflag)
      return;
    *V.t = tout;
    *rkf = 2;
    return;
  }

  /* prepare for step ...*/

  op = false;
  scale = 2 / reler;
  ae = scale * aber;

  do {
    hf = false;
    hmin = u26 * fabs(*V.t);
    dt = tout - *V.t;
    if (fabs(dt) < 2 * fabs(h)) {
      if (fabs(dt) > fabs(h))
	h = dt / 2;
      else {
	op = true;
	h = dt;
      }
    }

    /* test no. of fn evals ... */

    do {

      if (nf > maxnf) {
	*rkf = 4;
	kf = 4;
	return;
      }

      /* an approx step ... */

      _rkfehl(h, &V);
      nf += 5;

      /* test tolerance ...*/

      est = 0.0;
      FORLIM = V.n;
      for (k = 0; k < FORLIM; k++) {
	w2 = fabs(V.y0[k]) + fabs(V.f1[k]) + ae;
	if (w2 <= 0) {
	  *rkf = 5;
	  return;
	}
	w1 = fabs(21970 * V.f3[k] - 2090 * V.yy[k] - 15048 * V.f4[k] +
		  22528 * V.f2[k] - 27360 * V.f5[k]) / w2;
	if (w1 > est)
	  est = w1;
      }
      est = fabs(h) * est * scale / 752400.0;

      /* unsuccessful step - try again ...*/

      if (est > 1) {
	hf = true;
	op = false;
	s = 0.1;
	if (est < 59049.0)
	  s = 0.9 / exp(log(est) / 5);
	h = s * h;
	if (fabs(h) <= hmin) {
	  *rkf = 6;
	  kf = 6;
	  return;
	}
      }

    } while (est > 1);

    *V.t += h;
    memcpy(V.y0, V.f1, sizeof(odearray_));
    FORLIM = V.n;
    for (k = 0; k < FORLIM; k++)
      V.y[k] = V.y0[k];

    _inityp(&V);
    if (Errorflag)
      return;
    nf++;

    /* choose next h ...*/

    s = 5.0;
    if (est > 1.889568e-4)
      s = 0.9 / exp(log(est) / 5);
    if (hf && s > 1)
      s = 1.0;
    w1 = s * fabs(h);
    if (hmin > w1)
      w1 = hmin;
    h = w1 * _Sgn(h);

    /* another step ? ...*/

    if (op) {
      *V.t = tout;
      *rkf = 2;
      return;
    }

  } while (*rkf > 0);

  *rkf = -2;

}  /*of Rkf45*/

#undef MAXN

static void
Fft(double *f, int m, double sign)
{

  /* Fast Fourier Transform algorithm for Discrete Fourier Transforms.*/
  /*FFT works on n=2^m points and needs an array f[j,k].*/
  /*f[j,0] contains real part, f[j,1] contains imag part*/
  /*of the jth data value.  */
  double a, b, sn, cs, s1, c1;
  int i, i2, j, j2, k, l, l2, n, m1, n1;


  if (Errorflag)
    return;

  if (m < 1) {  /*'Power of 2 must be >= 1'*/
    _setmatherror(13);
    return;
  }

  /* force s=+/-1 */
  sign = _Sgn(_Sgn(sign) + 0.1);

  /* re-order f ....  */

  /* n:=2^m */
  n = 1;
  for (i = 1; i <= m; i++)
    n += n;
  n1 = n / 2;
  i = 1;

  for (j = 1; j < n; j++) {
    if (i > j) {
      for (l = -1; l <= 0; l++) {
	i2 = i + i + l;
	j2 = j + j + l;
	a = f[i2 - 1];
	f[i2 - 1] = f[j2 - 1];
	f[j2 - 1] = a;
      }
    }
    k = n1;
    while (i > k) {
      i -= k;
      k /= 2;
    }
    i += k;
  }

  /* now compute FT ... */
  n1 = 1;
  c1 = -1.0;
  s1 = 0.0;

  for (m1 = 1; m1 <= m; m1++) {
    k = n1;
    n1 += n1;
    cs = 1.0;
    sn = 0.0;

    for (i = 1; i <= k; i++) {
      j = i;
      do {
	j2 = j + j - 1;
	l2 = j2 + k + k;
	a = f[l2 - 1] * cs - f[l2] * sn;
	b = f[l2 - 1] * sn + f[l2] * cs;
	f[l2 - 1] = f[j2 - 1] - a;
	f[l2] = f[j2] - b;
	f[j2 - 1] += a;
	f[j2] += b;
	j += n1;
      } while (j <= n);

      a = cs * c1 - sn * s1;
      sn = sn * c1 + cs * s1;
      cs = a;
    }

    s1 = sign * sqrt((1 - c1) / 2);
    c1 = sqrt((1 + c1) / 2);
  }

  /* no errors */
  Errorflag = false;

}  /*fft*/


/*------------------------Fnj (Bessel function)--------------*/
double BesselCL(int n, double x)
{
  int i, k;
  double JJ, JJold, j, j0, j1, sum, tol, a;
  boolean odd;

  if (n < 0)
    return 0.0;   /*fnj*/
  else {
    tol = 1.0e-9;
    JJ = 0.0;

    /*Deal with x=0*/
    a = fabs(x);
    if (a < tol) {
      if (n == 0)
	return 1.0;
      else
	return 0.0;
    } else {
      /*Starting order K for recurrence*/
      if (a >= 5.0)
	i = (long)fabs(1.4 * x + 60.0 / x);
      else
	i = (long)(a + 6.0);
      k = (long)(n + 2.0 + a / 4.0);
      if (i > k)
	k = i;

      /*Loop until estimates agree*/
      do {
	JJold = JJ;
	j1 = 0.0;
	j0 = 1.0e-30;
	sum = 0.0;
	odd = ((k & 1) == 1);

	/*Apply recurrence down from current max order*/
	for (i = 1; i <= k - 2; i++) {
	  j = 2.0 * (k - i) * j0 / x - j1;
	  j1 = j0;
	  j0 = j;
	  if (k - i - 1 == n)
	    JJ = j;
	  odd = !odd;
	  if (odd)
	    sum += j * 2;
	}
	j = 2.0 * j / x - j1;
	if (n == 0)
	  JJ = j;
	sum += j;
	JJ /= sum;
	k += 3;
      } while (fabs(JJ - JJold) >= tol);
      return JJ;
    }
  }


}



void TridiagCL(int m, double *ap, double *bp, double *cp, double *rp, double *xp)
{
  /*--------------------------------------------------*/
  /*Solves set of tridiagonal equations*/
  /*m equations.           a[i], b[i], c[i]: LHS coefficients*/
  /* r[i]: RHS;         g[i], z[i]: temporary vectors              x[i]:  answers*/
  /*----------------------------------------------------*/
  int j, k, sz;
  unsigned short sw;
  double *gp, *zp;

  /*----------------------------------*/
  sz = m * sizeof(double);
  /*gp := arrayptr(newptr(sz * 2));*/
  sw = m * 2;
  gp = _New_array(sw);
  zp = (double *)((long)gp + sz);
  ap[0] = 0.0;
  cp[m - 1] = 0.0;
  zp[0] = rp[0] / bp[0];
  gp[0] = cp[0] / bp[0];
  /*---------------------------------forward recursion*/
  for (j = 2; j <= m; j++) {
    zp[j - 1] = (rp[j - 1] - ap[j - 1] * zp[j - 2]) /
		(bp[j - 1] - gp[j - 2] * ap[j - 1]);
    gp[j - 1] = cp[j - 1] / (bp[j - 1] - gp[j - 2] * ap[j - 1]);
  }
  /*----------------------------------backward recursion*/
  k = m;
  xp[k - 1] = zp[k - 1];
  do {
    k--;
    xp[k - 1] = zp[k - 1] - gp[k - 1] * xp[k];
  } while (k != 1);
  /*--------------------------------*/
  /*disposPtr(mptr(gp));*/
  _Free_array(&gp, sw);
}


typedef double fars[1][2];


/* Local variables for sinfft: */
struct LOC_sinfft {
  double *fvs;
  int n;
  double (*fp)[2];
} ;

/*----------------------extend array to give odd function for Fft*/
static void _oddfn(struct LOC_sinfft *LINK)
{
  int j, ny2, ny2j;
  double v;
  int FORLIM;

  ny2 = LINK->n * 2;
  FORLIM = LINK->n;
  for (j = 0; j <= FORLIM; j++) {
    ny2j = ny2 - j;
    v = LINK->fvs[j];
    LINK->fp[j][0] = v;
    LINK->fp[j][1] = 0.0;
    LINK->fp[ny2j][0] = -v;
    LINK->fp[ny2j][1] = 0.0;
  }
}


/*...........................................CTsinft...................................................*/
void FftSinCL(int pwr, double *fvs_, int v)
{
  /*----------------------------------------------------------*/
  /*Fourier sine transform - normal and inverse*/
  /*pwr:  2**pwr +1 array elements  (pwr=4  gives arrays 0..16)*/
  /*fvs:  input/output array.    v=1  normal.     v=-1  inverse.*/
  /*fvs doubled in temporary array for Fft, with pwa = pwr+1*/
  /*----------------------------------------------------------*/
  struct LOC_sinfft V;
  int k, pwa;
  unsigned short sw;
  double d;
  int FORLIM;
  double *TEMP;

  /*---------------------------------*/
  V.fvs = fvs_;
  V.n = 1;
  for (k = 1; k <= pwr; k++)
    V.n *= 2;
  d = V.n;
  sw = (V.n + 1) * 4;
  /*--------------------------------get memory for CSL Fft*/
  /*fp := farp(newptr(n1 * sizeof(citreal))); */
  V.fp = (double(*)[2])_New_array(sw);
  /*----------------------------------transforms in fp*/
  _oddfn(&V);
  pwa = pwr + 1;
  Fft(V.fp[0], pwa, 1.0);
  /*-------------------------update original array - normal and inverse*/
  if (v == 1) {
    FORLIM = V.n;
    for (k = 1; k < FORLIM; k++)   /*divide by 2n for normal transform */
      V.fvs[k] = V.fp[k][1] / 2.0;
  } else {
    FORLIM = V.n;
    for (k = 1; k < FORLIM; k++)
      V.fvs[k] = V.fp[k][1] / d;
  }
  TEMP = (double *)V.fp;
  /*---------------------------*/
  /*disposPtr(mptr(fp));*/
  _Free_array(&TEMP, sw);
}


/*--------------------------Poisson2d---------------------------------*/
static void
_Poisson2d(double *psp, double *ztp, int nnx, int nny, double dlx, double dly)
{
  /*---------------------------------------------------------*/
  /*Fast Poisson solver, given zeta matrix, returning psi matrix, using*/
  /*Fast Fourier transforms along y values, and solving tridiagonal*/
  /*equations along the x values*/
  /*psp, ztp point to calling arrays [0..nnx, 0..nny] - in effect, na0 = nny+1*/
  /*dlx,dly are x and y increments*/
  /*nny must be a power of 2*/
  /*---------------------------------------------------------*/
  int i, j, pw, k, kt, na0, sz, nn, na, nb;
  unsigned int sw1, sw2;
  double dlxx, dlyy, bi;
  double *ap, *bp, *cp, *rp, *ansp, *fztp;

  /*-------------------------------*/
  na0 = nny + 1;
  nn = 1;
  pw = 0;
  while (nn < nny) {
    nn += nn;
    pw++;
  }
  /*----------------------------*/
  if (nn != nny) {
    _CTerror("x dimension not a power of 2");
    Errorflag = true;
    return;
  }
  /*------------------------------memory pointers for tridiagonal coefficients*/
  sz = nnx * sizeof(double);
  /*ap := arrayptr(newptr(sz * 5));*/
  sw1 = nnx * 5;
  ap = _New_array(sw1);
  bp = (double *)((long)ap + sz);
  cp = (double *)((long)bp + sz);
  rp = (double *)((long)cp + sz);
  ansp = (double *)((long)rp + sz);
  /*----------------------------------memory pointer for zeta transform*/
  sw2 = (nnx + 1) * na0;
  kt = sw2;
  /*fztp := arrayptr(newptr(kt * sizeof(citreal)));*/
  fztp = _New_array(sw2);
  /*-----------------------------------------initialise psi and zeta transform*/
  dlxx = 1.0 / (dlx * dlx);
  dlyy = 1.0 / (dly * dly);
  for (k = 0; k < kt; k++)
    fztp[k] = ztp[k];
  /*--------------------------amend fzeta to allow for y-boundary conditions*/
  na = na0 + 1;
  nb = na0 * (nnx - 1) + 1;
  for (j = 0; j <= nny - 2; j++) {
    fztp[na + j] -= psp[j + 1] * dlxx;
    fztp[nb + j] -= psp[nb + na0 + j] * dlxx;
  }
  /*------------------------------x-boundary conditions*/
  na = 2;
  nb = na0 - 1;
  for (i = 1; i < nnx; i++) {
    na += na0;
    nb += na0;
    fztp[na - 1] -= psp[na - 2] * dlyy;
    fztp[nb - 1] -= psp[nb] * dlyy;
  }
  /*-------------------------------get all y fourier transforms*/
  for (i = 1; i < nnx; i++)
    FftSinCL(pw, &fztp[i * na0], 1);
  /*-------------------------------set up tridiagonal coefficients*/
  for (j = 1; j < nny; j++) {
    bi = 2 * ((cos(M_PI * j / nny) - 1.0) * dlyy - dlxx);
    for (i = 0; i <= nnx - 2; i++) {
      ap[i] = dlxx;
      bp[i] = bi;
      cp[i] = dlxx;
      rp[i] = fztp[(i + 1) * na0 + j];
    }
    /*-----------------------------------solve equations*/
    TridiagCL(nnx - 1, ap, bp, cp, rp, ansp);
    for (i = 1; i < nnx; i++)
      psp[i * na0 + j] = ansp[i - 1];
  }
  /*-----------------------------------inverse sine transform for psi*/
  for (i = 1; i < nnx; i++)
    FftSinCL(pw, &psp[i * na0], -1);
  /*-----------------------------------release memory*/
  /*disposPtr(mptr(ap));
  disposPtr(mptr(fztp));*/
  _Free_array(&ap, sw1);
  _Free_array(&fztp, sw2);
}


/*---------------------------normal distribution function*/
double _Npd(double x)
{
  /*dnh far;*/
  return (exp(-x * x * 0.5) / 2.5066282746);   /* / sqrt(2*PI*/
}


/*--------------------------cumulative p.d. for N(0,1)*/
double PhiCL(double x)
{
  double pbase, er;
  double (*TEMP)(double x);

  pbase = 0.0000317848;   /* pbase := _Romb(-50.0, -4.0, 7, Npd, err); */
  if (x <= -4.0)
    return pbase;
  else {
    TEMP = _Npd;
    return (_Romb(-4.0, x, 7, TEMP, &er) + pbase);
  }
}


/*----------------------------inverse Phi function*/
double InvPhiCL(double pa)
{
  double pbase, x, xi, pp, p1, p2;
  int k;

  pbase = 0.0000317848;
  if (pa < pbase)
    return -4.0;
  if (pa > 1.0 - pbase)
    return 4.0;
  /*----------------------*/
  x = -4.0;
  if (pa > 0.5)
    x = 0.0;
  xi = 1.0;
  for (k = 1; k <= 3; k++) {
    do {
      x += xi;
      pp = PhiCL(x);
    } while (pa >= pp);
    x -= xi;
    if (k != 3)
      xi *= 0.1;
  }
  p1 = PhiCL(x);
  p2 = PhiCL(x + xi);
  return (x + xi * (pa - p1) / (p2 - p1));
}


#define MAXN            30

#define MAXNSTRING      "30"

#define WORKN           450   /* = n(n-1)/2 */
#define MAXITER         100

#define MAXITERSTRING   "100"

static int pv[MAXN];

static void
_errormatrix(int code)
{
  char str[LINESTRINGSIZE + 1];

  switch (code) {

  case 1:
    strcpy(str, "Iteration for singular value failed");
    break;

  case 2:
    sprintf(str, "Accuracy not achieved after %s iterations in Power Method",
	    MAXITERSTRING);
    break;

  case 3:
    strcpy(str, "Unsuitable initial value for eigenvalue");
    break;

  case 4:
    sprintf(str,
      "Accuracy not achieved after %s iterations in Modified Power Method",
      MAXITERSTRING);
    break;

  case 5:
    strcpy(str, "Array size must be >= 1");
    break;

  case 6:
    sprintf(str, "Array size must be <= %s", MAXNSTRING);
    break;

  case 7:
    strcpy(str, "Accuracy must be > 0");
    break;
  }

  _Set_error(code + 21000, str);
}


/* check sizes and dimensions */
static void
_size_check(int n, int n0)
{
  if (n <= 0 || n0 <= 0)
    _errormatrix(5);
  else if (n > MAXN || n0 > MAXN)
    _errormatrix(6);
}

static double *_local_copy(double *a, int n)
{
  double *Result, *p;

  /* add this to prevent C compilation warnings*/
  Result = NULL;
  p = _New_array(n);
  /*2.41*/
  /*--- D   if Errorflag then*/
  if (p == NULL)
    return Result;
  /*------------------------------------££*/
  memmove(p, a, n * sizeof(double));
  /*commented in the original...*/
  /*BlockMove(Ptr(a), Ptr(p), n * SizeOf(Citreal));*/

  /*---------------------------------------*/
  return p;
}

static void
_CTsolve(int n, int n0, double *a, double *b)
{
  int i, k;
  double t;
  double *av[MAXN];

  /*2.41*/
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/

  _size_check(n, n0);
  if (Errorflag)
    return;

  for (i = 0; i < n; i++)
    av[i] = (void *)((long)a + i * n0 * sizeof(double));

  if (n > 1) {
    for (k = 0; k <= n - 2; k++) {   /* forward elimination */
      i = pv[k];
      t = b[i - 1];
      b[i - 1] = b[k];
      b[k] = t;
      for (i = k + 1; i < n; i++)
	b[i] += av[i][k] * t;
    }

    for (k = n; k >= 2; k--) {   /* back substitution */
      t = -(b[k - 1] / av[k - 1][k - 1]);
      b[k - 1] = -t;
      for (i = 0; i <= k - 2; i++)
	b[i] += av[i][k - 1] * t;
    }

  }


  b[0] /= av[0][0];
}  /*solve*/

static void
_CTdecomp(int n, int n0, double *a, double *det, double *cond)
{
  double ek, t, an, yn, zn;
  int i, j, k, m;
  double b[MAXN];
  double *av[MAXN];

  /*2.41*/
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/
  /*--- D   */
  _size_check(n, n0);
  if (Errorflag)
    return;

  pv[n - 1] = 1;
  *cond = 0.0;

  /* start permutation indicator */
  for (i = 0; i < n; i++)
    av[i] = (void *)((long)a + i * n0 * sizeof(double));

  /* trap 1x1 case */
  if (n == 1) {
    *cond = 1.0;
    *det = av[0][0];
    if (*det != 0)
      return;

    /* exact singularity */
    *cond = 1e32;
    return;

  }

  /* find norm A ... */

  an = 0.0;
  for (j = 0; j < n; j++) {
    t = 0.0;
    for (i = 0; i < n; i++)
      t += fabs(av[i][j]);
    if (t > an)
      an = t;
  }

  /* Gaussian elimination, partial pivot ...*/

  for (k = 1; k < n; k++) {  /* find pivot */
    m = k;
    for (i = k + 1; i <= n; i++) {
      if (fabs(av[i - 1][k - 1]) > fabs(av[m - 1][k - 1]))
	m = i;
    }

    pv[k - 1] = m;   /* update pivot vector */
    if (m != k)
      pv[n - 1] = -pv[n - 1];
    t = av[m - 1][k - 1];
    av[m - 1][k - 1] = av[k - 1][k - 1];
    av[k - 1][k - 1] = t;

    if (t == 0)
      *cond = -1.0;   /* skip elim if pivot=0 */
    else {
      for (i = k; i < n; i++)
	av[i][k - 1] = -(av[i][k - 1] / t);
      for (j = k; j < n; j++) {
	t = av[m - 1][j];
	av[m - 1][j] = av[k - 1][j];
	av[k - 1][j] = t;
	if (t != 0) {
	  for (i = k; i < n; i++)
	    av[i][j] += av[i][k - 1] * t;
	}
      }
    }
  }

  /* find condition of A ... */

  if (*cond == -1 || av[n - 1][n - 1] == 0) {
    *det = 0.0;
    *cond = 1e32;
    return;

  }

  for (k = 1; k <= n; k++) {   /* solve (A tr) Y = E */
    t = 0.0;
    if (k > 1) {
      for (i = 0; i <= k - 2; i++)
	t += av[i][k - 1] * b[i];
    }
    if (t < 0)
      ek = -1.0;
    else
      ek = 1.0;
    b[k - 1] = -((ek + t) / av[k - 1][k - 1]);
  }

  for (k = n - 2; k >= 0; k--) {
    t = 0.0;
    for (i = k + 1; i < n; i++)
      t += av[i][k] * b[k];
    b[k] = t;
    m = pv[k];
    if (m != k + 1) {
      t = b[m - 1];
      b[m - 1] = b[k];
      b[k] = t;
    }
  }

  yn = 0.0;
  for (i = 0; i < n; i++)
    yn += fabs(b[i]);

  _CTsolve(n, n0, a, b);

  /* solve A X = B */
  zn = 0.0;
  for (i = 0; i < n; i++)
    zn += fabs(b[i]);

  /* get cond from norms */
  *cond = an * zn / yn;
  if (*cond < 1)
    *cond = 1.0;

  /* find the determinant by multiplying */
  /* pv[n]*a[1,1]*a[2,2]*....*a[n,n]     */
  *det = pv[n - 1];
  for (i = 0; i < n; i++)
    *det *= av[i][i];

  /* no errors */
  Errorflag = false;

}  /*decomp*/


/* Local variables for CTinvert: */
struct LOC_CTinvert {
  int n;
  double *av[MAXN];
} ;

static double _fnorm(struct LOC_CTinvert *LINK)
{
  int i, j;
  double s, norm;
  int FORLIM, FORLIM1;

  norm = 0.0;
  FORLIM = LINK->n;
  for (i = 0; i < FORLIM; i++) {
    s = 0.0;
    FORLIM1 = LINK->n;
    for (j = 0; j < FORLIM1; j++)
      s += fabs(LINK->av[i][j]);
    if (s > norm)
      norm = s;
  }
  return norm;
}  /*fnorm*/


static void
_CTinvert(int n_, int n0, double *a, double *det, double *cond)
{
  struct LOC_CTinvert V;
  double aa, b, d, norm;
  int i, j, k, m;
  double a0[MAXN];
  int FORLIM, FORLIM1, FORLIM2;

  /*2.41*/
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/

  V.n = n_;
  _size_check(V.n, n0);
  if (Errorflag)
    return;

  *det = 1.0;
  *cond = 1e32;

  FORLIM = V.n;
  for (i = 0; i < FORLIM; i++)   /***/
    V.av[i] = (void *)((long)a + i * n0 * sizeof(double));

  /* trap n=1 ... */
  if (V.n == 1) {
    *det = V.av[0][0];
    if (*det == 0)
      return;
    else {
      V.av[0][0] = 1 / *det;
      *cond = 1.0;
      return;

    }

  }

  /* norm of initial matrix ... */
  norm = _fnorm(&V);

  FORLIM = V.n;
  /* initial pivot vector ... */
  for (i = 1; i <= FORLIM; i++)
    pv[i - 1] = i;

  /* main elimination loop ... */
  k = 0;
  do {
    k++;
    aa = 0.0;

    FORLIM = V.n;
    /* look for pivot ... */
    for (i = k; i <= FORLIM; i++) {
      if (fabs(V.av[i - 1][k - 1]) > aa) {
	m = i;
	aa = fabs(V.av[i - 1][k - 1]);
      }
    }

    if (aa == 0)
      *det = 0.0;   /* if det:=0 then matrix singular .... */
    else {
      /* swap rows ... */
      if (k != m) {
	FORLIM = V.n;
	for (j = 0; j < FORLIM; j++) {
	  aa = V.av[k - 1][j];
	  V.av[k - 1][j] = V.av[m - 1][j];
	  V.av[m - 1][j] = aa;
	}
	*det = -*det;
	i = pv[k - 1];
	pv[k - 1] = pv[m - 1];
	pv[m - 1] = i;
      }

      /* update determinant ... */
      d = V.av[k - 1][k - 1];
      *det *= d;

      FORLIM = V.n;
      /* eliminate terms below pivot ... */
      for (i = k; i < FORLIM; i++) {
	/* find factor ... */
	aa = V.av[i][k - 1] / d;

	FORLIM1 = V.n;
	/* work along row;  j<:=k in L, j>k in U ... */
	for (j = 1; j <= FORLIM1; j++) {
	  if (j == k)
	    V.av[i][j - 1] = -aa;
	  else {
	    V.av[i][j - 1] -= aa * V.av[k - 1][j - 1];

	  }
	}
      }

    }
  } while (k != V.n - 1 && *det != 0);

  if (*det == 0)
    return;


  /* normalise bottom row ... */
  aa = V.av[V.n - 1][V.n - 1];
  *det = aa * *det;
  if (aa == 0)
    return;


  d = 1 / aa;
  V.av[V.n - 1][V.n - 1] = d;
  FORLIM = V.n - 2;
  for (j = 0; j <= FORLIM; j++)
    V.av[V.n - 1][j] = d * V.av[V.n - 1][j];

  /* back substitute ... */
  for (i = V.n - 1; i >= 1; i--) {
    d = V.av[i - 1][i - 1];
    FORLIM1 = V.n;
    for (j = i; j < FORLIM1; j++)
      a0[j] = V.av[i - 1][j];

    FORLIM1 = V.n;
    /* subtract multiples of lower rows ... */
    for (j = 0; j < FORLIM1; j++) {
      aa = 0.0;
      FORLIM2 = V.n;
      for (m = i; m < FORLIM2; m++)
	aa += a0[m] * V.av[m][j];
      if (j + 1 < i)
	b = V.av[i - 1][j];
      else if (j + 1 == i)
	b = 1.0;
      else
	b = 0.0;
      V.av[i - 1][j] = (b - aa) / d;
    }
  }

  FORLIM = V.n;
  /* column permutations ... */
  for (j = 0; j < FORLIM; j++) {
    k = j;
    do {
      k++;
    } while (pv[k - 1] != j + 1);
    if (k > j + 1) {
      pv[k - 1] = pv[j];
      FORLIM1 = V.n;
      for (i = 0; i < FORLIM1; i++) {
	aa = V.av[i][j];
	V.av[i][j] = V.av[i][k - 1];
	V.av[i][k - 1] = aa;
      }
    }
  }

  *cond = norm * _fnorm(&V);

  /* no errors */
  Errorflag = false;




}  /*invert*/


/* Local variables for CTsvd: */
struct LOC_CTsvd {
  int n0;
} ;

static int _index_(int i, int j, struct LOC_CTsvd *LINK)
{
  return ((i - 1) * LINK->n0 + j);
}

static void
_CTsvd(int m, int n, int n0_, double *a, double *u, double *v, double *w, int *svdflag)
{
  struct LOC_CTsvd V;
  double c, f, g, h, s, x, y, z, scale, norm;
  int i, j, k, l, i1, k1, l1, mn, its;
  double rv1[MAXN];
  boolean rv1small, wsmall, fsmall, again;

  V.n0 = n0_;
  /*2.41*/
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/
  _size_check(m, V.n0);
  if (Errorflag)
    return;
  _size_check(n, V.n0);
  if (Errorflag)
    return;
  *svdflag = 0;
  g = 0.0;
  scale = 0.0;
  norm = 0.0;
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++)
      u[_index_(i, j, &V) - 1] = a[_index_(i, j, &V) - 1];
  }

  /* Householder reduction to bi-diagonal form ... */
  for (i = 1; i <= n; i++) {
    l = i + 1;
    rv1[i - 1] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i <= m) {
      for (k = i; k <= m; k++)
	scale += fabs(u[_index_(k, i, &V) - 1]);
      if (scale != 0) {
	for (k = i; k <= m; k++) {
	  u[_index_(k, i, &V) - 1] /= scale;
	  s += u[_index_(k, i, &V) - 1] * u[_index_(k, i, &V) - 1];
	}
	f = u[_index_(i, i, &V) - 1];
	g = -sqrt(s) * _Sgn(f);
	h = f * g - s;
	u[_index_(i, i, &V) - 1] = f - g;
	if (i < n) {
	  for (j = l; j <= n; j++) {
	    s = 0.0;
	    for (k = i; k <= m; k++)
	      s += u[_index_(k, i, &V) - 1] * u[_index_(k, j, &V) - 1];
	    f = s / h;
	    for (k = i; k <= m; k++)
	      u[_index_(k, j, &V) - 1] += f * u[_index_(k, i, &V) - 1];
	  }
	}
	for (k = i; k <= m; k++)
	  u[_index_(k, i, &V) - 1] = scale * u[_index_(k, i, &V) - 1];
      }
    }

    w[i - 1] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i <= m && i != n) {
      for (k = l; k <= n; k++)
	scale += fabs(u[_index_(i, k, &V) - 1]);
      if (scale != 0) {
	for (k = l; k <= n; k++) {
	  u[_index_(i, k, &V) - 1] /= scale;
	  s += u[_index_(i, k, &V) - 1] * u[_index_(i, k, &V) - 1];
	}
	f = u[_index_(i, l, &V) - 1];
	g = -sqrt(s) * _Sgn(f);
	h = f * g - s;
	u[_index_(i, l, &V) - 1] = f - g;
	for (k = l; k <= n; k++)
	  rv1[k - 1] = u[_index_(i, k, &V) - 1] / h;
	if (i != m) {
	  for (j = l; j <= m; j++) {
	    s = 0.0;
	    for (k = l; k <= n; k++)
	      s += u[_index_(j, k, &V) - 1] * u[_index_(i, k, &V) - 1];
	    for (k = l; k <= n; k++)
	      u[_index_(j, k, &V) - 1] += s * rv1[k - 1];
	  }
	}
	for (k = l; k <= n; k++)
	  u[_index_(i, k, &V) - 1] = scale * u[_index_(i, k, &V) - 1];
      }
    }

    c = fabs(w[i - 1]) + fabs(rv1[i - 1]);
    if (c > norm) {
      norm = c;

    }
  }

  /* Accumulation of r-hand transforms ... */
  for (i = n; i >= 1; i--) {
    if (i < n) {
      if (g != 0) {
	for (j = l; j <= n; j++)
	  v[_index_(j, i, &V) - 1] = u[_index_(i, j, &V) - 1] /
				    u[_index_(i, l, &V) - 1] / g;
	for (j = l; j <= n; j++) {
	  s = 0.0;
	  for (k = l; k <= n; k++)
	    s += u[_index_(i, k, &V) - 1] * v[_index_(k, j, &V) - 1];
	  for (k = l; k <= n; k++)
	    v[_index_(k, j, &V) - 1] += s * v[_index_(k, i, &V) - 1];
	}
      }
      for (j = l; j <= n; j++) {
	v[_index_(i, j, &V) - 1] = 0.0;
	v[_index_(j, i, &V) - 1] = 0.0;
      }
    }
    v[_index_(i, i, &V) - 1] = 1.0;
    g = rv1[i - 1];
    l = i;
  }

  /* Accumulation of l-hand transforms ... */
  mn = n;
  if (m < n)
    mn = m;
  for (i = mn; i >= 1; i--) {
    l = i + 1;
    g = w[i - 1];
    if (i < n) {
      for (j = l; j <= n; j++)
	u[_index_(i, j, &V) - 1] = 0.0;
    }
    if (g != 0) {
      if (i < mn) {
	for (j = l; j <= n; j++) {
	  s = 0.0;
	  for (k = l; k <= m; k++)
	    s += u[_index_(k, i, &V) - 1] * u[_index_(k, j, &V) - 1];
	  f = s / u[_index_(i, i, &V) - 1] / g;
	  for (k = i; k <= m; k++)
	    u[_index_(k, j, &V) - 1] += f * u[_index_(k, i, &V) - 1];
	}
      }
      for (j = i; j <= m; j++)
	u[_index_(j, i, &V) - 1] /= g;
    } else {
      for (j = i; j <= m; j++)
	u[_index_(j, i, &V) - 1] = 0.0;
    }
    u[_index_(i, i, &V) - 1]++;
  }

  /* Diagonalisation of bi-diag form ... */
  for (k = n; k >= 1; k--) {
    k1 = k - 1;
    its = 0;

    /* test for splitting ... */
    do {
      l = k + 1;
      do {
	l--;
	l1 = l - 1;
	rv1small = (fabs(rv1[l - 1]) + norm == norm);
	/*wsmall := (Abs(w^[l]) + norm = norm); ££ -- error here --££*/
	if (l1 > 0)
	  wsmall = (fabs(w[l1 - 1]) + norm == norm);
	else
	  wsmall = true;
      } while (!(rv1small || wsmall));
      /* rv1[1]=0, so loop will stop! */

      /* Cancellation of rv1[l] if l>1 .... */
      if (!rv1small) {
	c = 0.0;
	s = 1.0;
	i = l;
	do {
	  f = s * rv1[i - 1];
	  rv1[i - 1] = c * rv1[i - 1];
	  fsmall = (fabs(f) + norm == norm);
	  if (!fsmall) {
	    g = w[i - 1];
	    h = sqrt(f * f + g * g);
	    w[i - 1] = h;
	    c = g / h;
	    s = -(f / h);
	    for (j = 1; j <= m; j++) {
	      y = u[_index_(j, l1, &V) - 1];
	      z = u[_index_(j, i, &V) - 1];
	      u[_index_(j, l1, &V) - 1] = y * c + z * s;
	      u[_index_(j, i, &V) - 1] = z * c - y * s;
	    }
	    i++;
	  }
	} while (!(i > k || fsmall));
      }

      z = w[k - 1];

      again = (l != k);
      if (again) {
	/* shift from bottom 2x2 minor ... */
	if (its == 30) {
	      /*'Iteration for singular value failed to converge'*/
		_errormatrix(1);
	  *svdflag = k;
	  return;
	}
	its++;
	x = w[l - 1];
	y = w[k1 - 1];
	g = rv1[k1 - 1];
	h = rv1[k - 1];
	f = ((y - z) * (y + z) + (g - h) * (g + h)) / 2 / h / y;
	g = sqrt(f * f + 1);
	f = ((x - z) * (x + z) + h * (y / (f + fabs(g) * _Sgn(f)) - h)) / x;
	c = 1.0;
	s = 1.0;

	/* next QR transformation ... */
	for (i1 = l; i1 <= k1; i1++) {
	  i = i1 + 1;
	  g = rv1[i - 1];
	  y = w[i - 1];
	  h = s * g;
	  g = c * g;
	  z = sqrt(f * f + h * h);
	  rv1[i1 - 1] = z;
	  c = f / z;
	  s = h / z;
	  f = x * c + g * s;
	  g = g * c - x * s;
	  h = y * s;
	  y *= c;

	  for (j = 1; j <= n; j++) {
	    x = v[_index_(j, i1, &V) - 1];
	    z = v[_index_(j, i, &V) - 1];
	    v[_index_(j, i1, &V) - 1] = x * c + z * s;
	    v[_index_(j, i, &V) - 1] = z * c - x * s;
	  }

	  z = sqrt(f * f + h * h);
	  w[i1 - 1] = z;

	  /* rotation arbitrary if z:=0 ... */
	  if (z != 0) {
	    c = f / z;
	    s = h / z;
	  }
	  f = c * g + s * y;
	  x = c * y - s * g;
	  for (j = 1; j <= m; j++) {
	    y = u[_index_(j, i1, &V) - 1];
	    z = u[_index_(j, i, &V) - 1];
	    u[_index_(j, i1, &V) - 1] = y * c + z * s;
	    u[_index_(j, i, &V) - 1] = z * c - y * s;
	  }
	}


	rv1[l - 1] = 0.0;
	rv1[k - 1] = f;
	w[k - 1] = x;

      } else if (z < 0) {
	w[k - 1] = -z;   /* make w[k] >:=0 ... */
	for (j = 1; j <= n; j++)
	  v[_index_(j, k, &V) - 1] = -v[_index_(j, k, &V) - 1];
      }



    } while (again);
  }

  /* no errors */
  Errorflag = false;



}  /*svd*/


static void
_CTpower(int n, int n0, double acc, double *a, double *eigvec, double *eigval)
{
  int i, j, iter;
  double s, ss, s2, err;
  double d[MAXN];
  double *av[MAXN];

  _size_check(n, n0);
  if (Errorflag)
    return;
  if (acc <= 0.0) {
    _errormatrix(7);
    return;
  }

  for (i = 0; i < n; i++)
    av[i] = (void *)((long)a + i * n0 * sizeof(double));
  iter = 0;
  do {
    s2 = 0.0;
    ss = 0.0;
    for (i = 0; i < n; i++) {
      s = 0.0;
      for (j = 0; j < n; j++)
	s += av[i][j] * eigvec[j];
      d[i] = s;
      ss += eigvec[i] * s;
      s2 += s * s;
    }
    *eigval = ss;
    /* Rayleigh coefficient */
    for (i = 0; i < n; i++)
      d[i] /= sqrt(s2);
    err = 0.0;
    for (i = 0; i < n; i++)
      err += fabs(fabs(d[i]) - fabs(eigvec[i]));
    for (i = 0; i < n; i++)
      eigvec[i] = d[i];
    iter++;
  } while (iter <= MAXITER && err >= acc);

  /* no errors */
  if (err < acc)
    Errorflag = false;
  else  /*'Accuracy not achieved after n iteration'*/
    _errormatrix(2);
}

static void
_CTmodpower(int n, int n0, double acc, double *a, double *eigvec, double *eigval)
{
  int i, iter;
  double s, ss, s2, err, det, cond, lambda;
  double d[MAXN];
  double *av[MAXN];
  double *dptr;

  /*2.41*/
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/
  _size_check(n, n0);
  if (Errorflag)
    return;
  if (acc <= 0.0) {
    _errormatrix(7);
    return;
  }

  for (i = 0; i < n; i++)
    av[i] = (void *)((long)a + i * n0 * sizeof(double));
  iter = 0;
  dptr = d;
  lambda = *eigval;

  /* form A - lambda I */
  for (i = 0; i < n; i++)
    av[i][i] -= *eigval;

  /* decompose matrix */
  _CTdecomp(n, n0, a, &det, &cond);

  if (cond > 1e8) {  /*'Unsuitable initial value for eigenvalue'; */
    _errormatrix(3);
    return;
  }

  for (i = 0; i < n; i++)
    d[i] = eigvec[i];

  do {

    /* solve (A-LI)*X2=X1 */
    _CTsolve(n, n0, a, dptr);

    s2 = 0.0;
    ss = 0.0;
    for (i = 0; i < n; i++) {
      s = d[i];
      ss += eigvec[i] * s;
      s2 += s * s;
    }
    *eigval = ss;

    /* Rayleigh coefficient */
    for (i = 0; i < n; i++)
      d[i] /= sqrt(s2);
    err = 0.0;
    for (i = 0; i < n; i++)
      err += fabs(fabs(d[i]) - fabs(eigvec[i]));
    for (i = 0; i < n; i++)
      eigvec[i] = d[i];
    iter++;
  } while (iter <= MAXITER && err >= acc);

  /* no errors */
  if (err < acc && *eigval != 0.0) {
    Errorflag = false;
    *eigval = lambda + 1.0 / *eigval;
  } else  /*accuracy not achieved*/
    _errormatrix(4);
}


/* Local variables for CTjacobi: */
struct LOC_CTjacobi {
  int n;
  boolean vectors;
  double *eigval;
  double sq, sqend, ulim;
  int km;
  double u[WORKN];
  double *av[MAXN];
} ;

/* Local variables for jcrot: */
struct LOC_jcrot {
  struct LOC_CTjacobi *LINK;
  int k1;
  double cs, sn;
  int k2;
} ;

static void _jcdo(struct LOC_jcrot *LINK)
{
  double u1, u2;

  u1 = LINK->LINK->u[LINK->k1 - 1];
  u2 = LINK->LINK->u[LINK->k2 - 1];
  LINK->LINK->u[LINK->k1 - 1] = LINK->cs * u1 - LINK->sn * u2;
  LINK->LINK->u[LINK->k2 - 1] = LINK->sn * u1 + LINK->cs * u2;
}

static void _jcrot(int i1, int j1, int k1_, struct LOC_CTjacobi *LINK)
{
  /* one jacobi rotation */
  struct LOC_jcrot V;
  double th, di, dj, ai, aj, uk;
  int i, FORLIM;

  /*jcdo*/

  V.LINK = LINK;
  V.k1 = k1_;
  /* first find angle ... */
  di = LINK->eigval[i1 - 1];
  dj = LINK->eigval[j1 - 1];
  uk = LINK->u[V.k1 - 1];
  th = M_PI / 4;
  if (di != dj)
    th = atan(2 * uk / (dj - di)) / 2;
  V.cs = cos(th);
  V.sn = sin(th);

  /* do diag and pivot els ... */
  LINK->sq -= uk * uk;
  LINK->u[V.k1 - 1] = 0.0;
  LINK->ulim = sqrt(LINK->sq / LINK->km) / 2;
  uk = 2 * V.cs * V.sn * uk;
  LINK->eigval[i1 - 1] = di * V.cs * V.cs - uk + dj * V.sn * V.sn;
  LINK->eigval[j1 - 1] = di * V.sn * V.sn + uk + dj * V.cs * V.cs;

  /* do remaining els, start at top of cols ... */
  i = 1;
  V.k1 = i1 - 1;
  V.k2 = j1 - 1;

  /* while i<i1 do columns ... */
  while (i < i1) {
    _jcdo(&V);
    i++;
    V.k1 += LINK->n - i;
    V.k2 += LINK->n - i;
  }

  /* step over first diag element ... */
  V.k1++;
  i++;
  V.k2 += LINK->n - i;
  while (i < j1) {
    _jcdo(&V);
    V.k1++;
    i++;
    V.k2 += LINK->n - i;
  }

  /* step over second diag element ... */
  V.k1++;
  i++;
  V.k2++;

  /* while i<=n, do rows ... */
  while (i <= LINK->n) {
    _jcdo(&V);
    V.k1++;
    V.k2++;
    i++;
  }

  /* update eigenvectors ... */
  if (!LINK->vectors)
    return;
  FORLIM = LINK->n;
  for (i = 0; i < FORLIM; i++) {
    ai = LINK->av[i1 - 1][i];
    aj = LINK->av[j1 - 1][i];
    LINK->av[i1 - 1][i] = V.cs * ai - V.sn * aj;
    LINK->av[j1 - 1][i] = V.sn * ai + V.cs * aj;
  }

}  /*jcrot*/

static void _jcsq(struct LOC_CTjacobi *LINK)
{
  double a;
  int i, FORLIM;

  LINK->sq = 0.0;
  FORLIM = LINK->km;
  for (i = 0; i < FORLIM; i++)
    LINK->sq += LINK->u[i] * LINK->u[i];
  a = 0.0;
  FORLIM = LINK->n;
  for (i = 0; i < FORLIM; i++)
    a += fabs(LINK->eigval[i]);
  LINK->sqend = LINK->n * EPSILON * a;
}


/* of ModPower */


static void
_CTjacobi(int n_, int n0, boolean vectors_, double *a, double *eigval_)
{
  struct LOC_CTjacobi V;
  int i, j, k, FORLIM, FORLIM1;

  V.n = n_;
  V.vectors = vectors_;
  V.eigval = eigval_;
  /*2.41*/
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/
  _size_check(V.n, n0);
  if (Errorflag)
    return;

  FORLIM = V.n;
  for (i = 0; i < FORLIM; i++)
    V.av[i] = (void *)((long)a + i * n0 * sizeof(double));

  FORLIM = V.n;
  /* eigval initialised ... */
  for (i = 0; i < FORLIM; i++)
    V.eigval[i] = V.av[i][i];

  /* u[k] initialised ... */
  k = 0;
  FORLIM = V.n;
  for (i = 1; i < FORLIM; i++) {
    FORLIM1 = V.n;
    for (j = i; j < FORLIM1; j++) {
      k++;
      V.u[k - 1] = V.av[i - 1][j];
    }
  }

  V.km = k;

  /* a[i,j] initialised ... */
  if (V.vectors) {
    FORLIM = V.n;
    for (i = 0; i < FORLIM; i++) {
      FORLIM1 = V.n;
      for (j = 0; j < FORLIM1; j++)
	V.av[i][j] = 0.0;
      V.av[i][i] = 1.0;
    }
  }

  if (V.n <= 1)
    return;

  /* sum squares of off-diag elements & find thresholds ... */
  _jcsq(&V);
  V.sqend = sqrt(EPSILON) * V.sq;
  V.ulim = sqrt(V.sq / V.km) / 2;
  if (V.sq == 0)
    return;

  /* main loop ... */
  i = 1;
  j = 1;
  k = 0;
  do {
    j++;
    if (j > V.n) {
      i++;
      j = i + 1;
      if (i >= V.n) {
	i = 1;
	j = 2;
      }
    }

    k = k % V.km + 1;
/* p2c: nctmatrix.pas, line 1215:
 * Note: Using % for possibly-negative arguments [317] */

    /* test els over thresholds ... */
    if (fabs(V.u[k - 1]) > V.ulim)
      _jcrot(i, j, k, &V);

    /* test sum squares for reset ... */
    _jcsq(&V);
  } while (V.sq >= V.sqend);

  /* no errors */
  Errorflag = false;

}  /*jacobi*/


typedef double *arrayptrarray[1000];

/*intarray = array[1..1000] of integer;
intarrayptr = ^intarray;
*/



static void
_CTband(int n, int n0, int lbw, int ubw, double *a, double *bin, double *x, double *det)
{
  /*order of the matrix (=number of rows in the compact form) */
  /*second dimension of the passed matrix*/
  /*lower and upper bandwidths (there are n0>=lbw+ubw+1 cols in the comp. form)*/
  /*as in 'ax=b'*/

  /*the determinant*/
  double aa, d;
  int i, j, k, m;
  double **av;
  double *b;
  double temp;



  /*--- Ddnh getmem -> newptr*/
  /*--- D   getmem(av, sizeof(pointer)*n);*/
  /*--- D   getmem(b, sizeof(ctreal)*n);*/
  /*--- D*/
  av = (double **)_NewPtr(sizeof(void *) * n);
  b = (double *)_NewPtr(sizeof(double) * n);

  for (i = 0; i < n; i++) {
    /*--- Ddnh      getmem(av^[i],sizeof(ctreal)*(1+2*lbw+ubw));*/
    av[i] = (double *)_NewPtr(sizeof(double) * (lbw * 2 + ubw + 1));
    for (j = 0; j <= lbw + ubw; j++)
      av[i][j] = a[i * n0 + j];
    for (j = lbw + ubw + 1; j <= lbw * 2 + ubw; j++)
      av[i][j] = 0.0;
    b[i] = bin[i];
  }


  *det = 1.0;
  /* main elimination loop ... */
  k = 0;
  do {
    k++;
    aa = 0.0;
    /* look for pivot ... */
    m = k;
    while (m <= n && m <= k + lbw) {
      if (fabs(av[m - 1][lbw + k - m]) > aa) {
	j = m;
	aa = fabs(av[m - 1][lbw + k - m]);
      }
      m++;
    }
    if (aa == 0)
      *det = 0.0;   /* det:=0  if matrix singular */
    else {
      /* swap row indices to select pivot ...*/
      i = k;
      while (i <= n && i <= k + ubw + lbw) {
	temp = av[k - 1][lbw + i - k];
	av[k - 1][lbw + i - k] = av[j - 1][lbw + i - j];
	av[j - 1][lbw + i - j] = temp;
	i++;
      }
      temp = b[k - 1];
      b[k - 1] = b[j - 1];
      b[j - 1] = temp;
      *det = -*det;
    }
    /* update det ... */
    d = av[k - 1][lbw];
    *det *= d;
    /* eliminate terms below pivot ... */
    m = k + 1;
    while (m <= k + lbw && m <= n) {
      /* find factor ... */
      aa = av[m - 1][lbw + k - m] / d;
      /* work along this row -  a[m,k]:=0  but need not store ... */
      j = k + 1;
      while (j <= k + ubw + lbw && j <= n) {
	av[m - 1][lbw + j - m] -= aa * av[k - 1][lbw + j - k];
	j++;
      }
      /* and same for RHS vector ... */
      b[m - 1] -= aa * b[k - 1];

      /*                 end;                  */

      m++;
    }

  } while (k != n - 1 && *det != 0);
  /* back substitute ... */
  aa = av[n - 1][lbw];
  *det = aa * *det;
  if (*det != 0) {
    x[n - 1] = b[n - 1] / aa;
    for (k = n - 1; k >= 1; k--) {
      aa = 0.0;
      j = k + 1;
      while (j <= n && j <= k + lbw + ubw) {
	aa += av[k - 1][lbw + j - k] * x[j - 1];
	j++;
      }
      x[k - 1] = (b[k - 1] - aa) / av[k - 1][lbw];
    }
  }


  for (i = 0; i < n; i++)   /*dnh , sizeof(ctreal)*(1+2*lbw+ubw))*/
    Free(av[i]);
  /*dnh ,sizeof(pointer)*n)*/
  Free(av);   /*dnh , sizeof(ctreal)*n)*/
  Free(b);

}


typedef double *arrayptrarray_[1000];

/*intarray = array[1..1000] of integer;
intarrayptr = ^intarray;
*/


/*NEW*/

/* this routine should work with larger matrices (7500x7500 in theory is OK)
(but we still pass an arrayptr to a normal pascal-style 2 dimensional array
  to the routine) */

static void
_CTelim2(int n, int n0, double *a, double *bin, double *x, double *det)
{
  double aa, d;
  int i, j, k, m, p;
  double **av;
  int *pivot;
  double *b;

  /*2.41*/
  if (n < 1 || n > n0) {
    _CTerror("array dimension error");
    return;
  }
  /*--- D   if Errorflag then*/
  /*--- D      Exit;*/
  /*dnh - why didn't he use newptr in the original? */
  av = (double **)_NewPtr(sizeof(void *) * n);
  /*--- D    getmem(av, sizeof(pointer)*n);*/
  b = (double *)_NewPtr(sizeof(double) * n);
  /*--- D    getmem(b, sizeof(ctreal)*n);*/
  for (i = 0; i < n; i++) {
    av[i] = (double *)_NewPtr(sizeof(double) * n);
    /*--- D    getmem(av^[i], sizeof(ctreal)*n);*/
    for (j = 0; j < n; j++)
      av[i][j] = a[n0 * i + j];
    b[i] = bin[i];
  }
  pivot = (int *)_NewPtr(sizeof(double) * n);
  /*--- D    getmem(pivot, sizeof(integer)*n);*/

  /* now down to business ...
  The remained of this routine is a translation of 'elim', found at the start of the unit.  */
  *det = 1.0;
  /* initial order pivot vector ... */
  for (i = 1; i <= n; i++)
    pivot[i - 1] = i;
  /* main elimination loop ... */
  k = 0;
  do {
    k++;
    aa = 0.0;
    /* look for pivot ... */
    for (m = k; m <= n; m++) {
      i = pivot[m - 1];
      if (fabs(av[i - 1][k - 1]) > aa) {
	j = m;
	aa = fabs(av[i - 1][k - 1]);
      }
    }
    if (aa == 0)
      *det = 0.0;   /* det:=0  if matrix singular */
    else {
      /* swap row indices to select pivot ...*/
      m = pivot[k - 1];
      p = pivot[j - 1];
      if (k != j) {
	pivot[k - 1] = p;
	pivot[j - 1] = m;
	*det = -*det;
      }
      /* update det ... */
      d = av[p - 1][k - 1];
      *det *= d;
      /* eliminate terms below pivot ... */
      for (m = k; m < n; m++) {
	i = pivot[m];
	/* find factor ... */
	aa = av[i - 1][k - 1] / d;
	/* work along this row -  a[i,k]:=0  but need not store ... */
	for (j = k; j < n; j++)
	  av[i - 1][j] -= aa * av[p - 1][j];
	/* and same for RHS vector ... */
	b[i - 1] -= aa * b[p - 1];
      }
    }
  } while (k != n - 1 && *det != 0);
  /* back substitute ... */
  aa = av[pivot[n - 1] - 1][n - 1];
  *det = aa * *det;
  if (*det != 0) {
    x[n - 1] = b[pivot[n - 1] - 1] / aa;
    for (k = n - 1; k >= 1; k--) {
      i = pivot[k - 1];
      aa = 0.0;
      for (j = k; j < n; j++)
	aa += av[i - 1][j] * x[j];
      x[k - 1] = (b[i - 1] - aa) / av[i - 1][k - 1];
    }
  }

  /*Now get rid of the memory*/

  for (i = 0; i < n; i++)   /*DisposPtr(mptr(av^[i])); */
    Free(av[i]);
  /*dnh , sizeof(ctreal)*n)*/
  /* DisposPtr(mptr(av)); */
  Free(av);   /*dnh ,sizeof(pointer)*n)*/
  /*DisposPtr(mptr(pivot));*/
  Free(pivot);   /*dnh ,sizeof(integer)*n)*/
  /*    DisposPtr(mptr(b));   */
  /*dnh , sizeof(ctreal)*n)*/
  Free(b);

  /*2.41 does it like this*/
  /*--- D   for i:=1 to n do _DisposePtr(mptr(av^[i]), sizeof(ctreal)*n);*/
  /*--- D   Freemem(av^[i], sizeof(ctreal)*n);*/
  /*--- D   _DisposePtr(mptr(av), sizeof(pointer)*n);*/
  /*--- D   Freemem(av,sizeof(pointer)*n);*/
  /*--- D   _DisposePtr(mptr(pivot), sizeof(integer)*n);*/
  /*--- D   Freemem(pivot,sizeof(integer)*n);*/
  /*--- D   _DisposePtr(mptr(b), sizeof(ctreal)*n);*/
  /*--- D   Freemem(b, sizeof(ctreal)*n);*/
}


/* numore.c inserted 2000-05-20 */


/*  CATAM Software Library */
/*  Copyright  DAMTP  1994, 1998, 1999, 2000 */
/*  Department of Applied Mathematics and Theoretical Physics */
/*  University of Cambridge */
/*  C conversion 2.41 -> C2.0 - dnh - damtp - Spring 2000 */


/*
--------------------------------------------------------------------
  cclibn.c non-windows code
  code for Turbo Pascal library functions used by CATSL
  new functions for timing, plotting, 
--------------------------------------------------------------------
*/

long _Timer(void)
{
/*    time of day version (not mingw32) */
/*     struct timeval t; */
/*     struct timezone z; */
/*     gettimeofday(&t, &z); */
/*     timer1 = 100L * (long)t.tv_sec + t.tv_usec / 1000L; */

  /* clock version */
  return 100L * (long)clock() / (long)CLOCKS_PER_SEC;
}

char *_NewPtr(long nn)
{
  return (char*) malloc(nn);
}

void _DisposePtr(char *p, long nn)
{
/*    nn is not needed in C */
  free(p);
}

/******************************************************************/
/*                           new routines                         */
/******************************************************************/

static long elapsedtime;

/* reset timer */
void TickCL(void)
{
  elapsedtime = _Timer();
}

/* elapsed time in 100th of sec */
long TockCL(void)

{
  return _Timer() - elapsedtime;
}

/* Common string functions: */

/* from p2clib */

/* Delete the substring of length "len" at index "pos" from "s".
   Delete less if out-of-range. */


void _strdelete(char *s, int pos, int len)
{
    register int slen;

    if (--pos < 0)
        return;
    slen = strlen(s) - pos;
    if (slen <= 0)
        return;
    s += pos;
    if (slen <= len) {
        *s = 0;
        return;
    }
    while ((*s = s[len])) s++;
}


/* Insert string "src" at index "pos" of "dst". */

void _strinsert(char *src, char *dst, int pos)
{
    register int slen, dlen;

    if (--pos < 0)
        return;
    dlen = strlen(dst);
    dst += dlen;
    dlen -= pos;
    if (dlen <= 0) {
        strcpy(dst, src);
        return;
    }
    slen = strlen(src);
    do {
        dst[slen] = *dst;
        --dst;
    } while (--dlen >= 0);
    dst++;
    while (--slen >= 0)
        *dst++ = *src++;
}


/* Store in "ret" the substring of length "len" starting from "pos" (1-based).
   Store a shorter or null string if out-of-range.  Return "ret". */

char* _strsub(char *ret, char *s, int pos, int len)
{
  char *s2;
  
  if (--pos < 0 || len <= 0) {
    *ret = 0;
    return ret;
  }
  while (pos > 0) {
    if (!*s++) {
      *ret = 0;
      return ret;
    }
    pos--;
  }
  s2 = ret;
  while (--len >= 0) {
    if (!(*s2++ = *s++))
      return ret;
  }
  *s2 = 0;
  return ret;
}


/* Return the index of the first occurrence of "pat" as a substring of "s",
   starting at index "pos" (1-based).  Result is 1-based, 0 if not found. */

int _strpos2(char *s, char *pat, int pos)
{
    register char *cp, ch;
    register int slen;

    if (--pos < 0)
        return 0;
    slen = strlen(s) - pos;
    cp = s + pos;
    if (!(ch = *pat++))
        return 0;
    pos = strlen(pat);
    slen -= pos;
    while (--slen >= 0) {
        if (*cp++ == ch && !strncmp(cp, pat, pos))
            return cp - s;
    }
    return 0;
}



#define BAD_RND_ARG     22
static int seed1 = 13728;
static int seed2 =  8710;
static int seed3 = 14906;

double _rnd(double arg)
{
  double res;
  int rarg;

  /* of rnd */
  seed1 = seed1 % 177 * 171 - seed1 / 177 * 2;
/* p2c: wctcntrl.pas, line xxx:
 * Note: Using % for possibly-negative arguments [317] */
  if (seed1 < 0)
    seed1 += 30269;

  seed2 = seed2 % 176 * 172 - seed2 / 176 * 35;
/* p2c: wctcntrl.pas, line xxx:
 * Note: Using % for possibly-negative arguments [317] */
  if (seed2 < 0)
    seed2 += 30307;

  seed3 = seed3 % 178 * 170 - seed3 / 178 * 63;
/* p2c: wctcntrl.pas, line xxx:
 * Note: Using % for possibly-negative arguments [317] */
  if (seed3 < 0)
    seed3 += 30323;

  res = seed1 / 30269.0 + seed2 / 30307.0 + seed3 / 30323.0;

  res -= (long)res;

  if (arg >= 0.5 && arg <= 32766.0) {
    rarg = (long)floor(arg + 0.5);

    if (rarg == 1)
      return res;
    else
      return ((long)(rarg * res) + 1.0);
  }
#ifdef TODO
  /*  what to do about this - it's part of the evaluator */
  evaluate_error(BAD_RND_ARG);
#endif
}

void SetRandomCL(int seed)
{
  time_t t;
  int s;
  
  if (seed <= 0) {
    time(&t);
    s = abs((int) t);
  } else s = seed;
  seed1 = (s % 100 + 1) * 300;
  seed2 = (s % 6000 + 1) * 5;
  seed3 =  s % 30000 + 1;
}

double RandomCL(void) 
{
  return _rnd(1.);
}
int    RandomIntCL(int range)
{
  return _rnd((float) range);
}


double ZeroCL(double (*F)(double x), double a, double b, double tol)
{
  /*  	reorder parms Zeroin */
  return _Zeroin(a, b, tol, F);
}

double MinCL(double (*F)(double x), double a, double b, double tol)
{
   /* reorder parms Fmin */
  return _Fmin(a, b, tol, F);
}

double RombergCL(double (*F)(double x), double a, double b, double *error, int n)
{
  /* reorder F and n parms */
  return _Romb(a, b, n, F, error);
}

double QuadCL(double (*F)(double x), double a, double b, double aber,
	      double reler, double *errest, double *flag)
{
  /* F position, nf removed */
  int nf;
  return _Quanc8(a, b, aber, reler, F, errest, flag, &nf);
}

double GaussElimCL(int n, int n0, two_d_array a, double z[], double x[])
{
  /* return determinant */
  double det;
  _CTelim2(n, n0, (double *)a, z, x, &det);
  return det;
}

double DecomposeCL(int n, int n0, two_d_array a)
{
  /* return det, remove condition argument */
  double det, cond;
  _CTdecomp(n, n0, (double *)a, &det, &cond);
  return det;
}

double InvertCL(int n, int n0, two_d_array a)
{
  /* return det, remove condition argument */
  double det, cond;
  _CTinvert(n, n0, (double *)a, &det, &cond);
  return det;
}

int SvdCL(int m, int n, int n0, two_d_array a, two_d_array u, 
	  two_d_array v, double w[])
{
  /* return flag */
  int flag;
  _CTsvd(m, n, n0, (double *)a, (double *)u, (double *)v, w, &flag);
  return flag;
}

double EigenMaxCL(int n, int n0, double acc, two_d_array a, double *eigvec)
{
  /* return eigenvalue */
  double eigval;
  _CTpower(n, n0, acc, (double *)a, eigvec, &eigval);
  return eigval;
}

double EigenValCL(int n, int n0, double acc, two_d_array a, double eigvec[], double theta)
{
  /* return eigenvalue */
  double eigval = theta;
  _CTmodpower(n, n0, acc, (double *)a, eigvec, &eigval);
  return eigval;
}

void JacobiCL(int n, int n0, two_d_array a, double eigval[])
{
  /* remove boolean argument */
  _CTjacobi(n, n0, true, (double *)a, eigval);
}

boolean RkfCL(int n, double aber, double reler, double *t, double *dt, double dtmin, void (*F)(double t, double *y, double *ydash), double *y, double *yp, int *nleft) {

  boolean flag; /* return as function value */
  _Rkf(n, aber, reler, t, dt, dtmin, F, y, yp, nleft, &flag);
  return flag;
}

double BandCL(int n, int n0, int lbw, int ubw, two_d_array a, double z[], double x[])
{
  double det;
  _CTband(n, n0, lbw, ubw, (double *)a, z, x, &det);
}

void PoissonCL(two_d_array psp, two_d_array ztp, int nnx, int nny, double dlx, double dly)
{
  _Poisson2d((double *)psp, (double *)ztp, nnx, nny, dlx, dly);
}

void SolveCL(int n, int n0, two_d_array a, double z[])
{
  _CTsolve(n, n0, (double *)a, z);
}

void   
FftCL(double f[][2], int m, double sign)
{
Fft((double *)f, m, sign);
}
