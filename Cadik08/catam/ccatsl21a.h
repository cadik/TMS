#ifndef CATNUM_H
#define CATNUM_H

#ifndef  true
# define true    1
# define false   0
#endif
#ifndef  TRUE
# define TRUE    1		/* alternate upper-case spelling */
# define FALSE   0
#endif

typedef unsigned char   uchar;

# define _FNSIZE  120

/* File mode magic numbers */

#define fmClosed        0xd7b0L
#define fmInput         0xd7b1L
#define fmOutput        0xd7b2L
#define fmInOut         0xd7b3L


/* File attribute constants */

#define faReadOnly      0x1
#define faHidden        0x2
#define faSysFile       0x4
#define faVolumeID      0x8
#define faDirectory     0x10
#define faArchive       0x20
#define faAnyFile       0x3f


#define FileNotFound     10
#define LINESTRINGSIZE  80
typedef char spring[256];
typedef char linespring[LINESTRINGSIZE + 1];
typedef char Plinespring[LINESTRINGSIZE + 1];
typedef char Pspring[256];
typedef double CTreal;
typedef double citreal;

#define MAXREAL         1.7e+308
#define MINREAL         5.0e-324
#define EPSILON         1.2e-16

#define RTN             13
#define SPC             32
#define BS              8
#define ESC             27
#define BRK             3
#define BELL            7
#define ZERO            48
#define NINE            57
#define SPACE           ((char)SPC)
#define Blanks          "                                                                                "
typedef int Intarray[32000];
typedef long Larray[16000];
typedef float Rarray[10000];
typedef citreal Datarray[8000];
typedef char Linestring[LINESTRINGSIZE + 1];
typedef struct Textobj {
  char *str;
  struct Textobj *nxt;
} Textobj;
typedef void (*fn0)(void);
typedef double (*fn1)(double x);
typedef double (*fn2)(double x, double y);
typedef double (*fn3)(double x, double y, double z);
typedef double (*fn_1_real)(double x);
typedef double (*fn_2_real)(double x, double y);
typedef double (*fn_3_real)(double x, double y, double z);
typedef void (*fn_ode)(double t, double y[], double ydash[]);

#define FreeR(p)    (free((void *)(p)))    /* used if arg is an rvalue */
#define Free(p)     (free((void *)(p)), (p)=NULL)

/* this allows the ctcntrl startup code to be called from ctbase/ctinit2 */
#define Setup_dictionary Setup_dictionary_
#define Setup_random Setup_random_

/* external variables */
#undef  EVAR
#ifdef  NUM_C
#define EVAR
#else
#define EVAR extern
#endif

EVAR boolean Errorflag;
/* cancelflag, wmfflag -> windata.h */
EVAR int Errorcode;
static spring Errorstring;
#ifdef NUM_C
char *ErrorMessageCD = Errorstring;
#else
extern char *ErrorMessageCD;
#endif

_EXTERN_ char *_NewPtr(long nn);
_EXTERN_ double *_New_array(long nn);
_EXTERN_ double _Arcsin(double x);
_EXTERN_ void _CTerror(char *s);
_EXTERN_ void _DisposePtr(char *p, long nn);
_EXTERN_ double _Float(int i);
_EXTERN_ void _Free_array(double **p, long nn);
_EXTERN_ int _Imax(int a, int b);
_EXTERN_ int _Imin(int a, int b);
_EXTERN_ double _Max(double a, double b);
_EXTERN_ double _Min(double a, double b);
_EXTERN_ double _Npd(double x);
_EXTERN_ double _Nraise(double value, int np);
_EXTERN_ double _Raise(double value, double power);
_EXTERN_ void _Set_error(int errcode, char *errstr);
_EXTERN_ void _Setup_random_(void);
_EXTERN_ int _Sgn(double x);
_EXTERN_ double _Tan(double x);
_EXTERN_ long _Timer(void);
_EXTERN_ double _rnd(double arg);
_EXTERN_ void _strdelete(char *s, int pos, int len);
_EXTERN_ void _strinsert(char *src, char *dst, int pos);
_EXTERN_ int _strpos2(char *s, char *pat, int pos);
_EXTERN_ char* _strsub(char *ret, char *s, int pos, int len);

#endif /*CATNUM_H*/
