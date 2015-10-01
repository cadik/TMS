/* Radiance file format definitions
*/

#include  "RadianceC.h"
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <time.h>
#include  <ctype.h> 
 
#define  RED		0
#define  GRN		1
#define  BLU		2 
#define  EXP		3
#define  MINELEN	8	/* minimum scanline length for encoding */
#define  MAXELEN	0x7fff	/* maximum scanline length for encoding */
#define  MINRUN		4	/* minimum run length */ 
#define  copycolr(c1,c2)	(c1[0]=c2[0],c1[1]=c2[1],c1[2]=c2[2],c1[3]=c2[3]) 
#define	 MAXLINE	512 
#define  COLRFMT		"32-bit_rle_rgbe" 
#define  YMAJOR			4 
#define  XDECR			1
#define  YDECR			2 
#define  COLXS		128	/* excess used for exponent */ 

typedef unsigned char COLR[4];		/* red, green, blue (or X,Y,Z), exponent */
typedef double  COLOR[3];	/* red, green, blue (or X,Y,Z) */
 
class TMORadiance  
{
public:
	struct check 
	{
		FILE	*fp;
		char	fs[64];
		char	ex[64];
	}; 

	struct RESOLU 
	{
		int	rt;		/* orientation (from flags above) */
		int	xr, yr;		/* x and y resolution */
	}; 

	static int str2resolu(RESOLU  *rp, char  *buf);
	static int freadcolrs(COLR *scanline, int len, FILE *fp); 
	static int fwritecolrs(COLR *scanline, int len, FILE *fp);
	static int fgetresolu(int *sl, int *ns, FILE *fp);
	static int formatval(char  *r, char  *s);
	static int exposureval(char  *r, char  *s);
	static int mycheck(char *s, struct check  *cp);
	static int globmatch(char	*p, char *s);
	static int getheader(FILE  *fp, int  (*f)(char*, struct TMORadiance::check*), struct check *p);
	static int checkheader(FILE  *fin,char  * fmt, FILE  *fout);
	static int oldreadcolrs(COLR  *scanline, int  len, register FILE  *fp);
	static void colr_color(COLOR col,COLR clr);
	static void setcolr(COLR clr, double r, double g, double b);

	static double dExposure;

	TMORadiance();
	virtual ~TMORadiance();

	static char FMTSTR[];
	static char EXPSTR[];
};

