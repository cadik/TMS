/*
 *              This product includes Radiance software
 *              (http://radsite.lbl.gov/)
 *              developed by the Lawrence Berkeley National Laboratory
 *              (http://www.lbl.gov/).
 */

#include "TMORadiance.h"

char  TMORadiance::FMTSTR[] = "FORMAT=";	/* format identifier */
char  TMORadiance::EXPSTR[] = "EXPOSURE=";	/* exposure identifier */
double TMORadiance::dExposure = 1.0;

void TMORadiance::colr_color(COLOR col,COLR clr)		/* convert short to float color */
{
	double  f;
	
	if (clr[EXP] == 0)
		col[RED] = col[GRN] = col[BLU] = 0.0;
	else {
		f = ldexp(1.0, (int)clr[EXP]-(COLXS+8));
		col[RED] = (clr[RED] + 0.5)*f;
		col[GRN] = (clr[GRN] + 0.5)*f;
		col[BLU] = (clr[BLU] + 0.5)*f;
	}
} 

void TMORadiance::setcolr(COLR clr, double r, double g, double b)		/* assign a short color value */
{
	double  d;
	int  e;
	
	d = r > g ? r : g;
	if (b > d) d = b;

	if (d <= 1e-32) {
		clr[RED] = clr[GRN] = clr[BLU] = 0;
		clr[EXP] = 0;
		return;
	}

	d = frexp(d, &e) * 255.9999 / d;

	if (r > 0.0)
		clr[RED] = r * d;
	else
		clr[RED] = 0;
	if (g > 0.0)
		clr[GRN] = g * d;
	else
		clr[GRN] = 0;
	if (b > 0.0)
		clr[BLU] = b * d;
	else
		clr[BLU] = 0;

	clr[EXP] = e + COLXS;
}

int TMORadiance::str2resolu(RESOLU  *rp, char  *buf)		/* convert resolution line to struct */
{
	register char  *xndx, *yndx;
	register char  *cp;

	if (buf == NULL)
		return(0);
	xndx = yndx = NULL;
	for (cp = buf; *cp; cp++)
		if (*cp == 'X')
			xndx = cp;
		else if (*cp == 'Y')
			yndx = cp;
	if (xndx == NULL || yndx == NULL)
		return(0);
	rp->rt = 0;
	if (xndx > yndx) rp->rt |= YMAJOR;
	if (xndx[-1] == '-') rp->rt |= XDECR;
	if (yndx[-1] == '-') rp->rt |= YDECR;
	if ((rp->xr = atoi(xndx+1)) <= 0)
		return(0);
	if ((rp->yr = atoi(yndx+1)) <= 0)
		return(0);
	return(1);
} 

int TMORadiance::fgetresolu(int *sl, int *ns, FILE *fp)			/* get picture dimensions */
{
	RESOLU  rs;
	char  resolu_buf[32]; 

	if (!str2resolu(&rs, fgets(resolu_buf,32,fp))) throw -6;
		
	if (rs.rt & YMAJOR) {
		*sl = rs.xr;
		*ns = rs.yr;
	} else {
		*sl = rs.yr;
		*ns = rs.xr;
	}
	return(rs.rt);
} 
int TMORadiance::oldreadcolrs(COLR  *scanline, int  len, register FILE  *fp)		/* read in an old colr scanline */
{
	int  rshift;
	register int  i;
	
	rshift = 0;
	
	while (len > 0) {
		scanline[0][RED] = getc(fp);
		scanline[0][GRN] = getc(fp);
		scanline[0][BLU] = getc(fp);
		scanline[0][EXP] = getc(fp);
		if (feof(fp) || ferror(fp))
			return(-1);
		if (scanline[0][RED] == 1 &&
				scanline[0][GRN] == 1 &&
				scanline[0][BLU] == 1) {
			for (i = scanline[0][EXP] << rshift; i > 0; i--) {
				copycolr(scanline[0], scanline[-1]);
				scanline++;
				len--;
			}
			rshift += 8;
		} else {
			scanline++;
			len--;
			rshift = 0;
		}
	}
	return(0);
}
 

int TMORadiance::freadcolrs(COLR *scanline, int  len, FILE  *fp)		/* read in an encoded colr scanline */
{
	register int  i, j;
	int  code, val;
					/* determine scanline type */
	if ((len < MINELEN) | (len > MAXELEN))
		return(oldreadcolrs(scanline, len, fp));
	if ((i = getc(fp)) == EOF)
		return(-1);
	if (i != 2) {
		ungetc(i, fp);
		return(oldreadcolrs(scanline, len, fp));
	}
	scanline[0][GRN] = getc(fp);
	scanline[0][BLU] = getc(fp);
	if ((i = getc(fp)) == EOF)
		return(-1);
	if (scanline[0][GRN] != 2 || scanline[0][BLU] & 128) {
		scanline[0][RED] = 2;
		scanline[0][EXP] = i;
		return(oldreadcolrs(scanline+1, len-1, fp));
	}
	if ((scanline[0][BLU]<<8 | i) != len)
		return(-1);		/* length mismatch! */
					/* read each component */
	for (i = 0; i < 4; i++)
	    for (j = 0; j < len; ) {
		if ((code = getc(fp)) == EOF)
		    return(-1);
		if (code > 128) {	/* run */
		    code &= 127;
		    if ((val = getc(fp)) == EOF)
			return -1;
		    while (code--)
			scanline[j++][i] = val;
		} else			/* non-run */
		    while (code--) {
			if ((val = getc(fp)) == EOF)
			    return -1;
			scanline[j++][i] = val;
		    }
	    }
	return(0);
}

int TMORadiance::fwritecolrs(COLR  *scanline,  int len, FILE *fp)		/* write out a colr scanline */
{
	register int  i, j, beg, cnt = 1;
	int  c2;
	
	if ((len < MINELEN) | (len > MAXELEN))	/* OOBs, write out flat */
		return(fwrite((char *)scanline,sizeof(COLR),len,fp) - len);
					/* put magic header */
	putc(2, fp);
	putc(2, fp);
	putc(len>>8, fp);
	putc(len&255, fp);
					/* put components seperately */
	for (i = 0; i < 4; i++) {
	    for (j = 0; j < len; j += cnt) {	/* find next run */
		for (beg = j; beg < len; beg += cnt) {
		    for (cnt = 1; cnt < 127 && beg+cnt < len &&
			    scanline[beg+cnt][i] == scanline[beg][i]; cnt++)
			;
		    if (cnt >= MINRUN)
			break;			/* long enough */
		}
		if (beg-j > 1 && beg-j < MINRUN) {
		    c2 = j+1;
		    while (scanline[c2++][i] == scanline[j][i])
			if (c2 == beg) {	/* short run */
			    putc(128+beg-j, fp);
			    putc(scanline[j][i], fp);
			    j = beg;
			    break;
			}
		}
		while (j < beg) {		/* write out non-run */
		    if ((c2 = beg-j) > 128) c2 = 128;
		    putc(c2, fp);
		    while (c2--)
			putc(scanline[j++][i], fp);
		}
		if (cnt >= MINRUN) {		/* write out run */
		    putc(128+cnt, fp);
		    putc(scanline[beg][i], fp);
		} else
		    cnt = 0;
	    }
	}
	return(ferror(fp) ? -1 : 0);
}
 

int TMORadiance::formatval(char  *r, char  *s)			/* get format value (return true if format) */
{
	register char  *cp = FMTSTR;

	while (*cp) if (*cp++ != *s++) return(0);
	while (isspace(*s)) s++;
	if (!*s) return(0);
	if (r == NULL) return(1);
	do
		*r++ = *s++;
	while(*s && !isspace(*s));
	*r = '\0';
	return(1);
} 

int TMORadiance::exposureval(char  *r, char  *s)			/* get exposure value (return true if format) */
{
	register char  *cp = EXPSTR;

	while(*s && isspace(*s)) s++;
	while (*cp) if (*cp++ != *s++) return(0);
	while (isspace(*s)) s++;
	if (!*s) return(0);
	if (r == NULL) return(1);
	do
		*r++ = *s++;
	while(*s && !isspace(*s));
	*r = '\0';
	return(1);
} 

int TMORadiance::mycheck(char *s, struct check  *cp)			/* check a header line for format info. */
{
	if (!formatval(cp->fs, s))
	{
		if (exposureval(cp->ex, s)) dExposure = strtod(cp->ex, 0);
	}
	return(0);
} 

int TMORadiance::globmatch(char	*p, char *s)			/* check for match of s against pattern p */
{
	int	setmatch;

	do {
		switch (*p) {
		case '?':			/* match any character */
			if (!*s++)
				return(0);
			break;
		case '*':			/* match any string */
			while (p[1] == '*') p++;
			do
				if ( (p[1]=='?' || p[1]==*s) &&
						globmatch(p+1,s) )
					return(1);
			while (*s++);
			return(0);
		case '[':			/* character set */
			setmatch = *s == *++p;
			if (!*p)
				return(0);
			while (*++p != ']') {
				if (!*p)
					return(0);
				if (*p == '-') {
					setmatch += p[-1] <= *s && *s <= p[1];
					if (!*++p)
						break;
				} else
					setmatch += *p == *s;
			}
			if (!setmatch)
				return(0);
			s++;
			break;
		case '\\':			/* literal next */
			p++;
		/* fall through */
		default:			/* normal character */
			if (*p != *s)
				return(0);
			s++;
			break;
		}
	} while (*p++);
	return(1);
} 

int TMORadiance::getheader(FILE  *fp, int  (*f)(char*, struct TMORadiance::check*), struct TMORadiance::check *p)		/* get header from file */
{
	char  buf[MAXLINE];

	for ( ; ; ) {
		buf[MAXLINE-2] = '\n';
		if (fgets(buf, MAXLINE, fp) == NULL)
			return(-1);
		if (buf[0] == '\n')
			return(0);
#ifdef MSDOS
		if (buf[0] == '\r' && buf[1] == '\n')
			return(0);
#endif
		if (buf[MAXLINE-2] != '\n') {
			ungetc(buf[MAXLINE-2], fp);	/* prevent false end */
			buf[MAXLINE-2] = '\0';
		}
		if (f != NULL && (*f)(buf, p) < 0)
			return(-1);
	}
} 

int TMORadiance::checkheader(FILE  *fin,char  * fmt,FILE  *fout)
{
	struct check	cdat;
	register char	*cp;

	cdat.fp = fout;
	cdat.fs[0] = '\0';
	if (getheader(fin, mycheck, &cdat) < 0)
		return(-1);
	if (!cdat.fs[0])
		return(0);
	for (cp = fmt; *cp; cp++)		/* check for globbing */
		if ((*cp == '?') | (*cp == '*'))
			if (globmatch(fmt, cdat.fs)) {
				strcpy(fmt, cdat.fs);
				return(1);
			} else
				return(-1);
	return(strcmp(fmt, cdat.fs) ? -1 : 1);	/* literal match */
} 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TMORadiance::TMORadiance()
{

}

TMORadiance::~TMORadiance()
{

}
