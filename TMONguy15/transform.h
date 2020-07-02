/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* transform.h                   S. Paine rev. 2014 August 26
*
* Declarations for transform.c
************************************************************/

#ifndef AM_TRANSFORM_H
#define AM_TRANSFORM_H

void fft_dif(double*, unsigned long);
void ifft_dit(double*, unsigned long);
void fht_dif(double*, unsigned long);
void fht_dit(double*, unsigned long);
void hilbert(double*, unsigned long);
void bitrev_permute(double*, unsigned long);
void bitrev_permute_real(double*, unsigned long);
void ft_benchmarks(void);

#endif /* AM_TRANSFORM_H */
