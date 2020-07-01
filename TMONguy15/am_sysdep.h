/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* am_sysdep.h                   S. Paine rev. 2018 August 23
*
* Constants and declarations for am_sysdep.c
************************************************************/

#ifndef AM_AM_SYSDEP_H
#define AM_AM_SYSDEP_H

/*
 * For optimum performance, the default processor data cache size
 * settings here can be overridden by target-specific definitions
 * supplied at compile time.
 *
 * Line-by-line and CIA computations are blocked to fit in L1
 * cache.  If the cache size is set to 0, cache blocking is
 * turned off.
 *
 * For FFTs and FHTs, the L1 cache size setting controls the
 * point at which the computation switches over from recursive
 * to iterative.
 *
 * The L2 cache size is used for sizing internal benchmarks.
 */
#ifndef L1_CACHE_BYTES
    #define L1_CACHE_BYTES  0x8000
#endif

#ifndef L1_CACHE_WAYS
    #define L1_CACHE_WAYS   8
#endif

#ifndef L2_CACHE_BYTES
    #define L2_CACHE_BYTES  0x100000
#endif

void am_sleep(unsigned int);
double am_timer(double);

#endif /* AM_AM_SYSDEP_H */
