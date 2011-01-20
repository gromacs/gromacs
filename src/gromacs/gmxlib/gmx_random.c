/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gmx_matrix.c,v 1.4 2008/12/02 18:27:57 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gmx_random.h>

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <time.h>
#include <math.h>
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <process.h>
#endif

#include "maths.h"
#include "gmx_random_gausstable.h"

#define RNG_N 624
#define RNG_M 397
#define RNG_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define RNG_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define RNG_LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* Note that if you change the size of the Gaussian table you will also
 * have to generate new initialization data for the table in
 * gmx_random_gausstable.h
 *
 * We have included the routine print_gaussian_table() in this file
 * for convenience - use it if you need a different size of the table.
 */
#define GAUSS_TABLE 14 /* the size of the gauss table is 2^GAUSS_TABLE */
#define GAUSS_SHIFT (32 - GAUSS_TABLE)


struct gmx_rng {
  unsigned int  mt[624];  
  int           mti;
  int           has_saved;  
  double        gauss_saved;
};



int
gmx_rng_n(void)
{
  return RNG_N;
}


gmx_rng_t 
gmx_rng_init(unsigned int seed)
{
  struct gmx_rng *rng;
    
  if((rng=(struct gmx_rng *)malloc(sizeof(struct gmx_rng)))==NULL)
	return NULL;
  
  rng->has_saved=0; /* no saved gaussian number yet */

  rng->mt[0]= seed & 0xffffffffUL;
  for (rng->mti=1; rng->mti<RNG_N; rng->mti++) {
    rng->mt[rng->mti] = 
      (1812433253UL * (rng->mt[rng->mti-1] ^
		       (rng->mt[rng->mti-1] >> 30)) + rng->mti); 
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    rng->mt[rng->mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
  return rng;
}

gmx_rng_t 
gmx_rng_init_array(unsigned int seed[], int seed_length)
{
    int i, j, k;
    gmx_rng_t rng;

    if((rng=gmx_rng_init(19650218UL))==NULL)
	  return NULL;
	
    i=1; j=0;
    k = (RNG_N>seed_length ? RNG_N : seed_length);
    for (; k; k--) {
        rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^
				     (rng->mt[i-1] >> 30)) * 1664525UL))
          + seed[j] + j; /* non linear */
        rng->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=RNG_N) { rng->mt[0] = rng->mt[RNG_N-1]; i=1; }
        if (j>=seed_length) j=0;
    }
    for (k=RNG_N-1; k; k--) {
        rng->mt[i] = (rng->mt[i] ^ ((rng->mt[i-1] ^ 
				     (rng->mt[i-1] >> 30)) * 
				    1566083941UL))
          - i; /* non linear */
        rng->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=RNG_N) { rng->mt[0] = rng->mt[RNG_N-1]; i=1; }
    }

    rng->mt[0] = 0x80000000UL; 
    /* MSB is 1; assuring non-zero initial array */ 
    return rng;
}


void
gmx_rng_destroy(gmx_rng_t rng)
{
  if(rng)
    free(rng);
}


void
gmx_rng_get_state(gmx_rng_t rng, unsigned int *mt,int *mti)
{
  int i;

  for(i=0; i<RNG_N; i++) {
    mt[i] = rng->mt[i];
  }
  *mti = rng->mti;
}


void
gmx_rng_set_state(gmx_rng_t rng,  unsigned int *mt,int mti)
{
  int i;

  for(i=0; i<RNG_N; i++) {
    rng->mt[i] = mt[i];
  }
  rng->mti = mti;
}


unsigned int
gmx_rng_make_seed(void)
{
  FILE *fp;
  unsigned int data;
  int ret;
  long my_pid;

#ifdef HAVE_UNISTD_H
  fp=fopen("/dev/random","rb"); /* will return NULL if it is not present */
#else
  fp=NULL;
#endif
  if(fp!=NULL) {
    ret=fread(&data,sizeof(unsigned int),1,fp);
    fclose(fp);
  } else {
    /* No random device available, use time-of-day and process id */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    my_pid = (long)_getpid();
#else
    my_pid = (long)getpid();
#endif
    data=(unsigned int)(((long)time(NULL)+my_pid) % (long)1000000);
  }
  return data;
}


/* The random number state contains 624 entries that are returned one by
 * one as random numbers. When we run out of them, this routine is called to
 * regenerate 624 new entries.
 */
static void
gmx_rng_update(gmx_rng_t rng)
{
  unsigned int lastx,x1,x2,y,*mt;
  int mti,k;
  const unsigned int mag01[2]={0x0UL, RNG_MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    /* update random numbers */
  mt=rng->mt;   /* pointer to array - avoid repeated dereferencing */
  mti=rng->mti;
  
  x1=mt[0];
  for (k=0;k<RNG_N-RNG_M-3;k+=4) {
    x2=mt[k+1];
    y=(x1&RNG_UPPER_MASK)|(x2&RNG_LOWER_MASK);
    mt[k]= mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    x1=mt[k+2];
    y= (x2&RNG_UPPER_MASK)|(x1&RNG_LOWER_MASK);
    mt[k+1] = mt[k+RNG_M+1] ^ (y >> 1) ^ mag01[y & 0x1UL];
    x2=mt[k+3];
    y=(x1&RNG_UPPER_MASK)|(x2&RNG_LOWER_MASK);
    mt[k+2]= mt[k+RNG_M+2] ^ (y >> 1) ^ mag01[y & 0x1UL];
    x1=mt[k+4];
    y= (x2&RNG_UPPER_MASK)|(x1&RNG_LOWER_MASK);
    mt[k+3] = mt[k+RNG_M+3] ^ (y >> 1) ^ mag01[y & 0x1UL];
  }
  x2=mt[k+1];
  y=(x1&RNG_UPPER_MASK)|(x2&RNG_LOWER_MASK);
  mt[k]= mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
  k++;
  x1=mt[k+1];
  y=(x2&RNG_UPPER_MASK)|(x1&RNG_LOWER_MASK);
  mt[k]= mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
  k++;
  x2=mt[k+1];
  y=(x1&RNG_UPPER_MASK)|(x2&RNG_LOWER_MASK);
  mt[k]= mt[k+RNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
  k++;
  for (;k<RNG_N-1;k+=4) {
    x1=mt[k+1];
    y = (x2&RNG_UPPER_MASK)|(x1&RNG_LOWER_MASK);
    mt[k] = mt[k+(RNG_M-RNG_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    x2=mt[k+2];
    y = (x1&RNG_UPPER_MASK)|(x2&RNG_LOWER_MASK);
    mt[k+1] = mt[k+(RNG_M-RNG_N)+1] ^ (y >> 1) ^ mag01[y & 0x1UL];
    x1=mt[k+3];
    y = (x2&RNG_UPPER_MASK)|(x1&RNG_LOWER_MASK);
    mt[k+2] = mt[k+(RNG_M-RNG_N)+2] ^ (y >> 1) ^ mag01[y & 0x1UL];
    x2=mt[k+4];
    y = (x1&RNG_UPPER_MASK)|(x2&RNG_LOWER_MASK);
    mt[k+3] = mt[k+(RNG_M-RNG_N)+3] ^ (y >> 1) ^ mag01[y & 0x1UL];
  }
  y = (x2&RNG_UPPER_MASK)|(mt[0]&RNG_LOWER_MASK);
  mt[RNG_N-1] = mt[RNG_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
  
  rng->mti = 0;
}


real
gmx_rng_gaussian_real(gmx_rng_t rng)
{
  real x,y,r;

  if(rng->has_saved) {
    rng->has_saved=0;
    return rng->gauss_saved;
  } else {
    do {
      x=2.0*gmx_rng_uniform_real(rng)-1.0;
      y=2.0*gmx_rng_uniform_real(rng)-1.0;
      r=x*x+y*y;
    } while(r>1.0 || r==0.0);
    
    r=sqrt(-2.0*log(r)/r);
    rng->gauss_saved=y*r; /* save second random number */
    rng->has_saved=1;
    return x*r; /* return first random number */
  }
}




/* Return a random unsigned integer, i.e. 0..4294967295 
 * Provided in header file for performace reasons.
 * Unfortunately this function cannot be inlined, since
 * it needs to refer the internal-linkage gmx_rng_update().
 */
unsigned int
gmx_rng_uniform_uint32(gmx_rng_t rng)
{
  unsigned int y;
  
  if(rng->mti==624) 
    gmx_rng_update(rng);
  y=rng->mt[rng->mti++];
  
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return y;  
} 





/* Return a uniform floating point number on the interval 0<=x<1 */
real
gmx_rng_uniform_real(gmx_rng_t rng)
{
  if(sizeof(real)==sizeof(double))
    return ((double)gmx_rng_uniform_uint32(rng))*(1.0/4294967296.0); 
  else
    return ((float)gmx_rng_uniform_uint32(rng))*(1.0/4294967423.0); 
  /* divided by the smallest number that will generate a 
    * single precision real number on 0<=x<1.
    * This needs to be slightly larger than MAX_UNIT since
    * we are limited to an accuracy of 1e-7.
    */
}



real 
gmx_rng_gaussian_table(gmx_rng_t rng)
{
  unsigned int i;
  
  i = gmx_rng_uniform_uint32(rng);
  
  /* The Gaussian table is a static constant in this file */
  return gaussian_table[i >> GAUSS_SHIFT];
}


/*
 * Print a lookup table for Gaussian numbers with 4 entries on each
 * line, formatted for inclusion in this file. Size is 2^bits.
 */
void
print_gaussian_table(int bits)
{
  int n,nh,i,j;
  double invn,fac,x,invgauss,det,dx;
  real  *table;
  
  n = 1 << bits;
  table = (real *)malloc(n*sizeof(real));
  
  /* Fill a table of size n such that random draws from it
    * produce a Gaussian distribution.
    * We integrate the Gaussian distribution G approximating:
    *   integral(x->x+dx) G(y) dy
    * with:
    *   G(x) dx + G'(x) dx^2/2 = G(x) dx - G(x) x dx^2/2
    * Then we need to find dx such that the integral is 1/n.
    * The last step uses dx = 1/x as the approximation is not accurate enough.
    */
  invn = 1.0/n;
  fac = sqrt(2*M_PI);
  x = 0.5*fac*invn;
  nh = n/2;
  for(i=0; i<nh; i++) {
    if (i > 0) {
      if (i < nh-1) {
	invgauss = fac*exp(0.5*x*x);
	/* det is larger than 0 for all x, except for the last */
	det = 1 - 2*invn*x*invgauss;
	dx = (1 - sqrt(det))/x;
      } else {
	dx = 1/x;
      }
      x = x + dx;
    }
    table[nh-1-i] = -x;
    table[nh+i]   =  x;
  }
  printf("static const real *\ngaussian_table[%d] = {\n",n);
  for(i=0;i<n;i+=4) {
    printf("  ");
    for(j=0;j<4;j++) {
      printf("%14.7e",table[i+j]);
      if(i+j<(n-1))
	printf(",");
    }
    printf("\n");
  }
  printf("};\n");
  free(table);
}

