/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg and Magnus Lundborg
 * Copyright (c) 2010, 2013-2014 The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/compression/tng_compress.h"

/* 64 bit integers are not required in this part of the program. But
   they improve the running speed if present. If 64 bit integers are
   available define the symbol HAVE64BIT. It should get automatically
   defined by the defines in my64bit.h */
#include "../../include/compression/my64bit.h"

#include "../../include/compression/widemuldiv.h"

#ifndef TRAJNG_X86_GCC_INLINE_MULDIV
#if defined(__GNUC__) && defined(__i386__)
#define TRAJNG_X86_GCC_INLINE_MULDIV
#endif /* gcc & i386 */
#if defined(__GNUC__) && defined(__x86_64__)
#define TRAJNG_X86_GCC_INLINE_MULDIV
#endif /* gcc & x86_64 */
#endif /* TRAJNG X86 GCC INLINE MULDIV */

#ifdef USE_WINDOWS
#define TNG_INLINE __inline
#else
#define TNG_INLINE inline
#endif

/* Multiply two 32 bit unsigned integers returning a 64 bit unsigned value (in two integers) */
static TNG_INLINE void Ptngc_widemul(unsigned int i1, unsigned int i2, unsigned int *ohi, unsigned int *olo)
{
#if defined(TRAJNG_X86_GCC_INLINE_MULDIV)
  __asm__ __volatile__ ("mull %%edx\n\t"
                        : "=a" (i1), "=d" (i2)
                        : "a" (i1),"d" (i2)
                        : "cc");
  *ohi=i2;
  *olo=i1;
#else /* TRAJNG X86 GCC INLINE MULDIV */

#ifdef HAVE64BIT
  my_uint64_t res= ((my_uint64_t)i1) * ((my_uint64_t)i2);
  *olo=res & 0xFFFFFFFFU;
  *ohi=(res>>32) & 0xFFFFFFFFU;
#else /* HAVE64BIT */

  unsigned int bits=16;
  unsigned int L_m=(1<<bits)-1; /* Lower bits mask. */

  unsigned int a_U,a_L; /* upper and lower half of a */
  unsigned int b_U,b_L; /* upper and lower half of b */

  unsigned int x_UU,x_UL,x_LU,x_LL; /* temporary storage. */
  unsigned int x,x_U,x_L; /* temporary storage. */

  /* Extract partial ints */
  a_L=i1 & L_m;
  a_U=i1>>bits;
  b_L=i2 & L_m;
  b_U=i2>>bits;

  /* Do a*a=>2a multiply where a is half number of bits in an int */
  x=a_L*b_L;
  x_LL=x & L_m;
  x_LU=x>>bits;

  x=a_U*b_L;
  x_LU+=x & L_m;
  x_UL=x>>bits;

  x=a_L*b_U;
  x_LU+=x & L_m;
  x_UL+=x>>bits;

  x=a_U*b_U;
  x_UL+=x & L_m;
  x_UU=x>>bits;

  /* Combine results */
  x_UL+=x_LU>>bits;
  x_LU&=L_m;
  x_UU+=x_UL>>bits;
  x_UL&=L_m;

  x_U=(x_UU<<bits)|x_UL;
  x_L=(x_LU<<bits)|x_LL;

  /* Write back results */
  *ohi=x_U;
  *olo=x_L;
#endif /* HAVE64BIT */
#endif /* TRAJNG X86 GCC INLINE MULDIV */
}

/* Divide a 64 bit unsigned value in hi:lo with the 32 bit value i and
   return the result in result and the remainder in remainder */
static TNG_INLINE void Ptngc_widediv(unsigned int hi, unsigned int lo, const unsigned int i, unsigned int *result, unsigned int *remainder)
{
#if defined(TRAJNG_X86_GCC_INLINE_MULDIV)
  __asm__ __volatile__ ("divl %%ecx\n\t"
                        : "=a" (lo), "=d" (hi)
                        : "a" (lo),"d" (hi), "c" (i)
                        : "cc");
  *result=lo;
  *remainder=hi;
#else /* TRAJNG X86 GCC INLINE MULDIV */
#ifdef HAVE64BIT
  my_uint64_t v= (((my_uint64_t)hi)<<32) | ((my_uint64_t)lo);
  my_uint64_t res=v/((my_uint64_t)i);
  my_uint64_t rem=v-res*((my_uint64_t)i);
  *result=(unsigned int)res;
  *remainder=(unsigned int)rem;
#else /* HAVE64BIT */
  unsigned int res;
  unsigned int rmask;
  unsigned int s_U,s_L;
  unsigned int bits=16U;
  unsigned int bits2=bits*2U;
  unsigned int hibit=bits2-1U;
  unsigned int hibit_mask=1U<<hibit;
  unsigned int allbits=(hibit_mask-1U)|hibit_mask;

  /* Do division. */
  rmask=hibit_mask;
  res=0;
  s_U=i>>(bits2-hibit);
  s_L=(i<<hibit)&0xFFFFFFFFU;
  while (rmask)
    {
      if ((s_U<hi) || ((s_U==hi) && (s_L<=lo)))
        {
          /* Subtract */
          hi-=s_U;
          if (s_L>lo)
            {
              unsigned int t;
              hi--; /* Borrow */
              t=allbits-s_L;
              lo+=t+1;
            }
          else
            lo-=s_L;

          /* Set bit. */
          res|=rmask;
        }
      rmask>>=1;
      s_L>>=1;
      if (s_U & 1U)
        s_L|=hibit_mask;
      s_U>>=1;
    }
  *remainder=lo;
  *result=res;
#endif /* HAVE64BIT */
#endif /* TRAJNG X86 GCC INLINE MULDIV */
}

/* Add a unsigned int to a largeint. j determines which value in the
   largeint to add v1 to. */
static TNG_INLINE void largeint_add_gen(const unsigned int v1, unsigned int *largeint, const int n, int j)
{
  /* Add with carry. unsigned ints in C wrap modulo 2**bits when "overflowed". */
  unsigned int v2=(v1+largeint[j])&0xFFFFFFFFU; /* Add and cap at 32 bits */
  unsigned int carry=0;
  if ((((unsigned int)-1)&0xFFFFFFFFU) -v1<largeint[j])
    carry=1;
  largeint[j]=v2;
  j++;
  while ((j<n) && carry)
    {
      v2=(largeint[j]+carry)&0xFFFFFFFFU;
      carry=0;
      if ((((unsigned int)-1)&0xFFFFFFFFU) -1<largeint[j])
        carry=1;
      largeint[j]=v2;
      j++;
    }
}

/* Add a unsigned int to a largeint. */
void Ptngc_largeint_add(const unsigned int v1, unsigned int *largeint, const int n)
{
  largeint_add_gen(v1,largeint,n,0);
}

/* Multiply v1 with largeint_in and return result in largeint_out */
void Ptngc_largeint_mul(const unsigned int v1, unsigned int *largeint_in, unsigned int *largeint_out, const int n)
{
  int i;
  unsigned int lo,hi;

  memset(largeint_out, 0U, sizeof(unsigned int) * n);

  for (i=0; i<n-1; i++)
    {
      if (largeint_in[i]!=0U)
        {
          Ptngc_widemul(v1,largeint_in[i],&hi,&lo); /* 32x32->64 mul */
          largeint_add_gen(lo,largeint_out,n,i);
          largeint_add_gen(hi,largeint_out,n,i+1);
        }
    }
  if (largeint_in[i]!=0U)
    {
      Ptngc_widemul(v1,largeint_in[i],&hi,&lo); /* 32x32->64 mul */
      largeint_add_gen(lo,largeint_out,n,i);
    }
}

/* Return the remainder from dividing largeint_in with v1. Result of the division is returned in largeint_out */
unsigned int Ptngc_largeint_div(const unsigned int v1, unsigned int *largeint_in, unsigned int *largeint_out, const int n)
{
  unsigned int result,remainder=0;
  int i;
  unsigned int hi;
  /* Boot */
  hi=0U;
  i=n;
  while (i)
    {
      i--;
      Ptngc_widediv(hi,largeint_in[i],v1,&result,&remainder);
      largeint_out[i]=result;
      hi=remainder;
    }
  return remainder;
}
