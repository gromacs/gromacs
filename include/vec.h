/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _vec_h
#define _vec_h

static char *SRCID_vec_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) vec.h 1.8 12/16/92"
#endif /* HAVE_IDENT */

/*
  collection of in-line ready operations:
  
  lookup-table optimized scalar operations:
  real invsqrt(float x)
  void vecinvsqrt(real in[],real out[],int n)
  real recip(float x)
  void vecrecip(real in[],real out[],int n)
  real sqr(real x)
  
  vector operations:
  void rvec_add(const rvec a,const rvec b,rvec c)  c = a + b
  void rvec_inc(rvec a,rvec b)                     a += b
  void rvec_sub(const rvec a,const rvec b,rvec c)  c = a - b
  void rvec_dec(rvec a,rvec b)                     a -= b
  void copy_rvec(const rvec a,rvec b)              b = a (reals)
  void copy_ivec(const ivec a,ivec b)              b = a (integers)
  void ivec_sub(const ivec a,const ivec b,ivec c)  c = a - b
  void svmul(real a,rvec v1,rvec v2)               v2 = a * v1
  void clear_rvec(rvec a)                          a = 0
  void clear_rvecs(int n,rvec v[])
  real iprod(rvec a,rvec b)                        = a . b (inner product)
  real iiprod(ivec a,ivec b)                       = a . b (integers)
  real norm2(rvec a)                               = | a |^2 ( = x*y*z )
  real norm(rvec a)                                = | a |
  void oprod(rvec a,rvec b,rvec c)                 c = a x b (outer product)
  void dprod(rvec a,rvec b,rvec c)                 c = a * b (direct product)
  real cos_angle(rvec a,rvec b)
  real cos_angle_no_table(rvec a,rvec b)
  real distance2(rvec v1, rvec v2)                 = | v2 - v1 |^2
  void unitv(rvec src,rvec dest)                   dest = src / |src|
  void unitv_no_table(rvec src,rvec dest)          dest = src / |src|
  
  matrix (3x3) operations:
  void copy_mat(matrix a,matrix b)                 b = a
  void clear_mat(matrix a)			   a = 0
  void mmul(matrix a,matrix b,matrix dest)	   dest = a . b
  void transpose(matrix src,matrix dest)	   dest = src*
  void tmmul(matrix a,matrix b,matrix dest)	   dest = a* . b
  void mtmul(matrix a,matrix b,matrix dest)	   dest = a . b*
  real det(matrix a)				   = det(a)
  void m_add(matrix a,matrix b,matrix dest)	   dest = a + b
  void m_sub(matrix a,matrix b,matrix dest)	   dest = a - b
  void msmul(matrix m1,real r1,matrix dest)	   dest = r1 * m1
  void m_inv(matrix src,matrix dest)		   dest = src^-1
  void mvmul(matrix a,rvec src,rvec dest)	   dest = a . src
  real trace(matrix m)                             = trace(m)
*/

#ifdef USE_AXP_ASM
#include "axp_asm.h"
#endif
#include "maths.h"
#include "typedefs.h"
#include "sysstuff.h"
#include "macros.h"
#include "fatal.h"
#include "x86_cpu.h"

#define EXP_LSB         0x00800000
#define EXP_MASK        0x7f800000
#define EXP_SHIFT       23
#define FRACT_MASK      0x007fffff
#define FRACT_SIZE      11              /* significant part of fraction */
#define FRACT_SHIFT     (EXP_SHIFT-FRACT_SIZE)
#define EXP_ADDR(val)   (((val)&EXP_MASK)>>EXP_SHIFT)
#define FRACT_ADDR(val) (((val)&(FRACT_MASK|EXP_LSB))>>FRACT_SHIFT)

#define PR_VEC(a)       a[XX],a[YY],a[ZZ]

#ifdef SOFTWARE_SQRT
extern const unsigned int cinvsqrtexptab[];
extern const unsigned int cinvsqrtfracttab[];
#endif

#ifdef SOFTWARE_RECIP
extern const unsigned int crecipexptab[];
extern const unsigned int crecipfracttab[];
#endif

typedef union 
{
  unsigned int bval;
  float fval;
} t_convert;


#ifdef SOFTWARE_INVSQRT
static inline real invsqrt(float x)
{
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu;
  real        y;
#ifdef DOUBLE
  real        y2;
#endif
 
  bit_pattern.fval=x;
  exp   = EXP_ADDR(bit_pattern.bval);
  fract = FRACT_ADDR(bit_pattern.bval);
  result.bval=cinvsqrtexptab[exp] | cinvsqrtfracttab[fract];
  lu    = result.fval;
  
  y=(half*lu*(three-((x*lu)*lu)));
#ifdef DOUBLE
  y2=(half*y*(three-((x*y)*y)));
  
  return y2;                    /* 10 Flops */
#else
  return y;                     /* 5  Flops */
#endif
}
#define INVSQRT_DONE 
#endif /* gmx_invsqrt */

#ifndef INVSQRT_DONE
#define invsqrt(x) (1.0f/sqrt(x))
#endif



static inline void vecinvsqrt(real in[],real out[],int n)
{
#ifdef INVSQRT_DONE  
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu,x;
#ifdef DOUBLE
  real        y;
#endif
#endif /* INVSQRT_DONE */
  int i;
  
#ifndef DOUBLE
#ifdef USE_3DNOW
  if(cpu_capabilities & X86_3DNOW_SUPPORT)
    vecinvsqrt_3dnow(in,out,n);
  else
#endif
#ifdef USE_SSE
  if((cpu_capabilities & X86_SSE_SUPPORT) && !((unsigned long int)in & 0x1f) && !((unsigned long int)out & 0x1f)) /* SSE data must be cache aligned */
    vecinvsqrt_sse(in,out,n);
  else
#endif /* no x86 optimizations */
#endif /* not double */    
#ifdef INVSQRT_DONE
    for(i=0;i<n;i++) {
      x=in[i];
      bit_pattern.fval=x;
      exp   = EXP_ADDR(bit_pattern.bval);
      fract = FRACT_ADDR(bit_pattern.bval);
      result.bval=cinvsqrtexptab[exp] | cinvsqrtfracttab[fract];
      lu    = result.fval;
      
#ifdef DOUBLE
      y=(half*lu*(three-((x*lu)*lu)));
      out[i]=(half*y*(three-((x*y)*y)));
#else
      out[i]=(half*lu*(three-((x*lu)*lu)));
#endif
    }
#else  /* no gmx invsqrt */
    for(i=0;i<n;i++)
      out[i]=1.0f/sqrt(in[i]);
#endif /* INVSQRT_DONE */
}

#ifdef SOFTWARE_RECIP
static inline real recip(float x)
{
  const real  two=2.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu;
  real        y;
#ifdef DOUBLE
  real        y2;
#endif
 
  bit_pattern.fval=x;
  exp   = EXP_ADDR(bit_pattern.bval);
  fract = FRACT_ADDR(bit_pattern.bval);
  result.bval=crecipexptab[exp] | crecipfracttab[fract];
  lu    = result.fval;
  
  y=lu*(two-x*lu);
#ifdef DOUBLE
  y2=y*(two-x*y);
  
  return y2;                    /* 6 Flops */
#else
  return y;                     /* 3  Flops */
#endif
}
#endif


static inline void vecrecip(real in[],real out[],int n)
{
#ifdef SOFTWARE_RECIP
  const real  two=2.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu,x;
#ifdef DOUBLE
  real        y;
#endif
#endif /* SOFTWARE_RECIP */
  int i;

#ifndef DOUBLE  
#ifdef USE_3DNOW
  if(cpu_capabilities & X86_3DNOW_SUPPORT)
    vecrecip_3dnow(in,out,n);
  else
#endif
#ifdef USE_SSE
  if((cpu_capabilities & X86_SSE_SUPPORT) && !((unsigned long int)in & 0x1f) && !((unsigned long int)out & 0x1f)) /* SSE data must be cache aligned */
    vecrecip_sse(in,out,n);
  else
#endif /* no x86 optimizations */
#endif /* not double */    
#ifdef SOFTWARE_RECIP
    for(i=0;i<n;i++) {
      x=in[i];
      bit_pattern.fval=x;
      exp   = EXP_ADDR(bit_pattern.bval);
      fract = FRACT_ADDR(bit_pattern.bval);
      result.bval=crecipexptab[exp] | crecipfracttab[fract];
      lu    = result.fval;
      
#ifdef DOUBLE
      y=lu*(two-x*lu);
      out[i]=y*(two-x*y);
#else
      out[i]=lu*(two-x*lu);
#endif
    }
#else /* No gmx recip */ 
    for(i=0;i<n;i++)
      out[i]=1.0f/(in[i]);
#endif /* SOFTWARE_RECIP */
}

static inline real sqr(real x)
{
  return (x*x);
}

static inline void rvec_add(const rvec a,const rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void rvec_inc(rvec a,rvec b)
{
  real x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static inline void rvec_sub(const rvec a,const rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void rvec_dec(rvec a,rvec b)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static inline void copy_rvec(const rvec a,rvec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static inline void copy_ivec(const ivec a,ivec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static inline void ivec_sub(const ivec a,const ivec b,ivec c)
{
  int x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static inline void copy_mat(matrix a,matrix b)
{
  copy_rvec(a[XX],b[XX]);
  copy_rvec(a[YY],b[YY]);
  copy_rvec(a[ZZ],b[ZZ]);
}

static inline void svmul(real a,rvec v1,rvec v2)
{
  v2[XX]=a*v1[XX];
  v2[YY]=a*v1[YY];
  v2[ZZ]=a*v1[ZZ];
}

static inline real distance2(rvec v1, rvec v2)
{
  return sqr(v2[XX]-v1[XX]) + sqr(v2[YY]-v1[YY]) + sqr(v2[ZZ]-v1[ZZ]);
}

static inline void clear_rvec(rvec a)
{
  /* The ibm compiler has problems with inlining this 
   * when we use a const real variable
   */
  a[XX]=0.0;
  a[YY]=0.0;
  a[ZZ]=0.0;
}

static inline void clear_rvecs(int n,rvec v[])
{
  int i;
  
  for(i=0; (i<n); i++) 
    clear_rvec(v[i]);
}

static inline void clear_mat(matrix a)
{
  const real nul=0.0;
  
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=nul;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=nul;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=nul;
}

static inline real iprod(rvec a,rvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline real iiprod(ivec a,ivec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline real norm2(rvec a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

static inline real norm(rvec a)
{
  return sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}

static inline real cos_angle(rvec a,rvec b)
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  real   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; /* For accuracy these must be double! */
  
  ip=ipa=ipb=0;
  for(m=0; (m<DIM); m++) {		/* 18		*/
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip*invsqrt(ipa*ipb);		/*  7		*/
					/* 25 TOTAL	*/
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}

static inline real cos_angle_no_table(rvec a,rvec b)
{
  /* This version does not need the invsqrt lookup table */
  real   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; /* For accuracy these must be double! */
  
  ip=ipa=ipb=0;
  for(m=0; (m<DIM); m++) {		/* 18		*/
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb); 		/* 12		*/
					/* 30 TOTAL	*/
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}

static inline void oprod(rvec a,rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

static inline void mmul(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[XX][YY]*b[YY][XX]+a[XX][ZZ]*b[ZZ][XX];
  dest[YY][XX]=a[YY][XX]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[YY][ZZ]*b[ZZ][XX];
  dest[ZZ][XX]=a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
  dest[XX][YY]=a[XX][XX]*b[XX][YY]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[ZZ][YY];
  dest[YY][YY]=a[YY][XX]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[ZZ][YY];
  dest[ZZ][YY]=a[ZZ][XX]*b[XX][YY]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
  dest[XX][ZZ]=a[XX][XX]*b[XX][ZZ]+a[XX][YY]*b[YY][ZZ]+a[XX][ZZ]*b[ZZ][ZZ];
  dest[YY][ZZ]=a[YY][XX]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[YY][ZZ]*b[ZZ][ZZ];
  dest[ZZ][ZZ]=a[ZZ][XX]*b[XX][ZZ]+a[ZZ][YY]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void transpose(matrix src,matrix dest)
{
  dest[XX][XX]=src[XX][XX];
  dest[YY][XX]=src[XX][YY];
  dest[ZZ][XX]=src[XX][ZZ];
  dest[XX][YY]=src[YY][XX];
  dest[YY][YY]=src[YY][YY];
  dest[ZZ][YY]=src[YY][ZZ];
  dest[XX][ZZ]=src[ZZ][XX];
  dest[YY][ZZ]=src[ZZ][YY];
  dest[ZZ][ZZ]=src[ZZ][ZZ];
}

static inline void tmmul(matrix a,matrix b,matrix dest)
{
  /* Computes dest=mmul(transpose(a),b,dest) - used in do_pr_pcoupl */
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[YY][XX]*b[YY][XX]+a[ZZ][XX]*b[ZZ][XX];
  dest[XX][YY]=a[XX][XX]*b[XX][YY]+a[YY][XX]*b[YY][YY]+a[ZZ][XX]*b[ZZ][YY];
  dest[XX][ZZ]=a[XX][XX]*b[XX][ZZ]+a[YY][XX]*b[YY][ZZ]+a[ZZ][XX]*b[ZZ][ZZ];
  dest[YY][XX]=a[XX][YY]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[ZZ][YY]*b[ZZ][XX];
  dest[YY][YY]=a[XX][YY]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[ZZ][YY]*b[ZZ][YY];
  dest[YY][ZZ]=a[XX][YY]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[ZZ][YY]*b[ZZ][ZZ];
  dest[ZZ][XX]=a[XX][ZZ]*b[XX][XX]+a[YY][ZZ]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
  dest[ZZ][YY]=a[XX][ZZ]*b[XX][YY]+a[YY][ZZ]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
  dest[ZZ][ZZ]=a[XX][ZZ]*b[XX][ZZ]+a[YY][ZZ]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void mtmul(matrix a,matrix b,matrix dest)
{
  /* Computes dest=mmul(a,transpose(b),dest) - used in do_pr_pcoupl */
  dest[XX][XX]=a[XX][XX]*b[XX][XX]+a[XX][YY]*b[XX][YY]+a[XX][ZZ]*b[XX][ZZ];
  dest[XX][YY]=a[XX][XX]*b[YY][XX]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[YY][ZZ];
  dest[XX][ZZ]=a[XX][XX]*b[ZZ][XX]+a[XX][YY]*b[ZZ][YY]+a[XX][ZZ]*b[ZZ][ZZ];
  dest[YY][XX]=a[YY][XX]*b[XX][XX]+a[YY][YY]*b[XX][YY]+a[YY][ZZ]*b[XX][ZZ];
  dest[YY][YY]=a[YY][XX]*b[YY][XX]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[YY][ZZ];
  dest[YY][ZZ]=a[YY][XX]*b[ZZ][XX]+a[YY][YY]*b[ZZ][YY]+a[YY][ZZ]*b[ZZ][ZZ];
  dest[ZZ][XX]=a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[XX][YY]+a[ZZ][ZZ]*b[XX][ZZ];
  dest[ZZ][YY]=a[ZZ][XX]*b[YY][XX]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[YY][ZZ];
  dest[ZZ][ZZ]=a[ZZ][XX]*b[ZZ][XX]+a[ZZ][YY]*b[ZZ][YY]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline real det(matrix a)
{
  return ( a[XX][XX]*(a[YY][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[YY][ZZ])
	  -a[YY][XX]*(a[XX][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[XX][ZZ])
	  +a[ZZ][XX]*(a[XX][YY]*a[YY][ZZ]-a[YY][YY]*a[XX][ZZ]));
}

static inline void m_add(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]+b[XX][XX];
  dest[XX][YY]=a[XX][YY]+b[XX][YY];
  dest[XX][ZZ]=a[XX][ZZ]+b[XX][ZZ];
  dest[YY][XX]=a[YY][XX]+b[YY][XX];
  dest[YY][YY]=a[YY][YY]+b[YY][YY];
  dest[YY][ZZ]=a[YY][ZZ]+b[YY][ZZ];
  dest[ZZ][XX]=a[ZZ][XX]+b[ZZ][XX];
  dest[ZZ][YY]=a[ZZ][YY]+b[ZZ][YY];
  dest[ZZ][ZZ]=a[ZZ][ZZ]+b[ZZ][ZZ];
}

static inline void m_sub(matrix a,matrix b,matrix dest)
{
  dest[XX][XX]=a[XX][XX]-b[XX][XX];
  dest[XX][YY]=a[XX][YY]-b[XX][YY];
  dest[XX][ZZ]=a[XX][ZZ]-b[XX][ZZ];
  dest[YY][XX]=a[YY][XX]-b[YY][XX];
  dest[YY][YY]=a[YY][YY]-b[YY][YY];
  dest[YY][ZZ]=a[YY][ZZ]-b[YY][ZZ];
  dest[ZZ][XX]=a[ZZ][XX]-b[ZZ][XX];
  dest[ZZ][YY]=a[ZZ][YY]-b[ZZ][YY];
  dest[ZZ][ZZ]=a[ZZ][ZZ]-b[ZZ][ZZ];
}

static inline void msmul(matrix m1,real r1,matrix dest)
{
  dest[XX][XX]=r1*m1[XX][XX];
  dest[XX][YY]=r1*m1[XX][YY];
  dest[XX][ZZ]=r1*m1[XX][ZZ];
  dest[YY][XX]=r1*m1[YY][XX];
  dest[YY][YY]=r1*m1[YY][YY];
  dest[YY][ZZ]=r1*m1[YY][ZZ];
  dest[ZZ][XX]=r1*m1[ZZ][XX];
  dest[ZZ][YY]=r1*m1[ZZ][YY];
  dest[ZZ][ZZ]=r1*m1[ZZ][ZZ];
}

static inline void m_inv(matrix src,matrix dest)
{
  const real smallreal = 1.0e-18;
  const real largereal = 1.0e18;
  real  deter,c,fc;
  
  deter=det(src);
  c=1.0/deter;
  fc=fabs(c);
  
  if ((fc <= smallreal) || (fc >= largereal)) 
    fatal_error(0,"Determinant = %f",1.0/c);               
  
  dest[XX][XX]= c*(src[YY][YY]*src[ZZ][ZZ]-src[ZZ][YY]*src[YY][ZZ]);
  dest[XX][YY]=-c*(src[XX][YY]*src[ZZ][ZZ]-src[ZZ][YY]*src[XX][ZZ]);
  dest[XX][ZZ]= c*(src[XX][YY]*src[YY][ZZ]-src[YY][YY]*src[XX][ZZ]);
  dest[YY][XX]=-c*(src[YY][XX]*src[ZZ][ZZ]-src[ZZ][XX]*src[YY][ZZ]);
  dest[YY][YY]= c*(src[XX][XX]*src[ZZ][ZZ]-src[ZZ][XX]*src[XX][ZZ]);
  dest[YY][ZZ]=-c*(src[XX][XX]*src[YY][ZZ]-src[YY][XX]*src[XX][ZZ]);
  dest[ZZ][XX]= c*(src[YY][XX]*src[ZZ][YY]-src[ZZ][XX]*src[YY][YY]);
  dest[ZZ][YY]=-c*(src[XX][XX]*src[ZZ][YY]-src[ZZ][XX]*src[XX][YY]);
  dest[ZZ][ZZ]= c*(src[XX][XX]*src[YY][YY]-src[YY][XX]*src[XX][YY]);
}

static inline void mvmul(matrix a,rvec src,rvec dest)
{
  dest[XX]=a[XX][XX]*src[XX]+a[XX][YY]*src[YY]+a[XX][ZZ]*src[ZZ];
  dest[YY]=a[YY][XX]*src[XX]+a[YY][YY]*src[YY]+a[YY][ZZ]*src[ZZ];
  dest[ZZ]=a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
}

static inline void unitv(rvec src,rvec dest)
{
  real linv;
  
  linv=invsqrt(norm2(src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}

static inline void unitv_no_table(rvec src,rvec dest)
{
  real linv;
  
  linv=1.0/sqrt(norm2(src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}

static inline real trace(matrix m)
{
  return (m[XX][XX]+m[YY][YY]+m[ZZ][ZZ]);
}

static inline real _divide(real a,real b,char *file,int line)
{
  if (b == 0.0) 
    fatal_error(0,"Dividing by zero, file %s, line %d",file,line);
  return a/b;
}

static inline int _mod(int a,int b,char *file,int line)
{
  if (b == 0)
    fatal_error(0,"Modulo zero, file %s, line %d",file,line);
  return a % b;
}

#define divide(a,b) _divide((a),(b),__FILE__,__LINE__)
#define mod(a,b)    _mod((a),(b),__FILE__,__LINE__)

#endif	/* _vec_h */


