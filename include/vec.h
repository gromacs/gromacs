/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _vec_h
#define _vec_h

static char *SRCID_vec_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) vec.h 1.8 12/16/92"
#endif /* HAVE_IDENT */

#include "maths.h"
#include "typedefs.h"
#include "sysstuff.h"
#include "macros.h"
#include "fatal.h"
#include "lutab.h"

#ifdef _lnx_
#define gmx_inline inline
#else
#define gmx_inline
#endif

#ifdef CINVSQRT
static gmx_inline real invsqrt(float x)
{
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  word        exp,fract;
  float       lu;
  real        y;
#ifdef DOUBLE
  real        y2;
#endif

  bit_pattern.fval=x;
  exp   = EXP_ADDR(bit_pattern.bval);
  fract = FRACT_ADDR(bit_pattern.bval);
  result.bval=lookup_table.exp_seed[exp] | lookup_table.fract_seed[fract];
  lu    = result.fval;
  
  y=(half*lu*(three-((x*lu)*lu)));
#ifdef DOUBLE
  y2=(half*y*(three-((x*y)*y)));
  
  return y2;                    /* 10 Flops */
#else
  return y;			/* 5  Flops */
#endif
}
#define INVSQRT_DONE
#endif

#ifndef INVSQRT_DONE
#define invsqrt(x) (1.0f/sqrt(x))
#endif

static gmx_inline real sqr(real x)
{
  return (x*x);
}

static gmx_inline void rvec_add(rvec a,rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static gmx_inline void rvec_inc(rvec a,rvec b)
{
  real x,y,z;
  
  x=a[XX]+b[XX];
  y=a[YY]+b[YY];
  z=a[ZZ]+b[ZZ];
  
  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static gmx_inline void rvec_sub(rvec a,rvec b,rvec c)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static gmx_inline void rvec_dec(rvec a,rvec b)
{
  real x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  a[XX]=x;
  a[YY]=y;
  a[ZZ]=z;
}

static gmx_inline void copy_rvec(rvec a,rvec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static gmx_inline void copy_ivec(ivec a,ivec b)
{
  b[XX]=a[XX];
  b[YY]=a[YY];
  b[ZZ]=a[ZZ];
}

static gmx_inline void ivec_sub(ivec a,ivec b,ivec c)
{
  int x,y,z;
  
  x=a[XX]-b[XX];
  y=a[YY]-b[YY];
  z=a[ZZ]-b[ZZ];
  
  c[XX]=x;
  c[YY]=y;
  c[ZZ]=z;
}

static gmx_inline void copy_mat(matrix a,matrix b)
{
  copy_rvec(a[XX],b[XX]);
  copy_rvec(a[YY],b[YY]);
  copy_rvec(a[ZZ],b[ZZ]);
}

static gmx_inline void svmul(real a,rvec v1,rvec v2)
{
  v2[XX]=a*v1[XX];
  v2[YY]=a*v1[YY];
  v2[ZZ]=a*v1[ZZ];
}

static gmx_inline real distance2(rvec v1, rvec v2)
{
  return sqr(v2[XX]-v1[XX]) + sqr(v2[YY]-v1[YY]) + sqr(v2[ZZ]-v1[ZZ]);
}

static gmx_inline void clear_rvec(rvec a)
{
  const real nul=0.0;
  
  a[XX]=a[YY]=a[ZZ]=nul;
}

static gmx_inline void clear_rvecs(int n,rvec v[])
{
  int i;
  
  for(i=0; (i<n); i++) 
    clear_rvec(v[i]);
}

static gmx_inline void clear_mat(matrix a)
{
  const real nul=0.0;
  
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=nul;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=nul;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=nul;
}

static gmx_inline real iprod(rvec a,rvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static gmx_inline real iiprod(ivec a,ivec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static gmx_inline real norm2(rvec a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

static gmx_inline real norm(rvec a)
{
  return sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}

static gmx_inline real cos_angle(rvec a,rvec b)
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

static gmx_inline real cos_angle_no_table(rvec a,rvec b)
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

static gmx_inline void oprod(rvec a,rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

static gmx_inline void mmul(matrix a,matrix b,matrix dest)
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

static gmx_inline real det(matrix a)
{
  return ( a[XX][XX]*(a[YY][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[YY][ZZ])
	  -a[YY][XX]*(a[XX][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[XX][ZZ])
	  +a[ZZ][XX]*(a[XX][YY]*a[YY][ZZ]-a[YY][YY]*a[XX][ZZ]));
}

static gmx_inline void m_add(matrix a,matrix b,matrix dest)
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

static gmx_inline void m_sub(matrix a,matrix b,matrix dest)
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

static gmx_inline void msmul(matrix m1,real r1,matrix dest)
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

static gmx_inline void m_inv(matrix src,matrix dest)
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

static gmx_inline void mvmul(matrix a,rvec src,rvec dest)
{
  dest[XX]=a[XX][XX]*src[XX]+a[XX][YY]*src[YY]+a[XX][ZZ]*src[ZZ];
  dest[YY]=a[YY][XX]*src[XX]+a[YY][YY]*src[YY]+a[YY][ZZ]*src[ZZ];
  dest[ZZ]=a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
}

static gmx_inline void unitv(rvec src,rvec dest)
{
  real linv;
  
  linv=invsqrt(iprod(src,src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}

static gmx_inline void unitv_no_table(rvec src,rvec dest)
{
  real linv;
  
  linv=1.0/sqrt(iprod(src,src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}

static gmx_inline real trace(matrix m)
{
  return (m[XX][XX]+m[YY][YY]+m[ZZ][ZZ]);
}

static gmx_inline real _divide(real a,real b,char *file,int line)
{
  if (b == 0.0) 
    fatal_error(0,"Dividing by zero, file %s, line %d",file,line);
  return a/b;
}

static gmx_inline int _mod(int a,int b,char *file,int line)
{
  if (b == 0)
    fatal_error(0,"Modulo zero, file %s, line %d",file,line);
  return a % b;
}

#define divide(a,b) _divide((a),(b),__FILE__,__LINE__)
#define mod(a,b)    _mod((a),(b),__FILE__,__LINE__)

#endif	/* _vec_h */


