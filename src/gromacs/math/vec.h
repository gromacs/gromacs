/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_MATH_VEC_H
#define GMX_MATH_VEC_H

/*
   collection of in-line ready operations:

   swap:
   void real_swap(real *a, real *b)
   void int_swap(int *a, int *b)

   vector operations:
   void rvec_opp(rvec a)                            a = -a
   void rvec_add(const rvec a,const rvec b,rvec c)  c = a + b
   void dvec_add(const dvec a,const dvec b,dvec c)  c = a + b
   void ivec_add(const ivec a,const ivec b,ivec c)  c = a + b
   void rvec_inc(rvec a,const rvec b)               a += b
   void dvec_inc(dvec a,const dvec b)               a += b
   void ivec_inc(ivec a,const ivec b)               a += b
   void rvec_sub(const rvec a,const rvec b,rvec c)  c = a - b
   void dvec_sub(const dvec a,const dvec b,dvec c)  c = a - b
   void rvec_dec(rvec a,rvec b)                     a -= b
   void copy_rvec(const rvec a,rvec b)              b = a (reals)
   void copy_dvec(const dvec a,dvec b)              b = a (reals)
   void copy_rvec_to_dvec(const rvec a,dvec b)      b = a (reals)
   void copy_dvec_to_rvec(const dvec a,rvec b)      b = a (reals)
   void copy_ivec(const ivec a,ivec b)              b = a (integers)
   void ivec_sub(const ivec a,const ivec b,ivec c)  c = a - b
   void svmul(real a,rvec v1,rvec v2)               v2 = a * v1
   void dsvmul(double a,dvec v1,dvec v2)            v2 = a * v1
   void svdiv(real a, rvec v)                       v /= a
   void clear_rvec(rvec a)                          a = 0
   void clear_dvec(dvec a)                          a = 0
   void clear_ivec(rvec a)                          a = 0
   void clear_rvecs(int n,rvec v[])
   real iprod(rvec a,rvec b)                        = a . b (inner product)
   double diprod(dvec a,dvec b)                     = a . b (inner product)
   real iiprod(ivec a,ivec b)                       = a . b (integers)
   real norm2(rvec a)                               = | a |^2 ( = x*y*z )
   double dnorm2(dvec a)                            = | a |^2 ( = x*y*z )
   real norm(rvec a)                                = | a |
   double dnorm(dvec a)                             = | a |
   void cprod(rvec a,rvec b,rvec c)                 c = a x b (cross product)
   void dcprod(dvec a,dvec b,dvec c)                c = a x b (cross product)
   void dprod(rvec a,rvec b,rvec c)                 c = a * b (direct product)
   real cos_angle(rvec a,rvec b)
   real distance2(rvec v1, rvec v2)                 = | v2 - v1 |^2
   void unitv(rvec src,rvec dest)                   dest = src / |src|

   matrix (3x3) operations:
    ! indicates that dest should not be the same as a, b or src
    the _ur0 varieties work on matrices that have only zeros
    in the upper right part, such as box matrices, these varieties
    could produce less rounding errors, not due to the operations themselves,
    but because the compiler can easier recombine the operations
   void copy_mat(matrix a,matrix b)                 b = a
   void clear_mat(matrix a)                         a = 0
   void mmul(matrix a,matrix b,matrix dest)      !  dest = a . b
   void mmul_ur0(matrix a,matrix b,matrix dest)     dest = a . b
   void transpose(matrix src,matrix dest)        !  dest = src*
   void tmmul(matrix a,matrix b,matrix dest)     !  dest = a* . b
   void mtmul(matrix a,matrix b,matrix dest)     !  dest = a . b*
   real det(matrix a)                               = det(a)
   void m_add(matrix a,matrix b,matrix dest)        dest = a + b
   void m_sub(matrix a,matrix b,matrix dest)        dest = a - b
   void msmul(matrix m1,real r1,matrix dest)        dest = r1 * m1
   void mvmul(matrix a,rvec src,rvec dest)       !  dest = a . src
   void mvmul_ur0(matrix a,rvec src,rvec dest)      dest = a . src
   void tmvmul_ur0(matrix a,rvec src,rvec dest)     dest = a* . src
   real trace(matrix m)                             = trace(m)
 */

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

static inline void real_swap(real *a, real *b)
{
  real tmp;

  tmp = *a;
  *a = *b;
  *b = tmp;
}

static inline void int_swap(int *a, int *b)
{
  int tmp;

  tmp = *a;
  *a = *b;
  *b = tmp;
}

static inline void rvec_opp(rvec a)
{
  real x,y,z;

  x = -a[XX];
  y = -a[YY];
  z = -a[ZZ];

  a[XX] = x;
  a[YY] = y;
  a[ZZ] = z;
}

static inline void rvec_add(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void dvec_add(const dvec a, const dvec b, dvec c)
{
    double x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void ivec_add(const ivec a, const ivec b, ivec c)
{
    int x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void rvec_inc(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

static inline void dvec_inc(dvec a, const dvec b)
{
    double x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

static inline void rvec_sub(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void dvec_sub(const dvec a, const dvec b, dvec c)
{
    double x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void rvec_dec(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

static inline void copy_rvec(const rvec a, rvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

static inline void copy_rvec_to_dvec(const rvec a, dvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

static inline void copy_dvec_to_rvec(const dvec a, rvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

static inline void copy_rvecn(const rvec *a, rvec *b, int startn, int endn)
{
    int i;
    for (i = startn; i < endn; i++)
    {
        b[i][XX] = a[i][XX];
        b[i][YY] = a[i][YY];
        b[i][ZZ] = a[i][ZZ];
    }
}

static inline void copy_dvec(const dvec a, dvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

static inline void copy_ivec(const ivec a, ivec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

static inline void ivec_sub(const ivec a, const ivec b, ivec c)
{
    int x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void copy_mat(const matrix a, matrix b)
{
    copy_rvec(a[XX], b[XX]);
    copy_rvec(a[YY], b[YY]);
    copy_rvec(a[ZZ], b[ZZ]);
}

static inline void svmul(real a, const rvec v1, rvec v2)
{
    v2[XX] = a*v1[XX];
    v2[YY] = a*v1[YY];
    v2[ZZ] = a*v1[ZZ];
}

static inline void dsvmul(double a, const dvec v1, dvec v2)
{
    v2[XX] = a*v1[XX];
    v2[YY] = a*v1[YY];
    v2[ZZ] = a*v1[ZZ];
}

static inline void svdiv(real a, rvec v)
{
    v[XX] /= a;
    v[YY] /= a;
    v[ZZ] /= a;
}

static inline real distance2(const rvec v1, const rvec v2)
{
    return gmx::square(v2[XX]-v1[XX]) + gmx::square(v2[YY]-v1[YY]) + gmx::square(v2[ZZ]-v1[ZZ]);
}

static inline void clear_rvec(rvec a)
{
    /* The ibm compiler has problems with inlining this
     * when we use a const real variable
     */
    a[XX] = 0.0;
    a[YY] = 0.0;
    a[ZZ] = 0.0;
}

static inline void clear_dvec(dvec a)
{
    /* The ibm compiler has problems with inlining this
     * when we use a const real variable
     */
    a[XX] = 0.0;
    a[YY] = 0.0;
    a[ZZ] = 0.0;
}

static inline void clear_ivec(ivec a)
{
    a[XX] = 0;
    a[YY] = 0;
    a[ZZ] = 0;
}

static inline void clear_rvecs(int n, rvec v[])
{
    int i;

    for (i = 0; (i < n); i++)
    {
        clear_rvec(v[i]);
    }
}

static inline void clear_mat(matrix a)
{
    const real nul = 0.0;

    a[XX][XX] = a[XX][YY] = a[XX][ZZ] = nul;
    a[YY][XX] = a[YY][YY] = a[YY][ZZ] = nul;
    a[ZZ][XX] = a[ZZ][YY] = a[ZZ][ZZ] = nul;
}

static inline real iprod(const rvec a, const rvec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline double diprod(const dvec a, const dvec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline int iiprod(const ivec a, const ivec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline real norm2(const rvec a)
{
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

static inline double dnorm2(const dvec a)
{
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

/* WARNING:
 * As dnorm() uses sqrt() (which is slow) _only_ use it if you are sure you
 * don't need 1/dnorm(), otherwise use dnorm2()*dinvnorm(). */
static inline double dnorm(const dvec a)
{
    return std::sqrt(diprod(a, a));
}

/* WARNING:
 * As norm() uses sqrt() (which is slow) _only_ use it if you are sure you
 * don't need 1/norm(), otherwise use norm2()*invnorm(). */
static inline real norm(const rvec a)
{
    return std::sqrt(iprod(a, a));
}

static inline real invnorm(const rvec a)
{
    return gmx::invsqrt(norm2(a));
}

static inline real dinvnorm(const dvec a)
{
    return gmx::invsqrt(dnorm2(a));
}

/* WARNING:
 * Do _not_ use these routines to calculate the angle between two vectors
 * as acos(cos_angle(u,v)). While it might seem obvious, the acos function
 * is very flat close to -1 and 1, which will lead to accuracy-loss.
 * Instead, use the new gmx_angle() function directly.
 */
static inline real cos_angle(const rvec a, const rvec b)
{
    /*
     *                  ax*bx + ay*by + az*bz
     * cos-vec (a,b) =  ---------------------
     *                      ||a|| * ||b||
     */
    real   cosval;
    int    m;
    double aa, bb, ip, ipa, ipb, ipab; /* For accuracy these must be double! */

    ip = ipa = ipb = 0.0;
    for (m = 0; (m < DIM); m++) /* 18 */
    {
        aa   = a[m];
        bb   = b[m];
        ip  += aa*bb;
        ipa += aa*aa;
        ipb += bb*bb;
    }
    ipab = ipa*ipb;
    if (ipab > 0)
    {
        cosval = ip*gmx::invsqrt(ipab);  /*  7 */
    }
    else
    {
        cosval = 1;
    }
    /* 25 TOTAL */
    if (cosval > 1.0)
    {
        return 1.0;
    }
    if (cosval < -1.0)
    {
        return -1.0;
    }

    return cosval;
}

static inline void cprod(const rvec a, const rvec b, rvec c)
{
    c[XX] = a[YY]*b[ZZ]-a[ZZ]*b[YY];
    c[YY] = a[ZZ]*b[XX]-a[XX]*b[ZZ];
    c[ZZ] = a[XX]*b[YY]-a[YY]*b[XX];
}

static inline void dcprod(const dvec a, const dvec b, dvec c)
{
    c[XX] = a[YY]*b[ZZ]-a[ZZ]*b[YY];
    c[YY] = a[ZZ]*b[XX]-a[XX]*b[ZZ];
    c[ZZ] = a[XX]*b[YY]-a[YY]*b[XX];
}

/* This routine calculates the angle between a & b without any loss of accuracy close to 0/PI.
 * If you only need cos(theta), use the cos_angle() routines to save a few cycles.
 * This routine is faster than it might appear, since atan2 is accelerated on many CPUs (e.g. x86).
 */
static inline real gmx_angle(const rvec a, const rvec b)
{
    rvec w;
    real wlen, s;

    cprod(a, b, w);

    wlen  = norm(w);
    s     = iprod(a, b);

    return std::atan2(wlen, s);
}

static inline double gmx_angle_between_dvecs(const dvec a, const dvec b)
{
    dvec   w;
    double wlen, s;

    dcprod(a, b, w);

    wlen  = dnorm(w);
    s     = diprod(a, b);

    return std::atan2(wlen, s);
}

static inline void mmul_ur0(const matrix a, const matrix b, matrix dest)
{
    dest[XX][XX] = a[XX][XX]*b[XX][XX];
    dest[XX][YY] = 0.0;
    dest[XX][ZZ] = 0.0;
    dest[YY][XX] = a[YY][XX]*b[XX][XX]+a[YY][YY]*b[YY][XX];
    dest[YY][YY] =                     a[YY][YY]*b[YY][YY];
    dest[YY][ZZ] = 0.0;
    dest[ZZ][XX] = a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
    dest[ZZ][YY] =                     a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
    dest[ZZ][ZZ] =                                         a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void mmul(const matrix a, const matrix b, matrix dest)
{
    dest[XX][XX] = a[XX][XX]*b[XX][XX]+a[XX][YY]*b[YY][XX]+a[XX][ZZ]*b[ZZ][XX];
    dest[YY][XX] = a[YY][XX]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[YY][ZZ]*b[ZZ][XX];
    dest[ZZ][XX] = a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
    dest[XX][YY] = a[XX][XX]*b[XX][YY]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[ZZ][YY];
    dest[YY][YY] = a[YY][XX]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[ZZ][YY];
    dest[ZZ][YY] = a[ZZ][XX]*b[XX][YY]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
    dest[XX][ZZ] = a[XX][XX]*b[XX][ZZ]+a[XX][YY]*b[YY][ZZ]+a[XX][ZZ]*b[ZZ][ZZ];
    dest[YY][ZZ] = a[YY][XX]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[YY][ZZ]*b[ZZ][ZZ];
    dest[ZZ][ZZ] = a[ZZ][XX]*b[XX][ZZ]+a[ZZ][YY]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void transpose(const matrix src, matrix dest)
{
    dest[XX][XX] = src[XX][XX];
    dest[YY][XX] = src[XX][YY];
    dest[ZZ][XX] = src[XX][ZZ];
    dest[XX][YY] = src[YY][XX];
    dest[YY][YY] = src[YY][YY];
    dest[ZZ][YY] = src[YY][ZZ];
    dest[XX][ZZ] = src[ZZ][XX];
    dest[YY][ZZ] = src[ZZ][YY];
    dest[ZZ][ZZ] = src[ZZ][ZZ];
}

static inline void tmmul(const matrix a, const matrix b, matrix dest)
{
    /* Computes dest=mmul(transpose(a),b,dest) - used in do_pr_pcoupl */
    dest[XX][XX] = a[XX][XX]*b[XX][XX]+a[YY][XX]*b[YY][XX]+a[ZZ][XX]*b[ZZ][XX];
    dest[XX][YY] = a[XX][XX]*b[XX][YY]+a[YY][XX]*b[YY][YY]+a[ZZ][XX]*b[ZZ][YY];
    dest[XX][ZZ] = a[XX][XX]*b[XX][ZZ]+a[YY][XX]*b[YY][ZZ]+a[ZZ][XX]*b[ZZ][ZZ];
    dest[YY][XX] = a[XX][YY]*b[XX][XX]+a[YY][YY]*b[YY][XX]+a[ZZ][YY]*b[ZZ][XX];
    dest[YY][YY] = a[XX][YY]*b[XX][YY]+a[YY][YY]*b[YY][YY]+a[ZZ][YY]*b[ZZ][YY];
    dest[YY][ZZ] = a[XX][YY]*b[XX][ZZ]+a[YY][YY]*b[YY][ZZ]+a[ZZ][YY]*b[ZZ][ZZ];
    dest[ZZ][XX] = a[XX][ZZ]*b[XX][XX]+a[YY][ZZ]*b[YY][XX]+a[ZZ][ZZ]*b[ZZ][XX];
    dest[ZZ][YY] = a[XX][ZZ]*b[XX][YY]+a[YY][ZZ]*b[YY][YY]+a[ZZ][ZZ]*b[ZZ][YY];
    dest[ZZ][ZZ] = a[XX][ZZ]*b[XX][ZZ]+a[YY][ZZ]*b[YY][ZZ]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline void mtmul(const matrix a, const matrix b, matrix dest)
{
    /* Computes dest=mmul(a,transpose(b),dest) - used in do_pr_pcoupl */
    dest[XX][XX] = a[XX][XX]*b[XX][XX]+a[XX][YY]*b[XX][YY]+a[XX][ZZ]*b[XX][ZZ];
    dest[XX][YY] = a[XX][XX]*b[YY][XX]+a[XX][YY]*b[YY][YY]+a[XX][ZZ]*b[YY][ZZ];
    dest[XX][ZZ] = a[XX][XX]*b[ZZ][XX]+a[XX][YY]*b[ZZ][YY]+a[XX][ZZ]*b[ZZ][ZZ];
    dest[YY][XX] = a[YY][XX]*b[XX][XX]+a[YY][YY]*b[XX][YY]+a[YY][ZZ]*b[XX][ZZ];
    dest[YY][YY] = a[YY][XX]*b[YY][XX]+a[YY][YY]*b[YY][YY]+a[YY][ZZ]*b[YY][ZZ];
    dest[YY][ZZ] = a[YY][XX]*b[ZZ][XX]+a[YY][YY]*b[ZZ][YY]+a[YY][ZZ]*b[ZZ][ZZ];
    dest[ZZ][XX] = a[ZZ][XX]*b[XX][XX]+a[ZZ][YY]*b[XX][YY]+a[ZZ][ZZ]*b[XX][ZZ];
    dest[ZZ][YY] = a[ZZ][XX]*b[YY][XX]+a[ZZ][YY]*b[YY][YY]+a[ZZ][ZZ]*b[YY][ZZ];
    dest[ZZ][ZZ] = a[ZZ][XX]*b[ZZ][XX]+a[ZZ][YY]*b[ZZ][YY]+a[ZZ][ZZ]*b[ZZ][ZZ];
}

static inline real det(const matrix a)
{
    return ( a[XX][XX]*(a[YY][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[YY][ZZ])
             -a[YY][XX]*(a[XX][YY]*a[ZZ][ZZ]-a[ZZ][YY]*a[XX][ZZ])
             +a[ZZ][XX]*(a[XX][YY]*a[YY][ZZ]-a[YY][YY]*a[XX][ZZ]));
}


static inline void m_add(const matrix a, const matrix b, matrix dest)
{
    dest[XX][XX] = a[XX][XX]+b[XX][XX];
    dest[XX][YY] = a[XX][YY]+b[XX][YY];
    dest[XX][ZZ] = a[XX][ZZ]+b[XX][ZZ];
    dest[YY][XX] = a[YY][XX]+b[YY][XX];
    dest[YY][YY] = a[YY][YY]+b[YY][YY];
    dest[YY][ZZ] = a[YY][ZZ]+b[YY][ZZ];
    dest[ZZ][XX] = a[ZZ][XX]+b[ZZ][XX];
    dest[ZZ][YY] = a[ZZ][YY]+b[ZZ][YY];
    dest[ZZ][ZZ] = a[ZZ][ZZ]+b[ZZ][ZZ];
}

static inline void m_sub(const matrix a, const matrix b, matrix dest)
{
    dest[XX][XX] = a[XX][XX]-b[XX][XX];
    dest[XX][YY] = a[XX][YY]-b[XX][YY];
    dest[XX][ZZ] = a[XX][ZZ]-b[XX][ZZ];
    dest[YY][XX] = a[YY][XX]-b[YY][XX];
    dest[YY][YY] = a[YY][YY]-b[YY][YY];
    dest[YY][ZZ] = a[YY][ZZ]-b[YY][ZZ];
    dest[ZZ][XX] = a[ZZ][XX]-b[ZZ][XX];
    dest[ZZ][YY] = a[ZZ][YY]-b[ZZ][YY];
    dest[ZZ][ZZ] = a[ZZ][ZZ]-b[ZZ][ZZ];
}

static inline void msmul(const matrix m1, real r1, matrix dest)
{
    dest[XX][XX] = r1*m1[XX][XX];
    dest[XX][YY] = r1*m1[XX][YY];
    dest[XX][ZZ] = r1*m1[XX][ZZ];
    dest[YY][XX] = r1*m1[YY][XX];
    dest[YY][YY] = r1*m1[YY][YY];
    dest[YY][ZZ] = r1*m1[YY][ZZ];
    dest[ZZ][XX] = r1*m1[ZZ][XX];
    dest[ZZ][YY] = r1*m1[ZZ][YY];
    dest[ZZ][ZZ] = r1*m1[ZZ][ZZ];
}

static inline void mvmul(const matrix a, const rvec src, rvec dest)
{
    dest[XX] = a[XX][XX]*src[XX]+a[XX][YY]*src[YY]+a[XX][ZZ]*src[ZZ];
    dest[YY] = a[YY][XX]*src[XX]+a[YY][YY]*src[YY]+a[YY][ZZ]*src[ZZ];
    dest[ZZ] = a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
}


static inline void mvmul_ur0(const matrix a, const rvec src, rvec dest)
{
    dest[ZZ] = a[ZZ][XX]*src[XX]+a[ZZ][YY]*src[YY]+a[ZZ][ZZ]*src[ZZ];
    dest[YY] = a[YY][XX]*src[XX]+a[YY][YY]*src[YY];
    dest[XX] = a[XX][XX]*src[XX];
}

static inline void tmvmul_ur0(const matrix a, const rvec src, rvec dest)
{
    dest[XX] = a[XX][XX]*src[XX]+a[YY][XX]*src[YY]+a[ZZ][XX]*src[ZZ];
    dest[YY] =                   a[YY][YY]*src[YY]+a[ZZ][YY]*src[ZZ];
    dest[ZZ] =                                     a[ZZ][ZZ]*src[ZZ];
}

static inline void unitv(const rvec src, rvec dest)
{
    real linv;

    linv     = gmx::invsqrt(norm2(src));
    dest[XX] = linv*src[XX];
    dest[YY] = linv*src[YY];
    dest[ZZ] = linv*src[ZZ];
}

static inline real trace(const matrix m)
{
    return (m[XX][XX]+m[YY][YY]+m[ZZ][ZZ]);
}

#endif
