/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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


/* maths operations */
/* imported from cpu versions in math/vec.h */
__forceinline__ __device__
void svmul_gpu(real a, const rvec v1, rvec v2)
{
    v2[XX] = a*v1[XX];
    v2[YY] = a*v1[YY];
    v2[ZZ] = a*v1[ZZ];
}


__forceinline__ __device__
 void rvec_add_gpu(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__
void ivec_add_gpu(const ivec a, const ivec b, ivec c)
{
    int x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__
void rvec_inc_atomic(rvec a, const rvec b)
{
  atomicAdd(&a[XX],b[XX]);
  atomicAdd(&a[YY],b[YY]);
  atomicAdd(&a[ZZ],b[ZZ]);
}

__forceinline__ __device__
void rvec_inc_gpu(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

__forceinline__ __device__
void rvec_dec_atomic(rvec a, const rvec b)
{
/* real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;*/

  atomicAdd(&a[XX],-1.0f*b[XX]);
  atomicAdd(&a[YY],-1.0f*b[YY]);
  atomicAdd(&a[ZZ],-1.0f*b[ZZ]);
}

__forceinline__ __device__
void rvec_dec_gpu(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

__forceinline__ __device__
void cprod_gpu(const rvec a, const rvec b, rvec c)
{
    c[XX] = a[YY]*b[ZZ]-a[ZZ]*b[YY];
    c[YY] = a[ZZ]*b[XX]-a[XX]*b[ZZ];
    c[ZZ] = a[XX]*b[YY]-a[YY]*b[XX];
}

__forceinline__ __device__
real iprod_gpu(const rvec a, const rvec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

__forceinline__ __device__
real norm_gpu(const rvec a)
{
    return sqrt(iprod_gpu(a, a));
}

__forceinline__ __device__
real gmx_angle_gpu(const rvec a, const rvec b)
{
    rvec w;
    real wlen, s;

    cprod_gpu(a, b, w);

    wlen  = norm_gpu(w);
    s     = iprod_gpu(a, b);

    return atan2f(wlen, s); //requires float
}

__forceinline__ __device__
void clear_ivec_gpu(ivec a)
{
    a[XX] = 0;
    a[YY] = 0;
    a[ZZ] = 0;
}
__forceinline__ __device__
void rvec_sub_gpu(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__
real norm2_gpu(const rvec a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

__forceinline__ __device__
void copy_rvec_gpu(const rvec a, rvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

__forceinline__ __device__
void copy_ivec_gpu(const ivec a, ivec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

__forceinline__ __device__
real cos_angle_gpu(const rvec a, const rvec b)
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
        cosval = ip*rsqrt(ipab);  /*  7 */ //double precision
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



__device__ static inline float invsqrt(float x)
{
  return 1.0f/std::sqrt(x);
}


__device__ static inline void unitv_gpu(const rvec src, rvec dest)
{
  real linv;

  linv     = invsqrt(norm2_gpu(src));
  dest[XX] = linv*src[XX];
  dest[YY] = linv*src[YY];
  dest[ZZ] = linv*src[ZZ];
}




__device__ static inline int pbc_dx_aiuc_gpu(const t_pbc *pbc, const rvec x1, const rvec x2, rvec dx, const rvec pbc_hbox_diag, const matrix pbc_box, const rvec pbc_mhbox_diag, const rvec  pbc_fbox_diag,  const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift  ){


  int  i, j, is;
  rvec dx_start, trial;
  real d2min, d2trial;
  ivec ishift, ishift_start;

  rvec_sub_gpu(x1, x2, dx);
  clear_ivec_gpu(ishift);


  switch (pbc->ePBCDX)
  {
  case epbcdxRECTANGULAR:
      for (i = 0; i < DIM; i++)
      {
          if (dx[i] > pbc_hbox_diag[i])
          {
              dx[i] -=  pbc_fbox_diag[i];
              ishift[i]--;
          }
          else if (dx[i] <= pbc_mhbox_diag[i])
          {
              dx[i] +=  pbc_fbox_diag[i];
              ishift[i]++;
          }
      }
      break;

  case epbcdxTRICLINIC:
      
  for (i = DIM-1; i >= 1; i--)
    {
      if (dx[i] > pbc_hbox_diag[i])
	{
	  for (j = i; j >= 0; j--)
	    {
	      dx[j] -= pbc_box[i][j];
	    }
	  ishift[i]--;
	}
      else if (dx[i] <= pbc_mhbox_diag[i])
	{
	  for (j = i; j >= 0; j--)
	    {
	      dx[j] += pbc_box[i][j];
	    }
	  ishift[i]++;
	}
    }

  /* Allow 2 shifts in x */
  if (dx[XX] > pbc_hbox_diag[XX])
    {
      dx[XX] -= pbc_fbox_diag[XX];
      ishift[XX]--;
      if (dx[XX] > pbc_hbox_diag[XX])
	{
	  dx[XX] -= pbc_fbox_diag[XX];
	  ishift[XX]--;
	}
    }
  else if (dx[XX] <= pbc_mhbox_diag[XX])
    {
      dx[XX] += pbc_fbox_diag[XX];
      ishift[XX]++;
      if (dx[XX] <= pbc_mhbox_diag[XX])
	{
	  dx[XX] += pbc_fbox_diag[XX];
	  ishift[XX]++;
	}
    }
 
  /* dx is the distance in a rectangular box */
  d2min = norm2_gpu(dx);
  if (d2min > pbc->max_cutoff2)
    {
      copy_rvec_gpu(dx, dx_start);
      copy_ivec_gpu(ishift, ishift_start);
      d2min = norm2_gpu(dx);
      /* Now try all possible shifts, when the distance is within max_cutoff
       * it must be the shortest possible distance.
       */
      i = 0;
      while ((d2min > pbc->max_cutoff2) && (i < pbc->ntric_vec))
	{
	  rvec_add_gpu(dx_start, pbc_tric_vec[i], trial);
	  d2trial = norm2_gpu(trial);
	   if (d2trial < d2min)
	     {
	            copy_rvec_gpu(trial, dx);
	       ivec_add_gpu(ishift_start, pbc_tric_shift[i], ishift);
	       d2min = d2trial;
	     }
	  i++;
	}
    }
  
  break;
  case epbcdx2D_RECT:
      for (i = 0; i < DIM; i++)
      {
          if (i != pbc->dim)
          {
              if (dx[i] > pbc_hbox_diag[i])
              {
                  dx[i] -= pbc_fbox_diag[i];
                  ishift[i]--;
              }
              else if (dx[i] <= pbc_mhbox_diag[i])
              {
                  dx[i] += pbc_fbox_diag[i];
                  ishift[i]++;
              }
          }
      }
      break;
      
  }
  is = IVEC2IS(ishift);
  return is;
}

