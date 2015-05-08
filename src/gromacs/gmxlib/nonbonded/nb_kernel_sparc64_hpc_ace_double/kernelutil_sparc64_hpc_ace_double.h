/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
#ifndef _kernelutil_sparc64_hpc_ace_double_h_
#define _kernelutil_sparc64_hpc_ace_double_h_

/* Get gmx_simd_exp_d() */
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"

/* Fujitsu header borrows the name from SSE2, since some instructions have aliases.
 * Environment/compiler version GM-1.2.0-17 seems to be buggy; when -Xg is
 * defined to enable GNUC extensions, this sets _ISOC99_SOURCE, which in
 * turn causes all intrinsics to be declared inline _instead_ of static. This
 * leads to duplicate symbol errors at link time.
 * To work around this we unset this before including the HPC-ACE header, and
 * reset the value afterwards.
 */
#ifdef _ISOC99_SOURCE
#    undef _ISOC99_SOURCE
#    define SAVE_ISOC99_SOURCE
#endif

#include <emmintrin.h>

#ifdef SAVE_ISOC99_SOURCE
#    define _ISOC99_SOURCE
#    undef SAVE_ISOC99_SOURCE
#endif

#define GMX_FJSP_SHUFFLE2(x, y) (((x)<<1) | (y))

#define GMX_FJSP_TRANSPOSE2_V2R8(row0, row1) {           \
        _fjsp_v2r8 __gmx_t1 = row0;                          \
        row0           = _fjsp_unpacklo_v2r8(row0, row1);     \
        row1           = _fjsp_unpackhi_v2r8(__gmx_t1, row1); \
}


static void
gmx_fjsp_print_v2r8(const char *s, _fjsp_v2r8 a)
{
    double lo, hi;

    _fjsp_storel_v2r8(&lo, a);
    _fjsp_storeh_v2r8(&hi, a);
    printf("%s: %g %g\n", s, lo, hi);
}


static _fjsp_v2r8
gmx_fjsp_set1_v2r8(double d)
{
    return _fjsp_set_v2r8(d, d);
}

static _fjsp_v2r8
gmx_fjsp_load1_v2r8(const double * gmx_restrict ptr)
{
    return gmx_fjsp_set1_v2r8(*ptr);
}


static int
gmx_fjsp_any_lt_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b)
{
    union
    {
        double           d;
        long long int    i;
    }
    conv;

    a = _fjsp_cmplt_v2r8(a, b);
    a = _fjsp_or_v2r8(a, _fjsp_unpackhi_v2r8(a, a));
    _fjsp_storel_v2r8(&(conv.d), a);
    return (conv.i != 0);
}

/* 1.0/sqrt(x) */
static gmx_inline _fjsp_v2r8
gmx_fjsp_invsqrt_v2r8(_fjsp_v2r8 x)
{
    const _fjsp_v2r8 half  = gmx_fjsp_set1_v2r8(0.5);
    const _fjsp_v2r8 three = gmx_fjsp_set1_v2r8(3.0);
    _fjsp_v2r8       lu    = _fjsp_rsqrta_v2r8(x);

    lu = _fjsp_mul_v2r8(_fjsp_mul_v2r8(half, lu), _fjsp_nmsub_v2r8(_fjsp_mul_v2r8(lu, lu), x, three));
    /* The HPC-ACE instruction set is only available in double precision, while
     * single precision is typically sufficient for Gromacs. If you define
     * "GMX_RELAXED_DOUBLE_PRECISION" during compile, we stick to two Newton-Raphson
     * iterations and accept 32bits of accuracy in 1.0/sqrt(x) and 1.0/x, rather than full
     * double precision (53 bits). This is still clearly higher than single precision (24 bits).
     */
#ifndef GMX_RELAXED_DOUBLE_PRECISION
    lu = _fjsp_mul_v2r8(_fjsp_mul_v2r8(half, lu), _fjsp_nmsub_v2r8(_fjsp_mul_v2r8(lu, lu), x, three));
#endif
    return _fjsp_mul_v2r8(_fjsp_mul_v2r8(half, lu), _fjsp_nmsub_v2r8(_fjsp_mul_v2r8(lu, lu), x, three));
}


/* 1.0/x */
static gmx_inline _fjsp_v2r8
gmx_fjsp_inv_v2r8(_fjsp_v2r8 x)
{
    const _fjsp_v2r8 two  = gmx_fjsp_set1_v2r8(2.0);
    __m128d          lu   = _fjsp_rcpa_v2r8(x);

    /* Perform three N-R steps for double precision */
    lu         = _fjsp_mul_v2r8(lu, _fjsp_nmsub_v2r8(lu, x, two));
    /* The HPC-ACE instruction set is only available in double precision, while
     * single precision is typically sufficient for Gromacs. If you define
     * "GMX_RELAXED_DOUBLE_PRECISION" during compile, we stick to two Newton-Raphson
     * iterations and accept 32bits of accuracy in 1.0/sqrt(x) and 1.0/x, rather than full
     * double precision (53 bits). This is still clearly higher than single precision (24 bits).
     */
#ifndef GMX_RELAXED_DOUBLE_PRECISION
    lu         = _fjsp_mul_v2r8(lu, _fjsp_nmsub_v2r8(lu, x, two));
#endif
    return _fjsp_mul_v2r8(lu, _fjsp_nmsub_v2r8(lu, x, two));
}


static gmx_inline _fjsp_v2r8
gmx_fjsp_calc_rsq_v2r8(_fjsp_v2r8 dx, _fjsp_v2r8 dy, _fjsp_v2r8 dz)
{
    return _fjsp_madd_v2r8(dx, dx, _fjsp_madd_v2r8(dy, dy, _fjsp_mul_v2r8(dz, dz)));
}

/* Normal sum of four ymm registers */
#define gmx_fjsp_sum4_v2r8(t0, t1, t2, t3)  _fjsp_add_v2r8(_fjsp_add_v2r8(t0, t1), _fjsp_add_v2r8(t2, t3))





static _fjsp_v2r8
gmx_fjsp_load_2real_swizzle_v2r8(const double * gmx_restrict ptrA,
                                 const double * gmx_restrict ptrB)
{
    return _fjsp_unpacklo_v2r8(_fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA), _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB));
}

static _fjsp_v2r8
gmx_fjsp_load_1real_v2r8(const double * gmx_restrict ptrA)
{
    return _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA);
}


static void
gmx_fjsp_store_2real_swizzle_v2r8(double * gmx_restrict ptrA,
                                  double * gmx_restrict ptrB,
                                  _fjsp_v2r8            xmm1)
{
    _fjsp_v2r8 t2;

    t2       = _fjsp_unpackhi_v2r8(xmm1, xmm1);
    _fjsp_storel_v2r8(ptrA, xmm1);
    _fjsp_storel_v2r8(ptrB, t2);
}

static void
gmx_fjsp_store_1real_v2r8(double * gmx_restrict ptrA, _fjsp_v2r8 xmm1)
{
    _fjsp_storel_v2r8(ptrA, xmm1);
}


/* Similar to store, but increments value in memory */
static void
gmx_fjsp_increment_2real_swizzle_v2r8(double * gmx_restrict ptrA,
                                      double * gmx_restrict ptrB, _fjsp_v2r8 xmm1)
{
    _fjsp_v2r8 t1;

    t1   = _fjsp_unpackhi_v2r8(xmm1, xmm1);
    xmm1 = _fjsp_add_v2r8(xmm1, _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA));
    t1   = _fjsp_add_v2r8(t1, _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB));
    _fjsp_storel_v2r8(ptrA, xmm1);
    _fjsp_storel_v2r8(ptrB, t1);
}

static void
gmx_fjsp_increment_1real_v2r8(double * gmx_restrict ptrA, _fjsp_v2r8 xmm1)
{
    _fjsp_v2r8 tmp;

    tmp = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA);
    tmp = _fjsp_add_v2r8(tmp, xmm1);
    _fjsp_storel_v2r8(ptrA, tmp);
}



static gmx_inline void
gmx_fjsp_load_2pair_swizzle_v2r8(const double * gmx_restrict p1,
                                 const double * gmx_restrict p2,
                                 _fjsp_v2r8 * gmx_restrict   c6,
                                 _fjsp_v2r8 * gmx_restrict   c12)
{
    _fjsp_v2r8 t1, t2, t3;

    /* The c6/c12 array should be aligned */
    t1   = _fjsp_load_v2r8(p1);
    t2   = _fjsp_load_v2r8(p2);
    *c6  = _fjsp_unpacklo_v2r8(t1, t2);
    *c12 = _fjsp_unpackhi_v2r8(t1, t2);
}

static gmx_inline void
gmx_fjsp_load_1pair_swizzle_v2r8(const double * gmx_restrict p1,
                                 _fjsp_v2r8 * gmx_restrict   c6,
                                 _fjsp_v2r8 * gmx_restrict   c12)
{
    *c6     = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1);
    *c12    = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+1);
}


static gmx_inline void
gmx_fjsp_load_shift_and_1rvec_broadcast_v2r8(const double * gmx_restrict xyz_shift,
                                             const double * gmx_restrict xyz,
                                             _fjsp_v2r8 * gmx_restrict   x1,
                                             _fjsp_v2r8 * gmx_restrict   y1,
                                             _fjsp_v2r8 * gmx_restrict   z1)
{
    _fjsp_v2r8 mem_xy, mem_z, mem_sxy, mem_sz;

    mem_xy  = _fjsp_load_v2r8(xyz);
    mem_z   = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), xyz+2);
    mem_sxy = _fjsp_load_v2r8(xyz_shift);
    mem_sz  = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), xyz_shift+2);

    mem_xy  = _fjsp_add_v2r8(mem_xy, mem_sxy);
    mem_z   = _fjsp_add_v2r8(mem_z, mem_sz);

    *x1  = _fjsp_shuffle_v2r8(mem_xy, mem_xy, GMX_FJSP_SHUFFLE2(0, 0));
    *y1  = _fjsp_shuffle_v2r8(mem_xy, mem_xy, GMX_FJSP_SHUFFLE2(1, 1));
    *z1  = _fjsp_shuffle_v2r8(mem_z, mem_z, GMX_FJSP_SHUFFLE2(0, 0));
}


static gmx_inline void
gmx_fjsp_load_shift_and_3rvec_broadcast_v2r8(const double * gmx_restrict xyz_shift,
                                             const double * gmx_restrict xyz,
                                             _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1,
                                             _fjsp_v2r8 * gmx_restrict x2, _fjsp_v2r8 * gmx_restrict y2, _fjsp_v2r8 * gmx_restrict z2,
                                             _fjsp_v2r8 * gmx_restrict x3, _fjsp_v2r8 * gmx_restrict y3, _fjsp_v2r8 * gmx_restrict z3)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, sxy, sz, szx, syz;

    t1  = _fjsp_load_v2r8(xyz);
    t2  = _fjsp_load_v2r8(xyz+2);
    t3  = _fjsp_load_v2r8(xyz+4);
    t4  = _fjsp_load_v2r8(xyz+6);
    t5  = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), xyz+8);

    sxy = _fjsp_load_v2r8(xyz_shift);
    sz  = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), xyz_shift+2);
    szx = _fjsp_shuffle_v2r8(sz, sxy, GMX_FJSP_SHUFFLE2(0, 0));
    syz = _fjsp_shuffle_v2r8(sxy, sz, GMX_FJSP_SHUFFLE2(0, 1));

    t1  = _fjsp_add_v2r8(t1, sxy);
    t2  = _fjsp_add_v2r8(t2, szx);
    t3  = _fjsp_add_v2r8(t3, syz);
    t4  = _fjsp_add_v2r8(t4, sxy);
    t5  = _fjsp_add_v2r8(t5, sz);

    *x1  = _fjsp_shuffle_v2r8(t1, t1, GMX_FJSP_SHUFFLE2(0, 0));
    *y1  = _fjsp_shuffle_v2r8(t1, t1, GMX_FJSP_SHUFFLE2(1, 1));
    *z1  = _fjsp_shuffle_v2r8(t2, t2, GMX_FJSP_SHUFFLE2(0, 0));
    *x2  = _fjsp_shuffle_v2r8(t2, t2, GMX_FJSP_SHUFFLE2(1, 1));
    *y2  = _fjsp_shuffle_v2r8(t3, t3, GMX_FJSP_SHUFFLE2(0, 0));
    *z2  = _fjsp_shuffle_v2r8(t3, t3, GMX_FJSP_SHUFFLE2(1, 1));
    *x3  = _fjsp_shuffle_v2r8(t4, t4, GMX_FJSP_SHUFFLE2(0, 0));
    *y3  = _fjsp_shuffle_v2r8(t4, t4, GMX_FJSP_SHUFFLE2(1, 1));
    *z3  = _fjsp_shuffle_v2r8(t5, t5, GMX_FJSP_SHUFFLE2(0, 0));
}


static gmx_inline void
gmx_fjsp_load_shift_and_4rvec_broadcast_v2r8(const double * gmx_restrict xyz_shift,
                                             const double * gmx_restrict xyz,
                                             _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1,
                                             _fjsp_v2r8 * gmx_restrict x2, _fjsp_v2r8 * gmx_restrict y2, _fjsp_v2r8 * gmx_restrict z2,
                                             _fjsp_v2r8 * gmx_restrict x3, _fjsp_v2r8 * gmx_restrict y3, _fjsp_v2r8 * gmx_restrict z3,
                                             _fjsp_v2r8 * gmx_restrict x4, _fjsp_v2r8 * gmx_restrict y4, _fjsp_v2r8 * gmx_restrict z4)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, sxy, sz, szx, syz;

    t1  = _fjsp_load_v2r8(xyz);
    t2  = _fjsp_load_v2r8(xyz+2);
    t3  = _fjsp_load_v2r8(xyz+4);
    t4  = _fjsp_load_v2r8(xyz+6);
    t5  = _fjsp_load_v2r8(xyz+8);
    t6  = _fjsp_load_v2r8(xyz+10);

    sxy = _fjsp_load_v2r8(xyz_shift);
    sz  = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), xyz_shift+2);
    szx = _fjsp_shuffle_v2r8(sz, sxy, GMX_FJSP_SHUFFLE2(0, 0));
    syz = _fjsp_shuffle_v2r8(sxy, sz, GMX_FJSP_SHUFFLE2(0, 1));

    t1  = _fjsp_add_v2r8(t1, sxy);
    t2  = _fjsp_add_v2r8(t2, szx);
    t3  = _fjsp_add_v2r8(t3, syz);
    t4  = _fjsp_add_v2r8(t4, sxy);
    t5  = _fjsp_add_v2r8(t5, szx);
    t6  = _fjsp_add_v2r8(t6, syz);

    *x1  = _fjsp_shuffle_v2r8(t1, t1, GMX_FJSP_SHUFFLE2(0, 0));
    *y1  = _fjsp_shuffle_v2r8(t1, t1, GMX_FJSP_SHUFFLE2(1, 1));
    *z1  = _fjsp_shuffle_v2r8(t2, t2, GMX_FJSP_SHUFFLE2(0, 0));
    *x2  = _fjsp_shuffle_v2r8(t2, t2, GMX_FJSP_SHUFFLE2(1, 1));
    *y2  = _fjsp_shuffle_v2r8(t3, t3, GMX_FJSP_SHUFFLE2(0, 0));
    *z2  = _fjsp_shuffle_v2r8(t3, t3, GMX_FJSP_SHUFFLE2(1, 1));
    *x3  = _fjsp_shuffle_v2r8(t4, t4, GMX_FJSP_SHUFFLE2(0, 0));
    *y3  = _fjsp_shuffle_v2r8(t4, t4, GMX_FJSP_SHUFFLE2(1, 1));
    *z3  = _fjsp_shuffle_v2r8(t5, t5, GMX_FJSP_SHUFFLE2(0, 0));
    *x4  = _fjsp_shuffle_v2r8(t5, t5, GMX_FJSP_SHUFFLE2(1, 1));
    *y4  = _fjsp_shuffle_v2r8(t6, t6, GMX_FJSP_SHUFFLE2(0, 0));
    *z4  = _fjsp_shuffle_v2r8(t6, t6, GMX_FJSP_SHUFFLE2(1, 1));
}



static gmx_inline void
gmx_fjsp_load_1rvec_1ptr_swizzle_v2r8(const double * gmx_restrict p1,
                                      _fjsp_v2r8 * gmx_restrict x, _fjsp_v2r8 * gmx_restrict y, _fjsp_v2r8 * gmx_restrict z)
{
    *x            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1);
    *y            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+1);
    *z            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+2);
}

static gmx_inline void
gmx_fjsp_load_3rvec_1ptr_swizzle_v2r8(const double * gmx_restrict p1,
                                      _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1,
                                      _fjsp_v2r8 * gmx_restrict x2, _fjsp_v2r8 * gmx_restrict y2, _fjsp_v2r8 * gmx_restrict z2,
                                      _fjsp_v2r8 * gmx_restrict x3, _fjsp_v2r8 * gmx_restrict y3, _fjsp_v2r8 * gmx_restrict z3)
{
    *x1            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1);
    *y1            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+1);
    *z1            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+2);
    *x2            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+3);
    *y2            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+4);
    *z2            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+5);
    *x3            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+6);
    *y3            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+7);
    *z3            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+8);
}

static gmx_inline void
gmx_fjsp_load_4rvec_1ptr_swizzle_v2r8(const double * gmx_restrict p1,
                                      _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1,
                                      _fjsp_v2r8 * gmx_restrict x2, _fjsp_v2r8 * gmx_restrict y2, _fjsp_v2r8 * gmx_restrict z2,
                                      _fjsp_v2r8 * gmx_restrict x3, _fjsp_v2r8 * gmx_restrict y3, _fjsp_v2r8 * gmx_restrict z3,
                                      _fjsp_v2r8 * gmx_restrict x4, _fjsp_v2r8 * gmx_restrict y4, _fjsp_v2r8 * gmx_restrict z4)
{
    *x1            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1);
    *y1            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+1);
    *z1            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+2);
    *x2            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+3);
    *y2            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+4);
    *z2            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+5);
    *x3            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+6);
    *y3            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+7);
    *z3            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+8);
    *x4            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+9);
    *y4            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+10);
    *z4            = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), p1+11);
}


static gmx_inline void
gmx_fjsp_load_1rvec_2ptr_swizzle_v2r8(const double * gmx_restrict ptrA,
                                      const double * gmx_restrict ptrB,
                                      _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    t1           = _fjsp_load_v2r8(ptrA);
    t2           = _fjsp_load_v2r8(ptrB);
    t3           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+2);
    t4           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB+2);
    GMX_FJSP_TRANSPOSE2_V2R8(t1, t2);
    *x1          = t1;
    *y1          = t2;
    *z1          = _fjsp_unpacklo_v2r8(t3, t4);
}

static gmx_inline void
gmx_fjsp_load_3rvec_2ptr_swizzle_v2r8(const double * gmx_restrict ptrA, const double * gmx_restrict ptrB,
                                      _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1,
                                      _fjsp_v2r8 * gmx_restrict x2, _fjsp_v2r8 * gmx_restrict y2, _fjsp_v2r8 * gmx_restrict z2,
                                      _fjsp_v2r8 * gmx_restrict x3, _fjsp_v2r8 * gmx_restrict y3, _fjsp_v2r8 * gmx_restrict z3)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    t1           = _fjsp_load_v2r8(ptrA);
    t2           = _fjsp_load_v2r8(ptrB);
    t3           = _fjsp_load_v2r8(ptrA+2);
    t4           = _fjsp_load_v2r8(ptrB+2);
    t5           = _fjsp_load_v2r8(ptrA+4);
    t6           = _fjsp_load_v2r8(ptrB+4);
    t7           = _fjsp_load_v2r8(ptrA+6);
    t8           = _fjsp_load_v2r8(ptrB+6);
    t9           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+8);
    t10          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB+8);
    GMX_FJSP_TRANSPOSE2_V2R8(t1, t2);
    GMX_FJSP_TRANSPOSE2_V2R8(t3, t4);
    GMX_FJSP_TRANSPOSE2_V2R8(t5, t6);
    GMX_FJSP_TRANSPOSE2_V2R8(t7, t8);
    *x1          = t1;
    *y1          = t2;
    *z1          = t3;
    *x2          = t4;
    *y2          = t5;
    *z2          = t6;
    *x3          = t7;
    *y3          = t8;
    *z3          = _fjsp_unpacklo_v2r8(t9, t10);
}


static gmx_inline void
gmx_fjsp_load_4rvec_2ptr_swizzle_v2r8(const double * gmx_restrict ptrA, const double * gmx_restrict ptrB,
                                      _fjsp_v2r8 * gmx_restrict x1, _fjsp_v2r8 * gmx_restrict y1, _fjsp_v2r8 * gmx_restrict z1,
                                      _fjsp_v2r8 * gmx_restrict x2, _fjsp_v2r8 * gmx_restrict y2, _fjsp_v2r8 * gmx_restrict z2,
                                      _fjsp_v2r8 * gmx_restrict x3, _fjsp_v2r8 * gmx_restrict y3, _fjsp_v2r8 * gmx_restrict z3,
                                      _fjsp_v2r8 * gmx_restrict x4, _fjsp_v2r8 * gmx_restrict y4, _fjsp_v2r8 * gmx_restrict z4)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6;
    t1           = _fjsp_load_v2r8(ptrA);
    t2           = _fjsp_load_v2r8(ptrB);
    t3           = _fjsp_load_v2r8(ptrA+2);
    t4           = _fjsp_load_v2r8(ptrB+2);
    t5           = _fjsp_load_v2r8(ptrA+4);
    t6           = _fjsp_load_v2r8(ptrB+4);
    GMX_FJSP_TRANSPOSE2_V2R8(t1, t2);
    GMX_FJSP_TRANSPOSE2_V2R8(t3, t4);
    GMX_FJSP_TRANSPOSE2_V2R8(t5, t6);
    *x1          = t1;
    *y1          = t2;
    *z1          = t3;
    *x2          = t4;
    *y2          = t5;
    *z2          = t6;
    t1           = _fjsp_load_v2r8(ptrA+6);
    t2           = _fjsp_load_v2r8(ptrB+6);
    t3           = _fjsp_load_v2r8(ptrA+8);
    t4           = _fjsp_load_v2r8(ptrB+8);
    t5           = _fjsp_load_v2r8(ptrA+10);
    t6           = _fjsp_load_v2r8(ptrB+10);
    GMX_FJSP_TRANSPOSE2_V2R8(t1, t2);
    GMX_FJSP_TRANSPOSE2_V2R8(t3, t4);
    GMX_FJSP_TRANSPOSE2_V2R8(t5, t6);
    *x3          = t1;
    *y3          = t2;
    *z3          = t3;
    *x4          = t4;
    *y4          = t5;
    *z4          = t6;
}


static void
gmx_fjsp_decrement_1rvec_1ptr_swizzle_v2r8(double * gmx_restrict ptrA,
                                           _fjsp_v2r8 x1, _fjsp_v2r8 y1, _fjsp_v2r8 z1)
{
    _fjsp_v2r8 t1, t2, t3;

    t1           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA);
    t2           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+1);
    t3           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+2);

    t1           = _fjsp_sub_v2r8(t1, x1);
    t2           = _fjsp_sub_v2r8(t2, y1);
    t3           = _fjsp_sub_v2r8(t3, z1);
    _fjsp_storel_v2r8(ptrA, t1);
    _fjsp_storel_v2r8(ptrA+1, t2);
    _fjsp_storel_v2r8(ptrA+2, t3);
}

static void
gmx_fjsp_decrement_fma_1rvec_1ptr_swizzle_v2r8(double * gmx_restrict ptrA, _fjsp_v2r8 fscal,
                                               _fjsp_v2r8 dx1, _fjsp_v2r8 dy1, _fjsp_v2r8 dz1)
{
    _fjsp_v2r8 t1, t2, t3;

    t1           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA);
    t2           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+1);
    t3           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+2);

    t1           = _fjsp_nmsub_v2r8(fscal, dx1, t1);
    t2           = _fjsp_nmsub_v2r8(fscal, dy1, t2);
    t3           = _fjsp_nmsub_v2r8(fscal, dz1, t3);
    _fjsp_storel_v2r8(ptrA, t1);
    _fjsp_storel_v2r8(ptrA+1, t2);
    _fjsp_storel_v2r8(ptrA+2, t3);
}


static void
gmx_fjsp_decrement_3rvec_1ptr_swizzle_v2r8(double * gmx_restrict ptrA,
                                           _fjsp_v2r8 x1, _fjsp_v2r8 y1, _fjsp_v2r8 z1,
                                           _fjsp_v2r8 x2, _fjsp_v2r8 y2, _fjsp_v2r8 z2,
                                           _fjsp_v2r8 x3, _fjsp_v2r8 y3, _fjsp_v2r8 z3)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5;

    t1          = _fjsp_load_v2r8(ptrA);
    t2          = _fjsp_load_v2r8(ptrA+2);
    t3          = _fjsp_load_v2r8(ptrA+4);
    t4          = _fjsp_load_v2r8(ptrA+6);
    t5          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+8);

    x1          = _fjsp_unpacklo_v2r8(x1, y1);
    z1          = _fjsp_unpacklo_v2r8(z1, x2);
    y2          = _fjsp_unpacklo_v2r8(y2, z2);
    x3          = _fjsp_unpacklo_v2r8(x3, y3);
    /* nothing to be done for z3 */

    t1          = _fjsp_sub_v2r8(t1, x1);
    t2          = _fjsp_sub_v2r8(t2, z1);
    t3          = _fjsp_sub_v2r8(t3, y2);
    t4          = _fjsp_sub_v2r8(t4, x3);
    t5          = _fjsp_sub_v2r8(t5, z3);
    _fjsp_storel_v2r8(ptrA, t1);
    _fjsp_storeh_v2r8(ptrA+1, t1);
    _fjsp_storel_v2r8(ptrA+2, t2);
    _fjsp_storeh_v2r8(ptrA+3, t2);
    _fjsp_storel_v2r8(ptrA+4, t3);
    _fjsp_storeh_v2r8(ptrA+5, t3);
    _fjsp_storel_v2r8(ptrA+6, t4);
    _fjsp_storeh_v2r8(ptrA+7, t4);
    _fjsp_storel_v2r8(ptrA+8, t5);
}


static void
gmx_fjsp_decrement_4rvec_1ptr_swizzle_v2r8(double * gmx_restrict ptrA,
                                           _fjsp_v2r8 x1, _fjsp_v2r8 y1, _fjsp_v2r8 z1,
                                           _fjsp_v2r8 x2, _fjsp_v2r8 y2, _fjsp_v2r8 z2,
                                           _fjsp_v2r8 x3, _fjsp_v2r8 y3, _fjsp_v2r8 z3,
                                           _fjsp_v2r8 x4, _fjsp_v2r8 y4, _fjsp_v2r8 z4)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6;

    t1          = _fjsp_load_v2r8(ptrA);
    t2          = _fjsp_load_v2r8(ptrA+2);
    t3          = _fjsp_load_v2r8(ptrA+4);
    t4          = _fjsp_load_v2r8(ptrA+6);
    t5          = _fjsp_load_v2r8(ptrA+8);
    t6          = _fjsp_load_v2r8(ptrA+10);

    x1          = _fjsp_unpacklo_v2r8(x1, y1);
    z1          = _fjsp_unpacklo_v2r8(z1, x2);
    y2          = _fjsp_unpacklo_v2r8(y2, z2);
    x3          = _fjsp_unpacklo_v2r8(x3, y3);
    z3          = _fjsp_unpacklo_v2r8(z3, x4);
    y4          = _fjsp_unpacklo_v2r8(y4, z4);

    _fjsp_storel_v2r8(ptrA,    _fjsp_sub_v2r8( t1, x1 ));
    _fjsp_storeh_v2r8(ptrA+1,  _fjsp_sub_v2r8( t1, x1 ));
    _fjsp_storel_v2r8(ptrA+2,  _fjsp_sub_v2r8( t2, z1 ));
    _fjsp_storeh_v2r8(ptrA+3,  _fjsp_sub_v2r8( t2, z1 ));
    _fjsp_storel_v2r8(ptrA+4,  _fjsp_sub_v2r8( t3, y2 ));
    _fjsp_storeh_v2r8(ptrA+5,  _fjsp_sub_v2r8( t3, y2 ));
    _fjsp_storel_v2r8(ptrA+6,  _fjsp_sub_v2r8( t4, x3 ));
    _fjsp_storeh_v2r8(ptrA+7,  _fjsp_sub_v2r8( t4, x3 ));
    _fjsp_storel_v2r8(ptrA+8,  _fjsp_sub_v2r8( t5, z3 ));
    _fjsp_storeh_v2r8(ptrA+9,  _fjsp_sub_v2r8( t5, z3 ));
    _fjsp_storel_v2r8(ptrA+10, _fjsp_sub_v2r8( t6, y4 ));
    _fjsp_storeh_v2r8(ptrA+11, _fjsp_sub_v2r8( t6, y4 ));
}

static void
gmx_fjsp_decrement_1rvec_2ptr_swizzle_v2r8(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                           _fjsp_v2r8 x1, _fjsp_v2r8 y1, _fjsp_v2r8 z1)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7;

    t1          = _fjsp_load_v2r8(ptrA);
    t2          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+2);
    t3          = _fjsp_load_v2r8(ptrB);
    t4          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB+2);

    t5          = _fjsp_unpacklo_v2r8(x1, y1);
    t6          = _fjsp_unpackhi_v2r8(x1, y1);
    t7          = _fjsp_unpackhi_v2r8(z1, z1);

    t1          = _fjsp_sub_v2r8(t1, t5);
    t2          = _fjsp_sub_v2r8(t2, z1);

    t3          = _fjsp_sub_v2r8(t3, t6);
    t4          = _fjsp_sub_v2r8(t4, t7);

    _fjsp_storel_v2r8(ptrA, t1);
    _fjsp_storeh_v2r8(ptrA+1, t1);
    _fjsp_storel_v2r8(ptrA+2, t2);
    _fjsp_storel_v2r8(ptrB, t3);
    _fjsp_storeh_v2r8(ptrB+1, t3);
    _fjsp_storel_v2r8(ptrB+2, t4);
}


static void
gmx_fjsp_decrement_fma_1rvec_2ptr_swizzle_v2r8(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                               _fjsp_v2r8 fscal, _fjsp_v2r8 dx1, _fjsp_v2r8 dy1, _fjsp_v2r8 dz1)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7, fscalA, fscalB;

    t1          = _fjsp_load_v2r8(ptrA);
    t2          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+2);
    t3          = _fjsp_load_v2r8(ptrB);
    t4          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB+2);
    fscalA      = _fjsp_unpacklo_v2r8(fscal, fscal);
    fscalB      = _fjsp_unpackhi_v2r8(fscal, fscal);

    t5          = _fjsp_unpacklo_v2r8(dx1, dy1);
    t6          = _fjsp_unpackhi_v2r8(dx1, dy1);
    t7          = _fjsp_unpackhi_v2r8(dz1, dz1);

    t1          = _fjsp_nmsub_v2r8(fscalA, t5, t1);
    t2          = _fjsp_nmsub_v2r8(fscalA, dz1, t2);

    t3          = _fjsp_nmsub_v2r8(fscalB, t6, t3);
    t4          = _fjsp_nmsub_v2r8(fscalB, t7, t4);

    _fjsp_storel_v2r8(ptrA, t1);
    _fjsp_storeh_v2r8(ptrA+1, t1);
    _fjsp_storel_v2r8(ptrA+2, t2);
    _fjsp_storel_v2r8(ptrB, t3);
    _fjsp_storeh_v2r8(ptrB+1, t3);
    _fjsp_storel_v2r8(ptrB+2, t4);
}


static void
gmx_fjsp_decrement_3rvec_2ptr_swizzle_v2r8(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                           _fjsp_v2r8 x1, _fjsp_v2r8 y1, _fjsp_v2r8 z1,
                                           _fjsp_v2r8 x2, _fjsp_v2r8 y2, _fjsp_v2r8 z2,
                                           _fjsp_v2r8 x3, _fjsp_v2r8 y3, _fjsp_v2r8 z3)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    _fjsp_v2r8 tA, tB, tC, tD, tE, tF, tG, tH, tI;

    t1          = _fjsp_load_v2r8(ptrA);
    t2          = _fjsp_load_v2r8(ptrA+2);
    t3          = _fjsp_load_v2r8(ptrA+4);
    t4          = _fjsp_load_v2r8(ptrA+6);
    t5          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA+8);
    t6          = _fjsp_load_v2r8(ptrB);
    t7          = _fjsp_load_v2r8(ptrB+2);
    t8          = _fjsp_load_v2r8(ptrB+4);
    t9          = _fjsp_load_v2r8(ptrB+6);
    t10         = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB+8);

    tA          = _fjsp_unpacklo_v2r8(x1, y1);
    tB          = _fjsp_unpackhi_v2r8(x1, y1);
    tC          = _fjsp_unpacklo_v2r8(z1, x2);
    tD          = _fjsp_unpackhi_v2r8(z1, x2);
    tE          = _fjsp_unpacklo_v2r8(y2, z2);
    tF          = _fjsp_unpackhi_v2r8(y2, z2);
    tG          = _fjsp_unpacklo_v2r8(x3, y3);
    tH          = _fjsp_unpackhi_v2r8(x3, y3);
    tI          = _fjsp_unpackhi_v2r8(z3, z3);

    t1          = _fjsp_sub_v2r8(t1, tA);
    t2          = _fjsp_sub_v2r8(t2, tC);
    t3          = _fjsp_sub_v2r8(t3, tE);
    t4          = _fjsp_sub_v2r8(t4, tG);
    t5          = _fjsp_sub_v2r8(t5, z3);

    t6          = _fjsp_sub_v2r8(t6, tB);
    t7          = _fjsp_sub_v2r8(t7, tD);
    t8          = _fjsp_sub_v2r8(t8, tF);
    t9          = _fjsp_sub_v2r8(t9, tH);
    t10         = _fjsp_sub_v2r8(t10, tI);

    _fjsp_storel_v2r8(ptrA, t1);
    _fjsp_storeh_v2r8(ptrA+1, t1);
    _fjsp_storel_v2r8(ptrA+2, t2);
    _fjsp_storeh_v2r8(ptrA+3, t2);
    _fjsp_storel_v2r8(ptrA+4, t3);
    _fjsp_storeh_v2r8(ptrA+5, t3);
    _fjsp_storel_v2r8(ptrA+6, t4);
    _fjsp_storeh_v2r8(ptrA+7, t4);
    _fjsp_storel_v2r8(ptrA+8, t5);
    _fjsp_storel_v2r8(ptrB, t6);
    _fjsp_storeh_v2r8(ptrB+1, t6);
    _fjsp_storel_v2r8(ptrB+2, t7);
    _fjsp_storeh_v2r8(ptrB+3, t7);
    _fjsp_storel_v2r8(ptrB+4, t8);
    _fjsp_storeh_v2r8(ptrB+5, t8);
    _fjsp_storel_v2r8(ptrB+6, t9);
    _fjsp_storeh_v2r8(ptrB+7, t9);
    _fjsp_storel_v2r8(ptrB+8, t10);
}


static void
gmx_fjsp_decrement_4rvec_2ptr_swizzle_v2r8(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                           _fjsp_v2r8 x1, _fjsp_v2r8 y1, _fjsp_v2r8 z1,
                                           _fjsp_v2r8 x2, _fjsp_v2r8 y2, _fjsp_v2r8 z2,
                                           _fjsp_v2r8 x3, _fjsp_v2r8 y3, _fjsp_v2r8 z3,
                                           _fjsp_v2r8 x4, _fjsp_v2r8 y4, _fjsp_v2r8 z4)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
    _fjsp_v2r8 tA, tB, tC, tD, tE, tF, tG, tH, tI, tJ, tK, tL;

    t1          = _fjsp_load_v2r8(ptrA);
    t2          = _fjsp_load_v2r8(ptrA+2);
    t3          = _fjsp_load_v2r8(ptrA+4);
    t4          = _fjsp_load_v2r8(ptrA+6);
    t5          = _fjsp_load_v2r8(ptrA+8);
    t6          = _fjsp_load_v2r8(ptrA+10);
    t7          = _fjsp_load_v2r8(ptrB);
    t8          = _fjsp_load_v2r8(ptrB+2);
    t9          = _fjsp_load_v2r8(ptrB+4);
    t10         = _fjsp_load_v2r8(ptrB+6);
    t11         = _fjsp_load_v2r8(ptrB+8);
    t12         = _fjsp_load_v2r8(ptrB+10);

    tA          = _fjsp_unpacklo_v2r8(x1, y1);
    tB          = _fjsp_unpackhi_v2r8(x1, y1);
    tC          = _fjsp_unpacklo_v2r8(z1, x2);
    tD          = _fjsp_unpackhi_v2r8(z1, x2);
    tE          = _fjsp_unpacklo_v2r8(y2, z2);
    tF          = _fjsp_unpackhi_v2r8(y2, z2);
    tG          = _fjsp_unpacklo_v2r8(x3, y3);
    tH          = _fjsp_unpackhi_v2r8(x3, y3);
    tI          = _fjsp_unpacklo_v2r8(z3, x4);
    tJ          = _fjsp_unpackhi_v2r8(z3, x4);
    tK          = _fjsp_unpacklo_v2r8(y4, z4);
    tL          = _fjsp_unpackhi_v2r8(y4, z4);

    t1          = _fjsp_sub_v2r8(t1, tA);
    t2          = _fjsp_sub_v2r8(t2, tC);
    t3          = _fjsp_sub_v2r8(t3, tE);
    t4          = _fjsp_sub_v2r8(t4, tG);
    t5          = _fjsp_sub_v2r8(t5, tI);
    t6          = _fjsp_sub_v2r8(t6, tK);

    t7          = _fjsp_sub_v2r8(t7, tB);
    t8          = _fjsp_sub_v2r8(t8, tD);
    t9          = _fjsp_sub_v2r8(t9, tF);
    t10         = _fjsp_sub_v2r8(t10, tH);
    t11         = _fjsp_sub_v2r8(t11, tJ);
    t12         = _fjsp_sub_v2r8(t12, tL);

    _fjsp_storel_v2r8(ptrA,  t1);
    _fjsp_storeh_v2r8(ptrA+1, t1);
    _fjsp_storel_v2r8(ptrA+2, t2);
    _fjsp_storeh_v2r8(ptrA+3, t2);
    _fjsp_storel_v2r8(ptrA+4, t3);
    _fjsp_storeh_v2r8(ptrA+5, t3);
    _fjsp_storel_v2r8(ptrA+6, t4);
    _fjsp_storeh_v2r8(ptrA+7, t4);
    _fjsp_storel_v2r8(ptrA+8, t5);
    _fjsp_storeh_v2r8(ptrA+9, t5);
    _fjsp_storel_v2r8(ptrA+10, t6);
    _fjsp_storeh_v2r8(ptrA+11, t6);
    _fjsp_storel_v2r8(ptrB,  t7);
    _fjsp_storeh_v2r8(ptrB+1, t7);
    _fjsp_storel_v2r8(ptrB+2, t8);
    _fjsp_storeh_v2r8(ptrB+3, t8);
    _fjsp_storel_v2r8(ptrB+4, t9);
    _fjsp_storeh_v2r8(ptrB+5, t9);
    _fjsp_storel_v2r8(ptrB+6, t10);
    _fjsp_storeh_v2r8(ptrB+7, t10);
    _fjsp_storel_v2r8(ptrB+8, t11);
    _fjsp_storeh_v2r8(ptrB+9, t11);
    _fjsp_storel_v2r8(ptrB+10, t12);
    _fjsp_storeh_v2r8(ptrB+11, t12);
}



static gmx_inline void
gmx_fjsp_update_iforce_1atom_swizzle_v2r8(_fjsp_v2r8 fix1, _fjsp_v2r8 fiy1, _fjsp_v2r8 fiz1,
                                          double * gmx_restrict fptr,
                                          double * gmx_restrict fshiftptr)
{
    __m128d t1, t2, t3, t4;

    /* transpose data */
    t1   = fix1;
    fix1 = _fjsp_unpacklo_v2r8(fix1, fiy1); /* y0 x0 */
    fiy1 = _fjsp_unpackhi_v2r8(t1, fiy1);   /* y1 x1 */

    fix1 = _fjsp_add_v2r8(fix1, fiy1);
    fiz1 = _fjsp_add_v2r8( fiz1, _fjsp_unpackhi_v2r8(fiz1, fiz1 ));

    t4 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr), fix1 );
    _fjsp_storel_v2r8( fptr, t4 );
    _fjsp_storeh_v2r8( fptr+1, t4 );
    _fjsp_storel_v2r8( fptr+2, _fjsp_add_v2r8( _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), fptr+2), fiz1 ));

    t4 = _fjsp_add_v2r8( _fjsp_load_v2r8(fshiftptr), fix1 );
    _fjsp_storel_v2r8( fshiftptr, t4 );
    _fjsp_storeh_v2r8( fshiftptr+1, t4 );
    _fjsp_storel_v2r8( fshiftptr+2, _fjsp_add_v2r8( _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), fshiftptr+2), fiz1 ));
}

static gmx_inline void
gmx_fjsp_update_iforce_3atom_swizzle_v2r8(_fjsp_v2r8 fix1, _fjsp_v2r8 fiy1, _fjsp_v2r8 fiz1,
                                          _fjsp_v2r8 fix2, _fjsp_v2r8 fiy2, _fjsp_v2r8 fiz2,
                                          _fjsp_v2r8 fix3, _fjsp_v2r8 fiy3, _fjsp_v2r8 fiz3,
                                          double * gmx_restrict fptr,
                                          double * gmx_restrict fshiftptr)
{
    __m128d t1, t2, t3, t4, t5, t6;

    /* transpose data */
    GMX_FJSP_TRANSPOSE2_V2R8(fix1, fiy1);
    GMX_FJSP_TRANSPOSE2_V2R8(fiz1, fix2);
    GMX_FJSP_TRANSPOSE2_V2R8(fiy2, fiz2);
    t1   = fix3;
    fix3 = _fjsp_unpacklo_v2r8(fix3, fiy3); /* y0 x0 */
    fiy3 = _fjsp_unpackhi_v2r8(t1, fiy3);   /* y1 x1 */

    fix1 = _fjsp_add_v2r8(fix1, fiy1);
    fiz1 = _fjsp_add_v2r8(fiz1, fix2);
    fiy2 = _fjsp_add_v2r8(fiy2, fiz2);

    fix3 = _fjsp_add_v2r8(fix3, fiy3);
    fiz3 = _fjsp_add_v2r8( fiz3, _fjsp_unpackhi_v2r8(fiz3, fiz3));

    t3 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr), fix1 );
    t4 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+2), fiz1 );
    t5 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+4), fiy2 );
    t6 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+6), fix3 );

    _fjsp_storel_v2r8( fptr,   t3 );
    _fjsp_storeh_v2r8( fptr+1, t3 );
    _fjsp_storel_v2r8( fptr+2, t4 );
    _fjsp_storeh_v2r8( fptr+3, t4 );
    _fjsp_storel_v2r8( fptr+4, t5 );
    _fjsp_storeh_v2r8( fptr+5, t5 );
    _fjsp_storel_v2r8( fptr+6, t6 );
    _fjsp_storeh_v2r8( fptr+7, t6 );
    _fjsp_storel_v2r8( fptr+8, _fjsp_add_v2r8( _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), fptr+8), fiz3 ));

    fix1 = _fjsp_add_v2r8(fix1, fix3);
    t1   = _fjsp_shuffle_v2r8(fiz1, fiy2, GMX_FJSP_SHUFFLE2(0, 1));
    fix1 = _fjsp_add_v2r8(fix1, t1); /* x and y sums */

    t2   = _fjsp_shuffle_v2r8(fiy2, fiy2, GMX_FJSP_SHUFFLE2(1, 1));
    fiz1 = _fjsp_add_v2r8(fiz1, fiz3);
    fiz1 = _fjsp_add_v2r8(fiz1, t2); /* z sum */

    t3 = _fjsp_add_v2r8( _fjsp_load_v2r8(fshiftptr), fix1 );
    _fjsp_storel_v2r8( fshiftptr, t3 );
    _fjsp_storeh_v2r8( fshiftptr+1, t3 );
    _fjsp_storel_v2r8( fshiftptr+2, _fjsp_add_v2r8( _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), fshiftptr+2), fiz1 ));
}


static gmx_inline void
gmx_fjsp_update_iforce_4atom_swizzle_v2r8(_fjsp_v2r8 fix1, _fjsp_v2r8 fiy1, _fjsp_v2r8 fiz1,
                                          _fjsp_v2r8 fix2, _fjsp_v2r8 fiy2, _fjsp_v2r8 fiz2,
                                          _fjsp_v2r8 fix3, _fjsp_v2r8 fiy3, _fjsp_v2r8 fiz3,
                                          _fjsp_v2r8 fix4, _fjsp_v2r8 fiy4, _fjsp_v2r8 fiz4,
                                          double * gmx_restrict fptr,
                                          double * gmx_restrict fshiftptr)
{
    __m128d t1, t2, t3, t4, t5, t6, t7, t8;

    /* transpose data */
    GMX_FJSP_TRANSPOSE2_V2R8(fix1, fiy1);
    GMX_FJSP_TRANSPOSE2_V2R8(fiz1, fix2);
    GMX_FJSP_TRANSPOSE2_V2R8(fiy2, fiz2);
    GMX_FJSP_TRANSPOSE2_V2R8(fix3, fiy3);
    GMX_FJSP_TRANSPOSE2_V2R8(fiz3, fix4);
    GMX_FJSP_TRANSPOSE2_V2R8(fiy4, fiz4);

    fix1 = _fjsp_add_v2r8(fix1, fiy1);
    fiz1 = _fjsp_add_v2r8(fiz1, fix2);
    fiy2 = _fjsp_add_v2r8(fiy2, fiz2);
    fix3 = _fjsp_add_v2r8(fix3, fiy3);
    fiz3 = _fjsp_add_v2r8(fiz3, fix4);
    fiy4 = _fjsp_add_v2r8(fiy4, fiz4);

    t3 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr),    fix1 );
    t4 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+2),  fiz1 );
    t5 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+4),  fiy2 );
    t6 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+6),  fix3 );
    t7 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+8),  fiz3 );
    t8 = _fjsp_add_v2r8( _fjsp_load_v2r8(fptr+10), fiy4 );
    _fjsp_storel_v2r8( fptr,    t3 );
    _fjsp_storeh_v2r8( fptr+1,  t3 );
    _fjsp_storel_v2r8( fptr+2,  t4 );
    _fjsp_storeh_v2r8( fptr+3,  t4 );
    _fjsp_storel_v2r8( fptr+4,  t5 );
    _fjsp_storeh_v2r8( fptr+5,  t5 );
    _fjsp_storel_v2r8( fptr+6,  t6 );
    _fjsp_storeh_v2r8( fptr+7,  t6 );
    _fjsp_storel_v2r8( fptr+8,  t7 );
    _fjsp_storeh_v2r8( fptr+9,  t7 );
    _fjsp_storel_v2r8( fptr+10, t8 );
    _fjsp_storeh_v2r8( fptr+11, t8 );

    t1   = _fjsp_shuffle_v2r8(fiz1, fiy2, GMX_FJSP_SHUFFLE2(0, 1));
    fix1 = _fjsp_add_v2r8(fix1, t1);
    t2   = _fjsp_shuffle_v2r8(fiz3, fiy4, GMX_FJSP_SHUFFLE2(0, 1));
    fix3 = _fjsp_add_v2r8(fix3, t2);
    fix1 = _fjsp_add_v2r8(fix1, fix3); /* x and y sums */

    fiz1 = _fjsp_add_v2r8(fiz1, _fjsp_unpackhi_v2r8(fiy2, fiy2));
    fiz3 = _fjsp_add_v2r8(fiz3, _fjsp_unpackhi_v2r8(fiy4, fiy4));
    fiz1 = _fjsp_add_v2r8(fiz1, fiz3); /* z sum */

    t3 = _fjsp_add_v2r8( _fjsp_load_v2r8(fshiftptr), fix1 );
    _fjsp_storel_v2r8( fshiftptr, t3 );
    _fjsp_storeh_v2r8( fshiftptr+1, t3 );
    _fjsp_storel_v2r8( fshiftptr+2, _fjsp_add_v2r8( _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), fshiftptr+2), fiz1 ));
}



static gmx_inline void
gmx_fjsp_update_1pot_v2r8(_fjsp_v2r8 pot1, double * gmx_restrict ptrA)
{
    pot1 = _fjsp_add_v2r8(pot1, _fjsp_unpackhi_v2r8(pot1, pot1));
    _fjsp_storel_v2r8(ptrA, _fjsp_add_v2r8(pot1, _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA)));
}

static gmx_inline void
gmx_fjsp_update_2pot_v2r8(_fjsp_v2r8 pot1, double * gmx_restrict ptrA,
                          _fjsp_v2r8 pot2, double * gmx_restrict ptrB)
{
    GMX_FJSP_TRANSPOSE2_V2R8(pot1, pot2);
    pot1 = _fjsp_add_v2r8(pot1, pot2);
    pot2 = _fjsp_unpackhi_v2r8(pot1, pot1);

    _fjsp_storel_v2r8(ptrA, _fjsp_add_v2r8(pot1, _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrA)));
    _fjsp_storel_v2r8(ptrB, _fjsp_add_v2r8(pot2, _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), ptrB)));
}


#endif /* _kernelutil_sparc64_hpc_ace_double_h_ */
