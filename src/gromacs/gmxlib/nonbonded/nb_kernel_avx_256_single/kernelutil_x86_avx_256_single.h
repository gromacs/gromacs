/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
#ifndef _kernelutil_x86_avx_256_single_h_
#define _kernelutil_x86_avx_256_single_h_

#include "config.h"

#define gmx_mm_castsi128_ps(a) _mm_castsi128_ps(a)

static gmx_inline __m256 gmx_simdcall
gmx_mm256_unpack128lo_ps(__m256 xmm1, __m256 xmm2)
{
    return _mm256_permute2f128_ps(xmm1, xmm2, 0x20);
}

static gmx_inline __m256 gmx_simdcall
gmx_mm256_unpack128hi_ps(__m256 xmm1, __m256 xmm2)
{
    return _mm256_permute2f128_ps(xmm1, xmm2, 0x31);
}

static gmx_inline __m256 gmx_simdcall
gmx_mm256_set_m128(__m128 hi, __m128 lo)
{
    return _mm256_insertf128_ps(_mm256_castps128_ps256(lo), hi, 0x1);
}

/* Work around gcc bug with wrong type for mask formal parameter to maskload/maskstore */
#ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), _mm_castsi128_ps(mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)    _mm_maskstore_ps((mem), _mm_castsi128_ps(mask), (x))
#    define gmx_mm256_maskload_ps(mem, mask)    _mm256_maskload_ps((mem), _mm256_castsi256_ps(mask))
#    define gmx_mm256_maskstore_ps(mem, mask, x) _mm256_maskstore_ps((mem), _mm256_castsi256_ps(mask), (x))
#else
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), (mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)    _mm_maskstore_ps((mem), (mask), (x))
#    define gmx_mm256_maskload_ps(mem, mask)    _mm256_maskload_ps((mem), (mask))
#    define gmx_mm256_maskstore_ps(mem, mask, x) _mm256_maskstore_ps((mem), (mask), (x))
#endif

/* Transpose lower/upper half of 256-bit registers separately */
#define GMX_MM256_HALFTRANSPOSE4_PS(ymm0, ymm1, ymm2, ymm3) {            \
        __m256 __tmp0, __tmp1, __tmp2, __tmp3;                               \
                                                                      \
        __tmp0   = _mm256_unpacklo_ps((ymm0), (ymm1));                     \
        __tmp1   = _mm256_unpacklo_ps((ymm2), (ymm3));                     \
        __tmp2   = _mm256_unpackhi_ps((ymm0), (ymm1));                     \
        __tmp3   = _mm256_unpackhi_ps((ymm2), (ymm3));                     \
        ymm0     = _mm256_shuffle_ps(__tmp0, __tmp1, _MM_SHUFFLE(1, 0, 1, 0)); \
        ymm1     = _mm256_shuffle_ps(__tmp0, __tmp1, _MM_SHUFFLE(3, 2, 3, 2)); \
        ymm2     = _mm256_shuffle_ps(__tmp2, __tmp3, _MM_SHUFFLE(1, 0, 1, 0)); \
        ymm3     = _mm256_shuffle_ps(__tmp2, __tmp3, _MM_SHUFFLE(3, 2, 3, 2)); \
}


static gmx_inline __m256 gmx_simdcall
gmx_mm256_calc_rsq_ps(__m256 dx, __m256 dy, __m256 dz)
{
    return _mm256_add_ps( _mm256_add_ps( _mm256_mul_ps(dx, dx), _mm256_mul_ps(dy, dy) ), _mm256_mul_ps(dz, dz) );
}

/* Normal sum of four ymm registers */
#define gmx_mm256_sum4_ps(t0, t1, t2, t3)  _mm256_add_ps(_mm256_add_ps(t0, t1), _mm256_add_ps(t2, t3))


static gmx_inline int gmx_simdcall
gmx_mm256_any_lt(__m256 a, __m256 b)
{
    return _mm256_movemask_ps(_mm256_cmp_ps(a, b, _CMP_LT_OQ));
}


static gmx_inline __m256 gmx_simdcall
gmx_mm256_load_4real_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                const float * gmx_restrict ptrC, const float * gmx_restrict ptrD)
{
    __m128 t1, t2;

    t1 = _mm_unpacklo_ps(_mm_load_ss(ptrA), _mm_load_ss(ptrC));
    t2 = _mm_unpacklo_ps(_mm_load_ss(ptrB), _mm_load_ss(ptrD));
    return _mm256_castps128_ps256(_mm_unpacklo_ps(t1, t2));
}


static gmx_inline __m256 gmx_simdcall
gmx_mm256_load_8real_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                const float * gmx_restrict ptrE, const float * gmx_restrict ptrF,
                                const float * gmx_restrict ptrG, const float * gmx_restrict ptrH)
{
    __m256 t1, t2;

    t1 = gmx_mm256_load_4real_swizzle_ps(ptrA, ptrB, ptrC, ptrD);
    t2 = gmx_mm256_load_4real_swizzle_ps(ptrE, ptrF, ptrG, ptrH);

    return _mm256_permute2f128_ps(t1, t2, 0x20);
}



static gmx_inline void gmx_simdcall
gmx_mm256_store_4real_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                 float * gmx_restrict ptrC, float * gmx_restrict ptrD, __m256 xmm1)
{
    __m256 t2, t3, t4;

    t2       = _mm256_permute_ps(xmm1, _MM_SHUFFLE(1, 1, 1, 1));
    t3       = _mm256_permute_ps(xmm1, _MM_SHUFFLE(2, 2, 2, 2));
    t4       = _mm256_permute_ps(xmm1, _MM_SHUFFLE(3, 3, 3, 3));
    _mm_store_ss(ptrA, _mm256_castps256_ps128(xmm1));
    _mm_store_ss(ptrB, _mm256_castps256_ps128(t2));
    _mm_store_ss(ptrC, _mm256_castps256_ps128(t3));
    _mm_store_ss(ptrD, _mm256_castps256_ps128(t4));
}


static gmx_inline void gmx_simdcall
gmx_mm256_store_8real_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                 float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                 float * gmx_restrict ptrE, float * gmx_restrict ptrF,
                                 float * gmx_restrict ptrG, float * gmx_restrict ptrH, __m256 xmm1)
{
    __m256 t1;

    t1 = _mm256_permute2f128_ps(xmm1, xmm1, 0x11);

    gmx_mm256_store_4real_swizzle_ps(ptrA, ptrB, ptrC, ptrD, xmm1);
    gmx_mm256_store_4real_swizzle_ps(ptrE, ptrF, ptrG, ptrH, t1);
}


static gmx_inline void gmx_simdcall
gmx_mm256_increment_4real_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                     float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                     __m256 xmm1)
{
    __m128 t1, t2, t3, t4;

    t1   = _mm256_castps256_ps128(xmm1);
    t2   = _mm_permute_ps(t1, _MM_SHUFFLE(1, 1, 1, 1));
    t3   = _mm_permute_ps(t1, _MM_SHUFFLE(2, 2, 2, 2));
    t4   = _mm_permute_ps(t1, _MM_SHUFFLE(3, 3, 3, 3));

    t1   = _mm_add_ss(t1, _mm_load_ss(ptrA));
    t2   = _mm_add_ss(t2, _mm_load_ss(ptrB));
    t3   = _mm_add_ss(t3, _mm_load_ss(ptrC));
    t4   = _mm_add_ss(t4, _mm_load_ss(ptrD));

    _mm_store_ss(ptrA, t1);
    _mm_store_ss(ptrB, t2);
    _mm_store_ss(ptrC, t3);
    _mm_store_ss(ptrD, t4);
}

static gmx_inline void gmx_simdcall
gmx_mm256_increment_8real_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                     float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                     float * gmx_restrict ptrE, float * gmx_restrict ptrF,
                                     float * gmx_restrict ptrG, float * gmx_restrict ptrH,
                                     __m256 xmm1)
{
    __m256 t1;

    t1 = _mm256_permute2f128_ps(xmm1, xmm1, 0x11);

    gmx_mm256_increment_4real_swizzle_ps(ptrA, ptrB, ptrC, ptrD, xmm1);
    gmx_mm256_increment_4real_swizzle_ps(ptrE, ptrF, ptrG, ptrH, t1);
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_4pair_swizzle_ps(const float * gmx_restrict p1, const float * gmx_restrict p2,
                                const float * gmx_restrict p3, const float * gmx_restrict p4,
                                __m256 * gmx_restrict c6, __m256 * gmx_restrict c12)
{
    __m128 t1, t2, t3, t4;

    t1   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p1); /* - - c12a  c6a */
    t2   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p2); /* - - c12b  c6b */
    t3   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p3); /* - - c12c  c6c */
    t4   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)p4); /* - - c12d  c6d */

    t1   = _mm_unpacklo_ps(t1, t2);                     /* c12b c12a c6b c6a */
    t3   = _mm_unpacklo_ps(t3, t4);                     /* c12d c12c c6d c6c */

    *c6  = _mm256_castps128_ps256(_mm_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0)));
    *c12 = _mm256_castps128_ps256(_mm_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2)));
}

static gmx_inline void gmx_simdcall
gmx_mm256_load_8pair_swizzle_ps(const float * gmx_restrict p1, const float * gmx_restrict p2,
                                const float * gmx_restrict p3, const float * gmx_restrict p4,
                                const float * gmx_restrict p5, const float * gmx_restrict p6,
                                const float * gmx_restrict p7, const float * gmx_restrict p8,
                                __m256 * gmx_restrict c6, __m256 * gmx_restrict c12)
{
    __m256 c6l, c6h, c12l, c12h;

    gmx_mm256_load_4pair_swizzle_ps(p1, p2, p3, p4, &c6l, &c12l);
    gmx_mm256_load_4pair_swizzle_ps(p5, p6, p7, p8, &c6h, &c12h);

    *c6  = _mm256_permute2f128_ps(c6l, c6h, 0x20);
    *c12 = _mm256_permute2f128_ps(c12l, c12h, 0x20);
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_shift_and_1rvec_broadcast_ps(const float * gmx_restrict xyz_shift,
                                            const float * gmx_restrict xyz,
                                            __m256 * gmx_restrict      x1,
                                            __m256 * gmx_restrict      y1,
                                            __m256 * gmx_restrict      z1)
{
    __m128 t1, t2, t3, t4;

    t1   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)xyz_shift);
    t2   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)xyz);
    t3   = _mm_load_ss(xyz_shift+2);
    t4   = _mm_load_ss(xyz+2);
    t1   = _mm_add_ps(t1, t2);
    t3   = _mm_add_ss(t3, t4);

    t2   = _mm_permute_ps(t1, _MM_SHUFFLE(1, 1, 1, 1));
    t1   = _mm_permute_ps(t1, _MM_SHUFFLE(0, 0, 0, 0));
    t3   = _mm_permute_ps(t3, _MM_SHUFFLE(0, 0, 0, 0));

    *x1  = gmx_mm256_set_m128(t1, t1);
    *y1  = gmx_mm256_set_m128(t2, t2);
    *z1  = gmx_mm256_set_m128(t3, t3);
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_shift_and_3rvec_broadcast_ps(const float * gmx_restrict xyz_shift,
                                            const float * gmx_restrict xyz,
                                            __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1,
                                            __m256 * gmx_restrict x2, __m256 * gmx_restrict y2, __m256 * gmx_restrict z2,
                                            __m256 * gmx_restrict x3, __m256 * gmx_restrict y3, __m256 * gmx_restrict z3)
{
    __m128 tA, tB;
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9;

    tA   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)xyz_shift);
    tB   = _mm_load_ss(xyz_shift+2);

    t1   = _mm_loadu_ps(xyz);
    t2   = _mm_loadu_ps(xyz+4);
    t3   = _mm_load_ss(xyz+8);

    tA   = _mm_movelh_ps(tA, tB);
    t4   = _mm_permute_ps(tA, _MM_SHUFFLE(0, 2, 1, 0));
    t5   = _mm_permute_ps(tA, _MM_SHUFFLE(1, 0, 2, 1));
    t6   = _mm_permute_ps(tA, _MM_SHUFFLE(2, 1, 0, 2));

    t1   = _mm_add_ps(t1, t4);
    t2   = _mm_add_ps(t2, t5);
    t3   = _mm_add_ss(t3, t6);

    t9   = _mm_permute_ps(t3, _MM_SHUFFLE(0, 0, 0, 0));
    t8   = _mm_permute_ps(t2, _MM_SHUFFLE(3, 3, 3, 3));
    t7   = _mm_permute_ps(t2, _MM_SHUFFLE(2, 2, 2, 2));
    t6   = _mm_permute_ps(t2, _MM_SHUFFLE(1, 1, 1, 1));
    t5   = _mm_permute_ps(t2, _MM_SHUFFLE(0, 0, 0, 0));
    t4   = _mm_permute_ps(t1, _MM_SHUFFLE(3, 3, 3, 3));
    t3   = _mm_permute_ps(t1, _MM_SHUFFLE(2, 2, 2, 2));
    t2   = _mm_permute_ps(t1, _MM_SHUFFLE(1, 1, 1, 1));
    t1   = _mm_permute_ps(t1, _MM_SHUFFLE(0, 0, 0, 0));

    *x1  = gmx_mm256_set_m128(t1, t1);
    *y1  = gmx_mm256_set_m128(t2, t2);
    *z1  = gmx_mm256_set_m128(t3, t3);
    *x2  = gmx_mm256_set_m128(t4, t4);
    *y2  = gmx_mm256_set_m128(t5, t5);
    *z2  = gmx_mm256_set_m128(t6, t6);
    *x3  = gmx_mm256_set_m128(t7, t7);
    *y3  = gmx_mm256_set_m128(t8, t8);
    *z3  = gmx_mm256_set_m128(t9, t9);
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_shift_and_4rvec_broadcast_ps(const float * gmx_restrict xyz_shift,
                                            const float * gmx_restrict xyz,
                                            __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1,
                                            __m256 * gmx_restrict x2, __m256 * gmx_restrict y2, __m256 * gmx_restrict z2,
                                            __m256 * gmx_restrict x3, __m256 * gmx_restrict y3, __m256 * gmx_restrict z3,
                                            __m256 * gmx_restrict x4, __m256 * gmx_restrict y4, __m256 * gmx_restrict z4)
{
    __m128 tA, tB;
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;

    tA   = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)xyz_shift);
    tB   = _mm_load_ss(xyz_shift+2);

    t1   = _mm_loadu_ps(xyz);
    t2   = _mm_loadu_ps(xyz+4);
    t3   = _mm_loadu_ps(xyz+8);

    tA   = _mm_movelh_ps(tA, tB);
    t4   = _mm_permute_ps(tA, _MM_SHUFFLE(0, 2, 1, 0));
    t5   = _mm_permute_ps(tA, _MM_SHUFFLE(1, 0, 2, 1));
    t6   = _mm_permute_ps(tA, _MM_SHUFFLE(2, 1, 0, 2));

    t1   = _mm_add_ps(t1, t4);
    t2   = _mm_add_ps(t2, t5);
    t3   = _mm_add_ps(t3, t6);

    t12  = _mm_permute_ps(t3, _MM_SHUFFLE(3, 3, 3, 3));
    t11  = _mm_permute_ps(t3, _MM_SHUFFLE(2, 2, 2, 2));
    t10  = _mm_permute_ps(t3, _MM_SHUFFLE(1, 1, 1, 1));
    t9   = _mm_permute_ps(t3, _MM_SHUFFLE(0, 0, 0, 0));
    t8   = _mm_permute_ps(t2, _MM_SHUFFLE(3, 3, 3, 3));
    t7   = _mm_permute_ps(t2, _MM_SHUFFLE(2, 2, 2, 2));
    t6   = _mm_permute_ps(t2, _MM_SHUFFLE(1, 1, 1, 1));
    t5   = _mm_permute_ps(t2, _MM_SHUFFLE(0, 0, 0, 0));
    t4   = _mm_permute_ps(t1, _MM_SHUFFLE(3, 3, 3, 3));
    t3   = _mm_permute_ps(t1, _MM_SHUFFLE(2, 2, 2, 2));
    t2   = _mm_permute_ps(t1, _MM_SHUFFLE(1, 1, 1, 1));
    t1   = _mm_permute_ps(t1, _MM_SHUFFLE(0, 0, 0, 0));

    *x1  = gmx_mm256_set_m128(t1, t1);
    *y1  = gmx_mm256_set_m128(t2, t2);
    *z1  = gmx_mm256_set_m128(t3, t3);
    *x2  = gmx_mm256_set_m128(t4, t4);
    *y2  = gmx_mm256_set_m128(t5, t5);
    *z2  = gmx_mm256_set_m128(t6, t6);
    *x3  = gmx_mm256_set_m128(t7, t7);
    *y3  = gmx_mm256_set_m128(t8, t8);
    *z3  = gmx_mm256_set_m128(t9, t9);
    *x4  = gmx_mm256_set_m128(t10, t10);
    *y4  = gmx_mm256_set_m128(t11, t11);
    *z4  = gmx_mm256_set_m128(t12, t12);
}



static gmx_inline void gmx_simdcall
gmx_mm256_load_1rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                     const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                     __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1)
{
    __m128  t1, t2, t3, t4;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);
    t1             = gmx_mm_maskload_ps(ptrA, mask);
    t2             = gmx_mm_maskload_ps(ptrB, mask);
    t3             = gmx_mm_maskload_ps(ptrC, mask);
    t4             = gmx_mm_maskload_ps(ptrD, mask);
    _MM_TRANSPOSE4_PS(t1, t2, t3, t4);
    *x1           = _mm256_castps128_ps256(t1);
    *y1           = _mm256_castps128_ps256(t2);
    *z1           = _mm256_castps128_ps256(t3);
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_3rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                     const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                     __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1,
                                     __m256 * gmx_restrict x2, __m256 * gmx_restrict y2, __m256 * gmx_restrict z2,
                                     __m256 * gmx_restrict x3, __m256 * gmx_restrict y3, __m256 * gmx_restrict z3)
{
    __m128 t1, t2, t3, t4;
    t1            = _mm_loadu_ps(ptrA);
    t2            = _mm_loadu_ps(ptrB);
    t3            = _mm_loadu_ps(ptrC);
    t4            = _mm_loadu_ps(ptrD);
    _MM_TRANSPOSE4_PS(t1, t2, t3, t4);
    *x1           = _mm256_castps128_ps256(t1);
    *y1           = _mm256_castps128_ps256(t2);
    *z1           = _mm256_castps128_ps256(t3);
    *x2           = _mm256_castps128_ps256(t4);
    t1            = _mm_loadu_ps(ptrA+4);
    t2            = _mm_loadu_ps(ptrB+4);
    t3            = _mm_loadu_ps(ptrC+4);
    t4            = _mm_loadu_ps(ptrD+4);
    _MM_TRANSPOSE4_PS(t1, t2, t3, t4);
    *y2           = _mm256_castps128_ps256(t1);
    *z2           = _mm256_castps128_ps256(t2);
    *x3           = _mm256_castps128_ps256(t3);
    *y3           = _mm256_castps128_ps256(t4);
    t1            = _mm_load_ss(ptrA+8);
    t2            = _mm_load_ss(ptrB+8);
    t3            = _mm_load_ss(ptrC+8);
    t4            = _mm_load_ss(ptrD+8);
    t1            = _mm_unpacklo_ps(t1, t3);
    t3            = _mm_unpacklo_ps(t2, t4);
    *z3           = _mm256_castps128_ps256(_mm_unpacklo_ps(t1, t3));
}



static gmx_inline void gmx_simdcall
gmx_mm256_load_4rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                     const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                     __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1,
                                     __m256 * gmx_restrict x2, __m256 * gmx_restrict y2, __m256 * gmx_restrict z2,
                                     __m256 * gmx_restrict x3, __m256 * gmx_restrict y3, __m256 * gmx_restrict z3,
                                     __m256 * gmx_restrict x4, __m256 * gmx_restrict y4, __m256 * gmx_restrict z4)
{
    __m128 t1, t2, t3, t4;
    t1            = _mm_loadu_ps(ptrA);
    t2            = _mm_loadu_ps(ptrB);
    t3            = _mm_loadu_ps(ptrC);
    t4            = _mm_loadu_ps(ptrD);
    _MM_TRANSPOSE4_PS(t1, t2, t3, t4);
    *x1           = _mm256_castps128_ps256(t1);
    *y1           = _mm256_castps128_ps256(t2);
    *z1           = _mm256_castps128_ps256(t3);
    *x2           = _mm256_castps128_ps256(t4);
    t1            = _mm_loadu_ps(ptrA+4);
    t2            = _mm_loadu_ps(ptrB+4);
    t3            = _mm_loadu_ps(ptrC+4);
    t4            = _mm_loadu_ps(ptrD+4);
    _MM_TRANSPOSE4_PS(t1, t2, t3, t4);
    *y2           = _mm256_castps128_ps256(t1);
    *z2           = _mm256_castps128_ps256(t2);
    *x3           = _mm256_castps128_ps256(t3);
    *y3           = _mm256_castps128_ps256(t4);
    t1            = _mm_loadu_ps(ptrA+8);
    t2            = _mm_loadu_ps(ptrB+8);
    t3            = _mm_loadu_ps(ptrC+8);
    t4            = _mm_loadu_ps(ptrD+8);
    _MM_TRANSPOSE4_PS(t1, t2, t3, t4);
    *z3           = _mm256_castps128_ps256(t1);
    *x4           = _mm256_castps128_ps256(t2);
    *y4           = _mm256_castps128_ps256(t3);
    *z4           = _mm256_castps128_ps256(t4);
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_1rvec_8ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                     const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                     const float * gmx_restrict ptrE, const float * gmx_restrict ptrF,
                                     const float * gmx_restrict ptrG, const float * gmx_restrict ptrH,
                                     __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1)
{
    __m256  t1, t2, t3, t4, t5, t6, t7, t8;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    t1             = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrE, mask), gmx_mm_maskload_ps(ptrA, mask)); /*  - zE yE xE |  - zA yA xA */
    t2             = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrF, mask), gmx_mm_maskload_ps(ptrB, mask)); /*  - zF yF xF |  - zB yB xB */
    t3             = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrG, mask), gmx_mm_maskload_ps(ptrC, mask)); /*  - zG yG xG |  - zC yC xC */
    t4             = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrH, mask), gmx_mm_maskload_ps(ptrD, mask)); /*  - zH yH xH |  - zD yD xD */

    t5            = _mm256_unpacklo_ps(t1, t2);                                                          /* yF yE xF xE | yB yA xB xA */
    t6            = _mm256_unpacklo_ps(t3, t4);                                                          /* yH yG xH xG | yD yC xD xC */
    t7            = _mm256_unpackhi_ps(t1, t2);                                                          /*  -  - zF zE |  -  - zB zA */
    t8            = _mm256_unpackhi_ps(t3, t4);                                                          /*  -  - zH zG |  -  - zD zC */

    *x1           = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(1, 0, 1, 0));
    *y1           = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(3, 2, 3, 2));
    *z1           = _mm256_shuffle_ps(t7, t8, _MM_SHUFFLE(1, 0, 1, 0));
}


static gmx_inline void gmx_simdcall
gmx_mm256_load_3rvec_8ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                     const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                     const float * gmx_restrict ptrE, const float * gmx_restrict ptrF,
                                     const float * gmx_restrict ptrG, const float * gmx_restrict ptrH,
                                     __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1,
                                     __m256 * gmx_restrict x2, __m256 * gmx_restrict y2, __m256 * gmx_restrict z2,
                                     __m256 * gmx_restrict x3, __m256 * gmx_restrict y3, __m256 * gmx_restrict z3)
{
    __m256 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;

    t1           = _mm256_loadu_ps(ptrA);                                /* y3a x3a z2a y2a | x2a z1a y1a x1a */
    t2           = _mm256_loadu_ps(ptrB);                                /* y3b x3b z2b y2b | x2b z1b y1b x1b */
    t3           = _mm256_loadu_ps(ptrC);                                /* y3c x3c z2c y2c | x2c z1c y1c x1c */
    t4           = _mm256_loadu_ps(ptrD);                                /* y3d x3d z2d y2d | x2d z1d y1d x1d */
    t5           = _mm256_loadu_ps(ptrE);                                /* y3e x3e z2e y2e | x2e z1e y1e x1e */
    t6           = _mm256_loadu_ps(ptrF);                                /* y3f x3f z2f y2f | x2f z1f y1f x1f */
    t7           = _mm256_loadu_ps(ptrG);                                /* y3g x3g z2g y2g | x2g z1g y1g x1g */
    t8           = _mm256_loadu_ps(ptrH);                                /* y3h x3h z2h y2h | x2h z1h y1h x1h */

    t9           = _mm256_unpacklo_ps(t1, t2);                           /* z2b z2a y2b y2a | y1b y1a x1b x1a */
    t10          = _mm256_unpackhi_ps(t1, t2);                           /* y3b y3a x3b x3a | x2b x2a z1b z1a */
    t11          = _mm256_unpacklo_ps(t3, t4);                           /* z2d z2c y2d y2c | y1d y1c x1d x1c */
    t12          = _mm256_unpackhi_ps(t3, t4);                           /* y3d y3c x3d x3c | x2d x2c z1d z1c */
    t1           = _mm256_unpacklo_ps(t5, t6);                           /* z2f z2e y2f y2e | y1f y1e x1f x1e */
    t2           = _mm256_unpackhi_ps(t5, t6);                           /* y3f y3e x3f x3e | x2f x2e z1f z1e */
    t3           = _mm256_unpacklo_ps(t7, t8);                           /* z2h z2g y2h y2g | y1h y1g x1h x1g */
    t4           = _mm256_unpackhi_ps(t7, t8);                           /* y3h y3g x3h x3g | x2h x2g z1h z1g */

    t5           = _mm256_shuffle_ps(t9, t11, _MM_SHUFFLE(1, 0, 1, 0));  /* y2d y2c y2b y2a | x1d x1c x1b x1a */
    t6           = _mm256_shuffle_ps(t9, t11, _MM_SHUFFLE(3, 2, 3, 2));  /* z2d z2c z2b z2a | y1d y1c y1b y1a */
    t7           = _mm256_shuffle_ps(t10, t12, _MM_SHUFFLE(1, 0, 1, 0)); /* x3d x3c x3b x3a | z1d z1c z1b z1a */
    t8           = _mm256_shuffle_ps(t10, t12, _MM_SHUFFLE(3, 2, 3, 2)); /* y3d y3c y3b y3a | x2d x2c x2b x2a */

    t9           = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));   /* y2h y2g y2f y2e | x1h x1g x1f x1e */
    t10          = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));   /* z2h z2g z2f z2e | y1h y1g y1f y1e */
    t11          = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 0, 1, 0));   /* x3h x3g x3f x3e | z1h z1g z1f z1e */
    t12          = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(3, 2, 3, 2));   /* y3h y3g y3f y3e | x2h x2g x2f x2e */

    *x1          = _mm256_permute2f128_ps(t5, t9,  0x20);
    *y1          = _mm256_permute2f128_ps(t6, t10, 0x20);
    *z1          = _mm256_permute2f128_ps(t7, t11, 0x20);
    *x2          = _mm256_permute2f128_ps(t8, t12, 0x20);

    *y2          = _mm256_permute2f128_ps(t5, t9,  0x31);
    *z2          = _mm256_permute2f128_ps(t6, t10, 0x31);
    *x3          = _mm256_permute2f128_ps(t7, t11, 0x31);
    *y3          = _mm256_permute2f128_ps(t8, t12, 0x31);

    t1           = gmx_mm256_set_m128(_mm_load_ss(ptrE+8), _mm_load_ss(ptrA+8));
    t2           = gmx_mm256_set_m128(_mm_load_ss(ptrF+8), _mm_load_ss(ptrB+8));
    t3           = gmx_mm256_set_m128(_mm_load_ss(ptrG+8), _mm_load_ss(ptrC+8));
    t4           = gmx_mm256_set_m128(_mm_load_ss(ptrH+8), _mm_load_ss(ptrD+8));

    t1           = _mm256_unpacklo_ps(t1, t3);  /*  -   -  z3g z3e |  -   -  z3c z3a */
    t2           = _mm256_unpacklo_ps(t2, t4);  /*  -   -  z3h z3f |  -   -  z3d z3b */

    *z3          = _mm256_unpacklo_ps(t1, t2);
}



static gmx_inline void gmx_simdcall
gmx_mm256_load_4rvec_8ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                     const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                     const float * gmx_restrict ptrE, const float * gmx_restrict ptrF,
                                     const float * gmx_restrict ptrG, const float * gmx_restrict ptrH,
                                     __m256 * gmx_restrict x1, __m256 * gmx_restrict y1, __m256 * gmx_restrict z1,
                                     __m256 * gmx_restrict x2, __m256 * gmx_restrict y2, __m256 * gmx_restrict z2,
                                     __m256 * gmx_restrict x3, __m256 * gmx_restrict y3, __m256 * gmx_restrict z3,
                                     __m256 * gmx_restrict x4, __m256 * gmx_restrict y4, __m256 * gmx_restrict z4)
{
    __m256 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;

    t1           = _mm256_loadu_ps(ptrA);                                /* y3a x3a z2a y2a | x2a z1a y1a x1a */
    t2           = _mm256_loadu_ps(ptrB);                                /* y3b x3b z2b y2b | x2b z1b y1b x1b */
    t3           = _mm256_loadu_ps(ptrC);                                /* y3c x3c z2c y2c | x2c z1c y1c x1c */
    t4           = _mm256_loadu_ps(ptrD);                                /* y3d x3d z2d y2d | x2d z1d y1d x1d */
    t5           = _mm256_loadu_ps(ptrE);                                /* y3e x3e z2e y2e | x2e z1e y1e x1e */
    t6           = _mm256_loadu_ps(ptrF);                                /* y3f x3f z2f y2f | x2f z1f y1f x1f */
    t7           = _mm256_loadu_ps(ptrG);                                /* y3g x3g z2g y2g | x2g z1g y1g x1g */
    t8           = _mm256_loadu_ps(ptrH);                                /* y3h x3h z2h y2h | x2h z1h y1h x1h */

    t9           = _mm256_unpacklo_ps(t1, t2);                           /* z2b z2a y2b y2a | y1b y1a x1b x1a */
    t10          = _mm256_unpackhi_ps(t1, t2);                           /* y3b y3a x3b x3a | x2b x2a z1b z1a */
    t11          = _mm256_unpacklo_ps(t3, t4);                           /* z2d z2c y2d y2c | y1d y1c x1d x1c */
    t12          = _mm256_unpackhi_ps(t3, t4);                           /* y3d y3c x3d x3c | x2d x2c z1d z1c */
    t1           = _mm256_unpacklo_ps(t5, t6);                           /* z2f z2e y2f y2e | y1f y1e x1f x1e */
    t2           = _mm256_unpackhi_ps(t5, t6);                           /* y3f y3e x3f x3e | x2f x2e z1f z1e */
    t3           = _mm256_unpacklo_ps(t7, t8);                           /* z2h z2g y2h y2g | y1h y1g x1h x1g */
    t4           = _mm256_unpackhi_ps(t7, t8);                           /* y3h y3g x3h x3g | x2h x2g z1h z1g */

    t5           = _mm256_shuffle_ps(t9, t11, _MM_SHUFFLE(1, 0, 1, 0));  /* y2d y2c y2b y2a | x1d x1c x1b x1a */
    t6           = _mm256_shuffle_ps(t9, t11, _MM_SHUFFLE(3, 2, 3, 2));  /* z2d z2c z2b z2a | y1d y1c y1b y1a */
    t7           = _mm256_shuffle_ps(t10, t12, _MM_SHUFFLE(1, 0, 1, 0)); /* x3d x3c x3b x3a | z1d z1c z1b z1a */
    t8           = _mm256_shuffle_ps(t10, t12, _MM_SHUFFLE(3, 2, 3, 2)); /* y3d y3c y3b y3a | x2d x2c x2b x2a */
    t9           = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));   /* y2h y2g y2f y2e | x1h x1g x1f x1e */
    t10          = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));   /* z2h z2g z2f z2e | y1h y1g y1f y1e */
    t11          = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 0, 1, 0));   /* x3h x3g x3f x3e | z1h z1g z1f z1e */
    t12          = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(3, 2, 3, 2));   /* y3h y3g y3f y3e | x2h x2g x2f x2e */

    *x1          = _mm256_permute2f128_ps(t5, t9,  0x20);
    *y1          = _mm256_permute2f128_ps(t6, t10, 0x20);
    *z1          = _mm256_permute2f128_ps(t7, t11, 0x20);
    *x2          = _mm256_permute2f128_ps(t8, t12, 0x20);

    *y2          = _mm256_permute2f128_ps(t5, t9,  0x31);
    *z2          = _mm256_permute2f128_ps(t6, t10, 0x31);
    *x3          = _mm256_permute2f128_ps(t7, t11, 0x31);
    *y3          = _mm256_permute2f128_ps(t8, t12, 0x31);

    t1           = gmx_mm256_set_m128(_mm_loadu_ps(ptrE+8), _mm_loadu_ps(ptrA+8)); /* z4e y4e x4e z3e | z4a y4a x4a z3a */
    t2           = gmx_mm256_set_m128(_mm_loadu_ps(ptrF+8), _mm_loadu_ps(ptrB+8)); /* z4f y4f x4f z3f | z4b y4b x4b z3b */
    t3           = gmx_mm256_set_m128(_mm_loadu_ps(ptrG+8), _mm_loadu_ps(ptrC+8)); /* z4g y4g x4g z3g | z4c y4c x4c z3c */
    t4           = gmx_mm256_set_m128(_mm_loadu_ps(ptrH+8), _mm_loadu_ps(ptrD+8)); /* z4h y4h x4h z3h | z4d y4d x4d z3d */

    t5           = _mm256_unpacklo_ps(t1, t2);                                     /* x4f x4e z3f z3e | x4b x4a z3b z3a */
    t6           = _mm256_unpackhi_ps(t1, t2);                                     /* z4f z4e y4f y4e | z4b z4a y4b y4a */
    t7           = _mm256_unpacklo_ps(t3, t4);                                     /* x4h x4g z3h z3g | x4d x4c z3d z3c */
    t8           = _mm256_unpackhi_ps(t3, t4);                                     /* z4h z4g y4h y4g | z4d z4c y4d y4c */

    *z3          = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0));             /* z3h z3g z3f z3e | z3d z3c z3b z3a */
    *x4          = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2));             /* x4h x4g x4f x4e | x4d x4c x4b x4a */
    *y4          = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(1, 0, 1, 0));             /* y4h y4g y4f y4e | y4d y4c y4b y4a */
    *z4          = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(3, 2, 3, 2));             /* z4h z4g z4f z4e | z4d z4c z4b z4a */
}


static gmx_inline void gmx_simdcall
gmx_mm256_decrement_1rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                          float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                          __m256 x1, __m256 y1, __m256 z1)
{
    __m128  t1, t2, t3, t4, t5, t6, t7, t8;
    __m128i mask;

    /* Construct a mask without executing any data loads */
    mask        = _mm_blend_epi16(_mm_setzero_si128(), _mm_cmpeq_epi16(_mm_setzero_si128(), _mm_setzero_si128()), 0x3F);

    t3          = _mm_unpacklo_ps(_mm256_castps256_ps128(x1), _mm256_castps256_ps128(y1)); /* y1b x1b y1a x1a */
    t4          = _mm_unpackhi_ps(_mm256_castps256_ps128(x1), _mm256_castps256_ps128(y1)); /* y1d x1d y1c x1c */

    t1          = _mm_shuffle_ps(t3, _mm256_castps256_ps128(z1), _MM_SHUFFLE(0, 0, 1, 0)); /*  -  z1a y1a x1a */
    t2          = _mm_shuffle_ps(t3, _mm256_castps256_ps128(z1), _MM_SHUFFLE(0, 1, 3, 2)); /*  -  z1b y1b x1b */
    t3          = _mm_shuffle_ps(t4, _mm256_castps256_ps128(z1), _MM_SHUFFLE(0, 2, 1, 0)); /*  -  z1c y1c x1c */
    t4          = _mm_shuffle_ps(t4, _mm256_castps256_ps128(z1), _MM_SHUFFLE(0, 3, 3, 2)); /*  -  z1d y1d x1d */

    t5          = gmx_mm_maskload_ps(ptrA, mask);
    t6          = gmx_mm_maskload_ps(ptrB, mask);
    t7          = gmx_mm_maskload_ps(ptrC, mask);
    t8          = gmx_mm_maskload_ps(ptrD, mask);

    t5          = _mm_sub_ps(t5, t1);
    t6          = _mm_sub_ps(t6, t2);
    t7          = _mm_sub_ps(t7, t3);
    t8          = _mm_sub_ps(t8, t4);

    gmx_mm_maskstore_ps(ptrA, mask, t5);
    gmx_mm_maskstore_ps(ptrB, mask, t6);
    gmx_mm_maskstore_ps(ptrC, mask, t7);
    gmx_mm_maskstore_ps(ptrD, mask, t8);
}


static gmx_inline void gmx_simdcall
gmx_mm256_decrement_3rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                          float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                          __m256 x1, __m256 y1, __m256 z1,
                                          __m256 x2, __m256 y2, __m256 z2,
                                          __m256 x3, __m256 y3, __m256 z3)
{
    __m256 t1, t2, t3, t4, t5, t6;
    __m128 tA, tB, tC, tD;

    t1          = _mm256_loadu_ps(ptrA);
    t2          = _mm256_loadu_ps(ptrB);
    t3          = _mm256_loadu_ps(ptrC);
    t4          = _mm256_loadu_ps(ptrD);
    tA          = _mm_load_ss(ptrA+8);
    tB          = _mm_load_ss(ptrB+8);
    tC          = _mm_load_ss(ptrC+8);
    tD          = _mm_load_ss(ptrD+8);

    t5          = _mm256_unpacklo_ps(x1, y1);                                /* - - - - | y1b x1b y1a x1a */
    x1          = _mm256_unpackhi_ps(x1, y1);                                /* - - - - | y1d x1d y1c x1c */
    y1          = _mm256_unpacklo_ps(z1, x2);                                /* - - - - | x2b z1b x2a z1a */
    z1          = _mm256_unpackhi_ps(z1, x2);                                /* - - - - | x2d z1d x2c z1c */

    x2          = _mm256_unpacklo_ps(y2, z2);                                /* - - - - | z2b y2b z2a y2a */
    y2          = _mm256_unpackhi_ps(y2, z2);                                /* - - - - | z2d y2d z2c y2c */
    t6          = _mm256_unpacklo_ps(x3, y3);                                /* - - - - | y3b x3b y3a x3a */
    x3          = _mm256_unpackhi_ps(x3, y3);                                /* - - - - | y3d x3d y3c x3c */

    t5          = _mm256_insertf128_ps(t5, _mm256_castps256_ps128(x2), 0x1); /* z2b y2b z2a y2a | y1b x1b y1a x1a */
    x1          = _mm256_insertf128_ps(x1, _mm256_castps256_ps128(y2), 0x1); /* z2d y2d z2c y2c | y1d x1d y1c x1c */

    y1          = _mm256_insertf128_ps(y1, _mm256_castps256_ps128(t6), 0x1); /* y3b x3b y3a x3a | x2b z1b x2a z1a */
    z1          = _mm256_insertf128_ps(z1, _mm256_castps256_ps128(x3), 0x1); /* y3d x3d y3c x3c | x2d z1d x2c z1c */

    z2          = _mm256_shuffle_ps(t5, y1, _MM_SHUFFLE(1, 0, 1, 0));        /* y3a x3a z2a y2a | x2a z1a y1a x1a */
    t5          = _mm256_shuffle_ps(t5, y1, _MM_SHUFFLE(3, 2, 3, 2));        /* y3b x3b z2b y2b | x2b z1b y1b x1b */
    y1          = _mm256_shuffle_ps(x1, z1, _MM_SHUFFLE(1, 0, 1, 0));        /* y3c x3c z2c y2c | x2c z1c y1c x1c */
    x1          = _mm256_shuffle_ps(x1, z1, _MM_SHUFFLE(3, 2, 3, 2));        /* y3d x3d z2d y2d | x2d z1d y1d x1d */

    t1          = _mm256_sub_ps(t1, z2);
    t2          = _mm256_sub_ps(t2, t5);
    t3          = _mm256_sub_ps(t3, y1);
    t4          = _mm256_sub_ps(t4, x1);

    tA          = _mm_sub_ss(tA, _mm256_castps256_ps128(z3));
    tB          = _mm_sub_ss(tB, _mm_permute_ps(_mm256_castps256_ps128(z3), _MM_SHUFFLE(1, 1, 1, 1)));
    tC          = _mm_sub_ss(tC, _mm_permute_ps(_mm256_castps256_ps128(z3), _MM_SHUFFLE(2, 2, 2, 2)));
    tD          = _mm_sub_ss(tD, _mm_permute_ps(_mm256_castps256_ps128(z3), _MM_SHUFFLE(3, 3, 3, 3)));

    /* Here we store a full 256-bit value and a separate 32-bit one; no overlap can happen */
    _mm256_storeu_ps(ptrA, t1);
    _mm256_storeu_ps(ptrB, t2);
    _mm256_storeu_ps(ptrC, t3);
    _mm256_storeu_ps(ptrD, t4);
    _mm_store_ss(ptrA+8, tA);
    _mm_store_ss(ptrB+8, tB);
    _mm_store_ss(ptrC+8, tC);
    _mm_store_ss(ptrD+8, tD);
}


static gmx_inline void gmx_simdcall
gmx_mm256_decrement_4rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                          float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                          __m256 x1, __m256 y1, __m256 z1,
                                          __m256 x2, __m256 y2, __m256 z2,
                                          __m256 x3, __m256 y3, __m256 z3,
                                          __m256 x4, __m256 y4, __m256 z4)
{
    __m256 t1, t2, t3, t4, t5;
    __m128 tA, tB, tC, tD, tE, tF, tG, tH;

    t1          = _mm256_loadu_ps(ptrA);
    t2          = _mm256_loadu_ps(ptrB);
    t3          = _mm256_loadu_ps(ptrC);
    t4          = _mm256_loadu_ps(ptrD);
    tA          = _mm_loadu_ps(ptrA+8);
    tB          = _mm_loadu_ps(ptrB+8);
    tC          = _mm_loadu_ps(ptrC+8);
    tD          = _mm_loadu_ps(ptrD+8);

    t5          = _mm256_unpacklo_ps(x1, y1);                                                                      /* - - - - | y1b x1b y1a x1a */
    x1          = _mm256_unpackhi_ps(x1, y1);                                                                      /* - - - - | y1d x1d y1c x1c */
    y1          = _mm256_unpacklo_ps(z1, x2);                                                                      /* - - - - | x2b z1b x2a z1a */
    z1          = _mm256_unpackhi_ps(z1, x2);                                                                      /* - - - - | x2d z1d x2c z1c */

    x2          = _mm256_unpacklo_ps(y2, z2);                                                                      /* - - - - | z2b y2b z2a y2a */
    y2          = _mm256_unpackhi_ps(y2, z2);                                                                      /* - - - - | z2d y2d z2c y2c */
    z2          = _mm256_unpacklo_ps(x3, y3);                                                                      /* - - - - | y3b x3b y3a x3a */
    x3          = _mm256_unpackhi_ps(x3, y3);                                                                      /* - - - - | y3d x3d y3c x3c */

    y3          = _mm256_unpacklo_ps(z3, x4);                                                                      /* - - - - | x4b z3b x4a z3a */
    z3          = _mm256_unpackhi_ps(z3, x4);                                                                      /* - - - - | x4d z3d x4c z3c */
    x4          = _mm256_unpacklo_ps(y4, z4);                                                                      /* - - - - | z4b y4b z4a y4a */
    y4          = _mm256_unpackhi_ps(y4, z4);                                                                      /* - - - - | z4d y4d z4c y4c */

    x2          = _mm256_insertf128_ps(t5, _mm256_castps256_ps128(x2), 0x1);                                       /* z2b y2b z2a y2a | y1b x1b y1a x1a */
    x1          = _mm256_insertf128_ps(x1, _mm256_castps256_ps128(y2), 0x1);                                       /* z2d y2d z2c y2c | y1d x1d y1c x1c */
    y1          = _mm256_insertf128_ps(y1, _mm256_castps256_ps128(z2), 0x1);                                       /* y3b x3b y3a x3a | x2b z1b x2a z1a */
    z1          = _mm256_insertf128_ps(z1, _mm256_castps256_ps128(x3), 0x1);                                       /* y3d x3d y3c x3c | x2d z1d x2c z1c */

    z2          = _mm256_shuffle_ps(x2, y1, _MM_SHUFFLE(1, 0, 1, 0));                                              /* y3a x3a z2a y2a | x2a z1a y1a x1a */
    t5          = _mm256_shuffle_ps(x2, y1, _MM_SHUFFLE(3, 2, 3, 2));                                              /* y3b x3b z2b y2b | x2b z1b y1b x1b */
    y1          = _mm256_shuffle_ps(x1, z1, _MM_SHUFFLE(1, 0, 1, 0));                                              /* y3c x3c z2c y2c | x2c z1c y1c x1c */
    x1          = _mm256_shuffle_ps(x1, z1, _MM_SHUFFLE(3, 2, 3, 2));                                              /* y3d x3d z2d y2d | x2d z1d y1d x1d */

    tE          = _mm_shuffle_ps(_mm256_castps256_ps128(y3), _mm256_castps256_ps128(x4), _MM_SHUFFLE(1, 0, 1, 0)); /* z4a y4a x4a z3a */
    tF          = _mm_shuffle_ps(_mm256_castps256_ps128(y3), _mm256_castps256_ps128(x4), _MM_SHUFFLE(3, 2, 3, 2)); /* z4b y4b x4b z3b */

    tG          = _mm_shuffle_ps(_mm256_castps256_ps128(z3), _mm256_castps256_ps128(y4), _MM_SHUFFLE(1, 0, 1, 0)); /* z4c y4c x4c z3c */
    tH          = _mm_shuffle_ps(_mm256_castps256_ps128(z3), _mm256_castps256_ps128(y4), _MM_SHUFFLE(3, 2, 3, 2)); /* z4d y4d x4d z3d */

    t1          = _mm256_sub_ps(t1, z2);
    t2          = _mm256_sub_ps(t2, t5);
    t3          = _mm256_sub_ps(t3, y1);
    t4          = _mm256_sub_ps(t4, x1);

    tA          = _mm_sub_ps(tA, tE);
    tB          = _mm_sub_ps(tB, tF);
    tC          = _mm_sub_ps(tC, tG);
    tD          = _mm_sub_ps(tD, tH);

    /* Here we store a full 256-bit value and a separate 128-bit one; no overlap can happen */
    _mm256_storeu_ps(ptrA, t1);
    _mm256_storeu_ps(ptrB, t2);
    _mm256_storeu_ps(ptrC, t3);
    _mm256_storeu_ps(ptrD, t4);
    _mm_storeu_ps(ptrA+8, tA);
    _mm_storeu_ps(ptrB+8, tB);
    _mm_storeu_ps(ptrC+8, tC);
    _mm_storeu_ps(ptrD+8, tD);
}


static gmx_inline void gmx_simdcall
gmx_mm256_decrement_1rvec_8ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                          float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                          float * gmx_restrict ptrE, float * gmx_restrict ptrF,
                                          float * gmx_restrict ptrG, float * gmx_restrict ptrH,
                                          __m256 x1, __m256 y1, __m256 z1)
{
    __m256  t1, t2, t3, t4, t5, t6;
    __m256  tA, tB, tC, tD;
    __m128i mask;

    /* Construct a mask without executing any data loads */
    mask        = _mm_blend_epi16(_mm_setzero_si128(), _mm_cmpeq_epi16(_mm_setzero_si128(), _mm_setzero_si128()), 0x3F);

    tA          = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrE, mask), gmx_mm_maskload_ps(ptrA, mask));
    tB          = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrF, mask), gmx_mm_maskload_ps(ptrB, mask));
    tC          = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrG, mask), gmx_mm_maskload_ps(ptrC, mask));
    tD          = gmx_mm256_set_m128(gmx_mm_maskload_ps(ptrH, mask), gmx_mm_maskload_ps(ptrD, mask));
    t1          = _mm256_unpacklo_ps(x1, y1);                         /* y1f x1f y1e x1e | y1b x1b y1a x1a */
    t2          = _mm256_unpackhi_ps(x1, y1);                         /* y1h x1h y1g x1g | y1d x1d y1c x1c */

    t3          = _mm256_shuffle_ps(t1, z1, _MM_SHUFFLE(0, 0, 1, 0)); /*  -  z1e y1e x1e |  - z1a y1a x1a */
    t4          = _mm256_shuffle_ps(t1, z1, _MM_SHUFFLE(0, 1, 3, 2)); /*  -  z1f y1f x1f |  - z1b y1b x1b */
    t5          = _mm256_shuffle_ps(t2, z1, _MM_SHUFFLE(0, 2, 1, 0)); /*  -  z1g y1g x1g |  - z1c y1c x1c */
    t6          = _mm256_shuffle_ps(t2, z1, _MM_SHUFFLE(0, 3, 3, 2)); /*  -  z1h y1h x1h |  - z1d y1d x1d */

    tA          = _mm256_sub_ps(tA, t3);
    tB          = _mm256_sub_ps(tB, t4);
    tC          = _mm256_sub_ps(tC, t5);
    tD          = _mm256_sub_ps(tD, t6);

    gmx_mm_maskstore_ps(ptrA, mask, _mm256_castps256_ps128(tA));
    gmx_mm_maskstore_ps(ptrB, mask, _mm256_castps256_ps128(tB));
    gmx_mm_maskstore_ps(ptrC, mask, _mm256_castps256_ps128(tC));
    gmx_mm_maskstore_ps(ptrD, mask, _mm256_castps256_ps128(tD));
    gmx_mm_maskstore_ps(ptrE, mask, _mm256_extractf128_ps(tA, 0x1));
    gmx_mm_maskstore_ps(ptrF, mask, _mm256_extractf128_ps(tB, 0x1));
    gmx_mm_maskstore_ps(ptrG, mask, _mm256_extractf128_ps(tC, 0x1));
    gmx_mm_maskstore_ps(ptrH, mask, _mm256_extractf128_ps(tD, 0x1));
}



static gmx_inline void gmx_simdcall
gmx_mm256_decrement_3rvec_8ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                          float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                          float * gmx_restrict ptrE, float * gmx_restrict ptrF,
                                          float * gmx_restrict ptrG, float * gmx_restrict ptrH,
                                          __m256 x1, __m256 y1, __m256 z1,
                                          __m256 x2, __m256 y2, __m256 z2,
                                          __m256 x3, __m256 y3, __m256 z3)
{
    __m256 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
    __m256 tA, tB, tC, tD, tE, tF, tG, tH;
    __m256 tI, tJ, tK, tL;

    tA          = _mm256_loadu_ps(ptrA);
    tB          = _mm256_loadu_ps(ptrB);
    tC          = _mm256_loadu_ps(ptrC);
    tD          = _mm256_loadu_ps(ptrD);
    tE          = _mm256_loadu_ps(ptrE);
    tF          = _mm256_loadu_ps(ptrF);
    tG          = _mm256_loadu_ps(ptrG);
    tH          = _mm256_loadu_ps(ptrH);

    t1          = _mm256_unpacklo_ps(x1, y1);                         /* y1f x1f y1e x1e | y1b x1b y1a x1a */
    t2          = _mm256_unpackhi_ps(x1, y1);                         /* y1h x1h y1g x1g | y1d x1d y1c x1c */
    t3          = _mm256_unpacklo_ps(z1, x2);                         /* x2f z1f x2e z1e | x2b z1b x2a z1a */
    t4          = _mm256_unpackhi_ps(z1, x2);                         /* x2h z1h x2g z1g | x2d z1d x2c z1c */

    t5          = _mm256_unpacklo_ps(y2, z2);                         /* z2f y2f z2e y2e | z2b y2b z2a y2a */
    t6          = _mm256_unpackhi_ps(y2, z2);                         /* z2h y2h z2g y2g | z2d y2d z2c y2c */
    t7          = _mm256_unpacklo_ps(x3, y3);                         /* y3f x3f y3e x3e | y3b x3b y3a x3a */
    t8          = _mm256_unpackhi_ps(x3, y3);                         /* y3h x3h y3g x3g | y3d x3d y3c x3c */

    t9          = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0)); /* x2e z1e y1e x1e | x2a z1a y1a x1a */
    t10         = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2)); /* x2f z1f y1f x1f | x2b z1b y1b x1b */
    t11         = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 0, 1, 0)); /* x2g z1g y1g x1g | x2c z1c y1c x1c */
    t12         = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(3, 2, 3, 2)); /* x2h z1h y1h x1h | x2d z1d y1d x1d */

    t1          = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0)); /* y3e x3e z2e y2e | y3a x3a z2a y2a */
    t2          = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2)); /* y3f x3f z2f y2f | y3b x3b z2b y2b */
    t3          = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(1, 0, 1, 0)); /* y3g x3g z2g y2g | y3c x3c z2c y2c */
    t4          = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(3, 2, 3, 2)); /* y3h x3h z2h y2h | y3d x3d z2d y2d */

    t5          = gmx_mm256_unpack128lo_ps(t9, t1);                   /* y3a x3a z2a y2a | x2a z1a y1a x1a */
    t6          = gmx_mm256_unpack128hi_ps(t9, t1);                   /* y3e x3e z2e y2e | x2e z1e y1e x1e */
    t7          = gmx_mm256_unpack128lo_ps(t10, t2);                  /* y3b x3b z2b y2b | x2b z1b y1b x1b */
    t8          = gmx_mm256_unpack128hi_ps(t10, t2);                  /* y3f x3f z2f y2f | x2f z1f y1f x1f */
    t1          = gmx_mm256_unpack128lo_ps(t11, t3);                  /* y3c x3c z2c y2c | x2c z1c y1c x1c */
    t2          = gmx_mm256_unpack128hi_ps(t11, t3);                  /* y3g x3g z2g y2g | x2g z1g y1g x1g */
    t9          = gmx_mm256_unpack128lo_ps(t12, t4);                  /* y3d x3d z2d y2d | x2d z1d y1d x1d */
    t10         = gmx_mm256_unpack128hi_ps(t12, t4);                  /* y3h x3h z2h y2h | x2h z1h y1h x1h */

    tA          = _mm256_sub_ps(tA, t5);
    tB          = _mm256_sub_ps(tB, t7);
    tC          = _mm256_sub_ps(tC, t1);
    tD          = _mm256_sub_ps(tD, t9);
    tE          = _mm256_sub_ps(tE, t6);
    tF          = _mm256_sub_ps(tF, t8);
    tG          = _mm256_sub_ps(tG, t2);
    tH          = _mm256_sub_ps(tH, t10);

    _mm256_storeu_ps(ptrA, tA);
    _mm256_storeu_ps(ptrB, tB);
    _mm256_storeu_ps(ptrC, tC);
    _mm256_storeu_ps(ptrD, tD);
    _mm256_storeu_ps(ptrE, tE);
    _mm256_storeu_ps(ptrF, tF);
    _mm256_storeu_ps(ptrG, tG);
    _mm256_storeu_ps(ptrH, tH);

    tI          = gmx_mm256_set_m128(_mm_load_ss(ptrE+8), _mm_load_ss(ptrA+8));
    tJ          = gmx_mm256_set_m128(_mm_load_ss(ptrF+8), _mm_load_ss(ptrB+8));
    tK          = gmx_mm256_set_m128(_mm_load_ss(ptrG+8), _mm_load_ss(ptrC+8));
    tL          = gmx_mm256_set_m128(_mm_load_ss(ptrH+8), _mm_load_ss(ptrD+8));

    tI          = _mm256_unpacklo_ps(tI, tK);  /*  -  - zG zE |  -  - zC zA */
    tJ          = _mm256_unpacklo_ps(tJ, tL);  /*  -  - zH zF |  -  - zD zB */
    tI          = _mm256_unpacklo_ps(tI, tJ);  /* zH zG zF zE | zD zC zB zA */

    tI          = _mm256_sub_ps(tI, z3);
    tJ          = _mm256_permute_ps(tI, _MM_SHUFFLE(1, 1, 1, 1));
    tK          = _mm256_permute_ps(tI, _MM_SHUFFLE(2, 2, 2, 2));
    tL          = _mm256_permute_ps(tI, _MM_SHUFFLE(3, 3, 3, 3));

    _mm_store_ss(ptrA+8, _mm256_castps256_ps128(tI));
    _mm_store_ss(ptrB+8, _mm256_castps256_ps128(tJ));
    _mm_store_ss(ptrC+8, _mm256_castps256_ps128(tK));
    _mm_store_ss(ptrD+8, _mm256_castps256_ps128(tL));
    _mm_store_ss(ptrE+8, _mm256_extractf128_ps(tI, 0x1));
    _mm_store_ss(ptrF+8, _mm256_extractf128_ps(tJ, 0x1));
    _mm_store_ss(ptrG+8, _mm256_extractf128_ps(tK, 0x1));
    _mm_store_ss(ptrH+8, _mm256_extractf128_ps(tL, 0x1));
}


static gmx_inline void gmx_simdcall
gmx_mm256_decrement_4rvec_8ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                          float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                          float * gmx_restrict ptrE, float * gmx_restrict ptrF,
                                          float * gmx_restrict ptrG, float * gmx_restrict ptrH,
                                          __m256 x1, __m256 y1, __m256 z1,
                                          __m256 x2, __m256 y2, __m256 z2,
                                          __m256 x3, __m256 y3, __m256 z3,
                                          __m256 x4, __m256 y4, __m256 z4)
{
    __m256 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
    __m256 tA, tB, tC, tD, tE, tF, tG, tH;
    __m256 tI, tJ, tK, tL;

    tA          = _mm256_loadu_ps(ptrA);
    tB          = _mm256_loadu_ps(ptrB);
    tC          = _mm256_loadu_ps(ptrC);
    tD          = _mm256_loadu_ps(ptrD);
    tE          = _mm256_loadu_ps(ptrE);
    tF          = _mm256_loadu_ps(ptrF);
    tG          = _mm256_loadu_ps(ptrG);
    tH          = _mm256_loadu_ps(ptrH);

    t1          = _mm256_unpacklo_ps(x1, y1);                         /* y1f x1f y1e x1e | y1b x1b y1a x1a */
    t2          = _mm256_unpackhi_ps(x1, y1);                         /* y1h x1h y1g x1g | y1d x1d y1c x1c */
    t3          = _mm256_unpacklo_ps(z1, x2);                         /* x2f z1f x2e z1e | x2b z1b x2a z1a */
    t4          = _mm256_unpackhi_ps(z1, x2);                         /* x2h z1h x2g z1g | x2d z1d x2c z1c */

    t5          = _mm256_unpacklo_ps(y2, z2);                         /* z2f y2f z2e y2e | z2b y2b z2a y2a */
    t6          = _mm256_unpackhi_ps(y2, z2);                         /* z2h y2h z2g y2g | z2d y2d z2c y2c */
    t7          = _mm256_unpacklo_ps(x3, y3);                         /* y3f x3f y3e x3e | y3b x3b y3a x3a */
    t8          = _mm256_unpackhi_ps(x3, y3);                         /* y3h x3h y3g x3g | y3d x3d y3c x3c */

    t9          = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0)); /* x2e z1e y1e x1e | x2a z1a y1a x1a */
    t10         = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2)); /* x2f z1f y1f x1f | x2b z1b y1b x1b */
    t11         = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 0, 1, 0)); /* x2g z1g y1g x1g | x2c z1c y1c x1c */
    t12         = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(3, 2, 3, 2)); /* x2h z1h y1h x1h | x2d z1d y1d x1d */

    t1          = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0)); /* y3e x3e z2e y2e | y3a x3a z2a y2a */
    t2          = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2)); /* y3f x3f z2f y2f | y3b x3b z2b y2b */
    t3          = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(1, 0, 1, 0)); /* y3g x3g z2g y2g | y3c x3c z2c y2c */
    t4          = _mm256_shuffle_ps(t6, t8, _MM_SHUFFLE(3, 2, 3, 2)); /* y3h x3h z2h y2h | y3d x3d z2d y2d */

    t5          = gmx_mm256_unpack128lo_ps(t9, t1);                   /* y3a x3a z2a y2a | x2a z1a y1a x1a */
    t6          = gmx_mm256_unpack128hi_ps(t9, t1);                   /* y3e x3e z2e y2e | x2e z1e y1e x1e */
    t7          = gmx_mm256_unpack128lo_ps(t10, t2);                  /* y3b x3b z2b y2b | x2b z1b y1b x1b */
    t8          = gmx_mm256_unpack128hi_ps(t10, t2);                  /* y3f x3f z2f y2f | x2f z1f y1f x1f */
    t1          = gmx_mm256_unpack128lo_ps(t11, t3);                  /* y3c x3c z2c y2c | x2c z1c y1c x1c */
    t2          = gmx_mm256_unpack128hi_ps(t11, t3);                  /* y3g x3g z2g y2g | x2g z1g y1g x1g */
    t9          = gmx_mm256_unpack128lo_ps(t12, t4);                  /* y3d x3d z2d y2d | x2d z1d y1d x1d */
    t10         = gmx_mm256_unpack128hi_ps(t12, t4);                  /* y3h x3h z2h y2h | x2h z1h y1h x1h */

    tA          = _mm256_sub_ps(tA, t5);
    tB          = _mm256_sub_ps(tB, t7);
    tC          = _mm256_sub_ps(tC, t1);
    tD          = _mm256_sub_ps(tD, t9);
    tE          = _mm256_sub_ps(tE, t6);
    tF          = _mm256_sub_ps(tF, t8);
    tG          = _mm256_sub_ps(tG, t2);
    tH          = _mm256_sub_ps(tH, t10);

    _mm256_storeu_ps(ptrA, tA);
    _mm256_storeu_ps(ptrB, tB);
    _mm256_storeu_ps(ptrC, tC);
    _mm256_storeu_ps(ptrD, tD);
    _mm256_storeu_ps(ptrE, tE);
    _mm256_storeu_ps(ptrF, tF);
    _mm256_storeu_ps(ptrG, tG);
    _mm256_storeu_ps(ptrH, tH);

    tI          = gmx_mm256_set_m128(_mm_loadu_ps(ptrE+8), _mm_loadu_ps(ptrA+8));
    tJ          = gmx_mm256_set_m128(_mm_loadu_ps(ptrF+8), _mm_loadu_ps(ptrB+8));
    tK          = gmx_mm256_set_m128(_mm_loadu_ps(ptrG+8), _mm_loadu_ps(ptrC+8));
    tL          = gmx_mm256_set_m128(_mm_loadu_ps(ptrH+8), _mm_loadu_ps(ptrD+8));

    t1          = _mm256_unpacklo_ps(z3, x4);                         /* x4f z3f x4e z3e | x4b z3b x4a z3a */
    t2          = _mm256_unpackhi_ps(z3, x4);                         /* x4h z3h x4g z3g | x4d z3d x4c z3c */
    t3          = _mm256_unpacklo_ps(y4, z4);                         /* z4f y4f z4e y4e | z4b y4b z4a y4a */
    t4          = _mm256_unpackhi_ps(y4, z4);                         /* z4h y4h z4g y4g | z4d y4d z4c y4c */

    t5          = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0)); /* z4e y4e x4e z3e | z4a y4a x4a z3a */
    t6          = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2)); /* z4f y4f x4f z3f | z4b y4b x4b z3b */
    t7          = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(1, 0, 1, 0)); /* z4g y4g x4g z3g | z4c y4c x4c z3c */
    t8          = _mm256_shuffle_ps(t2, t4, _MM_SHUFFLE(3, 2, 3, 2)); /* z4h y4h x4h z3h | z4d y4d x4d z3d */

    tI          = _mm256_sub_ps(tI, t5);
    tJ          = _mm256_sub_ps(tJ, t6);
    tK          = _mm256_sub_ps(tK, t7);
    tL          = _mm256_sub_ps(tL, t8);

    _mm_storeu_ps(ptrA+8, _mm256_castps256_ps128(tI));
    _mm_storeu_ps(ptrB+8, _mm256_castps256_ps128(tJ));
    _mm_storeu_ps(ptrC+8, _mm256_castps256_ps128(tK));
    _mm_storeu_ps(ptrD+8, _mm256_castps256_ps128(tL));
    _mm_storeu_ps(ptrE+8, _mm256_extractf128_ps(tI, 0x1));
    _mm_storeu_ps(ptrF+8, _mm256_extractf128_ps(tJ, 0x1));
    _mm_storeu_ps(ptrG+8, _mm256_extractf128_ps(tK, 0x1));
    _mm_storeu_ps(ptrH+8, _mm256_extractf128_ps(tL, 0x1));
}


static gmx_inline void gmx_simdcall
gmx_mm256_update_iforce_1atom_swizzle_ps(__m256 fix1, __m256 fiy1, __m256 fiz1,
                                         float * gmx_restrict fptr,
                                         float * gmx_restrict fshiftptr)
{
    __m128 t1, t2, t3;

    fix1 = _mm256_hadd_ps(fix1, fix1);
    fiy1 = _mm256_hadd_ps(fiy1, fiz1);
    fix1 = _mm256_hadd_ps(fix1, fiy1); /* fiz1 fiy1 fix1 fix1 (in both lanes) */

    /* Add across the two lanes */
    t1   = _mm_add_ps(_mm256_castps256_ps128(fix1), _mm256_extractf128_ps(fix1, 0x1));

    t2 = _mm_load_ss(fptr);
    t2 = _mm_loadh_pi(t2, (__m64 *)(fptr+1));
    t3 = _mm_load_ss(fshiftptr);
    t3 = _mm_loadh_pi(t3, (__m64 *)(fshiftptr+1));

    t2 = _mm_add_ps(t2, t1);
    t3 = _mm_add_ps(t3, t1);

    _mm_store_ss(fptr, t2);
    _mm_storeh_pi((__m64 *)(fptr+1), t2);
    _mm_store_ss(fshiftptr, t3);
    _mm_storeh_pi((__m64 *)(fshiftptr+1), t3);
}


static gmx_inline void gmx_simdcall
gmx_mm256_update_iforce_3atom_swizzle_ps(__m256 fix1, __m256 fiy1, __m256 fiz1,
                                         __m256 fix2, __m256 fiy2, __m256 fiz2,
                                         __m256 fix3, __m256 fiy3, __m256 fiz3,
                                         float * gmx_restrict fptr,
                                         float * gmx_restrict fshiftptr)
{
    __m256 t1, t2, t3;
    __m128 tA, tB, tC;

    fix1 = _mm256_hadd_ps(fix1, fiy1);                /*  Y1g+Y1h Y1e+Y1f X1g+X1h X1e+X1f | Y1c+Y1d Y1a+Y1b X1c+X1d X1a+X1b */
    fiz1 = _mm256_hadd_ps(fiz1, fix2);                /*  X2g+X2h X2e+X2f Z1g+Z1h Z1e+Z1f | X2c+X2d X2a+X2b Z1c+Z1d Z1a+Z1b */
    fiy2 = _mm256_hadd_ps(fiy2, fiz2);                /*  Z2g+Z2h Z2e+Z2f Y2g+Y2h Y2e+Y2f | Z2c+Z2d Z2a+Z2b Y2c+Y2d Y2a+Y2b */
    fix3 = _mm256_hadd_ps(fix3, fiy3);                /*  Y3g+Y3h Y3e+Y3f X3g+X3h X3e+X3f | Y3c+Y3d Y3a+Y3b X3c+X3d X3a+X3b */
    fiz3 = _mm256_hadd_ps(fiz3, _mm256_setzero_ps()); /*  0       0       Z3g+Z3h Z3e+Z3f | 0       0       Z3c+Z3d Z3a+Z3b */

    fix1 = _mm256_hadd_ps(fix1, fiz1);                /*  X2e-h   Z1e-h   Y1e-h   X1e-h   | X2a-d   Z1a-d   Y1a-d   X1a-d   */
    fiy2 = _mm256_hadd_ps(fiy2, fix3);                /*  Y3e-h   X3e-h   Z2e-h   Y2e-h   | Y3a-d   X3a-d   Z2a-d   Y2a-d   */
    fiz3 = _mm256_hadd_ps(fiz3, _mm256_setzero_ps()); /*  0       0       0       Z3e-h   | 0       0       0       Z3a-d   */

    /* Add across the two lanes by swapping and adding back */
    t1   = gmx_mm256_unpack128lo_ps(fix1, fiy2);                                       /*  Y3a-d   X3a-d   Z2a-d   Y2a-d | X2a-d   Z1a-d   Y1a-d   X1a-d */
    t2   = gmx_mm256_unpack128hi_ps(fix1, fiy2);                                       /*  Y3e-h   X3e-h   Z2e-h   Y2e-h | X2e-h   Z1e-h   Y1e-h   X1e-h */
    t1   = _mm256_add_ps(t1, t2);                                                      /* y3 x3 z2 y2 | x2 z1 y1 x1 */

    tA   = _mm_add_ps(_mm256_castps256_ps128(fiz3), _mm256_extractf128_ps(fiz3, 0x1)); /* 0 0 0 z3 */

    t3   = _mm256_loadu_ps(fptr);
    t3   = _mm256_add_ps(t3, t1);
    _mm256_storeu_ps(fptr, t3);
    tB   = _mm_load_ss(fptr+8);
    tB   = _mm_add_ss(tB, tA);
    _mm_store_ss(fptr+8, tB);

    /* Add up shift force */
    tB   = _mm256_extractf128_ps(t1, 0x1);                                          /* y3 x3 z2 y2 */
    tC   = _mm_shuffle_ps(_mm256_castps256_ps128(t1), tB, _MM_SHUFFLE(1, 0, 3, 3)); /* z2 y2 x2 x2 */
    tB   = _mm_shuffle_ps(tB, tA, _MM_SHUFFLE(1, 0, 3, 2));                         /* 0 z3 y3 x3 */
    tC   = _mm_permute_ps(tC, _MM_SHUFFLE(3, 3, 2, 0));                             /*  - z2 y2 x2 */

    tB   = _mm_add_ps(tB, _mm256_castps256_ps128(t1));
    tA   = _mm_add_ps(tB, tC);                      /*  - z y x */

    tA   = _mm_blend_ps(_mm_setzero_ps(), tA, 0x7); /* 0 z y x */

    tC   = _mm_loadu_ps(fshiftptr);
    tC   = _mm_add_ps(tC, tA);
    _mm_storeu_ps(fshiftptr, tC);
}


static gmx_inline void gmx_simdcall
gmx_mm256_update_iforce_4atom_swizzle_ps(__m256 fix1, __m256 fiy1, __m256 fiz1,
                                         __m256 fix2, __m256 fiy2, __m256 fiz2,
                                         __m256 fix3, __m256 fiy3, __m256 fiz3,
                                         __m256 fix4, __m256 fiy4, __m256 fiz4,
                                         float * gmx_restrict fptr,
                                         float * gmx_restrict fshiftptr)
{
    __m256 t1, t2, t3;
    __m128 tA, tB, tC;

    fix1 = _mm256_hadd_ps(fix1, fiy1);                /*  Y1g+Y1h Y1e+Y1f X1g+X1h X1e+X1f | Y1c+Y1d Y1a+Y1b X1c+X1d X1a+X1b */
    fiz1 = _mm256_hadd_ps(fiz1, fix2);                /*  X2g+X2h X2e+X2f Z1g+Z1h Z1e+Z1f | X2c+X2d X2a+X2b Z1c+Z1d Z1a+Z1b */
    fiy2 = _mm256_hadd_ps(fiy2, fiz2);                /*  Z2g+Z2h Z2e+Z2f Y2g+Y2h Y2e+Y2f | Z2c+Z2d Z2a+Z2b Y2c+Y2d Y2a+Y2b */
    fix3 = _mm256_hadd_ps(fix3, fiy3);                /*  Y3g+Y3h Y3e+Y3f X3g+X3h X3e+X3f | Y3c+Y3d Y3a+Y3b X3c+X3d X3a+X3b */
    fiz3 = _mm256_hadd_ps(fiz3, fix4);                /*  X4g+X4h X4e+X4f Z3g+Z3h Z3e+Z3f | X4c+X4d X4a+X4b Z3c+Z3d Z3a+Z3b */
    fiy4 = _mm256_hadd_ps(fiy4, fiz4);                /*  Z4g+Z4h Z4e+Z4f Y4g+Y4h Y4e+Y4f | Z4c+Z4d Z4a+Z4b Y4c+Y4d Y4a+Y4b */

    fix1 = _mm256_hadd_ps(fix1, fiz1);                /*  X2e-h   Z1e-h   Y1e-h   X1e-h   | X2a-d   Z1a-d   Y1a-d   X1a-d   */
    fiy2 = _mm256_hadd_ps(fiy2, fix3);                /*  Y3e-h   X3e-h   Z2e-h   Y2e-h   | Y3a-d   X3a-d   Z2a-d   Y2a-d   */
    fiz3 = _mm256_hadd_ps(fiz3, fiy4);                /*  Z4e-h   Y4e-h   X4e-h   Z3e-h   | Z4a-d   Y4a-d   X4a-d   Z3a-d   */

    /* Add across the two lanes by swapping and adding back */
    t1   = gmx_mm256_unpack128lo_ps(fix1, fiy2);                                       /*  Y3a-d   X3a-d   Z2a-d   Y2a-d | X2a-d   Z1a-d   Y1a-d   X1a-d */
    t2   = gmx_mm256_unpack128hi_ps(fix1, fiy2);                                       /*  Y3e-h   X3e-h   Z2e-h   Y2e-h | X2e-h   Z1e-h   Y1e-h   X1e-h */
    t1   = _mm256_add_ps(t1, t2);                                                      /* y3 x3 z2 y2 | x2 z1 y1 x1 */

    tA   = _mm_add_ps(_mm256_castps256_ps128(fiz3), _mm256_extractf128_ps(fiz3, 0x1)); /* z4 y4 x4 z3 */

    t3   = _mm256_loadu_ps(fptr);
    t3   = _mm256_add_ps(t3, t1);
    _mm256_storeu_ps(fptr, t3);

    tB   = _mm_loadu_ps(fptr+8);
    tB   = _mm_add_ps(tB, tA);
    _mm_storeu_ps(fptr+8, tB);

    /* Add up shift force */
    tB   = _mm256_extractf128_ps(t1, 0x1);                                          /* y3 x3 z2 y2 */
    tC   = _mm_shuffle_ps(_mm256_castps256_ps128(t1), tB, _MM_SHUFFLE(1, 0, 3, 3)); /* z2 y2 x2 x2 */
    tB   = _mm_shuffle_ps(tB, tA, _MM_SHUFFLE(1, 0, 3, 2));                         /* 0 z3 y3 x3 */
    tC   = _mm_permute_ps(tC, _MM_SHUFFLE(3, 3, 2, 0));                             /*  - z2 y2 x2 */
    tA   = _mm_permute_ps(tA, _MM_SHUFFLE(0, 3, 2, 1));                             /* - z4 y4 x4 */

    tB   = _mm_add_ps(tB, _mm256_castps256_ps128(t1));
    tA   = _mm_add_ps(tA, tC);
    tA   = _mm_add_ps(tA, tB);

    tA   = _mm_blend_ps(_mm_setzero_ps(), tA, 0x7); /* 0 z y x */

    tC   = _mm_loadu_ps(fshiftptr);
    tC   = _mm_add_ps(tC, tA);
    _mm_storeu_ps(fshiftptr, tC);
}


static gmx_inline void gmx_simdcall
gmx_mm256_update_1pot_ps(__m256 pot1, float * gmx_restrict ptrA)
{
    __m128 t1;

    pot1 = _mm256_hadd_ps(pot1, pot1);
    pot1 = _mm256_hadd_ps(pot1, pot1);

    t1   = _mm_add_ps(_mm256_castps256_ps128(pot1), _mm256_extractf128_ps(pot1, 0x1));

    _mm_store_ss(ptrA, _mm_add_ss(_mm_load_ss(ptrA), t1));
}

static gmx_inline void gmx_simdcall
gmx_mm256_update_2pot_ps(__m256 pot1, float * gmx_restrict ptrA,
                         __m256 pot2, float * gmx_restrict ptrB)
{
    __m128 t1, t2;

    pot1 = _mm256_hadd_ps(pot1, pot2);
    pot1 = _mm256_hadd_ps(pot1, pot1);

    t1   = _mm_add_ps(_mm256_castps256_ps128(pot1), _mm256_extractf128_ps(pot1, 0x1));

    t2   = _mm_permute_ps(t1, _MM_SHUFFLE(1, 1, 1, 1));
    _mm_store_ss(ptrA, _mm_add_ss(_mm_load_ss(ptrA), t1));
    _mm_store_ss(ptrB, _mm_add_ss(_mm_load_ss(ptrB), t2));
}


#endif /* _kernelutil_x86_avx_256_single_h_ */
