/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifndef _gmx_math_x86_avx_256_single_h_
#define _gmx_math_x86_avx_256_single_h_

#include <math.h>

#include "gmx_x86_avx_256.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif



/************************
 *                      *
 * Simple math routines *
 *                      *
 ************************/

/* 1.0/sqrt(x), 256-bit wide version */
static gmx_inline __m256
gmx_mm256_invsqrt_ps(__m256 x)
{
    const __m256 half  = _mm256_set1_ps(0.5f);
    const __m256 three = _mm256_set1_ps(3.0f);

    __m256       lu = _mm256_rsqrt_ps(x);

    return _mm256_mul_ps(half, _mm256_mul_ps(_mm256_sub_ps(three, _mm256_mul_ps(_mm256_mul_ps(lu, lu), x)), lu));
}

/* 1.0/sqrt(x), 128-bit wide version */
static gmx_inline __m128
gmx_mm_invsqrt_ps(__m128 x)
{
    const __m128 half  = _mm_set_ps(0.5, 0.5, 0.5, 0.5);
    const __m128 three = _mm_set_ps(3.0, 3.0, 3.0, 3.0);

    __m128       lu = _mm_rsqrt_ps(x);

    return _mm_mul_ps(half, _mm_mul_ps(_mm_sub_ps(three, _mm_mul_ps(_mm_mul_ps(lu, lu), x)), lu));
}


/* sqrt(x) (256 bit) - Do NOT use this (but rather invsqrt) if you actually need 1.0/sqrt(x) */
static gmx_inline __m256
gmx_mm256_sqrt_ps(__m256 x)
{
    __m256 mask;
    __m256 res;

    mask = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_EQ_OQ);
    res  = _mm256_andnot_ps(mask, gmx_mm256_invsqrt_ps(x));

    res  = _mm256_mul_ps(x, res);

    return res;
}

/* sqrt(x) (128 bit) - Do NOT use this (but rather invsqrt) if you actually need 1.0/sqrt(x) */
static gmx_inline __m128
gmx_mm_sqrt_ps(__m128 x)
{
    __m128 mask;
    __m128 res;

    mask = _mm_cmpeq_ps(x, _mm_setzero_ps());
    res  = _mm_andnot_ps(mask, gmx_mm_invsqrt_ps(x));

    res  = _mm_mul_ps(x, res);

    return res;
}


/* 1.0/x, 256-bit wide */
static gmx_inline __m256
gmx_mm256_inv_ps(__m256 x)
{
    const __m256 two = _mm256_set1_ps(2.0f);

    __m256       lu = _mm256_rcp_ps(x);

    return _mm256_mul_ps(lu, _mm256_sub_ps(two, _mm256_mul_ps(lu, x)));
}

/* 1.0/x, 128-bit wide */
static gmx_inline __m128
gmx_mm_inv_ps(__m128 x)
{
    const __m128 two = _mm_set_ps(2.0f, 2.0f, 2.0f, 2.0f);

    __m128       lu = _mm_rcp_ps(x);

    return _mm_mul_ps(lu, _mm_sub_ps(two, _mm_mul_ps(lu, x)));
}


static gmx_inline __m256
gmx_mm256_abs_ps(__m256 x)
{
    const __m256 signmask  = _mm256_castsi256_ps( _mm256_set1_epi32(0x7FFFFFFF) );

    return _mm256_and_ps(x, signmask);
}

static gmx_inline __m128
gmx_mm_abs_ps(__m128 x)
{
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );

    return _mm_and_ps(x, signmask);
}


static __m256
gmx_mm256_log_ps(__m256 x)
{
    const __m256  expmask    = _mm256_castsi256_ps( _mm256_set1_epi32(0x7F800000) );
    const __m128i expbase_m1 = _mm_set1_epi32(127-1); /* We want non-IEEE format */
    const __m256  half       = _mm256_set1_ps(0.5f);
    const __m256  one        = _mm256_set1_ps(1.0f);
    const __m256  invsq2     = _mm256_set1_ps(1.0f/sqrt(2.0f));
    const __m256  corr1      = _mm256_set1_ps(-2.12194440e-4f);
    const __m256  corr2      = _mm256_set1_ps(0.693359375f);

    const __m256  CA_1        = _mm256_set1_ps(0.070376836292f);
    const __m256  CB_0        = _mm256_set1_ps(1.6714950086782716f);
    const __m256  CB_1        = _mm256_set1_ps(-2.452088066061482f);
    const __m256  CC_0        = _mm256_set1_ps(1.5220770854701728f);
    const __m256  CC_1        = _mm256_set1_ps(-1.3422238433233642f);
    const __m256  CD_0        = _mm256_set1_ps(1.386218787509749f);
    const __m256  CD_1        = _mm256_set1_ps(0.35075468953796346f);
    const __m256  CE_0        = _mm256_set1_ps(1.3429983063133937f);
    const __m256  CE_1        = _mm256_set1_ps(1.807420826584643f);

    __m256        fexp, fexp1;
    __m256i       iexp;
    __m128i       iexp128a, iexp128b;
    __m256        mask;
    __m256i       imask;
    __m128i       imask128a, imask128b;
    __m256        x1, x2;
    __m256        y;
    __m256        pA, pB, pC, pD, pE, tB, tC, tD, tE;

    /* Separate x into exponent and mantissa, with a mantissa in the range [0.5..1[ (not IEEE754 standard!) */
    fexp  = _mm256_and_ps(x, expmask);
    iexp  = _mm256_castps_si256(fexp);

    iexp128b = _mm256_extractf128_si256(iexp, 0x1);
    iexp128a = _mm256_castsi256_si128(iexp);

    iexp128a  = _mm_srli_epi32(iexp128a, 23);
    iexp128b  = _mm_srli_epi32(iexp128b, 23);
    iexp128a  = _mm_sub_epi32(iexp128a, expbase_m1);
    iexp128b  = _mm_sub_epi32(iexp128b, expbase_m1);

    x     = _mm256_andnot_ps(expmask, x);
    x     = _mm256_or_ps(x, one);
    x     = _mm256_mul_ps(x, half);

    mask  = _mm256_cmp_ps(x, invsq2, _CMP_LT_OQ);

    x     = _mm256_add_ps(x, _mm256_and_ps(mask, x));
    x     = _mm256_sub_ps(x, one);

    imask = _mm256_castps_si256(mask);

    imask128b = _mm256_extractf128_si256(imask, 0x1);
    imask128a = _mm256_castsi256_si128(imask);

    iexp128a  = _mm_add_epi32(iexp128a, imask128a);
    iexp128b  = _mm_add_epi32(iexp128b, imask128b);

    iexp  = _mm256_castsi128_si256(iexp128a);
    iexp  = _mm256_insertf128_si256(iexp, iexp128b, 0x1);

    x2    = _mm256_mul_ps(x, x);

    pA    = _mm256_mul_ps(CA_1, x);
    pB    = _mm256_mul_ps(CB_1, x);
    pC    = _mm256_mul_ps(CC_1, x);
    pD    = _mm256_mul_ps(CD_1, x);
    pE    = _mm256_mul_ps(CE_1, x);
    tB    = _mm256_add_ps(CB_0, x2);
    tC    = _mm256_add_ps(CC_0, x2);
    tD    = _mm256_add_ps(CD_0, x2);
    tE    = _mm256_add_ps(CE_0, x2);
    pB    = _mm256_add_ps(pB, tB);
    pC    = _mm256_add_ps(pC, tC);
    pD    = _mm256_add_ps(pD, tD);
    pE    = _mm256_add_ps(pE, tE);

    pA    = _mm256_mul_ps(pA, pB);
    pC    = _mm256_mul_ps(pC, pD);
    pE    = _mm256_mul_ps(pE, x2);
    pA    = _mm256_mul_ps(pA, pC);
    y     = _mm256_mul_ps(pA, pE);

    fexp  = _mm256_cvtepi32_ps(iexp);
    y     = _mm256_add_ps(y, _mm256_mul_ps(fexp, corr1));

    y     = _mm256_sub_ps(y, _mm256_mul_ps(half, x2));
    x2    = _mm256_add_ps(x, y);

    x2    = _mm256_add_ps(x2, _mm256_mul_ps(fexp, corr2));

    return x2;
}


static __m128
gmx_mm_log_ps(__m128 x)
{
    /* Same algorithm as cephes library */
    const __m128  expmask    = gmx_mm_castsi128_ps( _mm_set_epi32(0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000) );
    const __m128i expbase_m1 = _mm_set1_epi32(127-1); /* We want non-IEEE format */
    const __m128  half       = _mm_set1_ps(0.5f);
    const __m128  one        = _mm_set1_ps(1.0f);
    const __m128  invsq2     = _mm_set1_ps(1.0f/sqrt(2.0f));
    const __m128  corr1      = _mm_set1_ps(-2.12194440e-4f);
    const __m128  corr2      = _mm_set1_ps(0.693359375f);

    const __m128  CA_1        = _mm_set1_ps(0.070376836292f);
    const __m128  CB_0        = _mm_set1_ps(1.6714950086782716f);
    const __m128  CB_1        = _mm_set1_ps(-2.452088066061482f);
    const __m128  CC_0        = _mm_set1_ps(1.5220770854701728f);
    const __m128  CC_1        = _mm_set1_ps(-1.3422238433233642f);
    const __m128  CD_0        = _mm_set1_ps(1.386218787509749f);
    const __m128  CD_1        = _mm_set1_ps(0.35075468953796346f);
    const __m128  CE_0        = _mm_set1_ps(1.3429983063133937f);
    const __m128  CE_1        = _mm_set1_ps(1.807420826584643f);

    __m128        fexp, fexp1;
    __m128i       iexp;
    __m128        mask;
    __m128        x1, x2;
    __m128        y;
    __m128        pA, pB, pC, pD, pE, tB, tC, tD, tE;

    /* Separate x into exponent and mantissa, with a mantissa in the range [0.5..1[ (not IEEE754 standard!) */
    fexp  = _mm_and_ps(x, expmask);
    iexp  = gmx_mm_castps_si128(fexp);
    iexp  = _mm_srli_epi32(iexp, 23);
    iexp  = _mm_sub_epi32(iexp, expbase_m1);

    x     = _mm_andnot_ps(expmask, x);
    x     = _mm_or_ps(x, one);
    x     = _mm_mul_ps(x, half);

    mask  = _mm_cmplt_ps(x, invsq2);

    x     = _mm_add_ps(x, _mm_and_ps(mask, x));
    x     = _mm_sub_ps(x, one);
    iexp  = _mm_add_epi32(iexp, gmx_mm_castps_si128(mask)); /* 0xFFFFFFFF = -1 as int */

    x2    = _mm_mul_ps(x, x);

    pA    = _mm_mul_ps(CA_1, x);
    pB    = _mm_mul_ps(CB_1, x);
    pC    = _mm_mul_ps(CC_1, x);
    pD    = _mm_mul_ps(CD_1, x);
    pE    = _mm_mul_ps(CE_1, x);
    tB    = _mm_add_ps(CB_0, x2);
    tC    = _mm_add_ps(CC_0, x2);
    tD    = _mm_add_ps(CD_0, x2);
    tE    = _mm_add_ps(CE_0, x2);
    pB    = _mm_add_ps(pB, tB);
    pC    = _mm_add_ps(pC, tC);
    pD    = _mm_add_ps(pD, tD);
    pE    = _mm_add_ps(pE, tE);

    pA    = _mm_mul_ps(pA, pB);
    pC    = _mm_mul_ps(pC, pD);
    pE    = _mm_mul_ps(pE, x2);
    pA    = _mm_mul_ps(pA, pC);
    y     = _mm_mul_ps(pA, pE);

    fexp  = _mm_cvtepi32_ps(iexp);
    y     = _mm_add_ps(y, _mm_mul_ps(fexp, corr1));

    y     = _mm_sub_ps(y, _mm_mul_ps(half, x2));
    x2    = _mm_add_ps(x, y);

    x2    = _mm_add_ps(x2, _mm_mul_ps(fexp, corr2));

    return x2;
}


/*
 * 2^x function, 256-bit wide
 *
 * The 2^w term is calculated from a (6,0)-th order (no denominator) Minimax polynomia on the interval
 * [-0.5,0.5]. The coefficiencts of this was derived in Mathematica using the command:
 *
 * MiniMaxApproximation[(2^x), {x, {-0.5, 0.5}, 6, 0}, WorkingPrecision -> 15]
 *
 * The largest-magnitude exponent we can represent in IEEE single-precision binary format
 * is 2^-126 for small numbers and 2^127 for large ones. To avoid wrap-around problems, we set the
 * result to zero if the argument falls outside this range. For small numbers this is just fine, but
 * for large numbers you could be fancy and return the smallest/largest IEEE single-precision
 * number instead. That would take a few extra cycles and not really help, since something is
 * wrong if you are using single precision to work with numbers that cannot really be represented
 * in single precision.
 *
 * The accuracy is at least 23 bits.
 */
static __m256
gmx_mm256_exp2_ps(__m256 x)
{
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const __m256  arglimit = _mm256_set1_ps(126.0f);

    const __m128i expbase  = _mm_set1_epi32(127);
    const __m256  CC6      = _mm256_set1_ps(1.535336188319500E-004);
    const __m256  CC5      = _mm256_set1_ps(1.339887440266574E-003);
    const __m256  CC4      = _mm256_set1_ps(9.618437357674640E-003);
    const __m256  CC3      = _mm256_set1_ps(5.550332471162809E-002);
    const __m256  CC2      = _mm256_set1_ps(2.402264791363012E-001);
    const __m256  CC1      = _mm256_set1_ps(6.931472028550421E-001);
    const __m256  CC0      = _mm256_set1_ps(1.0f);

    __m256        p0, p1;
    __m256        valuemask;
    __m256i       iexppart;
    __m128i       iexppart128a, iexppart128b;
    __m256        fexppart;
    __m256        intpart;
    __m256        x2;


    iexppart  = _mm256_cvtps_epi32(x);
    intpart   = _mm256_round_ps(x, _MM_FROUND_TO_NEAREST_INT);

    iexppart128b = _mm256_extractf128_si256(iexppart, 0x1);
    iexppart128a = _mm256_castsi256_si128(iexppart);

    iexppart128a = _mm_slli_epi32(_mm_add_epi32(iexppart128a, expbase), 23);
    iexppart128b = _mm_slli_epi32(_mm_add_epi32(iexppart128b, expbase), 23);

    iexppart  = _mm256_castsi128_si256(iexppart128a);
    iexppart  = _mm256_insertf128_si256(iexppart, iexppart128b, 0x1);
    valuemask = _mm256_cmp_ps(arglimit, gmx_mm256_abs_ps(x), _CMP_GE_OQ);
    fexppart  = _mm256_and_ps(valuemask, _mm256_castsi256_ps(iexppart));

    x         = _mm256_sub_ps(x, intpart);
    x2        = _mm256_mul_ps(x, x);

    p0        = _mm256_mul_ps(CC6, x2);
    p1        = _mm256_mul_ps(CC5, x2);
    p0        = _mm256_add_ps(p0, CC4);
    p1        = _mm256_add_ps(p1, CC3);
    p0        = _mm256_mul_ps(p0, x2);
    p1        = _mm256_mul_ps(p1, x2);
    p0        = _mm256_add_ps(p0, CC2);
    p1        = _mm256_add_ps(p1, CC1);
    p0        = _mm256_mul_ps(p0, x2);
    p1        = _mm256_mul_ps(p1, x);
    p0        = _mm256_add_ps(p0, CC0);
    p0        = _mm256_add_ps(p0, p1);
    x         = _mm256_mul_ps(p0, fexppart);

    return x;
}


/* 2^x, 128 bit wide */
static __m128
gmx_mm_exp2_ps(__m128 x)
{
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128  arglimit = _mm_set1_ps(126.0f);

    const __m128i expbase  = _mm_set1_epi32(127);
    const __m128  CA6      = _mm_set1_ps(1.535336188319500E-004);
    const __m128  CA5      = _mm_set1_ps(1.339887440266574E-003);
    const __m128  CA4      = _mm_set1_ps(9.618437357674640E-003);
    const __m128  CA3      = _mm_set1_ps(5.550332471162809E-002);
    const __m128  CA2      = _mm_set1_ps(2.402264791363012E-001);
    const __m128  CA1      = _mm_set1_ps(6.931472028550421E-001);
    const __m128  CA0      = _mm_set1_ps(1.0f);

    __m128        valuemask;
    __m128i       iexppart;
    __m128        fexppart;
    __m128        intpart;
    __m128        x2;
    __m128        p0, p1;

    iexppart  = _mm_cvtps_epi32(x);
    intpart   = _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT);
    iexppart  = _mm_slli_epi32(_mm_add_epi32(iexppart, expbase), 23);
    valuemask = _mm_cmpge_ps(arglimit, gmx_mm_abs_ps(x));
    fexppart  = _mm_and_ps(valuemask, gmx_mm_castsi128_ps(iexppart));

    x         = _mm_sub_ps(x, intpart);
    x2        = _mm_mul_ps(x, x);

    p0        = _mm_mul_ps(CA6, x2);
    p1        = _mm_mul_ps(CA5, x2);
    p0        = _mm_add_ps(p0, CA4);
    p1        = _mm_add_ps(p1, CA3);
    p0        = _mm_mul_ps(p0, x2);
    p1        = _mm_mul_ps(p1, x2);
    p0        = _mm_add_ps(p0, CA2);
    p1        = _mm_add_ps(p1, CA1);
    p0        = _mm_mul_ps(p0, x2);
    p1        = _mm_mul_ps(p1, x);
    p0        = _mm_add_ps(p0, CA0);
    p0        = _mm_add_ps(p0, p1);
    x         = _mm_mul_ps(p0, fexppart);

    return x;
}


/* Exponential function, 256 bit wide. This could be calculated from 2^x as Exp(x)=2^(y),
 * where y=log2(e)*x, but there will then be a small rounding error since we lose some
 * precision due to the multiplication. This will then be magnified a lot by the exponential.
 *
 * Instead, we calculate the fractional part directly as a minimax approximation of
 * Exp(z) on [-0.5,0.5]. We use extended precision arithmetics to calculate the fraction
 * remaining after 2^y, which avoids the precision-loss.
 * The final result is correct to within 1 LSB over the entire argument range.
 */
static __m256
gmx_mm256_exp_ps(__m256 exparg)
{
    const __m256  argscale      = _mm256_set1_ps(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const __m256  arglimit      = _mm256_set1_ps(126.0f);
    const __m128i expbase       = _mm_set1_epi32(127);

    const __m256  invargscale0  = _mm256_set1_ps(0.693359375f);
    const __m256  invargscale1  = _mm256_set1_ps(-2.12194440e-4f);

    const __m256  CE5           = _mm256_set1_ps(1.9875691500e-4f);
    const __m256  CE4           = _mm256_set1_ps(1.3981999507e-3f);
    const __m256  CE3           = _mm256_set1_ps(8.3334519073e-3f);
    const __m256  CE2           = _mm256_set1_ps(4.1665795894e-2f);
    const __m256  CE1           = _mm256_set1_ps(1.6666665459e-1f);
    const __m256  CE0           = _mm256_set1_ps(5.0000001201e-1f);
    const __m256  one           = _mm256_set1_ps(1.0f);

    __m256        exparg2, exp2arg;
    __m256        pE0, pE1;
    __m256        valuemask;
    __m256i       iexppart;
    __m128i       iexppart128a, iexppart128b;
    __m256        fexppart;
    __m256        intpart;

    exp2arg = _mm256_mul_ps(exparg, argscale);

    iexppart  = _mm256_cvtps_epi32(exp2arg);
    intpart   = _mm256_round_ps(exp2arg, _MM_FROUND_TO_NEAREST_INT);

    iexppart128b = _mm256_extractf128_si256(iexppart, 0x1);
    iexppart128a = _mm256_castsi256_si128(iexppart);

    iexppart128a = _mm_slli_epi32(_mm_add_epi32(iexppart128a, expbase), 23);
    iexppart128b = _mm_slli_epi32(_mm_add_epi32(iexppart128b, expbase), 23);

    iexppart  = _mm256_castsi128_si256(iexppart128a);
    iexppart  = _mm256_insertf128_si256(iexppart, iexppart128b, 0x1);
    valuemask = _mm256_cmp_ps(arglimit, gmx_mm256_abs_ps(exp2arg), _CMP_GE_OQ);
    fexppart  = _mm256_and_ps(valuemask, _mm256_castsi256_ps(iexppart));

    /* Extended precision arithmetics */
    exparg    = _mm256_sub_ps(exparg, _mm256_mul_ps(invargscale0, intpart));
    exparg    = _mm256_sub_ps(exparg, _mm256_mul_ps(invargscale1, intpart));

    exparg2   = _mm256_mul_ps(exparg, exparg);

    pE1       = _mm256_mul_ps(CE5, exparg2);
    pE0       = _mm256_mul_ps(CE4, exparg2);
    pE1       = _mm256_add_ps(pE1, CE3);
    pE0       = _mm256_add_ps(pE0, CE2);
    pE1       = _mm256_mul_ps(pE1, exparg2);
    pE0       = _mm256_mul_ps(pE0, exparg2);
    pE1       = _mm256_add_ps(pE1, CE1);
    pE0       = _mm256_add_ps(pE0, CE0);
    pE1       = _mm256_mul_ps(pE1, exparg);
    pE0       = _mm256_add_ps(pE0, pE1);
    pE0       = _mm256_mul_ps(pE0, exparg2);
    exparg    = _mm256_add_ps(exparg, one);
    exparg    = _mm256_add_ps(exparg, pE0);

    exparg    = _mm256_mul_ps(exparg, fexppart);

    return exparg;
}


/* exp(), 128 bit wide. */
static __m128
gmx_mm_exp_ps(__m128 x)
{
    const __m128  argscale      = _mm_set1_ps(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const __m128  arglimit      = _mm_set1_ps(126.0f);
    const __m128i expbase       = _mm_set1_epi32(127);

    const __m128  invargscale0  = _mm_set1_ps(0.693359375f);
    const __m128  invargscale1  = _mm_set1_ps(-2.12194440e-4f);

    const __m128  CC5           = _mm_set1_ps(1.9875691500e-4f);
    const __m128  CC4           = _mm_set1_ps(1.3981999507e-3f);
    const __m128  CC3           = _mm_set1_ps(8.3334519073e-3f);
    const __m128  CC2           = _mm_set1_ps(4.1665795894e-2f);
    const __m128  CC1           = _mm_set1_ps(1.6666665459e-1f);
    const __m128  CC0           = _mm_set1_ps(5.0000001201e-1f);
    const __m128  one           = _mm_set1_ps(1.0f);

    __m128        y, x2;
    __m128        p0, p1;
    __m128        valuemask;
    __m128i       iexppart;
    __m128        fexppart;
    __m128        intpart;

    y = _mm_mul_ps(x, argscale);

    iexppart  = _mm_cvtps_epi32(y);
    intpart   = _mm_round_ps(y, _MM_FROUND_TO_NEAREST_INT);

    iexppart  = _mm_slli_epi32(_mm_add_epi32(iexppart, expbase), 23);
    valuemask = _mm_cmpge_ps(arglimit, gmx_mm_abs_ps(y));
    fexppart  = _mm_and_ps(valuemask, gmx_mm_castsi128_ps(iexppart));

    /* Extended precision arithmetics */
    x         = _mm_sub_ps(x, _mm_mul_ps(invargscale0, intpart));
    x         = _mm_sub_ps(x, _mm_mul_ps(invargscale1, intpart));

    x2        = _mm_mul_ps(x, x);

    p1        = _mm_mul_ps(CC5, x2);
    p0        = _mm_mul_ps(CC4, x2);
    p1        = _mm_add_ps(p1, CC3);
    p0        = _mm_add_ps(p0, CC2);
    p1        = _mm_mul_ps(p1, x2);
    p0        = _mm_mul_ps(p0, x2);
    p1        = _mm_add_ps(p1, CC1);
    p0        = _mm_add_ps(p0, CC0);
    p1        = _mm_mul_ps(p1, x);
    p0        = _mm_add_ps(p0, p1);
    p0        = _mm_mul_ps(p0, x2);
    x         = _mm_add_ps(x, one);
    x         = _mm_add_ps(x, p0);

    x         = _mm_mul_ps(x, fexppart);

    return x;
}



/* FULL precision erf(), 256-bit wide. Only errors in LSB */
static __m256
gmx_mm256_erf_ps(__m256 x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const __m256  CA6      = _mm256_set1_ps(7.853861353153693e-5f);
    const __m256  CA5      = _mm256_set1_ps(-8.010193625184903e-4f);
    const __m256  CA4      = _mm256_set1_ps(5.188327685732524e-3f);
    const __m256  CA3      = _mm256_set1_ps(-2.685381193529856e-2f);
    const __m256  CA2      = _mm256_set1_ps(1.128358514861418e-1f);
    const __m256  CA1      = _mm256_set1_ps(-3.761262582423300e-1f);
    const __m256  CA0      = _mm256_set1_ps(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const __m256  CB9      = _mm256_set1_ps(-0.0018629930017603923f);
    const __m256  CB8      = _mm256_set1_ps(0.003909821287598495f);
    const __m256  CB7      = _mm256_set1_ps(-0.0052094582210355615f);
    const __m256  CB6      = _mm256_set1_ps(0.005685614362160572f);
    const __m256  CB5      = _mm256_set1_ps(-0.0025367682853477272f);
    const __m256  CB4      = _mm256_set1_ps(-0.010199799682318782f);
    const __m256  CB3      = _mm256_set1_ps(0.04369575504816542f);
    const __m256  CB2      = _mm256_set1_ps(-0.11884063474674492f);
    const __m256  CB1      = _mm256_set1_ps(0.2732120154030589f);
    const __m256  CB0      = _mm256_set1_ps(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const __m256  CC10     = _mm256_set1_ps(-0.0445555913112064f);
    const __m256  CC9      = _mm256_set1_ps(0.21376355144663348f);
    const __m256  CC8      = _mm256_set1_ps(-0.3473187200259257f);
    const __m256  CC7      = _mm256_set1_ps(0.016690861551248114f);
    const __m256  CC6      = _mm256_set1_ps(0.7560973182491192f);
    const __m256  CC5      = _mm256_set1_ps(-1.2137903600145787f);
    const __m256  CC4      = _mm256_set1_ps(0.8411872321232948f);
    const __m256  CC3      = _mm256_set1_ps(-0.08670413896296343f);
    const __m256  CC2      = _mm256_set1_ps(-0.27124782687240334f);
    const __m256  CC1      = _mm256_set1_ps(-0.0007502488047806069f);
    const __m256  CC0      = _mm256_set1_ps(0.5642114853803148f);

    /* Coefficients for expansion of exp(x) in [0,0.1] */
    /* CD0 and CD1 are both 1.0, so no need to declare them separately */
    const __m256  CD2      = _mm256_set1_ps(0.5000066608081202f);
    const __m256  CD3      = _mm256_set1_ps(0.1664795422874624f);
    const __m256  CD4      = _mm256_set1_ps(0.04379839977652482f);

    const __m256  sieve    = _mm256_castsi256_ps( _mm256_set1_epi32(0xfffff000) );
    const __m256  signbit  = _mm256_castsi256_ps( _mm256_set1_epi32(0x80000000) );
    const __m256  one      = _mm256_set1_ps(1.0f);
    const __m256  two      = _mm256_set1_ps(2.0f);

    __m256        x2, x4, y;
    __m256        z, q, t, t2, w, w2;
    __m256        pA0, pA1, pB0, pB1, pC0, pC1;
    __m256        expmx2, corr;
    __m256        res_erf, res_erfc, res;
    __m256        mask;

    /* Calculate erf() */
    x2     = _mm256_mul_ps(x, x);
    x4     = _mm256_mul_ps(x2, x2);

    pA0  = _mm256_mul_ps(CA6, x4);
    pA1  = _mm256_mul_ps(CA5, x4);
    pA0  = _mm256_add_ps(pA0, CA4);
    pA1  = _mm256_add_ps(pA1, CA3);
    pA0  = _mm256_mul_ps(pA0, x4);
    pA1  = _mm256_mul_ps(pA1, x4);
    pA0  = _mm256_add_ps(pA0, CA2);
    pA1  = _mm256_add_ps(pA1, CA1);
    pA0  = _mm256_mul_ps(pA0, x4);
    pA1  = _mm256_mul_ps(pA1, x2);
    pA0  = _mm256_add_ps(pA0, pA1);
    pA0  = _mm256_add_ps(pA0, CA0);

    res_erf = _mm256_mul_ps(x, pA0);

    /* Calculate erfc */

    y       = gmx_mm256_abs_ps(x);
    t       = gmx_mm256_inv_ps(y);
    w       = _mm256_sub_ps(t, one);
    t2      = _mm256_mul_ps(t, t);
    w2      = _mm256_mul_ps(w, w);
    /*
     * We cannot simply calculate exp(-x2) directly in single precision, since
     * that will lose a couple of bits of precision due to the multiplication.
     * Instead, we introduce x=z+w, where the last 12 bits of precision are in w.
     * Then we get exp(-x2) = exp(-z2)*exp((z-x)*(z+x)).
     *
     * The only drawback with this is that it requires TWO separate exponential
     * evaluations, which would be horrible performance-wise. However, the argument
     * for the second exp() call is always small, so there we simply use a
     * low-order minimax expansion on [0,0.1].
     */

    z       = _mm256_and_ps(y, sieve);
    q       = _mm256_mul_ps( _mm256_sub_ps(z, y), _mm256_add_ps(z, y) );

    corr    = _mm256_mul_ps(CD4, q);
    corr    = _mm256_add_ps(corr, CD3);
    corr    = _mm256_mul_ps(corr, q);
    corr    = _mm256_add_ps(corr, CD2);
    corr    = _mm256_mul_ps(corr, q);
    corr    = _mm256_add_ps(corr, one);
    corr    = _mm256_mul_ps(corr, q);
    corr    = _mm256_add_ps(corr, one);

    expmx2  = gmx_mm256_exp_ps( _mm256_or_ps( signbit, _mm256_mul_ps(z, z) ) );
    expmx2  = _mm256_mul_ps(expmx2, corr);

    pB1  = _mm256_mul_ps(CB9, w2);
    pB0  = _mm256_mul_ps(CB8, w2);
    pB1  = _mm256_add_ps(pB1, CB7);
    pB0  = _mm256_add_ps(pB0, CB6);
    pB1  = _mm256_mul_ps(pB1, w2);
    pB0  = _mm256_mul_ps(pB0, w2);
    pB1  = _mm256_add_ps(pB1, CB5);
    pB0  = _mm256_add_ps(pB0, CB4);
    pB1  = _mm256_mul_ps(pB1, w2);
    pB0  = _mm256_mul_ps(pB0, w2);
    pB1  = _mm256_add_ps(pB1, CB3);
    pB0  = _mm256_add_ps(pB0, CB2);
    pB1  = _mm256_mul_ps(pB1, w2);
    pB0  = _mm256_mul_ps(pB0, w2);
    pB1  = _mm256_add_ps(pB1, CB1);
    pB1  = _mm256_mul_ps(pB1, w);
    pB0  = _mm256_add_ps(pB0, pB1);
    pB0  = _mm256_add_ps(pB0, CB0);

    pC0  = _mm256_mul_ps(CC10, t2);
    pC1  = _mm256_mul_ps(CC9, t2);
    pC0  = _mm256_add_ps(pC0, CC8);
    pC1  = _mm256_add_ps(pC1, CC7);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t2);
    pC0  = _mm256_add_ps(pC0, CC6);
    pC1  = _mm256_add_ps(pC1, CC5);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t2);
    pC0  = _mm256_add_ps(pC0, CC4);
    pC1  = _mm256_add_ps(pC1, CC3);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t2);
    pC0  = _mm256_add_ps(pC0, CC2);
    pC1  = _mm256_add_ps(pC1, CC1);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t);
    pC0  = _mm256_add_ps(pC0, pC1);
    pC0  = _mm256_add_ps(pC0, CC0);
    pC0  = _mm256_mul_ps(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = _mm256_cmp_ps(two, y, _CMP_LT_OQ);
    res_erfc = _mm256_blendv_ps(pB0, pC0, mask);
    res_erfc = _mm256_mul_ps(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OQ);
    res_erfc = _mm256_blendv_ps(res_erfc, _mm256_sub_ps(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = _mm256_cmp_ps(y, _mm256_set1_ps(0.75f), _CMP_LT_OQ);
    res  = _mm256_blendv_ps(_mm256_sub_ps(one, res_erfc), res_erf, mask);

    return res;
}


/* erf(), 128 bit wide */
static __m128
gmx_mm_erf_ps(__m128 x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const __m128  CA6      = _mm_set1_ps(7.853861353153693e-5f);
    const __m128  CA5      = _mm_set1_ps(-8.010193625184903e-4f);
    const __m128  CA4      = _mm_set1_ps(5.188327685732524e-3f);
    const __m128  CA3      = _mm_set1_ps(-2.685381193529856e-2f);
    const __m128  CA2      = _mm_set1_ps(1.128358514861418e-1f);
    const __m128  CA1      = _mm_set1_ps(-3.761262582423300e-1f);
    const __m128  CA0      = _mm_set1_ps(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const __m128  CB9      = _mm_set1_ps(-0.0018629930017603923f);
    const __m128  CB8      = _mm_set1_ps(0.003909821287598495f);
    const __m128  CB7      = _mm_set1_ps(-0.0052094582210355615f);
    const __m128  CB6      = _mm_set1_ps(0.005685614362160572f);
    const __m128  CB5      = _mm_set1_ps(-0.0025367682853477272f);
    const __m128  CB4      = _mm_set1_ps(-0.010199799682318782f);
    const __m128  CB3      = _mm_set1_ps(0.04369575504816542f);
    const __m128  CB2      = _mm_set1_ps(-0.11884063474674492f);
    const __m128  CB1      = _mm_set1_ps(0.2732120154030589f);
    const __m128  CB0      = _mm_set1_ps(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const __m128  CC10     = _mm_set1_ps(-0.0445555913112064f);
    const __m128  CC9      = _mm_set1_ps(0.21376355144663348f);
    const __m128  CC8      = _mm_set1_ps(-0.3473187200259257f);
    const __m128  CC7      = _mm_set1_ps(0.016690861551248114f);
    const __m128  CC6      = _mm_set1_ps(0.7560973182491192f);
    const __m128  CC5      = _mm_set1_ps(-1.2137903600145787f);
    const __m128  CC4      = _mm_set1_ps(0.8411872321232948f);
    const __m128  CC3      = _mm_set1_ps(-0.08670413896296343f);
    const __m128  CC2      = _mm_set1_ps(-0.27124782687240334f);
    const __m128  CC1      = _mm_set1_ps(-0.0007502488047806069f);
    const __m128  CC0      = _mm_set1_ps(0.5642114853803148f);

    /* Coefficients for expansion of exp(x) in [0,0.1] */
    /* CD0 and CD1 are both 1.0, so no need to declare them separately */
    const __m128  CD2      = _mm_set1_ps(0.5000066608081202f);
    const __m128  CD3      = _mm_set1_ps(0.1664795422874624f);
    const __m128  CD4      = _mm_set1_ps(0.04379839977652482f);

    const __m128  sieve    = gmx_mm_castsi128_ps( _mm_set1_epi32(0xfffff000) );
    const __m128  signbit  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x80000000) );
    const __m128  one      = _mm_set1_ps(1.0f);
    const __m128  two      = _mm_set1_ps(2.0f);

    __m128        x2, x4, y;
    __m128        z, q, t, t2, w, w2;
    __m128        pA0, pA1, pB0, pB1, pC0, pC1;
    __m128        expmx2, corr;
    __m128        res_erf, res_erfc, res;
    __m128        mask;

    /* Calculate erf() */
    x2     = _mm_mul_ps(x, x);
    x4     = _mm_mul_ps(x2, x2);

    pA0  = _mm_mul_ps(CA6, x4);
    pA1  = _mm_mul_ps(CA5, x4);
    pA0  = _mm_add_ps(pA0, CA4);
    pA1  = _mm_add_ps(pA1, CA3);
    pA0  = _mm_mul_ps(pA0, x4);
    pA1  = _mm_mul_ps(pA1, x4);
    pA0  = _mm_add_ps(pA0, CA2);
    pA1  = _mm_add_ps(pA1, CA1);
    pA0  = _mm_mul_ps(pA0, x4);
    pA1  = _mm_mul_ps(pA1, x2);
    pA0  = _mm_add_ps(pA0, pA1);
    pA0  = _mm_add_ps(pA0, CA0);

    res_erf = _mm_mul_ps(x, pA0);

    /* Calculate erfc */

    y       = gmx_mm_abs_ps(x);
    t       = gmx_mm_inv_ps(y);
    w       = _mm_sub_ps(t, one);
    t2      = _mm_mul_ps(t, t);
    w2      = _mm_mul_ps(w, w);
    /*
     * We cannot simply calculate exp(-x2) directly in single precision, since
     * that will lose a couple of bits of precision due to the multiplication.
     * Instead, we introduce x=z+w, where the last 12 bits of precision are in w.
     * Then we get exp(-x2) = exp(-z2)*exp((z-x)*(z+x)).
     *
     * The only drawback with this is that it requires TWO separate exponential
     * evaluations, which would be horrible performance-wise. However, the argument
     * for the second exp() call is always small, so there we simply use a
     * low-order minimax expansion on [0,0.1].
     */

    z       = _mm_and_ps(y, sieve);
    q       = _mm_mul_ps( _mm_sub_ps(z, y), _mm_add_ps(z, y) );

    corr    = _mm_mul_ps(CD4, q);
    corr    = _mm_add_ps(corr, CD3);
    corr    = _mm_mul_ps(corr, q);
    corr    = _mm_add_ps(corr, CD2);
    corr    = _mm_mul_ps(corr, q);
    corr    = _mm_add_ps(corr, one);
    corr    = _mm_mul_ps(corr, q);
    corr    = _mm_add_ps(corr, one);

    expmx2  = gmx_mm_exp_ps( _mm_or_ps( signbit, _mm_mul_ps(z, z) ) );
    expmx2  = _mm_mul_ps(expmx2, corr);

    pB1  = _mm_mul_ps(CB9, w2);
    pB0  = _mm_mul_ps(CB8, w2);
    pB1  = _mm_add_ps(pB1, CB7);
    pB0  = _mm_add_ps(pB0, CB6);
    pB1  = _mm_mul_ps(pB1, w2);
    pB0  = _mm_mul_ps(pB0, w2);
    pB1  = _mm_add_ps(pB1, CB5);
    pB0  = _mm_add_ps(pB0, CB4);
    pB1  = _mm_mul_ps(pB1, w2);
    pB0  = _mm_mul_ps(pB0, w2);
    pB1  = _mm_add_ps(pB1, CB3);
    pB0  = _mm_add_ps(pB0, CB2);
    pB1  = _mm_mul_ps(pB1, w2);
    pB0  = _mm_mul_ps(pB0, w2);
    pB1  = _mm_add_ps(pB1, CB1);
    pB1  = _mm_mul_ps(pB1, w);
    pB0  = _mm_add_ps(pB0, pB1);
    pB0  = _mm_add_ps(pB0, CB0);

    pC0  = _mm_mul_ps(CC10, t2);
    pC1  = _mm_mul_ps(CC9, t2);
    pC0  = _mm_add_ps(pC0, CC8);
    pC1  = _mm_add_ps(pC1, CC7);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t2);
    pC0  = _mm_add_ps(pC0, CC6);
    pC1  = _mm_add_ps(pC1, CC5);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t2);
    pC0  = _mm_add_ps(pC0, CC4);
    pC1  = _mm_add_ps(pC1, CC3);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t2);
    pC0  = _mm_add_ps(pC0, CC2);
    pC1  = _mm_add_ps(pC1, CC1);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t);
    pC0  = _mm_add_ps(pC0, pC1);
    pC0  = _mm_add_ps(pC0, CC0);
    pC0  = _mm_mul_ps(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = _mm_cmplt_ps(two, y);
    res_erfc = _mm_blendv_ps(pB0, pC0, mask);
    res_erfc = _mm_mul_ps(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = _mm_cmplt_ps(x, _mm_setzero_ps());
    res_erfc = _mm_blendv_ps(res_erfc, _mm_sub_ps(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = _mm_cmplt_ps(y, _mm_set1_ps(0.75f));
    res  = _mm_blendv_ps(_mm_sub_ps(one, res_erfc), res_erf, mask);

    return res;
}




/* FULL precision erfc(), 256 bit wide. Only errors in LSB */
static __m256
gmx_mm256_erfc_ps(__m256 x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const __m256  CA6      = _mm256_set1_ps(7.853861353153693e-5f);
    const __m256  CA5      = _mm256_set1_ps(-8.010193625184903e-4f);
    const __m256  CA4      = _mm256_set1_ps(5.188327685732524e-3f);
    const __m256  CA3      = _mm256_set1_ps(-2.685381193529856e-2f);
    const __m256  CA2      = _mm256_set1_ps(1.128358514861418e-1f);
    const __m256  CA1      = _mm256_set1_ps(-3.761262582423300e-1f);
    const __m256  CA0      = _mm256_set1_ps(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const __m256  CB9      = _mm256_set1_ps(-0.0018629930017603923f);
    const __m256  CB8      = _mm256_set1_ps(0.003909821287598495f);
    const __m256  CB7      = _mm256_set1_ps(-0.0052094582210355615f);
    const __m256  CB6      = _mm256_set1_ps(0.005685614362160572f);
    const __m256  CB5      = _mm256_set1_ps(-0.0025367682853477272f);
    const __m256  CB4      = _mm256_set1_ps(-0.010199799682318782f);
    const __m256  CB3      = _mm256_set1_ps(0.04369575504816542f);
    const __m256  CB2      = _mm256_set1_ps(-0.11884063474674492f);
    const __m256  CB1      = _mm256_set1_ps(0.2732120154030589f);
    const __m256  CB0      = _mm256_set1_ps(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const __m256  CC10     = _mm256_set1_ps(-0.0445555913112064f);
    const __m256  CC9      = _mm256_set1_ps(0.21376355144663348f);
    const __m256  CC8      = _mm256_set1_ps(-0.3473187200259257f);
    const __m256  CC7      = _mm256_set1_ps(0.016690861551248114f);
    const __m256  CC6      = _mm256_set1_ps(0.7560973182491192f);
    const __m256  CC5      = _mm256_set1_ps(-1.2137903600145787f);
    const __m256  CC4      = _mm256_set1_ps(0.8411872321232948f);
    const __m256  CC3      = _mm256_set1_ps(-0.08670413896296343f);
    const __m256  CC2      = _mm256_set1_ps(-0.27124782687240334f);
    const __m256  CC1      = _mm256_set1_ps(-0.0007502488047806069f);
    const __m256  CC0      = _mm256_set1_ps(0.5642114853803148f);

    /* Coefficients for expansion of exp(x) in [0,0.1] */
    /* CD0 and CD1 are both 1.0, so no need to declare them separately */
    const __m256  CD2      = _mm256_set1_ps(0.5000066608081202f);
    const __m256  CD3      = _mm256_set1_ps(0.1664795422874624f);
    const __m256  CD4      = _mm256_set1_ps(0.04379839977652482f);

    const __m256  sieve    = _mm256_castsi256_ps( _mm256_set1_epi32(0xfffff000) );
    const __m256  signbit  = _mm256_castsi256_ps( _mm256_set1_epi32(0x80000000) );
    const __m256  one      = _mm256_set1_ps(1.0f);
    const __m256  two      = _mm256_set1_ps(2.0f);

    __m256        x2, x4, y;
    __m256        z, q, t, t2, w, w2;
    __m256        pA0, pA1, pB0, pB1, pC0, pC1;
    __m256        expmx2, corr;
    __m256        res_erf, res_erfc, res;
    __m256        mask;

    /* Calculate erf() */
    x2     = _mm256_mul_ps(x, x);
    x4     = _mm256_mul_ps(x2, x2);

    pA0  = _mm256_mul_ps(CA6, x4);
    pA1  = _mm256_mul_ps(CA5, x4);
    pA0  = _mm256_add_ps(pA0, CA4);
    pA1  = _mm256_add_ps(pA1, CA3);
    pA0  = _mm256_mul_ps(pA0, x4);
    pA1  = _mm256_mul_ps(pA1, x4);
    pA0  = _mm256_add_ps(pA0, CA2);
    pA1  = _mm256_add_ps(pA1, CA1);
    pA0  = _mm256_mul_ps(pA0, x4);
    pA1  = _mm256_mul_ps(pA1, x2);
    pA0  = _mm256_add_ps(pA0, pA1);
    pA0  = _mm256_add_ps(pA0, CA0);

    res_erf = _mm256_mul_ps(x, pA0);

    /* Calculate erfc */
    y       = gmx_mm256_abs_ps(x);
    t       = gmx_mm256_inv_ps(y);
    w       = _mm256_sub_ps(t, one);
    t2      = _mm256_mul_ps(t, t);
    w2      = _mm256_mul_ps(w, w);
    /*
     * We cannot simply calculate exp(-x2) directly in single precision, since
     * that will lose a couple of bits of precision due to the multiplication.
     * Instead, we introduce x=z+w, where the last 12 bits of precision are in w.
     * Then we get exp(-x2) = exp(-z2)*exp((z-x)*(z+x)).
     *
     * The only drawback with this is that it requires TWO separate exponential
     * evaluations, which would be horrible performance-wise. However, the argument
     * for the second exp() call is always small, so there we simply use a
     * low-order minimax expansion on [0,0.1].
     */

    z       = _mm256_and_ps(y, sieve);
    q       = _mm256_mul_ps( _mm256_sub_ps(z, y), _mm256_add_ps(z, y) );

    corr    = _mm256_mul_ps(CD4, q);
    corr    = _mm256_add_ps(corr, CD3);
    corr    = _mm256_mul_ps(corr, q);
    corr    = _mm256_add_ps(corr, CD2);
    corr    = _mm256_mul_ps(corr, q);
    corr    = _mm256_add_ps(corr, one);
    corr    = _mm256_mul_ps(corr, q);
    corr    = _mm256_add_ps(corr, one);

    expmx2  = gmx_mm256_exp_ps( _mm256_or_ps( signbit, _mm256_mul_ps(z, z) ) );
    expmx2  = _mm256_mul_ps(expmx2, corr);

    pB1  = _mm256_mul_ps(CB9, w2);
    pB0  = _mm256_mul_ps(CB8, w2);
    pB1  = _mm256_add_ps(pB1, CB7);
    pB0  = _mm256_add_ps(pB0, CB6);
    pB1  = _mm256_mul_ps(pB1, w2);
    pB0  = _mm256_mul_ps(pB0, w2);
    pB1  = _mm256_add_ps(pB1, CB5);
    pB0  = _mm256_add_ps(pB0, CB4);
    pB1  = _mm256_mul_ps(pB1, w2);
    pB0  = _mm256_mul_ps(pB0, w2);
    pB1  = _mm256_add_ps(pB1, CB3);
    pB0  = _mm256_add_ps(pB0, CB2);
    pB1  = _mm256_mul_ps(pB1, w2);
    pB0  = _mm256_mul_ps(pB0, w2);
    pB1  = _mm256_add_ps(pB1, CB1);
    pB1  = _mm256_mul_ps(pB1, w);
    pB0  = _mm256_add_ps(pB0, pB1);
    pB0  = _mm256_add_ps(pB0, CB0);

    pC0  = _mm256_mul_ps(CC10, t2);
    pC1  = _mm256_mul_ps(CC9, t2);
    pC0  = _mm256_add_ps(pC0, CC8);
    pC1  = _mm256_add_ps(pC1, CC7);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t2);
    pC0  = _mm256_add_ps(pC0, CC6);
    pC1  = _mm256_add_ps(pC1, CC5);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t2);
    pC0  = _mm256_add_ps(pC0, CC4);
    pC1  = _mm256_add_ps(pC1, CC3);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t2);
    pC0  = _mm256_add_ps(pC0, CC2);
    pC1  = _mm256_add_ps(pC1, CC1);
    pC0  = _mm256_mul_ps(pC0, t2);
    pC1  = _mm256_mul_ps(pC1, t);
    pC0  = _mm256_add_ps(pC0, pC1);
    pC0  = _mm256_add_ps(pC0, CC0);
    pC0  = _mm256_mul_ps(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = _mm256_cmp_ps(two, y, _CMP_LT_OQ);
    res_erfc = _mm256_blendv_ps(pB0, pC0, mask);
    res_erfc = _mm256_mul_ps(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OQ);
    res_erfc = _mm256_blendv_ps(res_erfc, _mm256_sub_ps(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = _mm256_cmp_ps(y, _mm256_set1_ps(0.75f), _CMP_LT_OQ);
    res  = _mm256_blendv_ps(res_erfc, _mm256_sub_ps(one, res_erf), mask);

    return res;
}


/* erfc(), 128 bit wide */
static __m128
gmx_mm_erfc_ps(__m128 x)
{
    /* Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1] */
    const __m128  CA6      = _mm_set1_ps(7.853861353153693e-5f);
    const __m128  CA5      = _mm_set1_ps(-8.010193625184903e-4f);
    const __m128  CA4      = _mm_set1_ps(5.188327685732524e-3f);
    const __m128  CA3      = _mm_set1_ps(-2.685381193529856e-2f);
    const __m128  CA2      = _mm_set1_ps(1.128358514861418e-1f);
    const __m128  CA1      = _mm_set1_ps(-3.761262582423300e-1f);
    const __m128  CA0      = _mm_set1_ps(1.128379165726710f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2] */
    const __m128  CB9      = _mm_set1_ps(-0.0018629930017603923f);
    const __m128  CB8      = _mm_set1_ps(0.003909821287598495f);
    const __m128  CB7      = _mm_set1_ps(-0.0052094582210355615f);
    const __m128  CB6      = _mm_set1_ps(0.005685614362160572f);
    const __m128  CB5      = _mm_set1_ps(-0.0025367682853477272f);
    const __m128  CB4      = _mm_set1_ps(-0.010199799682318782f);
    const __m128  CB3      = _mm_set1_ps(0.04369575504816542f);
    const __m128  CB2      = _mm_set1_ps(-0.11884063474674492f);
    const __m128  CB1      = _mm_set1_ps(0.2732120154030589f);
    const __m128  CB0      = _mm_set1_ps(0.42758357702025784f);
    /* Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19] */
    const __m128  CC10     = _mm_set1_ps(-0.0445555913112064f);
    const __m128  CC9      = _mm_set1_ps(0.21376355144663348f);
    const __m128  CC8      = _mm_set1_ps(-0.3473187200259257f);
    const __m128  CC7      = _mm_set1_ps(0.016690861551248114f);
    const __m128  CC6      = _mm_set1_ps(0.7560973182491192f);
    const __m128  CC5      = _mm_set1_ps(-1.2137903600145787f);
    const __m128  CC4      = _mm_set1_ps(0.8411872321232948f);
    const __m128  CC3      = _mm_set1_ps(-0.08670413896296343f);
    const __m128  CC2      = _mm_set1_ps(-0.27124782687240334f);
    const __m128  CC1      = _mm_set1_ps(-0.0007502488047806069f);
    const __m128  CC0      = _mm_set1_ps(0.5642114853803148f);

    /* Coefficients for expansion of exp(x) in [0,0.1] */
    /* CD0 and CD1 are both 1.0, so no need to declare them separately */
    const __m128  CD2      = _mm_set1_ps(0.5000066608081202f);
    const __m128  CD3      = _mm_set1_ps(0.1664795422874624f);
    const __m128  CD4      = _mm_set1_ps(0.04379839977652482f);

    const __m128  sieve    = gmx_mm_castsi128_ps( _mm_set1_epi32(0xfffff000) );
    const __m128  signbit  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x80000000) );
    const __m128  one      = _mm_set1_ps(1.0f);
    const __m128  two      = _mm_set1_ps(2.0f);

    __m128        x2, x4, y;
    __m128        z, q, t, t2, w, w2;
    __m128        pA0, pA1, pB0, pB1, pC0, pC1;
    __m128        expmx2, corr;
    __m128        res_erf, res_erfc, res;
    __m128        mask;

    /* Calculate erf() */
    x2     = _mm_mul_ps(x, x);
    x4     = _mm_mul_ps(x2, x2);

    pA0  = _mm_mul_ps(CA6, x4);
    pA1  = _mm_mul_ps(CA5, x4);
    pA0  = _mm_add_ps(pA0, CA4);
    pA1  = _mm_add_ps(pA1, CA3);
    pA0  = _mm_mul_ps(pA0, x4);
    pA1  = _mm_mul_ps(pA1, x4);
    pA0  = _mm_add_ps(pA0, CA2);
    pA1  = _mm_add_ps(pA1, CA1);
    pA0  = _mm_mul_ps(pA0, x4);
    pA1  = _mm_mul_ps(pA1, x2);
    pA0  = _mm_add_ps(pA0, pA1);
    pA0  = _mm_add_ps(pA0, CA0);

    res_erf = _mm_mul_ps(x, pA0);

    /* Calculate erfc */
    y       = gmx_mm_abs_ps(x);
    t       = gmx_mm_inv_ps(y);
    w       = _mm_sub_ps(t, one);
    t2      = _mm_mul_ps(t, t);
    w2      = _mm_mul_ps(w, w);
    /*
     * We cannot simply calculate exp(-x2) directly in single precision, since
     * that will lose a couple of bits of precision due to the multiplication.
     * Instead, we introduce x=z+w, where the last 12 bits of precision are in w.
     * Then we get exp(-x2) = exp(-z2)*exp((z-x)*(z+x)).
     *
     * The only drawback with this is that it requires TWO separate exponential
     * evaluations, which would be horrible performance-wise. However, the argument
     * for the second exp() call is always small, so there we simply use a
     * low-order minimax expansion on [0,0.1].
     */

    z       = _mm_and_ps(y, sieve);
    q       = _mm_mul_ps( _mm_sub_ps(z, y), _mm_add_ps(z, y) );

    corr    = _mm_mul_ps(CD4, q);
    corr    = _mm_add_ps(corr, CD3);
    corr    = _mm_mul_ps(corr, q);
    corr    = _mm_add_ps(corr, CD2);
    corr    = _mm_mul_ps(corr, q);
    corr    = _mm_add_ps(corr, one);
    corr    = _mm_mul_ps(corr, q);
    corr    = _mm_add_ps(corr, one);

    expmx2  = gmx_mm_exp_ps( _mm_or_ps( signbit, _mm_mul_ps(z, z) ) );
    expmx2  = _mm_mul_ps(expmx2, corr);

    pB1  = _mm_mul_ps(CB9, w2);
    pB0  = _mm_mul_ps(CB8, w2);
    pB1  = _mm_add_ps(pB1, CB7);
    pB0  = _mm_add_ps(pB0, CB6);
    pB1  = _mm_mul_ps(pB1, w2);
    pB0  = _mm_mul_ps(pB0, w2);
    pB1  = _mm_add_ps(pB1, CB5);
    pB0  = _mm_add_ps(pB0, CB4);
    pB1  = _mm_mul_ps(pB1, w2);
    pB0  = _mm_mul_ps(pB0, w2);
    pB1  = _mm_add_ps(pB1, CB3);
    pB0  = _mm_add_ps(pB0, CB2);
    pB1  = _mm_mul_ps(pB1, w2);
    pB0  = _mm_mul_ps(pB0, w2);
    pB1  = _mm_add_ps(pB1, CB1);
    pB1  = _mm_mul_ps(pB1, w);
    pB0  = _mm_add_ps(pB0, pB1);
    pB0  = _mm_add_ps(pB0, CB0);

    pC0  = _mm_mul_ps(CC10, t2);
    pC1  = _mm_mul_ps(CC9, t2);
    pC0  = _mm_add_ps(pC0, CC8);
    pC1  = _mm_add_ps(pC1, CC7);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t2);
    pC0  = _mm_add_ps(pC0, CC6);
    pC1  = _mm_add_ps(pC1, CC5);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t2);
    pC0  = _mm_add_ps(pC0, CC4);
    pC1  = _mm_add_ps(pC1, CC3);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t2);
    pC0  = _mm_add_ps(pC0, CC2);
    pC1  = _mm_add_ps(pC1, CC1);
    pC0  = _mm_mul_ps(pC0, t2);
    pC1  = _mm_mul_ps(pC1, t);
    pC0  = _mm_add_ps(pC0, pC1);
    pC0  = _mm_add_ps(pC0, CC0);
    pC0  = _mm_mul_ps(pC0, t);

    /* SELECT pB0 or pC0 for erfc() */
    mask     = _mm_cmplt_ps(two, y);
    res_erfc = _mm_blendv_ps(pB0, pC0, mask);
    res_erfc = _mm_mul_ps(res_erfc, expmx2);

    /* erfc(x<0) = 2-erfc(|x|) */
    mask     = _mm_cmplt_ps(x, _mm_setzero_ps());
    res_erfc = _mm_blendv_ps(res_erfc, _mm_sub_ps(two, res_erfc), mask);

    /* Select erf() or erfc() */
    mask = _mm_cmplt_ps(y, _mm_set1_ps(0.75f));
    res  = _mm_blendv_ps(res_erfc, _mm_sub_ps(one, res_erf), mask);

    return res;
}



/* Calculate the force correction due to PME analytically.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables.
 *
 * The direct-space potential should be Erfc(beta*r)/r, but there
 * are some problems evaluating that:
 *
 * First, the error function is difficult (read: expensive) to
 * approxmiate accurately for intermediate to large arguments, and
 * this happens already in ranges of beta*r that occur in simulations.
 * Second, we now try to avoid calculating potentials in Gromacs but
 * use forces directly.
 *
 * We can simply things slight by noting that the PME part is really
 * a correction to the normal Coulomb force since Erfc(z)=1-Erf(z), i.e.
 *
 * V= 1/r - Erf(beta*r)/r
 *
 * The first term we already have from the inverse square root, so
 * that we can leave out of this routine.
 *
 * For pme tolerances of 1e-3 to 1e-8 and cutoffs of 0.5nm to 1.8nm,
 * the argument beta*r will be in the range 0.15 to ~4. Use your
 * favorite plotting program to realize how well-behaved Erf(z)/z is
 * in this range!
 *
 * We approximate f(z)=erf(z)/z with a rational minimax polynomial.
 * However, it turns out it is more efficient to approximate f(z)/z and
 * then only use even powers. This is another minor optimization, since
 * we actually WANT f(z)/z, because it is going to be multiplied by
 * the vector between the two atoms to get the vectorial force. The
 * fastest flops are the ones we can avoid calculating!
 *
 * So, here's how it should be used:
 *
 * 1. Calculate r^2.
 * 2. Multiply by beta^2, so you get z^2=beta^2*r^2.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *
 *       2*exp(-z^2)     erf(z)
 *       ------------ - --------
 *       sqrt(Pi)*z^2      z^3
 *
 * 5. Multiply the entire expression by beta^3. This will get you
 *
 *       beta^3*2*exp(-z^2)     beta^3*erf(z)
 *       ------------------  - ---------------
 *          sqrt(Pi)*z^2            z^3
 *
 *    or, switching back to r (z=r*beta):
 *
 *       2*beta*exp(-r^2*beta^2)   erf(r*beta)
 *       ----------------------- - -----------
 *            sqrt(Pi)*r^2            r^3
 *
 *
 *    With a bit of math exercise you should be able to confirm that
 *    this is exactly D[Erf[beta*r]/r,r] divided by r another time.
 *
 * 6. Add the result to 1/r^3, multiply by the product of the charges,
 *    and you have your force (divided by r). A final multiplication
 *    with the vector connecting the two particles and you have your
 *    vectorial force to add to the particles.
 *
 */
static __m256
gmx_mm256_pmecorrF_ps(__m256 z2)
{
    const __m256  FN6      = _mm256_set1_ps(-1.7357322914161492954e-8f);
    const __m256  FN5      = _mm256_set1_ps(1.4703624142580877519e-6f);
    const __m256  FN4      = _mm256_set1_ps(-0.000053401640219807709149f);
    const __m256  FN3      = _mm256_set1_ps(0.0010054721316683106153f);
    const __m256  FN2      = _mm256_set1_ps(-0.019278317264888380590f);
    const __m256  FN1      = _mm256_set1_ps(0.069670166153766424023f);
    const __m256  FN0      = _mm256_set1_ps(-0.75225204789749321333f);

    const __m256  FD4      = _mm256_set1_ps(0.0011193462567257629232f);
    const __m256  FD3      = _mm256_set1_ps(0.014866955030185295499f);
    const __m256  FD2      = _mm256_set1_ps(0.11583842382862377919f);
    const __m256  FD1      = _mm256_set1_ps(0.50736591960530292870f);
    const __m256  FD0      = _mm256_set1_ps(1.0f);

    __m256        z4;
    __m256        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = _mm256_mul_ps(z2, z2);

    polyFD0        = _mm256_mul_ps(FD4, z4);
    polyFD1        = _mm256_mul_ps(FD3, z4);
    polyFD0        = _mm256_add_ps(polyFD0, FD2);
    polyFD1        = _mm256_add_ps(polyFD1, FD1);
    polyFD0        = _mm256_mul_ps(polyFD0, z4);
    polyFD1        = _mm256_mul_ps(polyFD1, z2);
    polyFD0        = _mm256_add_ps(polyFD0, FD0);
    polyFD0        = _mm256_add_ps(polyFD0, polyFD1);

    polyFD0        = gmx_mm256_inv_ps(polyFD0);

    polyFN0        = _mm256_mul_ps(FN6, z4);
    polyFN1        = _mm256_mul_ps(FN5, z4);
    polyFN0        = _mm256_add_ps(polyFN0, FN4);
    polyFN1        = _mm256_add_ps(polyFN1, FN3);
    polyFN0        = _mm256_mul_ps(polyFN0, z4);
    polyFN1        = _mm256_mul_ps(polyFN1, z4);
    polyFN0        = _mm256_add_ps(polyFN0, FN2);
    polyFN1        = _mm256_add_ps(polyFN1, FN1);
    polyFN0        = _mm256_mul_ps(polyFN0, z4);
    polyFN1        = _mm256_mul_ps(polyFN1, z2);
    polyFN0        = _mm256_add_ps(polyFN0, FN0);
    polyFN0        = _mm256_add_ps(polyFN0, polyFN1);

    return _mm256_mul_ps(polyFN0, polyFD0);
}


static __m128
gmx_mm_pmecorrF_ps(__m128 z2)
{
    const __m128  FN6      = _mm_set1_ps(-1.7357322914161492954e-8f);
    const __m128  FN5      = _mm_set1_ps(1.4703624142580877519e-6f);
    const __m128  FN4      = _mm_set1_ps(-0.000053401640219807709149f);
    const __m128  FN3      = _mm_set1_ps(0.0010054721316683106153f);
    const __m128  FN2      = _mm_set1_ps(-0.019278317264888380590f);
    const __m128  FN1      = _mm_set1_ps(0.069670166153766424023f);
    const __m128  FN0      = _mm_set1_ps(-0.75225204789749321333f);

    const __m128  FD4      = _mm_set1_ps(0.0011193462567257629232f);
    const __m128  FD3      = _mm_set1_ps(0.014866955030185295499f);
    const __m128  FD2      = _mm_set1_ps(0.11583842382862377919f);
    const __m128  FD1      = _mm_set1_ps(0.50736591960530292870f);
    const __m128  FD0      = _mm_set1_ps(1.0f);

    __m128        z4;
    __m128        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = _mm_mul_ps(z2, z2);

    polyFD0        = _mm_mul_ps(FD4, z4);
    polyFD1        = _mm_mul_ps(FD3, z4);
    polyFD0        = _mm_add_ps(polyFD0, FD2);
    polyFD1        = _mm_add_ps(polyFD1, FD1);
    polyFD0        = _mm_mul_ps(polyFD0, z4);
    polyFD1        = _mm_mul_ps(polyFD1, z2);
    polyFD0        = _mm_add_ps(polyFD0, FD0);
    polyFD0        = _mm_add_ps(polyFD0, polyFD1);

    polyFD0        = gmx_mm_inv_ps(polyFD0);

    polyFN0        = _mm_mul_ps(FN6, z4);
    polyFN1        = _mm_mul_ps(FN5, z4);
    polyFN0        = _mm_add_ps(polyFN0, FN4);
    polyFN1        = _mm_add_ps(polyFN1, FN3);
    polyFN0        = _mm_mul_ps(polyFN0, z4);
    polyFN1        = _mm_mul_ps(polyFN1, z4);
    polyFN0        = _mm_add_ps(polyFN0, FN2);
    polyFN1        = _mm_add_ps(polyFN1, FN1);
    polyFN0        = _mm_mul_ps(polyFN0, z4);
    polyFN1        = _mm_mul_ps(polyFN1, z2);
    polyFN0        = _mm_add_ps(polyFN0, FN0);
    polyFN0        = _mm_add_ps(polyFN0, polyFN1);

    return _mm_mul_ps(polyFN0, polyFD0);
}



/* Calculate the potential correction due to PME analytically.
 *
 * See gmx_mm256_pmecorrF_ps() for details about the approximation.
 *
 * This routine calculates Erf(z)/z, although you should provide z^2
 * as the input argument.
 *
 * Here's how it should be used:
 *
 * 1. Calculate r^2.
 * 2. Multiply by beta^2, so you get z^2=beta^2*r^2.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *
 *        erf(z)
 *       --------
 *          z
 *
 * 5. Multiply the entire expression by beta and switching back to r (z=r*beta):
 *
 *       erf(r*beta)
 *       -----------
 *           r
 *
 * 6. Subtract the result from 1/r, multiply by the product of the charges,
 *    and you have your potential.
 */
static __m256
gmx_mm256_pmecorrV_ps(__m256 z2)
{
    const __m256  VN6      = _mm256_set1_ps(1.9296833005951166339e-8f);
    const __m256  VN5      = _mm256_set1_ps(-1.4213390571557850962e-6f);
    const __m256  VN4      = _mm256_set1_ps(0.000041603292906656984871f);
    const __m256  VN3      = _mm256_set1_ps(-0.00013134036773265025626f);
    const __m256  VN2      = _mm256_set1_ps(0.038657983986041781264f);
    const __m256  VN1      = _mm256_set1_ps(0.11285044772717598220f);
    const __m256  VN0      = _mm256_set1_ps(1.1283802385263030286f);

    const __m256  VD3      = _mm256_set1_ps(0.0066752224023576045451f);
    const __m256  VD2      = _mm256_set1_ps(0.078647795836373922256f);
    const __m256  VD1      = _mm256_set1_ps(0.43336185284710920150f);
    const __m256  VD0      = _mm256_set1_ps(1.0f);

    __m256        z4;
    __m256        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = _mm256_mul_ps(z2, z2);

    polyVD1        = _mm256_mul_ps(VD3, z4);
    polyVD0        = _mm256_mul_ps(VD2, z4);
    polyVD1        = _mm256_add_ps(polyVD1, VD1);
    polyVD0        = _mm256_add_ps(polyVD0, VD0);
    polyVD1        = _mm256_mul_ps(polyVD1, z2);
    polyVD0        = _mm256_add_ps(polyVD0, polyVD1);

    polyVD0        = gmx_mm256_inv_ps(polyVD0);

    polyVN0        = _mm256_mul_ps(VN6, z4);
    polyVN1        = _mm256_mul_ps(VN5, z4);
    polyVN0        = _mm256_add_ps(polyVN0, VN4);
    polyVN1        = _mm256_add_ps(polyVN1, VN3);
    polyVN0        = _mm256_mul_ps(polyVN0, z4);
    polyVN1        = _mm256_mul_ps(polyVN1, z4);
    polyVN0        = _mm256_add_ps(polyVN0, VN2);
    polyVN1        = _mm256_add_ps(polyVN1, VN1);
    polyVN0        = _mm256_mul_ps(polyVN0, z4);
    polyVN1        = _mm256_mul_ps(polyVN1, z2);
    polyVN0        = _mm256_add_ps(polyVN0, VN0);
    polyVN0        = _mm256_add_ps(polyVN0, polyVN1);

    return _mm256_mul_ps(polyVN0, polyVD0);
}


static __m128
gmx_mm_pmecorrV_ps(__m128 z2)
{
    const __m128  VN6      = _mm_set1_ps(1.9296833005951166339e-8f);
    const __m128  VN5      = _mm_set1_ps(-1.4213390571557850962e-6f);
    const __m128  VN4      = _mm_set1_ps(0.000041603292906656984871f);
    const __m128  VN3      = _mm_set1_ps(-0.00013134036773265025626f);
    const __m128  VN2      = _mm_set1_ps(0.038657983986041781264f);
    const __m128  VN1      = _mm_set1_ps(0.11285044772717598220f);
    const __m128  VN0      = _mm_set1_ps(1.1283802385263030286f);

    const __m128  VD3      = _mm_set1_ps(0.0066752224023576045451f);
    const __m128  VD2      = _mm_set1_ps(0.078647795836373922256f);
    const __m128  VD1      = _mm_set1_ps(0.43336185284710920150f);
    const __m128  VD0      = _mm_set1_ps(1.0f);

    __m128        z4;
    __m128        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = _mm_mul_ps(z2, z2);

    polyVD1        = _mm_mul_ps(VD3, z4);
    polyVD0        = _mm_mul_ps(VD2, z4);
    polyVD1        = _mm_add_ps(polyVD1, VD1);
    polyVD0        = _mm_add_ps(polyVD0, VD0);
    polyVD1        = _mm_mul_ps(polyVD1, z2);
    polyVD0        = _mm_add_ps(polyVD0, polyVD1);

    polyVD0        = gmx_mm_inv_ps(polyVD0);

    polyVN0        = _mm_mul_ps(VN6, z4);
    polyVN1        = _mm_mul_ps(VN5, z4);
    polyVN0        = _mm_add_ps(polyVN0, VN4);
    polyVN1        = _mm_add_ps(polyVN1, VN3);
    polyVN0        = _mm_mul_ps(polyVN0, z4);
    polyVN1        = _mm_mul_ps(polyVN1, z4);
    polyVN0        = _mm_add_ps(polyVN0, VN2);
    polyVN1        = _mm_add_ps(polyVN1, VN1);
    polyVN0        = _mm_mul_ps(polyVN0, z4);
    polyVN1        = _mm_mul_ps(polyVN1, z2);
    polyVN0        = _mm_add_ps(polyVN0, VN0);
    polyVN0        = _mm_add_ps(polyVN0, polyVN1);

    return _mm_mul_ps(polyVN0, polyVD0);
}


static int
gmx_mm256_sincos_ps(__m256  x,
                    __m256 *sinval,
                    __m256 *cosval)
{
    const __m256  two_over_pi = _mm256_set1_ps(2.0f/(float)M_PI);
    const __m256  half        = _mm256_set1_ps(0.5f);
    const __m256  one         = _mm256_set1_ps(1.0f);
    const __m256  zero        = _mm256_setzero_ps();

    const __m128i ione       = _mm_set1_epi32(1);

    const __m256  mask_one    = _mm256_castsi256_ps(_mm256_set1_epi32(1));
    const __m256  mask_two    = _mm256_castsi256_ps(_mm256_set1_epi32(2));
    const __m256  mask_three  = _mm256_castsi256_ps(_mm256_set1_epi32(3));

    const __m256  CA1         = _mm256_set1_ps(1.5703125f);
    const __m256  CA2         = _mm256_set1_ps(4.837512969970703125e-4f);
    const __m256  CA3         = _mm256_set1_ps(7.54978995489188216e-8f);

    const __m256  CC0         = _mm256_set1_ps(-0.0013602249f);
    const __m256  CC1         = _mm256_set1_ps(0.0416566950f);
    const __m256  CC2         = _mm256_set1_ps(-0.4999990225f);
    const __m256  CS0         = _mm256_set1_ps(-0.0001950727f);
    const __m256  CS1         = _mm256_set1_ps(0.0083320758f);
    const __m256  CS2         = _mm256_set1_ps(-0.1666665247f);

    const __m256  signbit    = _mm256_castsi256_ps( _mm256_set1_epi32(0x80000000) );

    __m256        y, y2;
    __m256        z;
    __m256i       iz;
    __m128i       iz_high, iz_low;
    __m256        offset_sin, offset_cos;
    __m256        mask_sin, mask_cos;
    __m256        tmp1, tmp2;
    __m256        tmp_sin, tmp_cos;

    y               = _mm256_mul_ps(x, two_over_pi);
    y               = _mm256_add_ps(y, _mm256_or_ps(_mm256_and_ps(y, signbit), half));

    iz              = _mm256_cvttps_epi32(y);
    z               = _mm256_round_ps(y, _MM_FROUND_TO_ZERO);

    offset_sin      = _mm256_and_ps(_mm256_castsi256_ps(iz), mask_three);

    iz_high         = _mm256_extractf128_si256(iz, 0x1);
    iz_low          = _mm256_castsi256_si128(iz);
    iz_low          = _mm_add_epi32(iz_low, ione);
    iz_high         = _mm_add_epi32(iz_high, ione);
    iz              = _mm256_castsi128_si256(iz_low);
    iz              = _mm256_insertf128_si256(iz, iz_high, 0x1);
    offset_cos      = _mm256_castsi256_ps(iz);

    /* Extended precision arithmethic to achieve full precision */
    y               = _mm256_mul_ps(z, CA1);
    tmp1            = _mm256_mul_ps(z, CA2);
    tmp2            = _mm256_mul_ps(z, CA3);
    y               = _mm256_sub_ps(x, y);
    y               = _mm256_sub_ps(y, tmp1);
    y               = _mm256_sub_ps(y, tmp2);

    y2              = _mm256_mul_ps(y, y);

    tmp1            = _mm256_mul_ps(CC0, y2);
    tmp1            = _mm256_add_ps(tmp1, CC1);
    tmp2            = _mm256_mul_ps(CS0, y2);
    tmp2            = _mm256_add_ps(tmp2, CS1);
    tmp1            = _mm256_mul_ps(tmp1, y2);
    tmp1            = _mm256_add_ps(tmp1, CC2);
    tmp2            = _mm256_mul_ps(tmp2, y2);
    tmp2            = _mm256_add_ps(tmp2, CS2);

    tmp1            = _mm256_mul_ps(tmp1, y2);
    tmp1            = _mm256_add_ps(tmp1, one);

    tmp2            = _mm256_mul_ps(tmp2, _mm256_mul_ps(y, y2));
    tmp2            = _mm256_add_ps(tmp2, y);

#ifdef __INTEL_COMPILER
    /* Intel Compiler version 12.1.3 20120130 is buggy if optimization is enabled unless we cast explicitly! */
    mask_sin        = _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(offset_sin, mask_one))), zero, _CMP_EQ_OQ);
    mask_cos        = _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(offset_cos, mask_one))), zero, _CMP_EQ_OQ);
#else
    mask_sin        = _mm256_cmp_ps( _mm256_and_ps(offset_sin, mask_one), zero, _CMP_EQ_OQ);
    mask_cos        = _mm256_cmp_ps( _mm256_and_ps(offset_cos, mask_one), zero, _CMP_EQ_OQ);
#endif
    tmp_sin         = _mm256_blendv_ps(tmp1, tmp2, mask_sin);
    tmp_cos         = _mm256_blendv_ps(tmp1, tmp2, mask_cos);

    tmp1            = _mm256_xor_ps(signbit, tmp_sin);
    tmp2            = _mm256_xor_ps(signbit, tmp_cos);

#ifdef __INTEL_COMPILER
    /* Intel Compiler version 12.1.3 20120130 is buggy if optimization is enabled unless we cast explicitly! */
    mask_sin        = _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(offset_sin, mask_two))), zero, _CMP_EQ_OQ);
    mask_cos        = _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(offset_cos, mask_two))), zero, _CMP_EQ_OQ);
#else
    mask_sin        = _mm256_cmp_ps( _mm256_and_ps(offset_sin, mask_two), zero, _CMP_EQ_OQ);
    mask_cos        = _mm256_cmp_ps( _mm256_and_ps(offset_cos, mask_two), zero, _CMP_EQ_OQ);

#endif
    *sinval         = _mm256_blendv_ps(tmp1, tmp_sin, mask_sin);
    *cosval         = _mm256_blendv_ps(tmp2, tmp_cos, mask_cos);

    return 0;
}

static int
gmx_mm_sincos_ps(__m128  x,
                 __m128 *sinval,
                 __m128 *cosval)
{
    const __m128  two_over_pi = _mm_set1_ps(2.0/M_PI);
    const __m128  half        = _mm_set1_ps(0.5);
    const __m128  one         = _mm_set1_ps(1.0);

    const __m128i izero      = _mm_set1_epi32(0);
    const __m128i ione       = _mm_set1_epi32(1);
    const __m128i itwo       = _mm_set1_epi32(2);
    const __m128i ithree     = _mm_set1_epi32(3);
    const __m128  signbit    = gmx_mm_castsi128_ps( _mm_set1_epi32(0x80000000) );

    const __m128  CA1         = _mm_set1_ps(1.5703125f);
    const __m128  CA2         = _mm_set1_ps(4.837512969970703125e-4f);
    const __m128  CA3         = _mm_set1_ps(7.54978995489188216e-8f);

    const __m128  CC0         = _mm_set1_ps(-0.0013602249f);
    const __m128  CC1         = _mm_set1_ps(0.0416566950f);
    const __m128  CC2         = _mm_set1_ps(-0.4999990225f);
    const __m128  CS0         = _mm_set1_ps(-0.0001950727f);
    const __m128  CS1         = _mm_set1_ps(0.0083320758f);
    const __m128  CS2         = _mm_set1_ps(-0.1666665247f);

    __m128        y, y2;
    __m128        z;
    __m128i       iz;
    __m128i       offset_sin, offset_cos;
    __m128        tmp1, tmp2;
    __m128        mask_sin, mask_cos;
    __m128        tmp_sin, tmp_cos;

    y          = _mm_mul_ps(x, two_over_pi);
    y          = _mm_add_ps(y, _mm_or_ps(_mm_and_ps(y, signbit), half));

    iz         = _mm_cvttps_epi32(y);
    z          = _mm_round_ps(y, _MM_FROUND_TO_ZERO);

    offset_sin = _mm_and_si128(iz, ithree);
    offset_cos = _mm_add_epi32(iz, ione);

    /* Extended precision arithmethic to achieve full precision */
    y               = _mm_mul_ps(z, CA1);
    tmp1            = _mm_mul_ps(z, CA2);
    tmp2            = _mm_mul_ps(z, CA3);
    y               = _mm_sub_ps(x, y);
    y               = _mm_sub_ps(y, tmp1);
    y               = _mm_sub_ps(y, tmp2);

    y2              = _mm_mul_ps(y, y);

    tmp1            = _mm_mul_ps(CC0, y2);
    tmp1            = _mm_add_ps(tmp1, CC1);
    tmp2            = _mm_mul_ps(CS0, y2);
    tmp2            = _mm_add_ps(tmp2, CS1);
    tmp1            = _mm_mul_ps(tmp1, y2);
    tmp1            = _mm_add_ps(tmp1, CC2);
    tmp2            = _mm_mul_ps(tmp2, y2);
    tmp2            = _mm_add_ps(tmp2, CS2);

    tmp1            = _mm_mul_ps(tmp1, y2);
    tmp1            = _mm_add_ps(tmp1, one);

    tmp2            = _mm_mul_ps(tmp2, _mm_mul_ps(y, y2));
    tmp2            = _mm_add_ps(tmp2, y);

    mask_sin        = gmx_mm_castsi128_ps(_mm_cmpeq_epi32( _mm_and_si128(offset_sin, ione), izero));
    mask_cos        = gmx_mm_castsi128_ps(_mm_cmpeq_epi32( _mm_and_si128(offset_cos, ione), izero));

    tmp_sin         = _mm_blendv_ps(tmp1, tmp2, mask_sin);
    tmp_cos         = _mm_blendv_ps(tmp1, tmp2, mask_cos);

    mask_sin        = gmx_mm_castsi128_ps(_mm_cmpeq_epi32( _mm_and_si128(offset_sin, itwo), izero));
    mask_cos        = gmx_mm_castsi128_ps(_mm_cmpeq_epi32( _mm_and_si128(offset_cos, itwo), izero));

    tmp1            = _mm_xor_ps(signbit, tmp_sin);
    tmp2            = _mm_xor_ps(signbit, tmp_cos);

    *sinval         = _mm_blendv_ps(tmp1, tmp_sin, mask_sin);
    *cosval         = _mm_blendv_ps(tmp2, tmp_cos, mask_cos);

    return 0;
}




/*
 * IMPORTANT: Do NOT call both sin & cos if you need both results, since each of them
 * will then call the sincos() routine and waste a factor 2 in performance!
 */
static __m256
gmx_mm256_sin_ps(__m256 x)
{
    __m256 s, c;
    gmx_mm256_sincos_ps(x, &s, &c);
    return s;
}

static __m128
gmx_mm_sin_ps(__m128 x)
{
    __m128 s, c;
    gmx_mm_sincos_ps(x, &s, &c);
    return s;
}


/*
 * IMPORTANT: Do NOT call both sin & cos if you need both results, since each of them
 * will then call the sincos() routine and waste a factor 2 in performance!
 */
static __m256
gmx_mm256_cos_ps(__m256 x)
{
    __m256 s, c;
    gmx_mm256_sincos_ps(x, &s, &c);
    return c;
}

static __m128
gmx_mm_cos_ps(__m128 x)
{
    __m128 s, c;
    gmx_mm_sincos_ps(x, &s, &c);
    return c;
}


static __m256
gmx_mm256_tan_ps(__m256 x)
{
    __m256 sinval, cosval;
    __m256 tanval;

    gmx_mm256_sincos_ps(x, &sinval, &cosval);

    tanval = _mm256_mul_ps(sinval, gmx_mm256_inv_ps(cosval));

    return tanval;
}

static __m128
gmx_mm_tan_ps(__m128 x)
{
    __m128 sinval, cosval;
    __m128 tanval;

    gmx_mm_sincos_ps(x, &sinval, &cosval);

    tanval = _mm_mul_ps(sinval, gmx_mm_inv_ps(cosval));

    return tanval;
}


static __m256
gmx_mm256_asin_ps(__m256 x)
{
    const __m256 signmask  = _mm256_castsi256_ps( _mm256_set1_epi32(0x7FFFFFFF) );
    const __m256 limitlow  = _mm256_set1_ps(1e-4f);
    const __m256 half      = _mm256_set1_ps(0.5f);
    const __m256 one       = _mm256_set1_ps(1.0f);
    const __m256 halfpi    = _mm256_set1_ps((float)M_PI/2.0f);

    const __m256 CC5        = _mm256_set1_ps(4.2163199048E-2f);
    const __m256 CC4        = _mm256_set1_ps(2.4181311049E-2f);
    const __m256 CC3        = _mm256_set1_ps(4.5470025998E-2f);
    const __m256 CC2        = _mm256_set1_ps(7.4953002686E-2f);
    const __m256 CC1        = _mm256_set1_ps(1.6666752422E-1f);

    __m256       sign;
    __m256       mask;
    __m256       xabs;
    __m256       z, z1, z2, q, q1, q2;
    __m256       pA, pB;

    sign  = _mm256_andnot_ps(signmask, x);
    xabs  = _mm256_and_ps(x, signmask);

    mask  = _mm256_cmp_ps(xabs, half, _CMP_GT_OQ);

    z1    = _mm256_mul_ps(half, _mm256_sub_ps(one, xabs));
    q1    = _mm256_mul_ps(z1, gmx_mm256_invsqrt_ps(z1));
    q1    = _mm256_andnot_ps(_mm256_cmp_ps(xabs, one, _CMP_EQ_OQ), q1);

    q2    = xabs;
    z2    = _mm256_mul_ps(q2, q2);

    z     = _mm256_blendv_ps(z2, z1, mask);
    q     = _mm256_blendv_ps(q2, q1, mask);

    z2    = _mm256_mul_ps(z, z);

    pA    = _mm256_mul_ps(CC5, z2);
    pB    = _mm256_mul_ps(CC4, z2);

    pA    = _mm256_add_ps(pA, CC3);
    pB    = _mm256_add_ps(pB, CC2);

    pA    = _mm256_mul_ps(pA, z2);
    pB    = _mm256_mul_ps(pB, z2);

    pA    = _mm256_add_ps(pA, CC1);
    pA    = _mm256_mul_ps(pA, z);

    z     = _mm256_add_ps(pA, pB);
    z     = _mm256_mul_ps(z, q);
    z     = _mm256_add_ps(z, q);

    q2    = _mm256_sub_ps(halfpi, z);
    q2    = _mm256_sub_ps(q2, z);

    z     = _mm256_blendv_ps(z, q2, mask);

    mask  = _mm256_cmp_ps(xabs, limitlow, _CMP_GT_OQ);
    z     = _mm256_blendv_ps(xabs, z, mask);

    z     = _mm256_xor_ps(z, sign);

    return z;
}

static __m128
gmx_mm_asin_ps(__m128 x)
{
    /* Same algorithm as cephes library */
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );
    const __m128 limitlow  = _mm_set1_ps(1e-4f);
    const __m128 half      = _mm_set1_ps(0.5f);
    const __m128 one       = _mm_set1_ps(1.0f);
    const __m128 halfpi    = _mm_set1_ps(M_PI/2.0f);

    const __m128 CC5        = _mm_set1_ps(4.2163199048E-2f);
    const __m128 CC4        = _mm_set1_ps(2.4181311049E-2f);
    const __m128 CC3        = _mm_set1_ps(4.5470025998E-2f);
    const __m128 CC2        = _mm_set1_ps(7.4953002686E-2f);
    const __m128 CC1        = _mm_set1_ps(1.6666752422E-1f);

    __m128       sign;
    __m128       mask;
    __m128       xabs;
    __m128       z, z1, z2, q, q1, q2;
    __m128       pA, pB;

    sign  = _mm_andnot_ps(signmask, x);
    xabs  = _mm_and_ps(x, signmask);

    mask  = _mm_cmpgt_ps(xabs, half);

    z1    = _mm_mul_ps(half, _mm_sub_ps(one, xabs));
    q1    = _mm_mul_ps(z1, gmx_mm_invsqrt_ps(z1));
    q1    = _mm_andnot_ps(_mm_cmpeq_ps(xabs, one), q1);

    q2    = xabs;
    z2    = _mm_mul_ps(q2, q2);

    z     = _mm_or_ps( _mm_and_ps(mask, z1), _mm_andnot_ps(mask, z2) );
    q     = _mm_or_ps( _mm_and_ps(mask, q1), _mm_andnot_ps(mask, q2) );

    z2    = _mm_mul_ps(z, z);

    pA    = _mm_mul_ps(CC5, z2);
    pB    = _mm_mul_ps(CC4, z2);

    pA    = _mm_add_ps(pA, CC3);
    pB    = _mm_add_ps(pB, CC2);

    pA    = _mm_mul_ps(pA, z2);
    pB    = _mm_mul_ps(pB, z2);

    pA    = _mm_add_ps(pA, CC1);
    pA    = _mm_mul_ps(pA, z);

    z     = _mm_add_ps(pA, pB);
    z     = _mm_mul_ps(z, q);
    z     = _mm_add_ps(z, q);

    q2    = _mm_sub_ps(halfpi, z);
    q2    = _mm_sub_ps(q2, z);

    z     = _mm_or_ps( _mm_and_ps(mask, q2), _mm_andnot_ps(mask, z) );

    mask  = _mm_cmpgt_ps(xabs, limitlow);
    z     = _mm_or_ps( _mm_and_ps(mask, z), _mm_andnot_ps(mask, xabs) );

    z = _mm_xor_ps(z, sign);

    return z;
}


static __m256
gmx_mm256_acos_ps(__m256 x)
{
    const __m256 signmask  = _mm256_castsi256_ps( _mm256_set1_epi32(0x7FFFFFFF) );
    const __m256 one_ps    = _mm256_set1_ps(1.0f);
    const __m256 half_ps   = _mm256_set1_ps(0.5f);
    const __m256 pi_ps     = _mm256_set1_ps((float)M_PI);
    const __m256 halfpi_ps = _mm256_set1_ps((float)M_PI/2.0f);

    __m256       mask1;
    __m256       mask2;
    __m256       xabs;
    __m256       z, z1, z2, z3;

    xabs  = _mm256_and_ps(x, signmask);
    mask1 = _mm256_cmp_ps(xabs, half_ps, _CMP_GT_OQ);
    mask2 = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_GT_OQ);

    z     = _mm256_mul_ps(half_ps, _mm256_sub_ps(one_ps, xabs));
    z     = _mm256_mul_ps(z, gmx_mm256_invsqrt_ps(z));
    z     = _mm256_andnot_ps(_mm256_cmp_ps(xabs, one_ps, _CMP_EQ_OQ), z);

    z     = _mm256_blendv_ps(x, z, mask1);
    z     = gmx_mm256_asin_ps(z);

    z2    = _mm256_add_ps(z, z);
    z1    = _mm256_sub_ps(pi_ps, z2);
    z3    = _mm256_sub_ps(halfpi_ps, z);

    z     = _mm256_blendv_ps(z1, z2, mask2);
    z     = _mm256_blendv_ps(z3, z, mask1);

    return z;
}

static __m128
gmx_mm_acos_ps(__m128 x)
{
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );
    const __m128 one_ps    = _mm_set1_ps(1.0f);
    const __m128 half_ps   = _mm_set1_ps(0.5f);
    const __m128 pi_ps     = _mm_set1_ps(M_PI);
    const __m128 halfpi_ps = _mm_set1_ps(M_PI/2.0f);

    __m128       mask1;
    __m128       mask2;
    __m128       xabs;
    __m128       z, z1, z2, z3;

    xabs  = _mm_and_ps(x, signmask);
    mask1 = _mm_cmpgt_ps(xabs, half_ps);
    mask2 = _mm_cmpgt_ps(x, _mm_setzero_ps());

    z     = _mm_mul_ps(half_ps, _mm_sub_ps(one_ps, xabs));
    z     = _mm_mul_ps(z, gmx_mm_invsqrt_ps(z));
    z     = _mm_andnot_ps(_mm_cmpeq_ps(xabs, one_ps), z);

    z     = _mm_blendv_ps(x, z, mask1);
    z     = gmx_mm_asin_ps(z);

    z2    = _mm_add_ps(z, z);
    z1    = _mm_sub_ps(pi_ps, z2);
    z3    = _mm_sub_ps(halfpi_ps, z);

    z     = _mm_blendv_ps(z1, z2, mask2);
    z     = _mm_blendv_ps(z3, z, mask1);

    return z;
}


static __m256
gmx_mm256_atan_ps(__m256 x)
{
    const __m256 signmask  = _mm256_castsi256_ps( _mm256_set1_epi32(0x7FFFFFFF) );
    const __m256 limit1    = _mm256_set1_ps(0.414213562373095f);
    const __m256 limit2    = _mm256_set1_ps(2.414213562373095f);
    const __m256 quarterpi = _mm256_set1_ps(0.785398163397448f);
    const __m256 halfpi    = _mm256_set1_ps(1.570796326794896f);
    const __m256 mone      = _mm256_set1_ps(-1.0f);
    const __m256 CC3       = _mm256_set1_ps(-3.33329491539E-1f);
    const __m256 CC5       = _mm256_set1_ps(1.99777106478E-1f);
    const __m256 CC7       = _mm256_set1_ps(-1.38776856032E-1);
    const __m256 CC9       = _mm256_set1_ps(8.05374449538e-2f);

    __m256       sign;
    __m256       mask1, mask2;
    __m256       y, z1, z2;
    __m256       x2, x4;
    __m256       sum1, sum2;

    sign  = _mm256_andnot_ps(signmask, x);
    x     = _mm256_and_ps(x, signmask);

    mask1 = _mm256_cmp_ps(x, limit1, _CMP_GT_OQ);
    mask2 = _mm256_cmp_ps(x, limit2, _CMP_GT_OQ);

    z1    = _mm256_mul_ps(_mm256_add_ps(x, mone), gmx_mm256_inv_ps(_mm256_sub_ps(x, mone)));
    z2    = _mm256_mul_ps(mone, gmx_mm256_inv_ps(x));

    y     = _mm256_and_ps(mask1, quarterpi);
    y     = _mm256_blendv_ps(y, halfpi, mask2);

    x     = _mm256_blendv_ps(x, z1, mask1);
    x     = _mm256_blendv_ps(x, z2, mask2);

    x2    = _mm256_mul_ps(x, x);
    x4    = _mm256_mul_ps(x2, x2);

    sum1  = _mm256_mul_ps(CC9, x4);
    sum2  = _mm256_mul_ps(CC7, x4);
    sum1  = _mm256_add_ps(sum1, CC5);
    sum2  = _mm256_add_ps(sum2, CC3);
    sum1  = _mm256_mul_ps(sum1, x4);
    sum2  = _mm256_mul_ps(sum2, x2);

    sum1  = _mm256_add_ps(sum1, sum2);
    sum1  = _mm256_sub_ps(sum1, mone);
    sum1  = _mm256_mul_ps(sum1, x);
    y     = _mm256_add_ps(y, sum1);

    y     = _mm256_xor_ps(y, sign);

    return y;
}

static __m128
gmx_mm_atan_ps(__m128 x)
{
    /* Same algorithm as cephes library */
    const __m128 signmask  = gmx_mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );
    const __m128 limit1    = _mm_set1_ps(0.414213562373095f);
    const __m128 limit2    = _mm_set1_ps(2.414213562373095f);
    const __m128 quarterpi = _mm_set1_ps(0.785398163397448f);
    const __m128 halfpi    = _mm_set1_ps(1.570796326794896f);
    const __m128 mone      = _mm_set1_ps(-1.0f);
    const __m128 CC3       = _mm_set1_ps(-3.33329491539E-1f);
    const __m128 CC5       = _mm_set1_ps(1.99777106478E-1f);
    const __m128 CC7       = _mm_set1_ps(-1.38776856032E-1);
    const __m128 CC9       = _mm_set1_ps(8.05374449538e-2f);

    __m128       sign;
    __m128       mask1, mask2;
    __m128       y, z1, z2;
    __m128       x2, x4;
    __m128       sum1, sum2;

    sign  = _mm_andnot_ps(signmask, x);
    x     = _mm_and_ps(x, signmask);

    mask1 = _mm_cmpgt_ps(x, limit1);
    mask2 = _mm_cmpgt_ps(x, limit2);

    z1    = _mm_mul_ps(_mm_add_ps(x, mone), gmx_mm_inv_ps(_mm_sub_ps(x, mone)));
    z2    = _mm_mul_ps(mone, gmx_mm_inv_ps(x));

    y     = _mm_and_ps(mask1, quarterpi);
    y     = _mm_blendv_ps(y, halfpi, mask2);

    x     = _mm_blendv_ps(x, z1, mask1);
    x     = _mm_blendv_ps(x, z2, mask2);

    x2    = _mm_mul_ps(x, x);
    x4    = _mm_mul_ps(x2, x2);

    sum1  = _mm_mul_ps(CC9, x4);
    sum2  = _mm_mul_ps(CC7, x4);
    sum1  = _mm_add_ps(sum1, CC5);
    sum2  = _mm_add_ps(sum2, CC3);
    sum1  = _mm_mul_ps(sum1, x4);
    sum2  = _mm_mul_ps(sum2, x2);

    sum1  = _mm_add_ps(sum1, sum2);
    sum1  = _mm_sub_ps(sum1, mone);
    sum1  = _mm_mul_ps(sum1, x);
    y     = _mm_add_ps(y, sum1);

    y     = _mm_xor_ps(y, sign);

    return y;
}


static __m256
gmx_mm256_atan2_ps(__m256 y, __m256 x)
{
    const __m256 pi          = _mm256_set1_ps( (float) M_PI);
    const __m256 minuspi     = _mm256_set1_ps( (float) -M_PI);
    const __m256 halfpi      = _mm256_set1_ps( (float) M_PI/2.0f);
    const __m256 minushalfpi = _mm256_set1_ps( (float) -M_PI/2.0f);

    __m256       z, z1, z3, z4;
    __m256       w;
    __m256       maskx_lt, maskx_eq;
    __m256       masky_lt, masky_eq;
    __m256       mask1, mask2, mask3, mask4, maskall;

    maskx_lt  = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OQ);
    masky_lt  = _mm256_cmp_ps(y, _mm256_setzero_ps(), _CMP_LT_OQ);
    maskx_eq  = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_EQ_OQ);
    masky_eq  = _mm256_cmp_ps(y, _mm256_setzero_ps(), _CMP_EQ_OQ);

    z         = _mm256_mul_ps(y, gmx_mm256_inv_ps(x));
    z         = gmx_mm256_atan_ps(z);

    mask1     = _mm256_and_ps(maskx_eq, masky_lt);
    mask2     = _mm256_andnot_ps(maskx_lt, masky_eq);
    mask3     = _mm256_andnot_ps( _mm256_or_ps(masky_lt, masky_eq), maskx_eq);
    mask4     = _mm256_and_ps(maskx_lt, masky_eq);
    maskall   = _mm256_or_ps( _mm256_or_ps(mask1, mask2), _mm256_or_ps(mask3, mask4) );

    z         = _mm256_andnot_ps(maskall, z);
    z1        = _mm256_and_ps(mask1, minushalfpi);
    z3        = _mm256_and_ps(mask3, halfpi);
    z4        = _mm256_and_ps(mask4, pi);

    z         = _mm256_or_ps( _mm256_or_ps(z, z1), _mm256_or_ps(z3, z4) );

    w         = _mm256_blendv_ps(pi, minuspi, masky_lt);
    w         = _mm256_and_ps(w, maskx_lt);

    w         = _mm256_andnot_ps(maskall, w);

    z         = _mm256_add_ps(z, w);

    return z;
}

static __m128
gmx_mm_atan2_ps(__m128 y, __m128 x)
{
    const __m128 pi          = _mm_set1_ps(M_PI);
    const __m128 minuspi     = _mm_set1_ps(-M_PI);
    const __m128 halfpi      = _mm_set1_ps(M_PI/2.0);
    const __m128 minushalfpi = _mm_set1_ps(-M_PI/2.0);

    __m128       z, z1, z3, z4;
    __m128       w;
    __m128       maskx_lt, maskx_eq;
    __m128       masky_lt, masky_eq;
    __m128       mask1, mask2, mask3, mask4, maskall;

    maskx_lt  = _mm_cmplt_ps(x, _mm_setzero_ps());
    masky_lt  = _mm_cmplt_ps(y, _mm_setzero_ps());
    maskx_eq  = _mm_cmpeq_ps(x, _mm_setzero_ps());
    masky_eq  = _mm_cmpeq_ps(y, _mm_setzero_ps());

    z         = _mm_mul_ps(y, gmx_mm_inv_ps(x));
    z         = gmx_mm_atan_ps(z);

    mask1     = _mm_and_ps(maskx_eq, masky_lt);
    mask2     = _mm_andnot_ps(maskx_lt, masky_eq);
    mask3     = _mm_andnot_ps( _mm_or_ps(masky_lt, masky_eq), maskx_eq);
    mask4     = _mm_and_ps(masky_eq, maskx_lt);

    maskall   = _mm_or_ps( _mm_or_ps(mask1, mask2), _mm_or_ps(mask3, mask4) );

    z         = _mm_andnot_ps(maskall, z);
    z1        = _mm_and_ps(mask1, minushalfpi);
    z3        = _mm_and_ps(mask3, halfpi);
    z4        = _mm_and_ps(mask4, pi);

    z         = _mm_or_ps( _mm_or_ps(z, z1), _mm_or_ps(z3, z4) );

    mask1     = _mm_andnot_ps(masky_lt, maskx_lt);
    mask2     = _mm_and_ps(maskx_lt, masky_lt);

    w         = _mm_or_ps( _mm_and_ps(mask1, pi), _mm_and_ps(mask2, minuspi) );
    w         = _mm_andnot_ps(maskall, w);

    z         = _mm_add_ps(z, w);

    return z;
}

#endif /* _gmx_math_x86_avx_256_single_h_ */
