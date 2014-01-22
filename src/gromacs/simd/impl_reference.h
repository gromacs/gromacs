/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team.
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

#ifndef GMX_SIMD_IMPL_REFERENCE_H
#define GMX_SIMD_IMPL_REFERENCE_H

#ifndef GMX_SIMD_SIMD_H
# error "Never include architecture-specific simd headers directly; use simd.h."
#endif

#include <math.h>

#include "gmx_fatal.h"

/* Float support present */
#define GMX_SIMD_HAVE_FLOAT

/* Double support present */
#define GMX_SIMD_HAVE_DOUBLE

/* Reference implementation uses software */
#undef GMX_SIMD_HAVE_HARDWARE

/* Reference implementation does not care about load alignment */
#define GMX_SIMD_HAVE_LOADU

/* Reference implementation does not care about store alignment */
#define GMX_SIMD_HAVE_STOREU

/* Support for logical operations on floating-point data */
#define GMX_SIMD_HAVE_LOGICAL

/* No guarantees about fused multiply-add in reference implementation */
#undef GMX_SIMD_HAVE_FMA

/* No guarantees about hardware fraction instruction in reference impl. */
#undef GMX_SIMD_HAVE_FRACTION

/* Integers corresponding to float SIMD. */
#define GMX_SIMD_HAVE_FINT32

/* We provide software routines to extract integers from SIMD types */
#define GMX_SIMD_HAVE_FINT32_EXTRACT

/* Logical ops supported for gmx_simd_fint32_t */
#define GMX_SIMD_HAVE_FINT32_LOGICAL

/* Arithmetics supported for gmx_simd_fint32_t */
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS

/* Integers corresponding to double SIMD. */
#define GMX_SIMD_HAVE_DINT32

/* We provide software routines to extract integers from SIMD types */
#define GMX_SIMD_HAVE_DINT32_EXTRACT

/* Logical ops supported for gmx_simd_dint32_t */
#define GMX_SIMD_HAVE_DINT32_LOGICAL

/* Arithmetics supported for gmx_simd_dint32_t */
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS

/* SIMD4 float enabled in reference implementation if float SIMD width==4 */
#define GMX_SIMD4_HAVE_FLOAT

/* SIMD4 double is enabled if normal double SIMD has width=4 */
#define GMX_SIMD4_HAVE_DOUBLE

/* Software reference is the only arch where we can change SIMD width. */
#define GMX_SIMD_FLOAT_WIDTH             4

/* Software reference is the only arch where we can change SIMD width. */
#define GMX_SIMD_DOUBLE_WIDTH            4

/* Software reference is the only arch where we can change SIMD width. */
#define GMX_SIMD_FINT32_WIDTH            GMX_SIMD_FLOAT_WIDTH

/* Software reference is the only arch where we can change SIMD width. */
#define GMX_SIMD_DINT32_WIDTH            GMX_SIMD_DOUBLE_WIDTH

/* The reference 1/sqrt(x) lookup provides full single precision */
#define GMX_SIMD_RSQRT_BITS             23

/* The reference 1/x lookup provides full single precision */
#define GMX_SIMD_RCP_BITS               23

/* SINGLE-PRECISION FLOATING-POINT OPERATIONS */

/* Float SIMD register.
 * Supported with GMX_SIMD_HAVE_FLOAT.
 */
typedef struct
{
    float r[GMX_SIMD_FLOAT_WIDTH];
}
gmx_simd_float_t;


/* Load GMX_SIMD_FLOAT_WIDTH numbers from aligned memory starting at m. */
static gmx_inline gmx_simd_float_t
gmx_simd_load_f(const float *m)
{
    gmx_simd_float_t  a;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/* Set all SIMD register elements to value pointed to by m (unaligned) */
static gmx_inline gmx_simd_float_t
gmx_simd_load1_f(float *m)
{
    gmx_simd_float_t  a;
    int               i;
    float             f = *m;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = f;
    }
    return a;
}

/* Set all SIMD register elements to the value r. */
static gmx_inline gmx_simd_float_t
gmx_simd_set1_f(float r)
{
    gmx_simd_float_t  a;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}


/* Store the contents of the SIMD register pr to aligned memory m. */
static gmx_inline void
gmx_simd_store_f(float *m, gmx_simd_float_t pr)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        m[i] = pr.r[i];
    }
}

/* Load from unaligned memory. Available with GMX_SIMD_HAVE_LOADU. */
#define gmx_simd_loadu_f gmx_simd_load_f


/* Store to unaligned memory. Available with GMX_SIMD_HAVE_STOREU. */
#define gmx_simd_storeu_f gmx_simd_store_f

/* Floating-point arithmetic */

/* Set all SIMD register elements to 0.0f */
static gmx_inline gmx_simd_float_t
gmx_simd_setzero_f()
{
    gmx_simd_float_t  a;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }
    return a;
}

/* Add two float SIMD registers */
static gmx_inline gmx_simd_float_t
gmx_simd_add_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }
    return c;
}

/* Subtract two SIMD registers */
static gmx_inline gmx_simd_float_t
gmx_simd_sub_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }
    return c;
}

/* Multiply two SIMD registers */
static gmx_inline gmx_simd_float_t
gmx_simd_mul_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i]*b.r[i];
    }
    return c;
}

/* Fused-multiply-add. Result is a*b+c.
 *
 * Supported in hardware if GMX_SIMD_HAVE_FMA is set, emulated otherwise.
 */
#define gmx_simd_fmadd_f(a, b, c) gmx_simd_add_f(gmx_simd_mul_f(a, b), c)


/* Fused-multiply-subtract. Result is a*b-c.
 *
 * Supported in hardware if GMX_SIMD_HAVE_FMA is set, emulated otherwise.
 */
#define gmx_simd_fmsub_f(a, b, c) gmx_simd_sub_f(gmx_simd_mul_f(a, b), c)


/* Fused-negated-multiply-add. Result is -a*b+c.
 *
 * Supported in hardware if GMX_SIMD_HAVE_FMA is set, emulated otherwise.
 */
#define gmx_simd_fnmadd_f(a, b, c) gmx_simd_sub_f(c, gmx_simd_mul_f(a, b))


/* Fused-negated-multiply-add. Result is -a*b-c.
 *
 * Supported in hardware if GMX_SIMD_HAVE_FMA is set, emulated otherwise.
 */
#define gmx_simd_fnmsub_f(a, b, c) gmx_simd_sub_f(gmx_simd_setzero_f(), gmx_simd_fmadd_f(a, b, c))


/* Bitwise and for two SIMD fp registers.
 * Supported with GMX_SIMD_HAVE_LOGICAL.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_and_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    int               val1, val2, res;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = val1 & val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        val1   = *((int *)&(a.r[i]));
        val2   = *((int *)&(b.r[i]));
        res    = val1 & val2;
        c.r[i] = *((float *)&res);
#endif
    }
    return c;
}

/* Bitwise andnot for two SIMD fp registers. c=(~a) & b
 * Supported with GMX_SIMD_HAVE_LOGICAL.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_andnot_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    int               val1, val2, res;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = (~val1) & val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        val1   = *((int *)&(a.r[i]));
        val2   = *((int *)&(b.r[i]));
        res    = (~val1) & val2;
        c.r[i] = *((float *)&res);
#endif
    }
    return c;
}

/* Bitwise or for two SIMD fp registers.
 * Supported with GMX_SIMD_HAVE_LOGICAL.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_or_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    int               val1, val2, res;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = val1 | val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        val1   = *((int *)&(a.r[i]));
        val2   = *((int *)&(b.r[i]));
        res    = val1 | val2;
        c.r[i] = *((float *)&res);
#endif
    }
    return c;
}

/* Bitwise xor for two SIMD fp registers.
 * Supported with GMX_SIMD_HAVE_LOGICAL.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_xor_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    int               val1, val2, res;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = val1 ^ val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        val1   = *((int *)&(a.r[i]));
        val2   = *((int *)&(b.r[i]));
        res    = val1 ^ val2;
        c.r[i] = *((float *)&res);
#endif
    }
    return c;
}

/* 1.0/sqrt(x) lookup. */
static gmx_inline gmx_simd_float_t
gmx_simd_rsqrt_f(gmx_simd_float_t x)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = (x.r[i] > 0.0f) ? 1.0f/sqrtf(x.r[i]) : 0.0f;
    }
    return b;
};

/* 1.0/x lookup */
static gmx_inline gmx_simd_float_t
gmx_simd_rcp_f(gmx_simd_float_t x)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = (x.r[i] != 0.0f) ? 1.0f/x.r[i] : 0.0f;
    }
    return b;
};

/* Floating-point abs */
static gmx_inline gmx_simd_float_t
gmx_simd_fabs_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = fabsf(a.r[i]);
    }
    return c;
}

/* Floating-point negate */
static gmx_inline gmx_simd_float_t
gmx_simd_fneg_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/* Set each SIMD element to the largest from two registers. */
static gmx_inline gmx_simd_float_t
gmx_simd_max_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/* Set each SIMD element to the smallest from two registers. */
static gmx_inline gmx_simd_float_t
gmx_simd_min_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = (a.r[i] <= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/* Round to nearest integer value (in floating-point format).
 *
 * Note that this reference implementation rounds exact half-way cases
 * away from zero, whereas most SIMD intrinsics will round to nearest even.
 * This could be fixed by using rint/rintf, but the bigger problem is that
 * MSVC does not support full C99, and none of the round or rint
 * functions are defined. It's much easier to approximately implement
 * round() than rint(), so we do that and hope we never get bitten in
 * testing. (Thanks, Microsoft.)
 */
static gmx_inline gmx_simd_float_t
gmx_simd_round_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef _MSC_VER
        int temp = (a.r[i] >= 0.0f) ? (a.r[i] + 0.5f) : (a.r[i] - 0.5f);
        b.r[i] = temp;
#else
        b.r[i] = roundf(a.r[i]);
#endif
    }
    return b;
}

/* Truncate, i.e. round towards zero - common hardware instruction.
 *
 * Note that this is truncation towards zero, not "floor". The reason for this
 * is that truncation is virtually always present as a dedicated hardware
 * instruction, but "floor" frequently isn't.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_trunc_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = truncf(a.r[i]);
    }
    return b;
}


/* Fraction of the floating point number.
 *
 * To maximize compatibility, we use the same definition of fractions as used
 * e.g. for the AMD64 hardware instructions. This relies on truncation towards
 * zero for the integer part, and the remaining fraction can thus be either
 * positive or negative. As an example, -1.42 would return the fraction -0.42.
 *
 * Hardware support with GMX_SIMD_HAVE_FRACTION, otherwise emulated.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_fraction_f(gmx_simd_float_t a)
{
    return gmx_simd_sub_f(a, gmx_simd_trunc_f(a));
}

/* Extract (integer) exponent from single precision floating-point number
 * and make it a float.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f(gmx_simd_float_t a)
{
    /* Mask with ones for the exponent field of single precision fp */
    const int          expmask = 0x7f800000;
    gmx_simd_float_t   b;
    int                i;
    union
    {
        float f;
        int   i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        /* Keep exponent, shift 23 right (float mantissa), remove bias (127) */
        b.r[i] = ((conv.i & expmask) >> 23) - 127;
    }
    return b;
}

/* Get the mantissa - similar to gmx_simd_get_exponent_f() */
static gmx_inline gmx_simd_float_t
gmx_simd_get_mantissa_f(gmx_simd_float_t a)
{
    const int          mantmask = 0x007fffff;
    const int          one      = 0x3f800000;
    gmx_simd_float_t   b;
    int                i;
    union
    {
        float f;
        int   i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        // remove current exponent, add a biased exponent for 1.0 (i.e., 2^0=1) */
        conv.i = (conv.i & (mantmask)) | one;
        b.r[i] = conv.f;
    }
    return b;
}

/* Set (integer) exponent from single precision floating-point number.
 * The argument will be _rounded_ to nearest integer since that is what we
 * need for the exponential functions, and this integer x will be set as the
 * exponent so the new fp number will be 2^x.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f(gmx_simd_float_t a)
{
    gmx_simd_float_t   b;
    int                i, iexp;
    union
    {
        float f;
        int   i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        // Critical to use same algorithm as for gmx_simd_round_f()
#ifdef _MSC_VER
        iexp = (a.r[i] >= 0.0f) ? (a.r[i] + 0.5f) : (a.r[i] - 0.5f);
#else
        iexp = roundf(a.r[i]);
#endif
        /* Add bias (127), and shift 23 bits left (mantissa size) */
        conv.i = (iexp + 127) << 23;
        b.r[i] = conv.f;
    }
    return b;
}


/* INTEGER SIMD datatypes */

/* Integer SIMD register type to use for conversions to/from float. This is
 * also the widest integer SIMD type.
 */
typedef struct
{
    int i[GMX_SIMD_FINT32_WIDTH];
}
gmx_simd_fint32_t;

#if (GMX_SIMD_FINT32_WIDTH < GMX_SIMD_FLOAT_WIDTH)
#    undef GMX_SIMD_HAVE_FINT32
#endif

/* Load aligned integer data */
static gmx_inline gmx_simd_fint32_t
gmx_simd_load_fi(int * m)
{
    gmx_simd_fint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/* set from integer */
static gmx_inline gmx_simd_fint32_t
gmx_simd_set1_fi(int b)
{
    gmx_simd_fint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/* Store aligned integer data */
static gmx_inline gmx_simd_fint32_t
gmx_simd_store_fi(int * m, gmx_simd_fint32_t a)
{
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
    return a;
};

#define gmx_simd_loadu_fi  gmx_simd_load_fi
#define gmx_simd_storeu_fi gmx_simd_store_fi

/* Set all SIMD register elements to 0. */
static gmx_inline gmx_simd_fint32_t
gmx_simd_setzero_fi()
{
    gmx_simd_fint32_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = 0;
    }
    return a;
}

/* Round single precision floating point to integer.
 *
 * Available with GMX_SIMD_HAVE_FINT.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_cvt_f2i(gmx_simd_float_t a)
{
    gmx_simd_fint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
#ifdef _MSC_VER
        b.i[i] = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
#else
        b.i[i] = roundf(a.r[i]);
#endif
    }
    return b;
};

/* Truncate single precision floating point to integer.
 *
 * Available with GMX_SIMD_HAVE_FINT.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_cvtt_f2i(gmx_simd_float_t a)
{
    gmx_simd_fint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.i[i] = a.r[i];
    }
    return b;
};

/* Convert integer to single precision floating-point.
 *
 * Available with GMX_SIMD_HAVE_FINT.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_cvt_i2f(gmx_simd_fint32_t a)
{
    gmx_simd_float_t   b;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.r[i] = a.i[i];
    }
    return b;
};

/* Extract element with index i from an integer SIMD register.
 * Available with GMX_SIMD_HAVE_FINT32_EXTRACT.
 */
static gmx_inline int
gmx_simd_extract_fi(gmx_simd_fint32_t a, int index)
{
    return a.i[index];
}

/* SIMD integer shift left, based on immediate value.
 * Available with GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_slli_fi(gmx_simd_fint32_t a, int n)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] << n;
    }
    return c;
}

/* SIMD integer shift right, based on immediate value.
 * Available with GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_srli_fi(gmx_simd_fint32_t a, int n)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] >> n;
    }
    return c;
}

/* Integer bitwise and.
 * Available with GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_and_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] & b.i[i];
    }
    return c;
}

/* Integer bitwise not-and.
 * Available with GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_andnot_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = (~a.i[i]) & b.i[i];
    }
    return c;
}

/* Integer bitwise or.
 * Available with GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_or_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] | b.i[i];
    }
    return c;
}

/* Integer bitwise xor.
 * Available with GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_xor_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] ^ b.i[i];
    }
    return c;
}

/* Add SIMD integers. Available with GMX_SIMD_HAVE_FINT32_ARIMTHETICS. */
static gmx_inline gmx_simd_fint32_t
gmx_simd_add_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] + b.i[i];
    }
    return c;
}

/* Subtract SIMD integers. Available with GMX_SIMD_HAVE_FINT32_ARIMTHETICS. */
static gmx_inline gmx_simd_fint32_t
gmx_simd_sub_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] - b.i[i];
    }
    return c;
}

/* Multiply SIMD integers. Available with GMX_SIMD_HAVE_FINT32_ARIMTHETICS. */
static gmx_inline gmx_simd_fint32_t
gmx_simd_mul_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i]*b.i[i];
    }
    return c;
}

/* BOOLEAN OPERATIONS */

/* Boolean type for float SIMD data. You should likely use gmx_simd_bool_t
 * (for gmx_simd_real_t) instead, unless you really know what you are doing.
 */
typedef struct
{
    int b[GMX_SIMD_FLOAT_WIDTH];
}
gmx_simd_fbool_t;

/* Equality comparison of two single precision fp values */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmpeq_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/* Less-than comparison of two single precision fp values */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmplt_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_fbool_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}

/* Less-than-or-equal comparison of two single precision fp values */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmple_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_fbool_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}

/* Logical AND on single precision SIMD booleans */
static gmx_inline gmx_simd_fbool_t
gmx_simd_and_fb(gmx_simd_fbool_t a, gmx_simd_fbool_t b)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/* Logical OR on single precision SIMD booleans */
static gmx_inline gmx_simd_fbool_t
gmx_simd_or_fb(gmx_simd_fbool_t a, gmx_simd_fbool_t b)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/* Returns non-zero if any of the boolean in x is True, otherwise 0.
 * The actual return value for "any true" will depend on the architecture!
 */
static gmx_inline int
gmx_simd_anytrue_fb(gmx_simd_fbool_t a)
{
    int             anytrue;
    int             i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}

/* Select from single precision fp SIMD register based on a boolean one. */
static gmx_inline gmx_simd_float_t
gmx_simd_blendzero_f(gmx_simd_float_t a, gmx_simd_fbool_t sel)
{
    gmx_simd_float_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? a.r[i] : 0.0;
    }
    return c;
}

/* Vector-blend instruction. For each element, select the value from register b
 * if the contents of the selector is true, otherwise from register a.
 * Emulated if not supported in hardware.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_blendv_f(gmx_simd_float_t a, gmx_simd_float_t b, gmx_simd_fbool_t sel)
{
    gmx_simd_float_t  d;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}


/* Boolean type for integer datatypes corresponding to float SIMD. */
typedef struct
{
    int b[GMX_SIMD_FINT32_WIDTH];
}
gmx_simd_fibool_t;

/* Equality comparison of two integers corresponding to float values.
 * Available with GMX_SIMD_HAVE_FINT32_ARIMTHETICS.
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_cmpeq_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] == b.i[i]);
    }
    return c;
}

/* Less-than comparison of two integers corresponding to float values.
 * Available with GMX_SIMD_HAVE_FINT32_ARIMTHETICS.
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_cmplt_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] < b.i[i]);
    }
    return c;
}

/* Logical AND on gmx_simd_fibool_t */
static gmx_inline gmx_simd_fibool_t
gmx_simd_and_fib(gmx_simd_fibool_t a, gmx_simd_fibool_t b)
{
    gmx_simd_fibool_t c;
    int               i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/* Logical OR on gmx_simd_fibool_t  */
static gmx_inline gmx_simd_fibool_t
gmx_simd_or_fib(gmx_simd_fibool_t a, gmx_simd_fibool_t b)
{
    gmx_simd_fibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/* Returns non-zero if any of the boolean in x is True, otherwise 0.
 * The actual return value for "any true" will depend on the architecture!
 */
static gmx_inline int
gmx_simd_anytrue_fib(gmx_simd_fibool_t a)
{
    int             anytrue;
    int             i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}

/* Select from single precision fp SIMD register based on a boolean. */
static gmx_inline gmx_simd_fint32_t
gmx_simd_blendzero_fi(gmx_simd_fint32_t a, gmx_simd_fibool_t sel)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? a.i[i] : 0.0;
    }
    return c;
}

/* Vector-blend instruction. For each element, select the value from register b
 * if the contents of the selector is true, otherwise from register a.
 * Emulated if not supported in hardware.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_blendv_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b, gmx_simd_fibool_t sel)
{
    gmx_simd_fint32_t d;
    int               i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        d.i[i] = sel.b[i] ? b.i[i] : a.i[i];
    }
    return d;
}

/* Convert from float boolean to corresponding integer boolean */
static gmx_inline gmx_simd_fibool_t
gmx_simd_cvt_fb2fib(gmx_simd_fbool_t a)
{
    gmx_simd_fibool_t  b;
    int                i;

    /* Integer width >= float width */
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/* Convert from integer boolean (corresponding to float) to float boolean */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cvt_fib2fb(gmx_simd_fibool_t a)
{
    gmx_simd_fbool_t  b;
    int               i;

    /* Integer width >= float width */
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}


/* Floating-point SIMD register type in double precision.
 * Supported with GMX_SIMD_HAVE_DOUBLE.
 */
typedef struct
{
    double r[GMX_SIMD_DOUBLE_WIDTH];
}
gmx_simd_double_t;

/* Double precision version of gmx_simd_load_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_load_d(const double *m)
{
    gmx_simd_double_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/* Double precision version of gmx_simd_load1_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_load1_d(double *m)
{
    gmx_simd_double_t  a;
    int                i;
    double             d = *m;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = d;
    }
    return a;
}

/* Double precision version of gmx_simd_set1_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_set1_d(double r)
{
    gmx_simd_double_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}

/* Double precision version of gmx_simd_store_f() */
static gmx_inline void
gmx_simd_store_d(double *m, gmx_simd_double_t pr)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        m[i] = pr.r[i];
    }
}

/* Double precision version of gmx_simd_loadu_f().
 *
 * Available with GMX_SIMD_HAVE_LOADU.
 */
#define gmx_simd_loadu_d gmx_simd_load_d

/* Double precision version of gmx_simd_storeu_f().
 *
 * Available with GMX_SIMD_HAVE_STOREU.
 */
#define gmx_simd_storeu_d gmx_simd_store_d

/* Double precision version of gmx_simd_setzero_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_setzero_d()
{
    gmx_simd_double_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }
    return a;
}

/* Double precision version of gmx_simd_add_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_add_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }
    return c;
}

/* Double precision version of gmx_simd_sub_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_sub_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }
    return c;
}

/* Double precision version of gmx_simd_mul_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_mul_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i]*b.r[i];
    }
    return c;
}

/* Double precision version of gmx_simd_madd_f() */
#define gmx_simd_fmadd_d(a, b, c) gmx_simd_add_d(gmx_simd_mul_d(a, b), c)

/* Double precision version of gmx_simd_msub_f() */
#define gmx_simd_fmsub_d(a, b, c) gmx_simd_sub_d(gmx_simd_mul_d(a, b), c)

/* Double precision version of gmx_simd_nmadd_f() */
#define gmx_simd_fnmadd_d(a, b, c) gmx_simd_sub_d(c, gmx_simd_mul_d(a, b))

/* Double precision version of gmx_simd_nmsub_f() */
#define gmx_simd_fnmsub_d(a, b, c) gmx_simd_sub_d(gmx_simd_setzero_d(), gmx_simd_fmadd_d(a, b, c))



/* Double precision version of gmx_simd_and_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_and_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    gmx_int64_t        val1, val2, res;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = val1 & val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        val1   = *((gmx_int64_t *)&(a.r[i]));
        val2   = *((gmx_int64_t *)&(b.r[i]));
        res    = val1 & val2;
        c.r[i] = *((double *)&res);
#endif
    }
    return c;
}

/* Double precision version of gmx_simd_andnot_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_andnot_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    gmx_int64_t        val1, val2, res;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = (~val1) & val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        val1   = *((gmx_int64_t *)&(a.r[i]));
        val2   = *((gmx_int64_t *)&(b.r[i]));
        res    = (~val1) & val2;
        c.r[i] = *((double *)&res);
#endif
    }
    return c;
}

/* Double precision version of gmx_simd_or_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_or_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    gmx_int64_t        val1, val2, res;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = val1 | val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        val1   = *((gmx_int64_t *)&(a.r[i]));
        val2   = *((gmx_int64_t *)&(b.r[i]));
        res    = val1 | val2;
        c.r[i] = *((double *)&res);
#endif
    }
    return c;
}

/* Double precision version of gmx_simd_xor_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_xor_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    gmx_int64_t        val1, val2, res;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = val1 ^ val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        val1   = *((gmx_int64_t *)&(a.r[i]));
        val2   = *((gmx_int64_t *)&(b.r[i]));
        res    = val1 ^ val2;
        c.r[i] = *((double *)&res);
#endif
    }
    return c;
}

/* 1.0/sqrt(x) lookup. */
static gmx_inline gmx_simd_double_t
gmx_simd_rsqrt_d(gmx_simd_double_t x)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        /* Sic - we only need single precision for the reference lookup, since
         * we have defined GMX_SIMD_RSQRT_BITS to 23.
         */
        b.r[i] = (x.r[i] > 0.0) ? 1.0f/sqrtf(x.r[i]) : 0.0;
    }
    return b;
};

/* 1.0/x lookup */
static gmx_inline gmx_simd_double_t
gmx_simd_rcp_d(gmx_simd_double_t x)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        /* Sic - we only need single precision for the reference lookup, since
         * we have defined GMX_SIMD_RCP_BITS to 23.
         */
        b.r[i] = (x.r[i] != 0.0) ? 1.0f/x.r[i] : 0.0;
    }
    return b;
};

/* Double precision version of gmx_simd_fabs_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_fabs_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = fabs(a.r[i]);
    }
    return c;
}

/* Double precision version of gmx_simd_fneg_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_fneg_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/* Double precision version of gmx_simd_max_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_max_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/* Double precision version of gmx_simd_min_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_min_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = (a.r[i] <= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/* Double precision version of gmx_simd_round_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_round_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef _MSC_VER
        int temp = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
        b.r[i] = temp;
#else
        b.r[i] = round(a.r[i]);
#endif
    }
    return b;
}

/* Double precision version of gmx_simd_trunc_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_trunc_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.r[i] = trunc(a.r[i]);
    }
    return b;
}

/* Double precision version of gmx_simd_fraction_f()
 *
 * Hardware with GMX_SIMD_HAVE_FRACTION, emulated as x-trunc(x) otherwise.
 */
static gmx_inline gmx_simd_double_t
gmx_simd_fraction_d(gmx_simd_double_t a)
{
    return gmx_simd_sub_d(a, gmx_simd_trunc_d(a));
}


/* Extract (integer) exponent from double precision floating-point number
 * and make it a float.
 */
static gmx_inline gmx_simd_double_t
gmx_simd_get_exponent_d(gmx_simd_double_t a)
{
    /* Mask with ones for the exponent field of double precision fp */
    const gmx_int64_t      expmask = 0x7ff0000000000000LL;
    gmx_simd_double_t      b;
    int                    i;
    union
    {
        double             d;
        gmx_int64_t        i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv.d = a.r[i];
        /* Zero everything but exponent field (remove sign),
         * shift 23 bits right (mantissa size), and remove exponent bias (1023).
         */
        b.r[i] = ((conv.i & expmask) >> 52) - 1023;
    }
    return b;
}

/* Get the mantissa - similar to gmx_simd_get_exponent_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_get_mantissa_d(gmx_simd_double_t a)
{
    const gmx_int64_t      mantmask = 0x000fffffffffffffLL;
    const gmx_int64_t      one      = 0x3ff0000000000000LL;
    gmx_simd_double_t      b;
    int                    i;
    union
    {
        double          d;
        gmx_int64_t     i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv.d = a.r[i];
        conv.i = (conv.i & (mantmask)) | one;
        b.r[i] = conv.d;
    }
    return b;
}

/* Set (integer) exponent from double precision floating-point number.
 * The argument will be _rounded_ to nearest integer since that is what we
 * need for the exponential functions, and this integer x will be set as the
 * exponent so the new fp number will be 2^x.
 */
static gmx_inline gmx_simd_double_t
gmx_simd_set_exponent_d(gmx_simd_double_t a)
{
    gmx_simd_double_t      b;
    int                    i;
    gmx_int64_t            iexp;
    union
    {
        double          d;
        gmx_int64_t     i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        // Critical to use same algorithm as for gmx_simd_round_d()
#ifdef _MSC_VER
        iexp = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
#else
        iexp = round(a.r[i]);
#endif
        /* Add bias (1023), and shift 52 bits left (mantissa size) */
        conv.i = (iexp + 1023) << 52;
        b.r[i] = conv.d;
    }
    return b;
}


/* Integer SIMD register type to use for conversions to/from double.
 * Available with GMX_SIMD_HAVE_DINT32.
 */
typedef struct
{
    int i[GMX_SIMD_DINT32_WIDTH];
}
gmx_simd_dint32_t;
#if (GMX_SIMD_DINT32_WIDTH < GMX_SIMD_DOUBLE_WIDTH)
#    undef GMX_SIMD_HAVE_DINT32
#endif

/* Load aligned integer data */
static gmx_inline gmx_simd_dint32_t
gmx_simd_load_di(int * m)
{
    gmx_simd_dint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/* set from integer */
static gmx_inline gmx_simd_dint32_t
gmx_simd_set1_di(int b)
{
    gmx_simd_dint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/* Store aligned integer data */
static gmx_inline gmx_simd_dint32_t
gmx_simd_store_di(int * m, gmx_simd_dint32_t a)
{
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
    return a;
};

#define gmx_simd_loadu_di  gmx_simd_load_di
#define gmx_simd_storeu_di gmx_simd_store_di

/* Set all SIMD register elements to 0. */
static gmx_inline gmx_simd_dint32_t
gmx_simd_setzero_di()
{
    gmx_simd_dint32_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = 0;
    }
    return a;
}

/* Round double precision floating point to integer.
 *
 * Available with GMX_SIMD_HAVE_DINT.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_cvt_d2i(gmx_simd_double_t a)
{
    gmx_simd_dint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
#ifdef _MSC_VER
        b.i[i] = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
#else
        b.i[i] = round(a.r[i]);
#endif
    }
    return b;
};

/* Truncate double precision floating point to integer.
 *
 * Available with GMX_SIMD_HAVE_DINT.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_cvtt_d2i(gmx_simd_double_t a)
{
    gmx_simd_dint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        b.i[i] = a.r[i];
    }
    return b;
};

/* Convert integer to single precision floating-point.
 *
 * Available with GMX_SIMD_HAVE_DINT.
 */
static gmx_inline gmx_simd_double_t
gmx_simd_cvt_i2d(gmx_simd_dint32_t a)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        b.r[i] = a.i[i];
    }
    return b;
};

/* Extract element with index i from an integer SIMD register.
 * Available with GMX_SIMD_HAVE_DINT32_EXTRACT.
 */
static gmx_inline int
gmx_simd_extract_di(gmx_simd_dint32_t a, int index)
{
    return a.i[index];
}

/* SIMD integer shift left, based on immediate value.
 * Available with GMX_SIMD_HAVE_DINT32_LOGICAL.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_slli_di(gmx_simd_dint32_t a, int n)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] << n;
    }
    return c;
}

/* SIMD integer shift right, based on immediate value.
 * Available with GMX_SIMD_HAVE_DINT32_LOGICAL.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_srli_di(gmx_simd_dint32_t a, int n)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] >> n;
    }
    return c;
}

/* Integer bitwise and.
 * Available with GMX_SIMD_HAVE_DINT32_LOGICAL.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_and_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] & b.i[i];
    }
    return c;
}

/* Integer bitwise not-and.
 * Available with GMX_SIMD_HAVE_DINT32_LOGICAL.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_andnot_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = (~a.i[i]) & b.i[i];
    }
    return c;
}

/* Integer bitwise or.
 * Available with GMX_SIMD_HAVE_DINT32_LOGICAL.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_or_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] | b.i[i];
    }
    return c;
}

/* Integer bitwise xor.
 * Available with GMX_SIMD_HAVE_DINT32_LOGICAL.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_xor_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] ^ b.i[i];
    }
    return c;
}

/* Add SIMD integers. Available with GMX_SIMD_HAVE_DINT32_ARIMTHETICS. */
static gmx_inline gmx_simd_dint32_t
gmx_simd_add_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] + b.i[i];
    }
    return c;
}

/* Subtract SIMD integers. Available with GMX_SIMD_HAVE_DINT32_ARIMTHETICS. */
static gmx_inline gmx_simd_dint32_t
gmx_simd_sub_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] - b.i[i];
    }
    return c;
}

/* Multiply SIMD integers. Available with GMX_SIMD_HAVE_DINT32_ARIMTHETICS. */
static gmx_inline gmx_simd_dint32_t
gmx_simd_mul_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i]*b.i[i];
    }
    return c;
}




/* Boolean type for double precision data. Use the generic gmx_simd_bool_t
 * (for gmx_simd_real_t) instead, unless you really know what you are doing.
 */
typedef struct
{
    int b[GMX_SIMD_DOUBLE_WIDTH];
}
gmx_simd_dbool_t;

/* Equality comparison of two double precision fp values.
 *
 * Available with GMX_SIMD_HAVE_DOUBLE.
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmpeq_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/* Less-than comparison of two double precision fp values
 *
 * Available with GMX_SIMD_HAVE_DOUBLE.
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmplt_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}

/* Less-than-or-equal comparison of two double precision fp values
 *
 * Available with GMX_SIMD_HAVE_DOUBLE.
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmple_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}


/* Logical AND on double precision SIMD booleans */
static gmx_inline gmx_simd_dbool_t
gmx_simd_and_db(gmx_simd_dbool_t a, gmx_simd_dbool_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/* Logical OR on double precision SIMD booleans */
static gmx_inline gmx_simd_dbool_t
gmx_simd_or_db(gmx_simd_dbool_t a, gmx_simd_dbool_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}


/* Returns non-zero if any of the boolean in x is True, otherwise 0.
 * The actual return value for "any true" will depend on the architecture!
 */
static gmx_inline int
gmx_simd_anytrue_db(gmx_simd_dbool_t a)
{
    int         anytrue;
    int         i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}


/* Select from a double precision fp SIMD register based on bool */
static gmx_inline gmx_simd_double_t
gmx_simd_blendzero_d(gmx_simd_double_t a, gmx_simd_dbool_t sel)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? a.r[i] : 0.0;
    }
    return c;
}

/* Double precision version of gmx_simd_blendv_f() */
static gmx_inline gmx_simd_double_t
gmx_simd_blendv_d(gmx_simd_double_t a, gmx_simd_double_t b, gmx_simd_dbool_t sel)
{
    gmx_simd_double_t  d;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}


/* Boolean type for integer datatypes corresponding to double SIMD. You should
 * likely use gmx_simd_ibool_t (for gmx_simd_int32_t) instead, unless you really
 * know what you are doing.
 */
typedef struct
{
    int b[GMX_SIMD_DINT32_WIDTH];
}
gmx_simd_dibool_t;

/* Equality comparison of two ints corresponding to double SIMD data.
 * Available with GMX_SIMD_HAVE_DINT32_ARIMTHETICS.
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_cmpeq_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] == b.i[i]);
    }
    return c;
}

/* Less-than comparison of two ints corresponding to double SIMD data
 * Available with GMX_SIMD_HAVE_DINT32_ARIMTHETICS.
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_cmplt_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] < b.i[i]);
    }
    return c;
}

/* Logical AND on gmx_simd_dibool_t */
static gmx_inline gmx_simd_dibool_t
gmx_simd_and_dib(gmx_simd_dibool_t a, gmx_simd_dibool_t b)
{
    gmx_simd_dibool_t c;
    int               i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/* Logical OR on gmx_simd_dibool_t */
static gmx_inline gmx_simd_dibool_t
gmx_simd_or_dib(gmx_simd_dibool_t a, gmx_simd_dibool_t b)
{
    gmx_simd_dibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/* Returns non-zero if any of the boolean in x is True, otherwise 0.
 * The actual return value for "any true" will depend on the architecture!
 */
static gmx_inline int
gmx_simd_anytrue_dib(gmx_simd_dibool_t a)
{
    int             anytrue;
    int             i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}

/* Select from SIMD ints (corresponding to double) based on a boolean one. */
static gmx_inline gmx_simd_dint32_t
gmx_simd_blendzero_di(gmx_simd_dint32_t a, gmx_simd_dibool_t sel)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? a.i[i] : 0.0;
    }
    return c;
}

/* Vector-blend instruction. For each element, select the value from register b
 * if the contents of the selector is true, otherwise from register a.
 * Emulated if not supported in hardware.
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_blendv_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b, gmx_simd_dibool_t sel)
{
    gmx_simd_dint32_t  d;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        d.i[i] = sel.b[i] ? b.i[i] : a.i[i];
    }
    return d;
}

/* Convert from double boolean to corresponding integer boolean */
static gmx_inline gmx_simd_dibool_t
gmx_simd_cvt_db2dib(gmx_simd_dbool_t a)
{
    gmx_simd_dibool_t  b;
    int                i;

    /* Integer width >= double width */
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/* Convert from integer boolean (corresponding to double) to double boolean */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cvt_dib2db(gmx_simd_dibool_t a)
{
    gmx_simd_dbool_t  b;
    int               i;

    /* Integer width >= double width */
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/* Convert float to double.
 *
 * This version is available if GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 */
static gmx_inline gmx_simd_double_t
gmx_simd_cvt_f2d(gmx_simd_float_t f)
{
    gmx_simd_double_t d;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int               i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = f.r[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_f2d() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings */
    d = gmx_simd_setzero_d();
#endif
    return d;
}

/* Convert double to float.
 *
 * This version is available if GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH.
 *
 * See gmx_simd_cvt_d2f().
 */
static gmx_inline gmx_simd_float_t
gmx_simd_cvt_d2f(gmx_simd_double_t d)
{
    gmx_simd_float_t f;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.r[i] = d.r[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_d2f() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings */
    f = gmx_simd_setzero_f();
#endif
    return f;
}

/* Convert float to double.
 *
 * This version is available if GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH.
 */
static gmx_inline void
gmx_simd_cvt_f2dd(gmx_simd_float_t f, gmx_simd_double_t *d0, gmx_simd_double_t *d1)
{
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d0->r[i] = f.r[i];
        d1->r[i] = f.r[GMX_SIMD_DOUBLE_WIDTH+i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_f2dd() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments */
    d0->r[0] = f.r[0];
    d1->r[0] = f.r[0];
#endif
}

/* Convert double to float.
 *
 * This version is available if GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_cvt_dd2f(gmx_simd_double_t d0, gmx_simd_double_t d1)
{
    gmx_simd_float_t f;
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.r[i]                       = d0.r[i];
        f.r[GMX_SIMD_DOUBLE_WIDTH+i] = d1.r[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_dd2f() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments & uninitialized f */
    f.r[0] = d0.r[0] + d1.r[0];
#endif
    return f;
}




#if (GMX_SIMD_FLOAT_WIDTH == 4)

/* SIMD4 float type. Available with GMX_SIMD_HAVE_SIMD4_FLOAT, but you should
 * typically use GMX_SIMD_HAVE_SIMD4_REAL to check default Gromacs precision.
 */
#    define gmx_simd4_float_t    gmx_simd_float_t

/* SIMD4 version of gmx_simd_load_f() */
#    define gmx_simd4_load_f     gmx_simd_load_f

/* Set all SIMD4 register elements to value pointed to by m (unaligned) */
#    define gmx_simd4_load1_f    gmx_simd_load1_f

/* Set all SIMD4 register elements to the value r. */
#    define gmx_simd4_set1_f     gmx_simd_set1_f

/* Store the contents of the SIMD4 register pr to aligned memory m. */
#    define gmx_simd4_store_f    gmx_simd_store_f

/* Load from unaligned memory. Available with GMX_SIMD_HAVE_LOADU. */
#    define gmx_simd4_loadu_f    gmx_simd_loadu_f

/* Store to unaligned memory. Available with GMX_SIMD_HAVE_STOREU. */
#    define gmx_simd4_storeu_f   gmx_simd_storeu_f

/* Set all SIMD4 register elements to 0. */
#    define gmx_simd4_setzero_f  gmx_simd_setzero_f

/* Add two SIMD4 registers */
#    define gmx_simd4_add_f      gmx_simd_add_f

/* Subtract two SIMD4 registers */
#    define gmx_simd4_sub_f      gmx_simd_sub_f

/* Multiply two SIMD4 registers */
#    define gmx_simd4_mul_f      gmx_simd_mul_f

/* Fused-multiply-add for SIMD4. Result is a*b+c. */
#    define gmx_simd4_fmadd_f    gmx_simd_fmadd_f

/* Fused-multiply-subtract for SIMD4. Result is a*b-c. */
#    define gmx_simd4_fmsub_f    gmx_simd_fmsub_f

/* Fused-negated-multiply-add. Result is -a*b+c. */
#    define gmx_simd4_fnmadd_f   gmx_simd_fnmadd_f

/* Fused-negated-multiply-add. Result is -a*b-c. */
#    define gmx_simd4_fnmsub_f   gmx_simd_fnmsub_f

/* Bitwise and for two SIMD4 fp registers. */
#    define gmx_simd4_and_f      gmx_simd_and_f

/* Bitwise andnot for two SIMD4 fp registers. c=(~a) & b */
#    define gmx_simd4_andnot_f   gmx_simd_andnot_f

/* Bitwise or for two SIMD4 fp registers. */
#    define gmx_simd4_or_f       gmx_simd_or_f

/* Bitwise xor for two SIMD4 fp registers. */
#    define gmx_simd4_xor_f      gmx_simd_xor_f

/* Lookup for approximate 1/sqrt(x). Accuracy (bits) is GMX_SIMD_RSQRT_BITS */
#    define gmx_simd4_rsqrt_f    gmx_simd_rsqrt_f

/* Floating-point absolute value */
#    define gmx_simd4_fabs_f     gmx_simd_fabs_f

/* Floating-point negate */
#    define gmx_simd4_fneg_f     gmx_simd_fneg_f

/* Set each SIMD4 element to the largest from two registers. */
#    define gmx_simd4_max_f      gmx_simd_max_f

/* Set each SIMD4 element to the smallest from two registers. */
#    define gmx_simd4_min_f      gmx_simd_min_f

/* Round to nearest integer value (in floating-point format). */
#    define gmx_simd4_round_f    gmx_simd_round_f

/* Round to largest integral value (fp format) not greater than x. */
#    define gmx_simd4_trunc_f    gmx_simd_trunc_f

static gmx_inline float
gmx_simd4_dotproduct3_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    return a.r[0]*b.r[0]+a.r[1]*b.r[1]+a.r[2]*b.r[2];
}

/* SIMD4 register type to use for logical comparisons on floats */
#    define gmx_simd4_fbool_t   gmx_simd_fbool_t

/* Equality comparison of two single precision fp values */
#    define gmx_simd4_cmpeq_f   gmx_simd_cmpeq_f

/* Less-than comparison of two single precision fp values */
#    define gmx_simd4_cmplt_f   gmx_simd_cmplt_f

/* Less-than comparison of two single precision fp values */
#    define gmx_simd4_cmple_f   gmx_simd_cmple_f

/* Logical AND on float SIMD booleans */
#    define gmx_simd4_and_fb gmx_simd_and_fb

/* Logical OR on float SIMD booleans */
#    define gmx_simd4_or_fb gmx_simd_or_fb

/* Returns non-zero if any of the boolean in x is True. */
#    define gmx_simd4_anytrue_fb gmx_simd_anytrue_fb

/* Select from single precision fp SIMD register based on a boolean one. */
#    define gmx_simd4_blendzero_f gmx_simd_blendzero_f

/* Vector-blend instruction. */
#    define gmx_simd4_blendv_f  gmx_simd_blendv_f

#else /* GMX_SIMD_FLOAT_WIDTH!=4 */
#    undef GMX_SIMD4_HAVE_FLOAT
#endif


#if (GMX_SIMD_DOUBLE_WIDTH == 4)

/* Floating-point SIMD4 register type in double precision. */
#    define gmx_simd4_double_t   gmx_simd_double_t

/* Double precision version of gmx_simd_load_f() */
#    define gmx_simd4_load_d     gmx_simd_load_d

/* Double precision version of gmx_simd4_load1_f() */
#    define gmx_simd4_load1_d    gmx_simd_load1_d

/* Double precision version of gmx_simd_set1_f() */
#    define gmx_simd4_set1_d     gmx_simd_set1_d

/* Double precision version of gmx_simd_store_f() */
#    define gmx_simd4_store_d   gmx_simd_store_d

/* Double precision version of gmx_simd_loadu_f(). */
#    define gmx_simd4_loadu_d   gmx_simd_loadu_d

/* Double precision version of gmx_simd_storeu_f(). */
#    define gmx_simd4_storeu_d  gmx_simd_storeu_d

/* Double precision version of gmx_simd_setzero_f() */
#    define gmx_simd4_setzero_d gmx_simd_setzero_d

/* Double precision version of gmx_simd_add_f() */
#    define gmx_simd4_add_d     gmx_simd_add_d

/* Double precision version of gmx_simd_sub_f() */
#    define gmx_simd4_sub_d     gmx_simd_sub_d

/* Double precision version of gmx_simd_mul_f() */
#    define gmx_simd4_mul_d     gmx_simd_mul_d

/* Double precision version of gmx_simd_madd_f() */
#    define gmx_simd4_fmadd_d   gmx_simd_fmadd_d

/* Double precision version of gmx_simd_msub_f() */
#    define gmx_simd4_fmsub_d   gmx_simd_fmsub_d

/* Double precision version of gmx_simd_nmadd_f() */
#    define gmx_simd4_fnmadd_d  gmx_simd_fnmadd_d

/* Double precision version of gmx_simd_nmsub_f() */
#    define gmx_simd4_fnmsub_d  gmx_simd_fnmsub_d

/* Double precision version of gmx_simd_and_f() */
#    define gmx_simd4_and_d     gmx_simd_and_d

/* Double precision version of gmx_simd_andnot_f() */
#    define gmx_simd4_andnot_d  gmx_simd_andnot_d

/* Double precision version of gmx_simd_or_f() */
#    define gmx_simd4_or_d      gmx_simd_or_d

/* Double precision version of gmx_simd_xor_f() */
#    define gmx_simd4_xor_d     gmx_simd_xor_d

/* Double precision version of gmx_simd_rsqrt_f() */
#    define gmx_simd4_rsqrt_d   gmx_simd_rsqrt_d

/* Double precision version of gmx_simd_fabs_f() */
#    define gmx_simd4_fabs_d    gmx_simd_fabs_d

/* Double precision version of gmx_simd_fneg_f() */
#    define gmx_simd4_fneg_d    gmx_simd_fneg_d

/* Double precision version of gmx_simd_max_f() */
#    define gmx_simd4_max_d     gmx_simd_max_d

/* Double precision version of gmx_simd_min_f() */
#    define gmx_simd4_min_d     gmx_simd_min_d

/* Double precision version of gmx_simd_round_f() */
#    define gmx_simd4_round_d   gmx_simd_round_d

/* Double precision version of gmx_simd_trunc_f() */
#    define gmx_simd4_trunc_d   gmx_simd_trunc_d

static gmx_inline double
gmx_simd4_dotproduct3_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    return a.r[0]*b.r[0]+a.r[1]*b.r[1]+a.r[2]*b.r[2];
}

/* SIMD4 register type to use for logical comparisons on doubles */
#    define gmx_simd4_dbool_t   gmx_simd_dbool_t

/* Equality comparison of two double precision fp values */
#    define gmx_simd4_cmpeq_d   gmx_simd_cmpeq_d

/* Less-than comparison of two double precision fp values */
#    define gmx_simd4_cmplt_d   gmx_simd_cmplt_d

/* Less-than comparison of two double precision fp values */
#    define gmx_simd4_cmple_d   gmx_simd_cmple_d

/* Logical AND on double SIMD booleans */
#    define gmx_simd4_and_db gmx_simd_and_db

/* Logical OR on double SIMD booleans */
#    define gmx_simd4_or_db gmx_simd_or_db

/* Returns non-zero if any of the boolean in x is True. */
#    define gmx_simd4_anytrue_db gmx_simd_anytrue_db

/* Select from double precision fp SIMD register based on a boolean one. */
#    define gmx_simd4_blendzero_d gmx_simd_blendzero_d

/* Vector-blend instruction. */
#    define gmx_simd4_blendv_d  gmx_simd_blendv_d

#else /* GMX_SIMD4_DOUBLE_WIDTH!=4 */
#    undef GMX_SIMD4_HAVE_DOUBLE
#endif

#endif /* GMX_SIMD_IMPL_REFERENCE_H */
