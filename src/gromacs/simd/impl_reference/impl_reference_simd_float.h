/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H

/*! \libinternal \file
 *
 * \brief Reference implementation, SIMD single precision.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include <math.h>

#include "impl_reference_common.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* \name SIMD implementation data types
 * \{
 */
/*! \libinternal \brief Float SIMD variable. Supported if GMX_SIMD_HAVE_FLOAT.
 */
typedef struct
{
    float r[GMX_SIMD_FLOAT_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_float_t;

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from float.
 *
 * This is also the widest integer SIMD type.
 */
typedef struct
{
    gmx_int32_t i[GMX_SIMD_FINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_fint32_t;

/*! \libinternal \brief Boolean type for float SIMD data.
 *
 * You should likely use gmx_simd_bool_t
 * (for gmx_simd_real_t) instead, unless you really know what you are doing.
 */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_FLOAT_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_fbool_t;

/*! \libinternal \brief Boolean type for integer datatypes corresponding to float SIMD. */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_FINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_fibool_t;

/*! \}
 *
 * \name SIMD implementation load/store operations for single precision floating point
 * \{
 */

/*! \brief Load \ref GMX_SIMD_FLOAT_WIDTH numbers from aligned memory.
 *
 * \param m Pointer to memory aligned to the SIMD width.
 * \return SIMD variable with data loaded.
 */
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

/*! \brief Set all SIMD variable elements to float pointed to by m (unaligned).
 *
 * \param m Pointer to single value in memory.
 * \return SIMD variable with all elements set to *m.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_load1_f(const float *m)
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

/*! \brief Set all SIMD float variable elements to the value r.
 *
 *  \param r floating-point constant
 *  \return SIMD variable with all elements set to r.
 */
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

/*! \brief Set all SIMD float variable elements to 0.0f.
 *
 *  \return The value 0.0 in all elements of a SIMD variable.
 */
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

/*! \brief Store the contents of the SIMD float variable pr to aligned memory m.
 *
 * \param[out] m Pointer to memory, aligned to SIMD width.
 * \param a SIMD variable to store
 */
static gmx_inline void
gmx_simd_store_f(float *m, gmx_simd_float_t a)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \brief Load SIMD float from unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \param m Pointer to memory, no alignment requirement.
 * \return SIMD variable with data loaded.
 */
#define gmx_simd_loadu_f gmx_simd_load_f

/*! \brief Store SIMD float to unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a SIMD variable to store.
 */
#define gmx_simd_storeu_f gmx_simd_store_f

/*! \}
 *
 * \name SIMD implementation load/store operations for integers (corresponding to float)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_load_i.
 *
 * \param m Pointer to memory, aligned to integer SIMD width.
 * \return SIMD integer variable.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_load_fi(const gmx_int32_t * m)
{
    gmx_simd_fint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/*! \brief Set SIMD from integer, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_set1_i.
 *
 *  \param b integer value to set variable to.
 *  \return SIMD variable with all elements set to b.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_set1_fi(gmx_int32_t b)
{
    gmx_simd_fint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/*! \brief Set all SIMD variable elements to 0, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_setzero_i.
 *
 * \return SIMD integer variable with all bits set to zero.
 */
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

/*! \brief Store aligned SIMD integer data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_store_i.
 *
 * \param m Memory aligned to integer SIMD width.
 * \param a SIMD variable to store.
 */
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

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_loadu_i.
 *
 * Supported with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \param m Pointer to memory, no alignment requirements.
 * \return SIMD integer variable.
 */
#define gmx_simd_loadu_fi  gmx_simd_load_fi

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_storeu_i.
 *
 * Supported with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \param m Memory pointer, no alignment requirements.
 * \param a SIMD variable to store.
 */
#define gmx_simd_storeu_fi gmx_simd_store_fi

/*! \brief Extract element with index i from \ref gmx_simd_fint32_t.
 *
 * You should typically call the real-precision \ref gmx_simd_extract_i.
 *
 * Available with \ref GMX_SIMD_HAVE_FINT32_EXTRACT.
 *
 * \param a SIMD variable
 * \param index Position to extract integer from
 * \return Single integer from position index in SIMD variable.
 */
static gmx_inline gmx_int32_t
gmx_simd_extract_fi(gmx_simd_fint32_t a, int index)
{
    return a.i[index];
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point bitwise logical operations
 * \{
 */

/*! \brief Bitwise and for two SIMD float variables. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_and_r.
 *
 * \param a data1
 * \param b data2
 * \return data1 & data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_and_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise andnot for SIMD float. c=(~a) & b. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_andnot_r.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_andnot_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = (~conv1.i) & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise or for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_or_r.
 *
 * \param a data1
 * \param b data2
 * \return data1 | data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_or_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i | conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise xor for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_r.
 *
 * \param a data1
 * \param b data2
 * \return data1 ^ data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_xor_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i ^ conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point arithmetics
 * \{
 */
/*! \brief Add two float SIMD variables.
 *
 * You should typically call the real-precision \ref gmx_simd_add_r.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
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

/*! \brief Subtract two SIMD variables.
 *
 * You should typically call the real-precision \ref gmx_simd_sub_r.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
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

/*! \brief Multiply two SIMD variables.
 *
 * You should typically call the real-precision \ref gmx_simd_mul_r.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 */
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

/*! \brief Fused-multiply-add. Result is a*b+c.
 *
 * You should typically call the real-precision \ref gmx_simd_fmadd_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return a*b+c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fmadd_f(a, b, c) gmx_simd_add_f(gmx_simd_mul_f(a, b), c)


/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * You should typically call the real-precision \ref gmx_simd_fmsub_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return a*b-c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fmsub_f(a, b, c) gmx_simd_sub_f(gmx_simd_mul_f(a, b), c)


/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * You should typically call the real-precision \ref gmx_simd_fnmadd_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return -a*b+c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fnmadd_f(a, b, c) gmx_simd_sub_f(c, gmx_simd_mul_f(a, b))


/*! \brief Fused-negated-multiply-sub. Result is -a*b-c.
 *
 * You should typically call the real-precision \ref gmx_simd_fnmsub_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return -a*b-c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fnmsub_f(a, b, c) gmx_simd_sub_f(gmx_simd_setzero_f(), gmx_simd_fmadd_f(a, b, c))

/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * You should typically call the real-precision \ref gmx_simd_rsqrt_r.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 */
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

/*! \brief SIMD 1.0/x lookup.
 *
 * You should typically call the real-precision \ref gmx_simd_rcp_r.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the reciprocal in simd_math.h.
 *
 * \param x Argument, x!=0
 * \return Approximation of 1/x, accuracy is \ref GMX_SIMD_RCP_BITS.
 */
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

/*! \brief SIMD Floating-point fabs().
 *
 * You should typically call the real-precision \ref gmx_simd_fabs_r.
 *
 * \param a any floating point values
 * \return fabs(a) for each element.
 */
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

/*! \brief SIMD floating-point negate.
 *
 * You should typically call the real-precision \ref gmx_simd_fneg_r.
 *
 * \param a Any floating-point value
 * \return -a
 */
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

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * You should typically call the real-precision \ref gmx_simd_max_r.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 */
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

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * You should typically call the real-precision \ref gmx_simd_min_r.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return min(a,b) for each element.
 */
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

/*! \brief Round to nearest integer value (in floating-point format).
 *
 * You should typically call the real-precision \ref gmx_simd_round_r.
 *
 * \param a Any floating-point value
 * \return The nearest integer, represented in floating-point format.
 *
 * \note The reference implementation rounds exact half-way cases
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

/*! \brief Truncate SIMD, i.e. round towards zero - common hardware instruction.
 *
 * You should typically call the real-precision \ref gmx_simd_trunc_r.
 *
 * \param a Any floating-point value
 * \return Integer rounded towards zero, represented in floating-point format.
 *
 * \note This is truncation towards zero, not floor(). The reason for this
 * is that truncation is virtually always present as a dedicated hardware
 * instruction, but floor() frequently isn't.
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


/*! \brief Fraction of the SIMD floating point number.
 *
 * You should typically call the real-precision \ref gmx_simd_fraction_r.
 *
 * \param a Any floating-point value
 * \return a-trunc(r)
 *
 * To maximize compatibility, we use the same definition of fractions as used
 * e.g. for the AMD64 hardware instructions. This relies on truncation towards
 * zero for the integer part, and the remaining fraction can thus be either
 * positive or negative. As an example, -1.42 would return the fraction -0.42.
 *
 * Hardware support with \ref GMX_SIMD_HAVE_FRACTION, otherwise emulated.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_fraction_f(gmx_simd_float_t a)
{
    return gmx_simd_sub_f(a, gmx_simd_trunc_f(a));
}

/*! \brief Extract (integer) exponent from single precision SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_get_exponent_r.
 *
 * \param a Any floating-point value
 * \return Exponent value, represented in floating-point format.
 *
 * The IEEE754 exponent field is selected, the bias removed, and it is converted to
 * a normal floating-point SIMD.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f(gmx_simd_float_t a)
{
    /* Mask with ones for the exponent field of single precision fp */
    const gmx_int32_t  expmask = 0x7f800000;
    gmx_simd_float_t   b;
    int                i;
    union
    {
        float        f;
        gmx_int32_t  i;
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

/*! \brief Get SIMD mantissa.
 *
 * You should typically call the real-precision \ref gmx_simd_get_mantissa_r.
 *
 * \param a Any floating-point value
 * \return Mantissa, represented in floating-point format.
 *
 * The mantissa field is selected, and a new neutral exponent created.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_get_mantissa_f(gmx_simd_float_t a)
{
    const gmx_int32_t  mantmask = 0x007fffff;
    const gmx_int32_t  one      = 0x3f800000;
    gmx_simd_float_t   b;
    int                i;
    union
    {
        float        f;
        gmx_int32_t  i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        /* remove current exponent, add a biased exponent for 1.0 (i.e., 2^0=1) */
        conv.i = (conv.i & (mantmask)) | one;
        b.r[i] = conv.f;
    }
    return b;
}

/*! \brief Set (integer) exponent from single precision floating-point SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_set_exponent_r.
 *
 * \param a A floating point value that will not overflow as 2^a.
 * \return 2^(round(a)).
 *
 * The input is \a rounded to the nearest integer, the exponent bias is added
 * to this integer, and the bits are shifted to the IEEE754 exponent part of the number.
 *
 * \note The argument will be \a rounded to nearest integer since that is what
 * we need for the exponential functions, and this integer x will be set as the
 * exponent so the new floating-point number will be 2^x.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f(gmx_simd_float_t a)
{
    gmx_simd_float_t   b;
    gmx_int32_t        iexp;
    int                i;
    union
    {
        float        f;
        gmx_int32_t  i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        /* Critical to use same algorithm as for gmx_simd_round_f() */
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

/*! \}
 *
 * \name SIMD implementation single precision floating-point comparisons, boolean, selection.
 * \{
 */
/*! \brief SIMD a==b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmpeq_r.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a==b.
 *
 * Beware that exact floating-point comparisons are difficult.
 */
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

/*! \brief SIMD a<b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmplt_r.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<b.
 */
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

/*! \brief SIMD a<=b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmple_r.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<=b.
 */
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

/*! \brief Logical \a and on single precision SIMD booleans.
 *
 * You should typically call the real-precision \ref gmx_simd_and_r.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a \& b are true.
 *
 * \note This is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 *
 * \sa gmx_simd_and_ib
 */
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

/*! \brief Logical \a or on single precision SIMD booleans.
 *
 * You should typically call the real-precision \ref gmx_simd_or_r.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a or b is true.
 *
 * Note that this is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 *
 * \sa gmx_simd_or_ib
 */
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

/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * You should typically call the real-precision \ref gmx_simd_anytrue_b.
 *
 * \param a Logical variable.
 * \return non-zero if any element in a is true, otherwise 0.
 *
 * The actual return value for truth will depend on the architecture,
 * so any non-zero value is considered truth.
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

/*! \brief Select from single precision SIMD variable where boolean is true.
 *
 * You should typically call the real-precision \ref gmx_simd_blendzero_r.
 *
 * \param a Floating-point variable to select from
 * \param sel Boolean selector
 * \return  For each element, a is selected for true, 0 for false.
 */
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

/*! \brief Select from single precision SIMD variable where boolean is false.
 *
 * You should typically call the real-precision \ref gmx_simd_blendnotzero_r.
 *
 * \param a Floating-point variable to select from
 * \param sel Boolean selector
 * \return  For each element, a is selected for false, 0 for true (sic).
 */
static gmx_inline gmx_simd_float_t
gmx_simd_blendnotzero_f(gmx_simd_float_t a, gmx_simd_fbool_t sel)
{
    gmx_simd_float_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? 0.0 : a.r[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection.
 *
 * You should typically call the real-precision \ref gmx_simd_blendv_r.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
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

/*! \brief Return sum of all elements in SIMD float variable.
 *
 * You should typically call the real-precision \ref gmx_simd_reduce_r.
 *
 * \param a SIMD variable to reduce/sum.
 * \return The sum of all elements in the argument variable.
 *
 */
static gmx_inline float
gmx_simd_reduce_f(gmx_simd_float_t a)
{
    float     sum = 0.0;
    int       i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        sum += a.r[i];
    }
    return sum;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) bitwise logical operations
 * \{
 */

/*! \brief SIMD integer shift left logical, based on immediate value.
 *
 * You should typically call the real-precision \ref gmx_simd_slli_i.
 *
 *  Logical shift. Each element is shifted (independently) up to 32 positions
 *  left, while zeros are shifted in from the right. Only available if
 * \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single) or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
 *  (double) is 1.
 *
 * \param a integer data to shift
 * \param n number of positions to shift left. n<=32.
 * \return shifted values
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

/*! \brief SIMD integer shift right logical, based on immediate value.
 *
 * You should typically call the real-precision \ref gmx_simd_srli_i.
 *
 *  Logical shift. Each element is shifted (independently) up to 32 positions
 *  right, while zeros are shifted in from the left. Only available if
 * \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single) or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
 *  (double) is 1.
 *
 * \param a integer data to shift
 * \param n number of positions to shift right. n<=32.
 * \return shifted values
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

/*! \brief Integer SIMD bitwise and.
 *
 * You should typically call the real-precision \ref gmx_simd_and_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \note You can \a not use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx_simd_blendzero_i instead.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \& b (bitwise and)
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

/*! \brief Integer SIMD bitwise not-and.
 *
 * You should typically call the real-precision \ref gmx_simd_andnot_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * Note that you can NOT use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx_simd_blendnotzero_i instead.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return (~a) \& b (bitwise andnot)
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

/*! \brief Integer SIMD bitwise or.
 *
 * You should typically call the real-precision \ref gmx_simd_or_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \| b (bitwise or)
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

/*! \brief Integer SIMD bitwise xor.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a ^ b (bitwise xor)
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

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) arithmetics
 * \{
 */
/*! \brief Add SIMD integers.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
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

/*! \brief Subtract SIMD integers.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
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

/*! \brief Multiply SIMD integers.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 *
 * \note Only the low 32 bits are retained, so this can overflow.
 */
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

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) comparisons, boolean, selection
 * \{
 */

/*! \brief Equality comparison of two integers corresponding to float values.
 *
 * You should typically call the real-precision \ref gmx_simd_cmpeq_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a==b
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

/*! \brief Less-than comparison of two SIMD integers corresponding to float values.
 *
 * You should typically call the real-precision \ref gmx_simd_cmplt_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a<b
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

/*! \brief Logical AND on gmx_simd_fibool_t.
 *
 * You should typically call the real-precision \ref gmx_simd_and_ib.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
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

/*! \brief Logical OR on gmx_simd_fibool_t.
 *
 * You should typically call the real-precision \ref gmx_simd_or_ib.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
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

/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * You should typically call the real-precision \ref gmx_simd_anytrue_ib.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * The actual return value for "any true" will depend on the architecture.
 * Any non-zero value should be considered truth.
 *
 * \param a SIMD boolean
 * \return Nonzero integer if any of the elements in a is true, otherwise 0.
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

/*! \brief Select from \ref gmx_simd_fint32_t variable where boolean is true.
 *
 * You should typically call the real-precision \ref gmx_simd_blendzero_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer to select from
 * \param sel Boolean selector
 * \return Elements from a where sel is true, 0 otherwise.
 */
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

/*! \brief Select from \ref gmx_simd_fint32_t variable where boolean is false.
 *
 * You should typically call the real-precision \ref gmx_simd_blendnotzero_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer to select from
 * \param sel Boolean selector
 * \return Elements from a where sel is false, 0 otherwise (sic).
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_blendnotzero_fi(gmx_simd_fint32_t a, gmx_simd_fibool_t sel)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? 0.0 : a.i[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection.
 *
 * You should typically call the real-precision \ref gmx_simd_blendv_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
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

/*! \}
 *
 * \name SIMD implementation conversion operations
 * \{
 */

/*! \brief Round single precision floating point to integer.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_r2i.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, rounded to nearest integer.
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

/*! \brief Truncate single precision floating point to integer.
 *
 * You should typically call the real-precision \ref gmx_simd_cvtt_r2i.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, truncated towards zero.
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

/*! \brief Convert integer to single precision floating-point.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_i2r.
 *
 * \param a SIMD integer
 * \return SIMD floating-pint
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

/*! \brief Convert from float boolean to corresponding integer boolean.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_b2ib.
 *
 * \param a Boolean corresponding to SIMD floating-point
 * \return Boolean that can be applied to SIMD integer operations.
 */
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

/*! \brief Convert from integer boolean (corresponding to float) to float boolean.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_ib2b.
 *
 * \param a Boolean corresponding to SIMD integer
 * \return Boolean that can be applied to SIMD floating-point.
 */
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

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H */
