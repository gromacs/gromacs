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

#ifndef GMX_SIMD_IMPL_REFERENCE_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_REFERENCE_SIMD_DOUBLE_H

/*! \libinternal \file
 *
 * \brief Reference implementation, SIMD double precision.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include <math.h>

#include "gromacs/utility/fatalerror.h"

#include "impl_reference_common.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* \name SIMD implementation data types
 * \{
 */

/*! \libinternal \brief Floating-point SIMD variable type in double precision.
 *
 * Supported with GMX_SIMD_HAVE_DOUBLE.
 */
typedef struct
{
    double r[GMX_SIMD_DOUBLE_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_double_t;

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from double.
 *
 * Available with GMX_SIMD_HAVE_DINT32.
 */
typedef struct
{
    gmx_int32_t i[GMX_SIMD_DINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_dint32_t;

/*! \libinternal \brief Boolean type for double precision SIMD data.
 *
 * Use the generic gmx_simd_bool_t
 * (for gmx_simd_real_t) instead, unless you really know what you are doing.
 */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_DOUBLE_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_dbool_t;

/*! \libinternal \brief Boolean type for integer datatypes corresponding to double SIMD.
 *
 * You should likely use gmx_simd_ibool_t (for gmx_simd_int32_t) instead,
 * unless you really know what you are doing.
 */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_DINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_dibool_t;

/*! \}
 *
 * \name SIMD implementation load/store operations for double precision floating point
 * \{
 */

/*! \brief Load \ref GMX_SIMD_DOUBLE_WIDTH numbers from aligned memory.
 *
 * \copydetails gmx_simd_load_f
 */
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

/*! \brief Set all SIMD variable elements to double pointed to by m (unaligned).
 *
 * \copydetails gmx_simd_load1_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load1_d(const double *m)
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

/*! \brief Set all SIMD double variable elements to the value r.
 *
 * \copydetails gmx_simd_set1_f
 */
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

/*! \brief Set all SIMD double variable elements to 0.0.
 *
 * \copydetails gmx_simd_setzero_f
 */
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

/*! \brief Store the contents of the SIMD double variable pr to aligned memory m.
 *
 * \copydetails gmx_simd_store_f
 */
static gmx_inline void
gmx_simd_store_d(double *m, gmx_simd_double_t a)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \brief Load SIMD double from unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \copydetails gmx_simd_loadu_f
 */
#define gmx_simd_loadu_d gmx_simd_load_d

/*! \brief Store SIMD double to unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \copydetails gmx_simd_storeu_f
 */
#define gmx_simd_storeu_d gmx_simd_store_d

/*! \}
 *
 * \name SIMD implementation load/store operations for integers (corresponding to double)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_load_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_load_di(const gmx_int32_t * m)
{
    gmx_simd_dint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/*! \brief Set SIMD from integer, width corresponds to \ref gmx_simd_double_t.
 *
 *  \copydetails gmx_simd_set1_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_set1_di(gmx_int32_t b)
{
    gmx_simd_dint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/*! \brief Set all SIMD variable elements to 0, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_setzero_fi
 */
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

/*! \brief Store aligned SIMD integer data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_store_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_store_di(gmx_int32_t * m, gmx_simd_dint32_t a)
{
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
    return a;
};

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_loadu_fi
 */
#define gmx_simd_loadu_di  gmx_simd_load_di

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_storeu_fi
 */
#define gmx_simd_storeu_di gmx_simd_store_di

/*! \brief Extract element with index i from \ref gmx_simd_dint32_t.
 *
 * \copydetails gmx_simd_extract_fi
 */
static gmx_inline gmx_int32_t
gmx_simd_extract_di(gmx_simd_dint32_t a, int index)
{
    return a.i[index];
}

/*! \}
 *
 * \name SIMD implementation double precision floating-point bitwise logical operations
 * \{
 */
/*! \brief Bitwise and for two SIMD double variables. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_and_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_and_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise andnot for SIMD double. c=(~a) & b. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_andnot_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_andnot_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = (~conv1.i) & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise or for SIMD double. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_or_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_or_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i | conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise xor for SIMD double. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_xor_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_xor_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
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
 * \name SIMD implementation double precision floating-point arithmetics
 * \{
 */
/*! \brief Add two double SIMD variables.
 *
 * \copydetails gmx_simd_add_f
 */
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

/*! \brief Add two float SIMD variables.
 *
 * \copydetails gmx_simd_sub_f
 */
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

/*! \brief Multiply two SIMD variables.
 *
 * \copydetails gmx_simd_mul_f
 */
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

/*! \brief Fused-multiply-add. Result is a*b+c.
 *
 * \copydetails gmx_simd_fmadd_f
 */
#define gmx_simd_fmadd_d(a, b, c) gmx_simd_add_d(gmx_simd_mul_d(a, b), c)

/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * \copydetails gmx_simd_fmsub_f
 */
#define gmx_simd_fmsub_d(a, b, c) gmx_simd_sub_d(gmx_simd_mul_d(a, b), c)

/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * \copydetails gmx_simd_fnmadd_f
 */
#define gmx_simd_fnmadd_d(a, b, c) gmx_simd_sub_d(c, gmx_simd_mul_d(a, b))

/*! \brief Fused-negated-multiply-add. Result is -a*b-c.
 *
 * \copydetails gmx_simd_fnmsub_f
 */
#define gmx_simd_fnmsub_d(a, b, c) gmx_simd_sub_d(gmx_simd_setzero_d(), gmx_simd_fmadd_d(a, b, c))

/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * \copydetails gmx_simd_rsqrt_f
 */
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

/*! \brief 1.0/x lookup.
 *
 * \copydetails gmx_simd_rcp_f
 */
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

/*! \brief SIMD Floating-point fabs().
 *
 * \copydetails gmx_simd_fabs_f
 */
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

/*! \brief SIMD floating-point negate.
 *
 * \copydetails gmx_simd_fneg_f
 */
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

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * \copydetails gmx_simd_max_f
 */
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

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * \copydetails gmx_simd_min_f
 */
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

/*! \brief Round to nearest integer value (in double floating-point format).
 *
 * \copydetails gmx_simd_round_f
 */
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

/*! \brief Truncate SIMD, i.e. round towards zero - common hardware instruction.
 *
 * \copydetails gmx_simd_trunc_f
 */
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

/*! \brief Fraction of the SIMD floating point number.
 *
 * \copydetails gmx_simd_fraction_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_fraction_d(gmx_simd_double_t a)
{
    return gmx_simd_sub_d(a, gmx_simd_trunc_d(a));
}


/*! \brief Extract (integer) exponent from double precision SIMD.
 *
 * \copydetails gmx_simd_get_exponent_f
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

/*! \brief Get SIMD doublemantissa.
 *
 * \copydetails gmx_simd_get_mantissa_f
 */
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

/*! \brief Set (integer) exponent from double precision floating-point SIMD.
 *
 * \copydetails gmx_simd_set_exponent_f
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
        /* Critical to use same algorithm as for gmx_simd_round_d() */
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

/*! \}
 *
 * \name SIMD implementation double precision floating-point comparison, boolean, selection.
 * \{
 */
/*! \brief SIMD a==b for double SIMD.
 *
 * \copydetails gmx_simd_cmpeq_f
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

/*! \brief SIMD a<b for double SIMD.
 *
 * \copydetails gmx_simd_cmplt_f
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

/*! \brief SIMD a<=b for double SIMD.
 *
 * \copydetails gmx_simd_cmple_f
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


/*! \brief Logical \a and on double precision SIMD booleans.
 *
 * \copydetails gmx_simd_and_fb
 */
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

/*! \brief Logical \a or on double precision SIMD booleans.
 *
 * \copydetails gmx_simd_or_fb
 */
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


/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * \copydetails gmx_simd_anytrue_fb
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


/*! \brief Select from double SIMD variable where boolean is true.
 *
 * \copydetails gmx_simd_blendzero_f
 */
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

/*! \brief Select from double SIMD variable where boolean is false.
 *
 * \copydetails gmx_simd_blendnotzero_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_blendnotzero_d(gmx_simd_double_t a, gmx_simd_dbool_t sel)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? 0.0 : a.r[i];
    }
    return c;
}

/*! \brief Vector-blend double SIMD selection.
 *
 * \copydetails gmx_simd_blendv_f
 */
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

/*! \brief Return sum of all elements in SIMD double variable.
 *
 * \copydetails gmx_simd_reduce_f
 *
 */
static gmx_inline double
gmx_simd_reduce_d(gmx_simd_double_t a)
{
    double    sum = 0.0;
    int       i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        sum += a.r[i];
    }
    return sum;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to double) bitwise logical operations
 * \{
 */

/*! \brief SIMD integer shift left, based on immediate value.
 *
 * \copydetails gmx_simd_slli_fi
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

/*! \brief SIMD integer shift right, based on immediate value.
 *
 * \copydetails gmx_simd_srli_fi
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

/*! \brief Integer bitwise and for SIMD variables.
 *
 * \copydetails gmx_simd_and_fi
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

/*! \brief Integer bitwise not-and for SIMD variables.
 *
 * \copydetails gmx_simd_andnot_fi
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

/*! \brief Integer bitwise or for SIMD variables.
 *
 * \copydetails gmx_simd_or_fi
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

/*! \brief Integer bitwise xor for SIMD variables.
 *
 * \copydetails gmx_simd_xor_fi
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

/*! \}
 *
 * \name SIMD implementation integer (corresponding to double) arithmetics
 * \{
 */
/*! \brief Add SIMD integers, corresponding to double precision.
 *
 * \copydetails gmx_simd_add_fi
 */
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

/*! \brief Subtract SIMD integers, corresponding to double precision.
 *
 * \copydetails gmx_simd_sub_fi
 */
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

/*! \brief Multiply SIMD integers, corresponding to double precision.
 *
 * \copydetails gmx_simd_mul_fi
 */
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

/*! \}
 *
 * \name SIMD implementation integer (corresponding to double) comparisons, boolean selection
 * \{
 */

/*! \brief Equality comparison of two ints corresponding to double SIMD data.
 *
 * \copydetails gmx_simd_cmpeq_fi
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

/*! \brief Less-than comparison of two ints corresponding to double SIMD data.
 *
 * \copydetails gmx_simd_cmplt_fi
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

/*! \brief Logical AND on gmx_simd_dibool_t.
 *
 * \copydetails gmx_simd_and_fib
 */
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

/*! \brief Logical OR on gmx_simd_dibool_t.
 *
 * \copydetails gmx_simd_or_fib
 */
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

/*! \brief Returns non-zero if any of the double-int SIMD booleans in x is True, otherwise 0.
 *
 * \copydetails gmx_simd_anytrue_fib
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

/*! \brief Select from SIMD ints (corresponding to double) where boolean is true.
 *
 * \copydetails gmx_simd_blendzero_fi
 */
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

/*! \brief Select from SIMD ints (corresponding to double) where boolean is false.
 *
 * \copydetails gmx_simd_blendnotzero_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_blendnotzero_di(gmx_simd_dint32_t a, gmx_simd_dibool_t sel)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? 0.0 : a.i[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection for double-int SIMD.
 *
 * \copydetails gmx_simd_blendv_fi
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

/*! \}
 *
 * \name SIMD implementation conversion operations
 * \{
 */

/*! \brief Round double precision floating point to integer.
 *
 * \copydetails gmx_simd_cvt_f2i
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

/*! \brief Truncate double precision floating point to integer.
 *
 * \copydetails gmx_simd_cvtt_f2i
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

/*! \brief Convert integer to double precision floating-point.
 *
 * \copydetails gmx_simd_cvt_i2f
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

/*! \brief Convert from double boolean to corresponding integer boolean.
 *
 * \copydetails gmx_simd_cvt_fb2fib
 */
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

/*! \brief Convert from integer boolean (corresponding to double) to double boolean.
 *
 * \copydetails gmx_simd_cvt_fib2fb
 */
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

/*! \brief Convert SIMD float to double.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is identical to
 * \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param f Single-precision SIMD variable
 * \return Double-precision SIMD variable of the same width
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
    d.r[0] = f.r[0];
#endif
    return d;
}

/*! \brief Convert SIMD double to float.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is identical to
 * \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param d Double-precision SIMD variable
 * \return Single-precision SIMD variable of the same width
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
    f.r[0] = d.r[0];
#endif
    return f;
}

/*! \brief Convert SIMD float to double.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is twice as large
 * as \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param f Single-precision SIMD variable
 * \param[out] d0 Double-precision SIMD variable, first half of values from f.
 * \param[out] d1 Double-precision SIMD variable, second half of values from f.
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

/*! \brief Convert SIMD double to float.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is twice as large
 * as \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param d0 Double-precision SIMD variable, first half of values to put in f.
 * \param d1 Double-precision SIMD variable, second half of values to put in f.
 * \return Single-precision SIMD variable with all values.
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

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD_DOUBLE_H */
