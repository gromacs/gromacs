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

#include <cmath>
#include <cstdint>

#include "gromacs/utility/fatalerror.h"

#include "impl_reference_common.h"
#include "impl_reference_simd_float.h"

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
SimdDouble;

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from double.
 *
 * Available with GMX_SIMD_HAVE_DINT32.
 */
typedef struct
{
    std::int32_t i[GMX_SIMD_DINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
SimdDInt32;

/*! \libinternal \brief Boolean type for double precision SIMD data.
 *
 * Use the generic SimdBool
 * (for SimdReal) instead, unless you really know what you are doing.
 */
typedef struct
{
    std::int32_t b[GMX_SIMD_DOUBLE_WIDTH]; /**< Implementation dependent. Don't touch. */
}
SimdDBool;

/*! \libinternal \brief Boolean type for integer datatypes corresponding to double SIMD.
 *
 * You should likely use SimdIBool (for SimdInt32) instead,
 * unless you really know what you are doing.
 */
typedef struct
{
    std::int32_t b[GMX_SIMD_DINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
SimdDIBool;

/*! \}
 *
 * \name SIMD implementation load/store operations for double precision floating point
 * \{
 */

/*! \brief Load \ref GMX_SIMD_DOUBLE_WIDTH numbers from aligned memory.
 *
 * \copydetails simdLoadF
 */
static inline SimdDouble
simdLoadD(const double *m)
{
    SimdDouble         a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/*! \brief Set all SIMD variable elements to double pointed to by m (unaligned).
 *
 * \copydetails simdLoad1F
 */
static inline SimdDouble
simdLoad1D(const double *m)
{
    SimdDouble         a;
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
 * \copydetails simdSet1F
 */
static inline SimdDouble
simdSet1D(double r)
{
    SimdDouble         a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}

/*! \brief Set all SIMD double variable elements to 0.0.
 *
 * \copydetails simdSetZeroF
 */
static inline SimdDouble
simdSetZeroD()
{
    SimdDouble         a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }
    return a;
}

/*! \brief Store the contents of the SIMD double variable pr to aligned memory m.
 *
 * \copydetails simdStoreF
 */
static inline void
simdStoreD(double *m, SimdDouble a)
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
 * \copydetails simdLoadUF
 */
#define simdLoadUD simdLoadD

/*! \brief Store SIMD double to unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \copydetails simdStoreUF
 */
#define simdStoreUD simdStoreD

/*! \}
 *
 * \name SIMD implementation load/store operations for integers (corresponding to double)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref SimdDouble.
 *
 * \copydetails simdLoadFI
 */
static inline SimdDInt32
simdLoadDI(const std::int32_t * m)
{
    SimdDInt32         a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/*! \brief Set SIMD from integer, width corresponds to \ref SimdDouble.
 *
 *  \copydetails simdSet1FI
 */
static inline SimdDInt32
simdSet1DI(std::int32_t b)
{
    SimdDInt32         a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/*! \brief Set all SIMD variable elements to 0, width corresponds to \ref SimdDouble.
 *
 * \copydetails simdSetZeroFI
 */
static inline SimdDInt32
simdSetZeroDI()
{
    SimdDInt32         a;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = 0;
    }
    return a;
}

/*! \brief Store aligned SIMD integer data, width corresponds to \ref SimdDouble.
 *
 * \copydetails simdStoreFI
 */
static inline SimdDInt32
simdStoreDI(std::int32_t * m, SimdDInt32 a)
{
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
    return a;
};

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref SimdDouble.
 *
 * \copydetails simdLoadUFI
 */
#define simdLoadUDI  simdLoadDI

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref SimdDouble.
 *
 * \copydetails simdStoreUFI
 */
#define simdStoreUDI simdStoreDI

/*! \brief Extract element with index i from \ref SimdDInt32.
 *
 * \copydetails simdExtractFI
 */
static inline std::int32_t
simdExtractDI(SimdDInt32 a, int index)
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
 * \copydetails simdAndF
 */
static inline SimdDouble
simdAndD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;
    union
    {
        double        r;
        std::int64_t  i;
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
 * \copydetails simdAndNotF
 */
static inline SimdDouble
simdAndNotD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;
    union
    {
        double        r;
        std::int64_t  i;
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
 * \copydetails simdOrF
 */
static inline SimdDouble
simdOrD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;
    union
    {
        double        r;
        std::int64_t  i;
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
 * \copydetails simdXorF
 */
static inline SimdDouble
simdXorD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;
    union
    {
        double        r;
        std::int64_t  i;
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
 * \copydetails simdAddF
 */
static inline SimdDouble
simdAddD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }
    return c;
}

/*! \brief Add two float SIMD variables.
 *
 * \copydetails simdSubF
 */
static inline SimdDouble
simdSubD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }
    return c;
}

/*! \brief Multiply two SIMD variables.
 *
 * \copydetails simdMulF
 */
static inline SimdDouble
simdMulD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i]*b.r[i];
    }
    return c;
}

/*! \brief Fused-multiply-add. Result is a*b+c.
 *
 * \copydetails simdFmaddF
 */
#define simdFmaddD(a, b, c) simdAddD(simdMulD(a, b), c)

/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * \copydetails simdFmsubF
 */
#define simdFmsubD(a, b, c) simdSubD(simdMulD(a, b), c)

/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * \copydetails simdFnmaddF
 */
#define simdFnmaddD(a, b, c) simdSubD(c, simdMulD(a, b))

/*! \brief Fused-negated-multiply-add. Result is -a*b-c.
 *
 * \copydetails simdFnmsubF
 */
#define simdFnmsubD(a, b, c) simdSubD(simdSetZeroD(), simdFmaddD(a, b, c))

/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * \copydetails simdRsqrtF
 */
static inline SimdDouble
simdRsqrtD(SimdDouble x)
{
    SimdDouble         b;
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
 * \copydetails simdRcpF
 */
static inline SimdDouble
simdRcpD(SimdDouble x)
{
    SimdDouble         b;
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
 * \copydetails simdAbsF
 */
static inline SimdDouble
simdAbsD(SimdDouble a)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = fabs(a.r[i]);
    }
    return c;
}

/*! \brief SIMD floating-point negate.
 *
 * \copydetails simdNegF
 */
static inline SimdDouble
simdNegD(SimdDouble a)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * \copydetails simdMaxF
 */
static inline SimdDouble
simdMaxD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * \copydetails simdMinF
 */
static inline SimdDouble
simdMinD(SimdDouble a, SimdDouble b)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = (a.r[i] <= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/*! \brief Round to nearest integer value (in double floating-point format).
 *
 * \copydetails simdRoundF
 */
static inline SimdDouble
simdRoundD(SimdDouble a)
{
    SimdDouble         b;
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
 * \copydetails simdTruncF
 */
static inline SimdDouble
simdTruncD(SimdDouble a)
{
    SimdDouble         b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.r[i] = trunc(a.r[i]);
    }
    return b;
}

/*! \brief Fraction of the SIMD floating point number.
 *
 * \copydetails simdFractionF
 */
static inline SimdDouble
simdFractionD(SimdDouble a)
{
    return simdSubD(a, simdTruncD(a));
}


/*! \brief Extract (integer) exponent from double precision SIMD.
 *
 * \copydetails simdGetExponentF
 */
static inline SimdDouble
simdGetExponentD(SimdDouble a)
{
    /* Mask with ones for the exponent field of double precision fp */
    const std::int64_t      expmask = 0x7ff0000000000000LL;
    SimdDouble              b;
    int                     i;
    union
    {
        double              d;
        std::int64_t        i;
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
 * \copydetails simdGetMantissaF
 */
static inline SimdDouble
simdGetMantissaD(SimdDouble a)
{
    const std::int64_t      mantmask = 0x000fffffffffffffLL;
    const std::int64_t      one      = 0x3ff0000000000000LL;
    SimdDouble              b;
    int                     i;
    union
    {
        double           d;
        std::int64_t     i;
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
 * \copydetails simdSetExponentF
 */
static inline SimdDouble
simdSetExponentD(SimdDouble a)
{
    SimdDouble              b;
    int                     i;
    std::int64_t            iexp;
    union
    {
        double           d;
        std::int64_t     i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        /* Critical to use same algorithm as for simdRoundD() */
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
 * \copydetails simdCmpEqF
 */
static inline SimdDBool
simdCmpEqD(SimdDouble a, SimdDouble b)
{
    SimdDBool         c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/*! \brief SIMD a<b for double SIMD.
 *
 * \copydetails simdCmpLtF
 */
static inline SimdDBool
simdCmpLtD(SimdDouble a, SimdDouble b)
{
    SimdDBool         c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}

/*! \brief SIMD a<=b for double SIMD.
 *
 * \copydetails simdCmpLeF
 */
static inline SimdDBool
simdCmpLeD(SimdDouble a, SimdDouble b)
{
    SimdDBool         c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}


/*! \brief Logical \a and on double precision SIMD booleans.
 *
 * \copydetails simdAndFB
 */
static inline SimdDBool
simdAndDB(SimdDBool a, SimdDBool b)
{
    SimdDBool         c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical \a or on double precision SIMD booleans.
 *
 * \copydetails simdOrFB
 */
static inline SimdDBool
simdOrDB(SimdDBool a, SimdDBool b)
{
    SimdDBool         c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}


/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * \copydetails simdAnyTrueFB
 */
static inline int
simdAnyTrueDB(SimdDBool a)
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
 * \copydetails simdMaskF
 */
static inline SimdDouble
simdMaskD(SimdDouble a, SimdDBool mask)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = mask.b[i] ? a.r[i] : 0.0;
    }
    return c;
}

/*! \brief Select from double SIMD variable where boolean is false.
 *
 * \copydetails simdMaskNotF
 */
static inline SimdDouble
simdMaskNotD(SimdDouble a, SimdDBool mask)
{
    SimdDouble         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = mask.b[i] ? 0.0 : a.r[i];
    }
    return c;
}

/*! \brief Vector-blend double SIMD selection.
 *
 * \copydetails simdBlendF
 */
static inline SimdDouble
simdBlendD(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    SimdDouble         d;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}

/*! \brief Return sum of all elements in SIMD double variable.
 *
 * \copydetails simdReduceF
 *
 */
static inline double
simdReduceD(SimdDouble a)
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
 * \copydetails simdSlliFI
 */
static inline SimdDInt32
simdSlliDI(SimdDInt32 a, int n)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] << n;
    }
    return c;
}

/*! \brief SIMD integer shift right, based on immediate value.
 *
 * \copydetails simdSrliFI
 */
static inline SimdDInt32
simdSrliDI(SimdDInt32 a, int n)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] >> n;
    }
    return c;
}

/*! \brief Integer bitwise and for SIMD variables.
 *
 * \copydetails simdAndFI
 */
static inline SimdDInt32
simdAndDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] & b.i[i];
    }
    return c;
}

/*! \brief Integer bitwise not-and for SIMD variables.
 *
 * \copydetails simdAndNotFI
 */
static inline SimdDInt32
simdAndNotDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = (~a.i[i]) & b.i[i];
    }
    return c;
}

/*! \brief Integer bitwise or for SIMD variables.
 *
 * \copydetails simdOrFI
 */
static inline SimdDInt32
simdOrDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] | b.i[i];
    }
    return c;
}

/*! \brief Integer bitwise xor for SIMD variables.
 *
 * \copydetails simdXorFI
 */
static inline SimdDInt32
simdXorDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
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
 * \copydetails simdAddFI
 */
static inline SimdDInt32
simdAddDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] + b.i[i];
    }
    return c;
}

/*! \brief Subtract SIMD integers, corresponding to double precision.
 *
 * \copydetails simdSubFI
 */
static inline SimdDInt32
simdSubDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] - b.i[i];
    }
    return c;
}

/*! \brief Multiply SIMD integers, corresponding to double precision.
 *
 * \copydetails simdMulFI
 */
static inline SimdDInt32
simdMulDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDInt32         c;
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
 * \copydetails simdCmpEqFI
 */
static inline SimdDIBool
simdCmpEqDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDIBool         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] == b.i[i]);
    }
    return c;
}

/*! \brief Less-than comparison of two ints corresponding to double SIMD data.
 *
 * \copydetails simdCmpLtFI
 */
static inline SimdDIBool
simdCmpLtDI(SimdDInt32 a, SimdDInt32 b)
{
    SimdDIBool         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] < b.i[i]);
    }
    return c;
}

/*! \brief Logical AND on SimdDIBool.
 *
 * \copydetails simdAndFIB
 */
static inline SimdDIBool
simdAndDIB(SimdDIBool a, SimdDIBool b)
{
    SimdDIBool        c;
    int               i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical OR on SimdDIBool.
 *
 * \copydetails simdOrFIB
 */
static inline SimdDIBool
simdOrDIB(SimdDIBool a, SimdDIBool b)
{
    SimdDIBool         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/*! \brief Returns non-zero if any of the double-int SIMD booleans in x is True, otherwise 0.
 *
 * \copydetails simdAnyTrueFIB
 */
static inline int
simdAnyTrueDIB(SimdDIBool a)
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
 * \copydetails simdMaskFI
 */
static inline SimdDInt32
simdMaskDI(SimdDInt32 a, SimdDIBool mask)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = mask.b[i] ? a.i[i] : 0.0;
    }
    return c;
}

/*! \brief Select from SIMD ints (corresponding to double) where boolean is false.
 *
 * \copydetails simdMaskNotFI
 */
static inline SimdDInt32
simdMaskNotDI(SimdDInt32 a, SimdDIBool mask)
{
    SimdDInt32         c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = mask.b[i] ? 0.0 : a.i[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection for double-int SIMD.
 *
 * \copydetails simdBlendFI
 */
static inline SimdDInt32
simdBlendDI(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    SimdDInt32         d;
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
 * \copydetails simdCvtF2I
 */
static inline SimdDInt32
simdCvtD2I(SimdDouble a)
{
    SimdDInt32         b;
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
 * \copydetails simdCvttF2I
 */
static inline SimdDInt32
simdCvttD2I(SimdDouble a)
{
    SimdDInt32         b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        b.i[i] = a.r[i];
    }
    return b;
};

/*! \brief Convert integer to double precision floating-point.
 *
 * \copydetails simdCvtI2F
 */
static inline SimdDouble
simdCvtI2D(SimdDInt32 a)
{
    SimdDouble         b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        b.r[i] = a.i[i];
    }
    return b;
};

/*! \brief Convert from double boolean to corresponding integer boolean.
 *
 * \copydetails simdCvtFB2FIB
 */
static inline SimdDIBool
simdCvtDB2DIB(SimdDBool a)
{
    SimdDIBool         b;
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
 * \copydetails simdCvtFIB2FB
 */
static inline SimdDBool
simdCvtDIB2DB(SimdDIBool a)
{
    SimdDBool         b;
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
static inline SimdDouble
simdCvtF2D(SimdFloat f)
{
    SimdDouble        d;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int               i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = f.r[i];
    }
#else
    gmx_fatal(FARGS, "simdCvtF2D() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
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
static inline SimdFloat
simdCvtD2F(SimdDouble d)
{
    SimdFloat        f;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.r[i] = d.r[i];
    }
#else
    gmx_fatal(FARGS, "simdCvtD2F() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
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
static inline void
simdCvtF2DD(SimdFloat f, SimdDouble *d0, SimdDouble *d1)
{
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d0->r[i] = f.r[i];
        d1->r[i] = f.r[GMX_SIMD_DOUBLE_WIDTH+i];
    }
#else
    gmx_fatal(FARGS, "simdCvtF2DD() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
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
static inline SimdFloat
simdCvtDD2F(SimdDouble d0, SimdDouble d1)
{
    SimdFloat        f;
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.r[i]                       = d0.r[i];
        f.r[GMX_SIMD_DOUBLE_WIDTH+i] = d1.r[i];
    }
#else
    gmx_fatal(FARGS, "simdCvtDD2F() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments & uninitialized f */
    f.r[0] = d0.r[0] + d1.r[0];
#endif
    return f;
}

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD_DOUBLE_H */
