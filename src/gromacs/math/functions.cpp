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
/*! \internal \file
 * \brief
 * Implements simple math functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_math
 */

#include "gmxpre.h"

#include "functions.h"

#include "config.h"

#include <cstdint>

#include <array>
#include <limits>

#if GMX_NATIVE_WINDOWS
#    include <intrin.h> // _BitScanReverse, _BitScanReverse64
#endif

#include "gromacs/math/utilities.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

unsigned int
log2I(std::uint32_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(0) is undefined");
#if HAVE_BUILTIN_CLZ
    // gcc, clang. xor with sign bit should be optimized out
    return __builtin_clz(n) ^ 31U;
#elif HAVE_BITSCANREVERSE
    // icc, MSVC
    {
        unsigned long res;
        _BitScanReverse(&res, static_cast<unsigned long>(n));
        return static_cast<unsigned int>(res);
    }
#elif HAVE_CNTLZ4
    return 31 - __cntlz4(n);
#else
    // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup

    static const std::array<char, 256>
    log2TableByte =
    {{
         0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
     }};

    unsigned int result;
    unsigned int tmp1, tmp2;

    if ((tmp1 = n >> 16) != 0)
    {
        result = ((tmp2 = tmp1 >> 8) != 0) ? 24 + log2TableByte[tmp2] : 16 + log2TableByte[tmp1];
    }
    else
    {
        result = ((tmp2 = n >> 8) != 0) ? 8 + log2TableByte[tmp2] : log2TableByte[n];
    }
    return result;
#endif
}


unsigned int
log2I(std::uint64_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(0) is undefined");
#if HAVE_BUILTIN_CLZLL
    // gcc, icc, clang. xor with sign bit should be optimized out
    return __builtin_clzll(n) ^ 63U;
#elif HAVE_BITSCANREVERSE64
    unsigned long res;
    _BitScanReverse64(&res, static_cast<unsigned __int64>(n));
    return static_cast<unsigned int>(res);
#elif HAVE_CNTLZ8
    return 63 - __cntlz8(n);
#else

    // No 64-bit log2 instrinsic available. Solve it by calling our internal
    // 32-bit version (which in turn might defer to a software solution)

    std::uint32_t high32Bits = static_cast<std::uint32_t>(n>>32);
    std::uint32_t result;

    if (high32Bits)
    {
        result = log2I(high32Bits) + 32;
    }
    else
    {
        result = log2I(static_cast<std::uint32_t>(n));
    }

    return result;
#endif
}

unsigned int
log2I(std::int32_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(n) for n<=0 is undefined");
    return log2I(static_cast<std::uint32_t>(n));
}

unsigned int
log2I(std::int64_t n)
{
    GMX_ASSERT(n > 0, "The behavior of log(n) for n<=0 is undefined");
    return log2I(static_cast<std::uint64_t>(n));
}

std::int64_t
greatestCommonDivisor(std::int64_t   p,
                      std::int64_t   q)
{
    while (q != 0)
    {
        std::int64_t tmp = q;
        q                = p % q;
        p                = tmp;
    }
    return p;
}

double
erfinv(double x)
{
    double xabs = std::abs(x);

    if (xabs > 1.0)
    {
        return std::nan("");
    }

    if (x == 1.0)
    {
        return std::numeric_limits<double>::infinity();
    }

    if (x == -1.0)
    {
        return -std::numeric_limits<double>::infinity();
    }

    double res;

    if (xabs <= 0.7)
    {
        // Rational approximation in range [0,0.7]
        double z = x*x;
        double P = (((-0.140543331 * z + 0.914624893) * z - 1.645349621) * z + 0.886226899);
        double Q = ((((0.012229801 * z - 0.329097515) * z + 1.442710462) * z - 2.118377725) * z + 1.0);
        res = x * P/Q;
    }
    else
    {
        // Rational approximation in range [0.7,1)
        double z = std::sqrt(-std::log((1.0 - std::abs(x))/2.0));
        double P = ((1.641345311 * z + 3.429567803) * z - 1.624906493) * z - 1.970840454;
        double Q = (1.637067800 * z + 3.543889200) * z + 1.0;
        res = std::copysign(1.0, x) * P/Q;
    }

    // Double precision requires two N-R iterations
    res = res - (std::erf(res) - x)/( (2.0/std::sqrt(M_PI))*std::exp(-res*res));
    res = res - (std::erf(res) - x)/( (2.0/std::sqrt(M_PI))*std::exp(-res*res));

    return res;
}

float
erfinv(float x)
{
    float xabs = std::abs(x);

    if (xabs > 1.0f)
    {
        return std::nan("");
    }

    if (x == 1.0f)
    {
        return std::numeric_limits<float>::infinity();
    }

    if (x == -1.0f)
    {
        return -std::numeric_limits<float>::infinity();
    }

    float res;

    if (xabs <= 0.7f)
    {
        // Rational approximation in range [0,0.7]
        float z = x*x;
        float P = (((-0.140543331f * z + 0.914624893f) * z - 1.645349621f) * z + 0.886226899f);
        float Q = ((((0.012229801f * z - 0.329097515f) * z + 1.442710462f) * z - 2.118377725f) * z + 1.0f);
        res = x * P/Q;
    }
    else
    {
        // Rational approximation in range [0.7,1)
        float z = std::sqrt(-std::log((1.0 - std::abs(x))/2.0f));
        float P = ((1.641345311f * z + 3.429567803f) * z - 1.624906493f) * z - 1.970840454f;
        float Q = (1.637067800f * z + 3.543889200f) * z + 1.0f;
        res = std::copysign(1.0f, x) * P/Q;
    }

    // Single N-R iteration sufficient for single precision
    res = res - (std::erf(res) - x)/( (2.0f/std::sqrt(M_PI))*std::exp(-res*res));

    return res;
}

} // namespace gmx
