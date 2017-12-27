/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "simd4.h"

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"

#include "data.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#if GMX_SIMD4_HAVE_REAL

const Simd4Real rSimd4_c0c1c2   = setSimd4RealFrom3R( c0, c1, c2);
const Simd4Real rSimd4_c3c4c5   = setSimd4RealFrom3R( c3, c4, c5);
const Simd4Real rSimd4_c6c7c8   = setSimd4RealFrom3R( c6, c7, c8);
const Simd4Real rSimd4_c3c0c4   = setSimd4RealFrom3R( c3, c0, c4);
const Simd4Real rSimd4_c4c6c8   = setSimd4RealFrom3R( c4, c6, c8);
const Simd4Real rSimd4_c7c2c3   = setSimd4RealFrom3R( c7, c2, c3);
const Simd4Real rSimd4_m0m1m2   = setSimd4RealFrom3R(-c0, -c1, -c2);
const Simd4Real rSimd4_m3m0m4   = setSimd4RealFrom3R(-c3, -c0, -c4);
const Simd4Real rSimd4_2p25     = setSimd4RealFrom1R(2.25);
const Simd4Real rSimd4_3p75     = setSimd4RealFrom1R(3.75);
const Simd4Real rSimd4_m2p25    = setSimd4RealFrom1R(-2.25);
const Simd4Real rSimd4_m3p75    = setSimd4RealFrom1R(-3.75);

#if GMX_SIMD_HAVE_LOGICAL
// The numbers below all have exponent (2^0), which will not change with AND/OR operations.
// We also leave the last part of the mantissa as zeros, to avoid rounding issues in the compiler
#if GMX_DOUBLE
const Simd4Real rSimd4_logicalA         = setSimd4RealFrom1R(1.3333333332557231188); // mantissa 01010101010101010101010101010101
const Simd4Real rSimd4_logicalB         = setSimd4RealFrom1R(1.7999999998137354851); // mantissa 11001100110011001100110011001100
const Simd4Real rSimd4_logicalResultAnd = setSimd4RealFrom1R(1.266666666604578495);  // mantissa 01000100010001000100010001000100
const Simd4Real rSimd4_logicalResultOr  = setSimd4RealFrom1R(1.8666666664648801088); // mantissa 11011101110111011101110111011101
#else                                                                                // GMX_DOUBLE
const Simd4Real rSimd4_logicalA         = setSimd4RealFrom1R(1.3333282470703125);    // mantissa 0101010101010101
const Simd4Real rSimd4_logicalB         = setSimd4RealFrom1R(1.79998779296875);      // mantissa 1100110011001100
const Simd4Real rSimd4_logicalResultAnd = setSimd4RealFrom1R(1.26666259765625);      // mantissa 0100010001000100
const Simd4Real rSimd4_logicalResultOr  = setSimd4RealFrom1R(1.8666534423828125);    // mantissa 1101110111011101
#endif                                                                               // GMX_DOUBLE
#endif                                                                               // GMX_SIMD_HAVE_LOGICAL

::std::vector<real>
simd4Real2Vector(const Simd4Real simd4)
{
    alignas(GMX_SIMD_ALIGNMENT) real  mem[GMX_SIMD4_WIDTH];

    store4(mem, simd4);
    std::vector<real>   v(mem, mem+GMX_SIMD4_WIDTH);

    return v;
}

Simd4Real
vector2Simd4Real(const std::vector<real> &v)
{
    alignas(GMX_SIMD_ALIGNMENT) real  mem[GMX_SIMD4_WIDTH];

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        mem[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return load4(mem);
}

Simd4Real
setSimd4RealFrom3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2Simd4Real(v);
}

Simd4Real
setSimd4RealFrom1R(real value)
{
    std::vector<real> v(GMX_SIMD4_WIDTH);
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2Simd4Real(v);
}

testing::AssertionResult
Simd4Test::compareSimd4RealUlp(const char *  refExpr,     const char *  tstExpr,
                               const Simd4Real ref, const Simd4Real tst)
{
    return compareVectorRealUlp(refExpr, tstExpr, simd4Real2Vector(ref), simd4Real2Vector(tst));
}

testing::AssertionResult
Simd4Test::compareSimd4RealEq(const char * refExpr, const char * tstExpr,
                              const Simd4Real ref, const Simd4Real tst)
{
    return compareVectorEq(refExpr, tstExpr, simd4Real2Vector(ref), simd4Real2Vector(tst));
}

#endif  // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace

#endif // GMX_SIMD
