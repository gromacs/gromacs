/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests for functionality that does boolean logic in SIMD.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"

namespace SIMDTests
{

#ifdef GMX_SIMD_HAVE_ANYTRUE

namespace ReferenceFunctions
{
/* gmx_anytrue_pb reads a gmx_mm_pb, but all the x86 implementations
 * just do an implicit conversion of SIMD real to SIMD bool. To test,
 * we have to do the same thing, which we do by wrapping the reference
 * code to map a possible SIMD bool to a reference-SIMD bool.
 */
static gmx_inline int
gmx_simd_ref_anytrue_pr(gmx_simd_ref_pr a)
{
    int             i;
    gmx_simd_ref_pb b;

    /* /todo This might need a change to support QPX SIMD, where all
       less than zero are interpreted as FALSE. Or it might be
       redundant? */
    real simdFalse = 0;
    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        b.r[i] = (simdFalse == a.r[i]) ? 0 : 1;
    }

    return gmx_simd_ref_anytrue_pb(b);
}

}

namespace TestFunctions
{

/* On x86, gmx_anytrue_pb is implemented with *movemask*, which
 * returns a bitmask that specifies which fields were true. However,
 * the return value is only used as a simple boolean and that is how
 * the reference code is implemented. So we wrap the real SIMD
 * function with the logic to reduce it to a simple boolean.
 */
static gmx_inline int
gmx_anytrue_pr(gmx_mm_pr a)
{
    return (0 != gmx_anytrue_pb(a));
}

}

typedef SimdFunctionTest<real> SimdFunctionAnytrue;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionAnytrue::callFunction(SimdFunctionSet &simdFunctionSet,
                                  FunctionType     function,
                                  real            *_result)
{
    typename SimdFunctionSet::realType a;
    int result;

    a          = simdFunctionSet.load_pr(inputs_[0]);
    result     = function(a);
    _result[0] = (real) (0 != result);
}

TEST_F(SimdFunctionAnytrue, gmx_anytrue_pb_Works)
{
    prepare(1, 1);
    Tester(ReferenceFunctions::gmx_simd_ref_anytrue_pr,
           TestFunctions::gmx_anytrue_pr, booleanReal(0));
}

#endif /* GMX_SIMD_HAVE_ANYTRUE */

}      // namespace
