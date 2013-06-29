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
 * Tests for SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"

namespace SIMDTests
{

/* Ideally tests for SIMD logical operations would load values into
 * SIMD registers and test them. However, there is no implementation
 * of a SIMD load of a vector of bool, because there is no need for
 * that in the kernels. There is also the issue that all current
 * (i.e. x86) implementations of these functions are actually done as
 * bitwise operations on register types that are identical to the
 * floating-point ones. */

/********************************************************************
 * Tests for SIMD wrapper functions
 */

namespace ReferenceFunctions
{
/* Bit-wise AND on SIMD reals, to mimic x86. The official
   gmx_simd_ref_and_pb does a logical AND, so we can't use that
   for testing the actual SIMD implementations. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_bitwise_and_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr result;

    BitManipulater  manipulater_a, manipulater_b, manipulater_c;
    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        manipulater_a.r = a.r[i];
        manipulater_b.r = b.r[i];
        manipulater_c.i = manipulater_a.i & manipulater_b.i;
        result.r[i]     = manipulater_c.r;
    }

    return result;
}

/* Bit-wise OR on SIMD reals, to mimic x86. The official
   gmx_simd_ref_and_pb does a logical OR, so we can't use that for
   testing the actual SIMD implementations. */
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_bitwise_or_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr result;

    BitManipulater  manipulater_a, manipulater_b, manipulater_c;
    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        manipulater_a.r = a.r[i];
        manipulater_b.r = b.r[i];
        manipulater_c.i = manipulater_a.i | manipulater_b.i;
        result.r[i]     = manipulater_c.r;
    }

    return result;
}

}

typedef SimdFunctionTest<2,1,real> SimdFunction_logical;

/* Specialization of the class for testing logical
 * functions.
 *
 * Irritatingly, this duplicates the content of
 * SimdFunctionWithTwoRealInputs::callFunction. The only difference is in the
 * InputKind template parameter of SimdFunctionTest<>. I don't think
 * this can be avoided, because a template function (ie. callFunction)
 * cannot be specialized without also fully specializing the class
 * type (ie. SimdFunctionTest). Neither does it seem worth trying to
 * actually call SimdFunctionWithTwoRealInputs::callFunction(). */
template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunction_logical::callFunction(SimdFunctionSet &simdFunctionSet,
                                       FunctionType function,
                                       real *_result)
{
    typename SimdFunctionSet::realType a, b, result;
    
    a      = simdFunctionSet.load_pr(inputs[0]);
    b      = simdFunctionSet.load_pr(inputs[1]);
    result = function(a, b);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunction_logical, gmx_and_pb_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_bitwise_and_pr,
           TestFunctions::gmx_and_pb, booleanReal(0));
}

TEST_F(SimdFunction_logical, gmx_or_pb_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_bitwise_or_pr,
           TestFunctions::gmx_or_pb, booleanReal(0));
}

// ===

#ifdef GMX_SIMD_HAVE_ANYTRUE

namespace ReferenceFunctions
{
/* gmx_anytrue_pb reads a gmx_mm_pb, but all the implementations
 * (i.e. x86) just do an implicit conversion of SIMD real to SIMD
 * bool. To test, we have to do the same thing, which we do by
 * wrapping the reference code.
 */
static gmx_inline int
gmx_simd_ref_anytrue_pr(gmx_simd_ref_pr a)
{
    int             i;
    gmx_simd_ref_pb b;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        /* This might need a change to support QPX SIMD, where all
           less than zero are interpreted as FALSE. */
        b.r[i] = (0.0 != a.r[i]);
    }

    return gmx_simd_ref_anytrue_pb(b);
}

}

namespace TestFunctions
{
/* gmx_anytrue_pb reads a gmx_mm_pb, but all the implementations
 * (i.e. x86) just do an implicit conversion of SIMD real to SIMD
 * bool. To test, we have to do the same thing, which we do by
 * wrapping the reference code.
 */
static gmx_inline int
gmx_anytrue_pr(gmx_mm_pr a)
{
    return (real) (0 != gmx_anytrue_pb(a));
}

}

typedef SimdFunctionTest<1,1,real> SimdFunctionAnytrue;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionAnytrue::callFunction(SimdFunctionSet &simdFunctionSet,
                                  FunctionType function,
                                  real *_result)
{
    typename SimdFunctionSet::realType a;
    int result;

    a      = simdFunctionSet.load_pr(inputs[0]);
    result = function(a);
    _result[0] = (real) (0 != result);
}

TEST_F(SimdFunctionAnytrue, gmx_anytrue_pb_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_anytrue_pr,
           TestFunctions::gmx_anytrue_pr, booleanReal(0));
}

#endif /* GMX_SIMD_HAVE_ANYTRUE */

} // namespace
