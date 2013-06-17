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

namespace ReferenceFunctions
{
//! Helper typedef for reference versions of SIMD logical-operation functions
typedef gmx_simd_ref_pb (*gmx_simd_ref_logical_pb)(gmx_simd_ref_pb, gmx_simd_ref_pb);

/*! \brief Helper function for doing SIMD-logical operations on SIMD
    reals.

    \tparam logical_function A SIMD reference implementation of a
    logical function (like "and" or "or)

    This is complex because true and false vary with the SIMD
    implementation, and because x86 true is a Nan (so that you cannot
    compare with it, because comparison with a NaN is always
    false).

    So, the hardware SIMD wrapper classes implement static functions
    to test for and set SIMD boolean values in real variables, so that
    this wrapper for the reference implementation can correctly mimic
    the hardware behaviour.

    Finally, we template the wrapper function so that we can test
    multiple SIMD logical functions without needing to duplicate
    code. */

template <gmx_simd_ref_logical_pb logical_function>
gmx_inline gmx_simd_ref_pr
gmx_simd_ref_boolean_wrapper_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    gmx_simd_ref_pr result;
    gmx_simd_ref_pb a_pb, b_pb, result_pb;

    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a_pb.r[i] = TestFunctions::HardwareSimdFunctionSet::isSimdTrue(a.r[i]);
        b_pb.r[i] = TestFunctions::HardwareSimdFunctionSet::isSimdTrue(b.r[i]);
    }

    result_pb = logical_function(a_pb, b_pb);

    for (int i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        result.r[i] = TestFunctions::HardwareSimdFunctionSet::getSimdBool(result_pb.r[i]);
    }

    return result;
}

/* TODO remove this?
static gmx_inline gmx_simd_ref_pr
gmx_simd_ref_and_wrapper_pr(gmx_simd_ref_pr a, gmx_simd_ref_pr b)
{
    return gmx_simd_ref_boolean_wrapper_pr<gmx_simd_ref_and_pb>(a, b);
}
*/

}

//! Typedef for the test fixture
typedef SimdFunctionTest<real> SimdFunction_logical;

/* Irritatingly, this duplicates the content of
 * SimdFunctionWithTwoRealInputs::callFunction. The only difference is in the
 * InputKind template parameter of SimdFunctionTest<>. I don't think
 * this can be avoided, because a template function (ie. callFunction)
 * cannot be specialized without also fully specializing the class
 * type (ie. SimdFunctionTest). Neither does it seem worth trying to
 * actually call SimdFunctionWithTwoRealInputs::callFunction(). */
template<class SimdFunctionSet,
         typename FunctionType>
void
callFunction(SimdFunctionSet     &simdFunctionSet,
             FunctionType         function,
             std::vector<real*>   inputs,
             real                *result)
{
    typename SimdFunctionSet::realType a, b, result_pr;

    a         = simdFunctionSet.load_pr(inputs[0]);
    b         = simdFunctionSet.load_pr(inputs[1]);
    result_pr = function(a, b);
    simdFunctionSet.store_pr(result, result_pr);
}

/* Wrapper function because MSVC can't cope otherwise. */
static gmx_inline gmx_mm_pb gmx_and_pb_wrapper(gmx_mm_pb a, gmx_mm_pb b)
{
    return gmx_and_pb(a, b);
}

ReferenceFunctions::gmx_simd_ref_boolean_wrapper_pr<ReferenceFunctions::gmx_simd_ref_and_pb>;

TEST_F(SimdFunction_logical, gmx_and_pb_Works)
{
    prepare(2, 1);
    RunTest(ReferenceFunctions::gmx_simd_ref_boolean_wrapper_pr<ReferenceFunctions::gmx_simd_ref_and_pb>,
            // TODO remove this
            //    RunTest(ReferenceFunctions::gmx_simd_ref_and_wrapper_pr,
            gmx_and_pb_wrapper, booleanReal(0));
}

/* Wrapper function because MSVC can't cope otherwise. */
static gmx_inline gmx_mm_pb gmx_or_pb_wrapper(gmx_mm_pb a, gmx_mm_pb b)
{
    return gmx_or_pb(a, b);
}

ReferenceFunctions::gmx_simd_ref_boolean_wrapper_pr<ReferenceFunctions::gmx_simd_ref_or_pb>;

TEST_F(SimdFunction_logical, gmx_or_pb_Works)
{
    prepare(2, 1);
    RunTest(ReferenceFunctions::gmx_simd_ref_boolean_wrapper_pr<ReferenceFunctions::gmx_simd_ref_or_pb>,
            gmx_or_pb_wrapper, booleanReal(0));
}

}      // namespace
