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
 * Tests for functionality that sets a SIMD vector to specified values
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"
#include <algorithm>

namespace SIMDTests
{

/* It's hard to test set and store independently, since we need one of
   them to be in a position to test the other. If a test here fails
   because of store_pr, so will almost all the others. If a test here
   fails because of the set function under test, then most of the
   other tests will pass, because they do not use a set. */

/* Test "set" functions */

typedef SimdFunctionTest<1,1,real> SimdFunctionSet1;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionSet1::callFunction(SimdFunctionSet &simdFunctionSet,
                               FunctionType function,
                               real *_result)
{
    typename SimdFunctionSet::realType result;
    
    result = function(inputs[0][0]);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionSet1, gmx_set1_pr_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_set1_pr,
           TestFunctions::gmx_set1_pr, real(0));
}

/*! \brief Test fixture for SIMD "set" functions
 *
 * \tparam NOutputs The number of output vectors produced by the SIMD function
 * \tparam InputType The type of the input vectors (real or unsigned int)
 * 
 * Ideally we could use
 *
 * typedef SimdFunctionTest<0,1,real> SimdFunctionSetzero;
 *
 * but lots of compilers complain about things that follow from the
 * (templated) number of inputs being zero. So we derive a new class
 * instead. Fortunately, that is simple.
 */
template <size_t NOutputs,
          typename InputType>
class SimdFunctionWithNoInputs :
            public SimdFunctionTestBase<NOutputs>
{
    public:
        typedef SimdFunctionTestBase<NOutputs> Parent;

        SimdFunctionWithNoInputs() : Parent(), maxUlps(1)
        {
        }

        //! This function gets defined differently in different
        // specializations.
        template<class SimdFunctionSet,
                 typename FunctionType>
        void callFunction(SimdFunctionSet &simdFunctionSet,
                          FunctionType function,
                          real *_result);

        template<typename ReferenceFunctionType,
                 typename TestFunctionType>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType testFunction)
        {
            for (int k = 0; k < numberOfRepeats; ++k)
            {
                InputType magicValue = -999;

                /* Ensure the result vectors are initially non-zero,
                   so that we test setzero properly. */
                std::fill(referenceResult, referenceResult + NOutputs * GMX_SIMD_WIDTH_HERE, magicValue);
                std::fill(testResult, testResult + NOutputs * GMX_SIMD_WIDTH_HERE, magicValue);

                callFunction(referenceSimd, referenceFunction, referenceResult);
                callFunction(hardwareSimd, testFunction, testResult);

                EXPECT_TRUE(RealArraysAreEqual(referenceResult,
                                               testResult,
                                               GMX_SIMD_WIDTH_HERE * NOutputs,
                                               maxUlps));
            }
        }

        // Help the compiler find symbols in templated parent class
        using Parent::referenceSimd;
        using Parent::hardwareSimd;
        using Parent::rng;
        using Parent::referenceResult;
        using Parent::testResult;
        using Parent::numberOfRepeats;
        
        // Maximum acceptable difference in "units in last place"
        // (ulps) for equality when comparing reals. Default is 1, and
        // produces standard Google Test behaviour. Larger values are
        // necessary for (e.g.) reduced-precision SIMD followed by
        // Newton-Raphson iterates, particularly in double precision.
        real maxUlps;
};

typedef SimdFunctionWithNoInputs<1,real> SimdFunctionSetzero;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionSetzero::callFunction(SimdFunctionSet &simdFunctionSet,
                               FunctionType function,
                               real *_result)
{
    typename SimdFunctionSet::realType result;
    
    result = function();
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionSetzero, gmx_setzero_pr_Works)
{
    Tester(ReferenceFunctions::gmx_simd_ref_setzero_pr,
           TestFunctions::gmx_setzero_pr);
}

} // namespace
