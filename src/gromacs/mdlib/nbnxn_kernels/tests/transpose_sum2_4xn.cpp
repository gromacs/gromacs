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
 * Tests for transpose SIMD functionality for 4xn kernels
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general_kernel.h"

namespace SIMDTests
{

#ifdef GMX_NBNXN_SIMD_4XN

#if GMX_SIMD_WIDTH_HERE == 2

//! Typedef for the test fixture
typedef SimdFunctionTest_kernel<real> SimdFunction4xnWithTwoRealInputs;

template<class SimdFunctionSet,
         typename FunctionType>
void
callFunctionLocal(SimdFunctionSet   &simdFunctionSet,
                  FunctionType       function,
                  std::vector<real*> inputs,
                  real              *result)
{
    typename SimdFunctionSet::realType a, b;
    typename SimdFunctionSet::realType result_pr;

    a         = simdFunctionSet.load_pr(inputs[0]);
    b         = simdFunctionSet.load_pr(inputs[1]);
    result_pr = function(a, b);
    simdFunctionSet.store_pr(result, result_pr);
}

template<>
template<typename ReferenceFunctionType,
         typename TestFunctionType,
         typename InputKind>
void
SimdFunction4xnWithTwoRealInputs::RunTest(ReferenceFunctionType referenceFunction,
                                          TestFunctionType      testFunction,
                                          int                   outputSimdWidth,
                                          InputKind             inputKind)
{
    GMX_ASSERT(bDonePrepare_, "Must call prepare() before RunTest()");
    for (int k = 0; k < g_numberOfRepeats; ++k)
    {
        generateTheInputs(inputKind);
        callFunctionLocal(referenceSimd_, referenceFunction, inputs_, referenceResult_);
        callFunctionLocal(hardwareSimd_, testFunction, inputs_, testResult_);
        testTheOutputs(outputSimdWidth);
    }
}

// TODO does not pass in SSE2 double, i.e. SIMD width 2, which is the
// only time the test is useful!
#if 2 != GMX_SIMD_WIDTH_HERE

/* The order of the addition in the reduction can be significant when
   testing. So the tolerance is quite loose. */
TEST_F(SimdFunction4xnWithTwoRealInputs, gmx_mm_transpose_sum2_pr_WorksFor4xn)
{
#ifndef GMX_DOUBLE
    maxUlps_ = 1e7; // empirically determined to be enough on x86
#endif
    prepare(2, 1, -0.5, 5);
    RunTest(gmx_simd_ref_transpose_sum2_pr,
            gmx_mm_transpose_sum2_pr,
            NBNXN_CPU_CLUSTER_I_SIZE,
            real(0));
}

#endif /* 2 != GMX_SIMD_WIDTH_HERE */

#endif

#endif

} // namespace
