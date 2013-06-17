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
 * Tests for functionality that loads memory value(s) to a SIMD vector.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"

namespace SIMDTests
{

/* It's hard to test load and store independently, since you need one
   of them to be in a position to test the other. If a test here fails
   because of store_pr, so will almost all the others. If a test here
   fails because of the load function under test, then those other
   tests that don't load a real may still pass. */

//! Typedef for the test fixture
typedef SimdFunctionTest<real> SimdFunctionLoad;

template<class SimdFunctionSet,
         typename FunctionType> void
callFunction(SimdFunctionSet     &simdFunctionSet,
             FunctionType         function,
             std::vector<real*>   inputs,
             real                *result)
{
    typename SimdFunctionSet::realType result_pr;

    result_pr = function(inputs[0]);
    simdFunctionSet.store_pr(result, result_pr);
}

/* Wrapper function because the macro wrapped is defined with a
   parameter. */
static gmx_inline gmx_mm_pr gmx_load_pr_wrapper(const real *a)
{
    return gmx_load_pr(a);
}

TEST_F(SimdFunctionLoad, gmx_load_pr_Works)
{
    prepare(1, 1);
    RunTest(ReferenceFunctions::gmx_simd_ref_load_pr,
            gmx_load_pr_wrapper, real(0));
}

/* Wrapper function because the macro wrapped is defined with a
   parameter. */
static gmx_inline gmx_mm_pr gmx_load1_pr_wrapper(const real *a)
{
    return gmx_load1_pr(a);
}

TEST_F(SimdFunctionLoad, gmx_load1_pr_Works)
{
    prepare(1, 1);
    RunTest(ReferenceFunctions::gmx_simd_ref_load1_pr,
            gmx_load1_pr_wrapper, real(0));
}

} // namespace
