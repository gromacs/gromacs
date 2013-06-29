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
 * Tests for functionality that does masking-style operations in SIMD.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "general.h"

namespace SIMDTests
{

/* We can't actually do a SIMD load of a vector of bool from memory,
 * because it is not implemented since there is no need for that in
 * the kernels. Instead, a vector of bool is always generated from a
 * comparison, so we do that in testing also. gmx_cmplt_pr is tested
 * separately. */

//! Typedef for the test fixture
typedef SimdFunctionTest<real> SimdFunctionMasknotadd;

template<>
template<class SimdFunctionSet,
         typename FunctionType> void
SimdFunctionMasknotadd::callFunction(SimdFunctionSet &simdFunctionSet,
                                     FunctionType     function,
                                     real            *_result)
{
    typename SimdFunctionSet::realType a, b, c, d, result;
    typename SimdFunctionSet::boolType mask;

    a      = simdFunctionSet.load_pr(inputs_[0]);
    b      = simdFunctionSet.load_pr(inputs_[1]);
    c      = simdFunctionSet.load_pr(inputs_[2]);
    d      = simdFunctionSet.load_pr(inputs_[3]);
    mask   = simdFunctionSet.cmplt_pr(c, d);
    result = function(mask, a, b);
    simdFunctionSet.store_pr(_result, result);
}

TEST_F(SimdFunctionMasknotadd, gmx_masknot_add_pr_Works)
{
    prepare(4, 1, 5, 0);
    Tester(ReferenceFunctions::gmx_simd_ref_masknot_add_pr,
           TestFunctions::gmx_masknot_add_pr, real(0));
}
} // namespace
