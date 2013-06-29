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

/********************************************************************
 * Tests for SIMD wrapper functions
 */

/* It's hard to test set and store independently, since you need one
   of them to be in a position to test the other. If a test here fails
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
           TestFunctions::gmx_set1_pr);
}

typedef SimdFunctionTest<0,1,real> SimdFunctionSetzero;

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
