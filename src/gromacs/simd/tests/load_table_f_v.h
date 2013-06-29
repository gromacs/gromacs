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
 * \brief Tests for functionality for table loads common to both 2xnn
 * and 4xn kernels. TAB_FDV0 can change between the two, so we compile
 * the two test fixtures in separate object files, and the common code
 * goes here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

/* TODO this makes more sense as a nbnxn-specific file, rather than a
   general SIMD test file. */

typedef SimdFunctionDoingTableLoad SimdFunctionDoingTableLoad_f_v;

template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunctionDoingTableLoad_f_v::callFunction(SimdFunctionSet &simdFunctionSet,
                                                  FunctionType     function,
                                                  real            *_result)
{
    typename SimdFunctionSet::realType indices;
    indices   = simdFunctionSet.load_pr(index);

    typename SimdFunctionSet::epiType index_epi;
    index_epi = simdFunctionSet.cvttpr_epi32(indices);

    typename SimdFunctionSet::realType *result;
    //    snew_aligned(result, numOutputs_, GMX_SIMD_WIDTH_HERE * sizeof(real));
    snew(result, numOutputs_);

#ifdef TAB_FDV0
    function(table, index_epi, ti, &result[0], &result[1], &result[2]);
#else
    function(table, table, index_epi, ti, &result[0], &result[1], &result[2]);
#endif

    for (size_t i = 0; i != numOutputs_; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}
