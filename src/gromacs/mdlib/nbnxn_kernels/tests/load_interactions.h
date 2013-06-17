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
 * Function definitions for 4xn kernel functionality for loading interaction masks.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */
#ifndef _simd_test_load_interactions_h_
#define _simd_test_load_interactions_h_

#include "general_kernel.h"

namespace SIMDTests
{

template<class SimdReal4FunctionSet>
static inline void
prepareToReadExclusions(SimdReal4FunctionSet                          &simdReal4FunctionSet,
                        int                                            filter_stride,
                        int                                            unrollj,
                        int                                            interaction_mask,
                        nbnxn_atomdata_t *                            &nbat,
                        typename SimdReal4FunctionSet::exclfilterType *masks,
                        char                                           interaction_mask_indices[GMX_SIMD_WIDTH_HERE])
{
    snew(nbat, 1);
    nbnxn_atomdata_init_simple_exclusion_masks(nbat);

    // Prepare support array for x86 for filtering out individual (i,j) interactions
    for (size_t i = 0; i < UNROLLI; ++i)
    {
        unsigned *simd_exclusion_filter;
        if (1 == filter_stride)
        {
            simd_exclusion_filter = nbat->simd_exclusion_filter1;
        }
        else
        {
            simd_exclusion_filter = nbat->simd_exclusion_filter2;
        }
        masks[i] = simdReal4FunctionSet.load_exclusion_filter(simd_exclusion_filter + i*unrollj*filter_stride);
    }

#ifdef GMX_CPU_ACCELERATION_IBM_QPX
    // Prepare support array for QPX. x86 uses 4x4 bits in an int to
    // encode as a bit mask which (i,j) pairs will interact. QPX uses
    // 4 bits in each of 4 chars to encode the same information. The
    // value of each char is used to look up simd_interaction_array
    // the values that can be loaded as SIMD booleans to do the actual
    // masking in the kernels.
    for (int j = 0; j < GMX_SIMD_WIDTH_HERE; j++)
    {
        int jj = j % UNROLLI;
        /* Get the four bits of the input for this j atom */
        interaction_mask_indices[j] = (interaction_mask >> (jj*4)) & 0xF;
    }
#endif
}

template<class SimdFunctionSet>
static inline void
cleanUpFromReadingExclusions(typename SimdFunctionSet::boolType       *bResults,
                             nbnxn_atomdata_t                         *nbat,
                             typename SimdFunctionSet::exclfilterType *masks)
{
    sfree(bResults);
    sfree(masks);
    sfree(nbat->simd_4xn_diagonal_j_minus_i);
    sfree(nbat->simd_2xnn_diagonal_j_minus_i);
    sfree(nbat->simd_exclusion_filter1);
    sfree(nbat->simd_exclusion_filter2);
#ifdef GMX_CPU_ACCELERATION_IBM_QPX
    sfree(nbat->simd_interaction_array);
#endif
    sfree(nbat);
}

/*! This function needs to be defined differently for 2xnn vs 4xn exclusion masking */
template<class SimdFunctionSet,
         typename FunctionType>
void
callFunctionLocal(SimdFunctionSet           &simdFunctionSet,
                  FunctionType               function,
                  int                        filter_stride,
                  int                        unrollj,
                  std::vector<unsigned int*> inputs,
                  real                      *result);

template <typename InputType>
template<typename ReferenceFunctionType,
         typename TestFunctionType,
         typename ReferenceSimdFunctionSet,
         typename ReferenceSimdReal4FunctionSet,
         typename TestSimdFunctionSet,
         typename TestSimdReal4FunctionSet,
         typename InputKind>
void
SimdFunctionTest_kernel<InputType>::
    RunTest(ReferenceFunctionType         referenceFunction,
            TestFunctionType              testFunction,
            int                           filter_stride,
            int                           unrollj,
            ReferenceSimdFunctionSet      referenceSimdFunctionSet,
            ReferenceSimdReal4FunctionSet referenceSimdReal4FunctionSet,
            TestSimdFunctionSet           testSimdFunctionSet,
            TestSimdReal4FunctionSet      testSimdReal4FunctionSet,
            InputKind                     inputKind)
{
    GMX_ASSERT(Parent::bDonePrepare_, "Must call prepare() before RunTest()");
    for (int k = 0; k < g_numberOfRepeats; ++k)
    {
        Parent::generateTheInputs(inputKind);
        callFunctionLocal(referenceSimdFunctionSet, referenceSimdReal4FunctionSet, referenceFunction, filter_stride, unrollj, Parent::inputs_, Parent::referenceResult_);
        callFunctionLocal(testSimdFunctionSet, testSimdReal4FunctionSet, testFunction, filter_stride, unrollj, Parent::inputs_, Parent::testResult_);
        Parent::testTheOutputs(GMX_SIMD_WIDTH_HERE);
    }
}

}      // namespace

#endif /* _simd_test_load_interactions_h_ */
