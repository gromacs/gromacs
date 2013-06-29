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

#include "base.h"

namespace SIMDTests
{

#if (defined GMX_NBNXN_SIMD_2XNN) || (defined GMX_NBNXN_SIMD_4XN)

#ifdef GMX_NBNXN_SIMD_4XN
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_4xn_outer_header.h"
#endif

#ifdef GMX_NBNXN_SIMD_2XNN
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_2xnn_outer_header.h"
#endif

namespace TestFunctions
{
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils.h"
}

namespace ReferenceFunctions
{
#define GMX_SIMD_REFERENCE_PLAIN_C
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_ref.h"
#undef GMX_SIMD_REFERENCE_PLAIN_C
}

//! TODO this makes more sense as a nbnxn-specific test, rather than a
//! general SIMD test

template <size_t NOutputs>
class SimdFunctionDoingTableLoad :
            public SimdFunctionTestBase<NOutputs>
{
    public:
        typedef SimdFunctionTestBase<NOutputs> Parent;

        SimdFunctionDoingTableLoad() : Parent()
        {
            snew_aligned(table, table_size, GMX_SIMD_WIDTH_HERE * sizeof(real));
            for(size_t i = 0; i != table_size; ++i)
            {
                /* Generate random floating-point numbers. It is not clear
                 * what range makes any particular sense, so pick 0.5 >= x >
                 * -0.5. */
                table[i] = gmx_rng_uniform_real(rng) - 0.5;
            }

            snew_aligned(index, GMX_SIMD_WIDTH_HERE, GMX_SIMD_WIDTH_HERE * sizeof(real));

            ti = TestFunctions::gmx_simd_align_int(ti_array);
        }

        ~SimdFunctionDoingTableLoad()
        {
            sfree(table);
            sfree(index);
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
            for (int i = 0; i < SIMD_TEST_NUM_REAL_VARIABLES; ++i)
            {
                // make random non-overflowing indices into the table
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    index[j] = gmx_rng_uniform_real(rng) * (table_size / GMX_SIMD_WIDTH_HERE);
                }

                callFunction(referenceSimd, referenceFunction, referenceResult);
                callFunction(hardwareSimd, testFunction, testResult);

                EXPECT_TRUE(RealArraysAreEqual(referenceResult,
                                               testResult,
                                               GMX_SIMD_WIDTH_HERE * NOutputs,
                                               1));
            }
        }

        using Parent::referenceSimd;
        using Parent::hardwareSimd;
        using Parent::rng;
        using Parent::referenceResult;
        using Parent::testResult;

        // Picked an arbitrary table size comparable with the size of
        // tables that are actually used.
        static const size_t table_size = 1024 * GMX_SIMD_WIDTH_HERE;

        // Table of random floats to look up in tests
        real *table;

        // SIMD vector of indices to look up in the table. In the kernels,
        // this is computed from reals, so we do the same here.
        real *index;

        // Temp array required by some implementations of load_table*
        // Allocate twice the necessary size to be sure alignment can
        // be achieved.
        int ti_array[GMX_SIMD_WIDTH_HERE+GMX_SIMD_WIDTH_HERE], *ti;
};

// ===

/* Specialization of the class for testing load_table_f */
typedef SimdFunctionDoingTableLoad<2> SimdFunctionDoingTableLoad_f;

template<>
template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunctionDoingTableLoad_f::callFunction(SimdFunctionSet &simdFunctionSet,
                                                FunctionType function,
                                                real *_result)
{
    typename SimdFunctionSet::epiType index_epi;
    typename SimdFunctionSet::realType result[2];
    typename SimdFunctionSet::realType indices;

    indices = simdFunctionSet.load_pr(index);
    index_epi = simdFunctionSet.cvttpr_epi32(indices);
    
    function(table, index_epi, ti, &result[0], &result[1]);
    
    for (size_t i = 0; i != numOutputs; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}

TEST_F(SimdFunctionDoingTableLoad_f, load_table_f_Works)
{
    
    Tester(ReferenceFunctions::load_table_f,
           TestFunctions::load_table_f);
}

// ===

/* Specialization of the class for testing load_table_f_v */
typedef SimdFunctionDoingTableLoad<3> SimdFunctionDoingTableLoad_f_v;

template<>
template<class SimdFunctionSet,
         typename FunctionType>
void SimdFunctionDoingTableLoad_f_v::callFunction(SimdFunctionSet &simdFunctionSet,
                                                  FunctionType function,
                                                  real *_result)
{
    typename SimdFunctionSet::epiType index_epi;
    typename SimdFunctionSet::realType result[numOutputs];
    typename SimdFunctionSet::realType indices;

    indices = simdFunctionSet.load_pr(index);
    index_epi = simdFunctionSet.cvttpr_epi32(indices);
    
#ifdef TAB_FDV0
    function(table, index_epi, ti, &result[0], &result[1], &result[2]);
#else
    function(table, table, index_epi, ti, &result[0], &result[1], &result[2]);
#endif

    for (size_t i = 0; i != numOutputs; ++i)
    {
        simdFunctionSet.store_pr(_result + i * GMX_SIMD_WIDTH_HERE, result[i]);
    }
}

TEST_F(SimdFunctionDoingTableLoad_f_v, load_table_f_v_Works)
{
    
    Tester(ReferenceFunctions::load_table_f_v,
           TestFunctions::load_table_f_v);
}

#endif

} // namespace
