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

/*! \brief Test fixture for table-load functions for 2xnn and 4xn kernels.
 *
 * \tparam NOutputs The number of output vectors produced by the SIMD function
 */
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

        /*! When a test fails, provide information about the inputs to
          help diagnose the problem. Any failure is a serious
          problem. */
        std::string const generateTraceString(size_t outputVectorIndex)
        {
            std::stringstream ss;
            ss << "The above failure occurred while testing " 
               << ::testing::UnitTest::GetInstance()->current_test_info()->name()
               << " with indices \n";
            for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
            {
                ss << index[j] << " ";
            }
            ss << "\nError occured in the above element of output vector "
               << (outputVectorIndex+1);
            return ss.str();
        }

        /*! \brief Do repetitive testing of the SIMD and reference versions of the function.
         *
         * \tparam ReferenceFunctionType The type signature of the reference SIMD function (only the compiler needs to worry about this)
         * \tparam TestFunctionType The type signature of the reference SIMD function (only the compiler needs to worry about this)
         * \tparam InputKind The kind of data present in the input vectors, which is distinct from the type of data. 
         *
         * Writes detailed output in failing cases. */
        template<typename ReferenceFunctionType,
                 typename TestFunctionType>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType testFunction)
        {
            for (int k = 0; k < numberOfRepeats; ++k)
            {
                // Make random non-overflowing indices into the table
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    index[j] = gmx_rng_uniform_real(rng) * numTableEntries;
                }

                callFunction(referenceSimd, referenceFunction, referenceResult);
                callFunction(hardwareSimd, testFunction, testResult);

                ::testing::AssertionResult outputRealsAreEqual
                      = RealArraysAreEqual(referenceResult,
                                           testResult,
                                           GMX_SIMD_WIDTH_HERE * NOutputs,
                                           1);
                EXPECT_TRUE(outputRealsAreEqual);
            }
        }

        // Help the compiler find symbols in templated parent class
        using Parent::referenceSimd;
        using Parent::hardwareSimd;
        using Parent::rng;
        using Parent::referenceResult;
        using Parent::testResult;
        using Parent::numberOfRepeats;

        // Picked an arbitrary table size comparable with the size of
        // tables that are actually used.
        static const size_t numTableEntries = 1024;
        static const size_t table_size = numTableEntries * GMX_SIMD_WIDTH_HERE;

        //! Table of random floats to look up in tests
        real *table;

        /*! SIMD vector of indices to look up in the table. In the kernels,
          this is computed from reals, so we do the same here. */
        real *index;

        /*! Temp array and pointer required by some
            implementations of load_table*. Allocate twice the
            necessary size to be sure SIMD alignment can be
            achieved within it. */
        int ti_array[2*GMX_SIMD_WIDTH_HERE], *ti;
};

// ===

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

// ===

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
