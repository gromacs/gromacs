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
class SimdFunctionDoingTableLoad :
    public SimdFunctionTestBase
{
    public:
        typedef SimdFunctionTestBase Parent;
        using Parent::rng_;

        SimdFunctionDoingTableLoad() : Parent(), index(0),
            ti(TestFunctions::gmx_simd_align_int(ti_array)),
            bDonePrepare(false)
        {
            snew_aligned(index, GMX_SIMD_WIDTH_HERE, GMX_SIMD_WIDTH_HERE * sizeof(real));
        }

        /*! Prepare the table and output vectors. Ideally, this would
            take place in the constructor or SetUp(), but neither of
            those can take any arguments. One could also template the
            whole test fixture class to avoid needing arguments to a
            constructor or function call at all, but you'd still need
            a separate typedef for most of the TEST_F fixtures to make
            things work. */
        void prepare(size_t numTableEntries, size_t numOutputs)
        {
            numTableEntries_ = numTableEntries;
            numOutputs_ = numOutputs;

            OutputAllocater outputAllocater(numOutputs_);
            outputAllocater(referenceResult_);
            outputAllocater(testResult_);

            snew_aligned(table, GMX_SIMD_WIDTH_HERE * numTableEntries_, GMX_SIMD_WIDTH_HERE * sizeof(real));
            for (size_t i = 0; i != GMX_SIMD_WIDTH_HERE * numTableEntries_; ++i)
            {
                /* Generate random floating-point numbers. It is not clear
                 * what range makes any particular sense, so pick 0.5 >= x >
                 * -0.5. */
                table[i] = gmx_rng_uniform_real(rng_) - 0.5;
            }
            bDonePrepare = true;
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
                          FunctionType     function,
                          real            *_result);

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
         * \tparam ReferenceFunctionType The type signature of the reference SIMD function
         * \tparam TestFunctionType The type signature of the reference SIMD function
         * \tparam InputKind The kind of data present in the input vectors, which is distinct from the type of data.
         *
         * Only the compiler needs to worry about the actual type that
         * is ReferenceFunctionType and TestFunctionType. Tests can
         * just pass in the reference and test versions of the
         * function and forget about it.
         *
         * Writes detailed output in failing cases. */
        template<typename ReferenceFunctionType,
                 typename TestFunctionType>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType      testFunction)
        {
            GMX_ASSERT(bDonePrepare, "Must call prepare() before Tester()");
            for (int k = 0; k < numberOfRepeats_; ++k)
            {
                // Make random non-overflowing indices into the table
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    index[j] = gmx_rng_uniform_real(rng_) * numTableEntries_;
                }

                callFunction(referenceSimd_, referenceFunction, referenceResult_);
                callFunction(hardwareSimd_, testFunction, testResult_);

                ::testing::AssertionResult outputRealsAreEqual
                    = RealArraysAreEqual(referenceResult_,
                                         testResult_,
                                         GMX_SIMD_WIDTH_HERE * numOutputs_,
                                         1);
                EXPECT_TRUE(outputRealsAreEqual);
            }
        }

        using Parent::hardwareSimd_;
        using Parent::referenceSimd_;

        //! Number of entries in the table that will be looked up
        //! during tests. The table will actually contain
        //! numTableEntries_ * GMX_SIMD_WIDTH_HERE reals.
        size_t numTableEntries_;

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

    private:
        //! Support preventing incorrect use of Tester() without prepare()
        bool bDonePrepare;
};
