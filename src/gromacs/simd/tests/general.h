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
/*! \internal \file \brief Defines general test fixture and supporting
 * classes for testing SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#include "base.h"
#include <string>
#include <sstream>

namespace SIMDTests
{

/*! We need some types to use so that we can distinguish the type of
 * the test input data and the kind of data it contains. This lets us
 * call the correct function to generate the right kind of data. */

/*! \brief Describes unsigned integers with only a single bit set. */
class TestWithSingleBitsSet
{
    public:
        TestWithSingleBitsSet(unsigned int) {};
};

/*! \brief Describes reals that contain booleans that the current SIMD
 * hardware will interpret as either true or false. */
class booleanReal
{
    public:
        booleanReal(unsigned int) {};
};

/*! \brief Test fixture for testing SIMD functions
 *
 * \tparam NInputs The number of input vectors expected by the SIMD function
 * \tparam NOutputs The number of output vectors produced by the SIMD function
 * \tparam InputType The type of the input vectors (real or unsigned int)
 *
 * The methods support the ability to loops over the number of repeats
 * specified on the command line with -numrepeats, generating random
 * input values, calling the reference and hardware SIMD functions on
 * the same inputs, and verifing that the output vector(s) are within
 * the specified tolerance of each other.
 *
 */

template <size_t NInputs,
          size_t NOutputs,
          typename InputType>
class SimdFunctionTest :
    public SimdFunctionTestBase<NOutputs>
{
    public:
        typedef SimdFunctionTestBase<NOutputs> Parent;

        SimdFunctionTest() : Parent(), shift(NInputs), scale(NInputs),
                             maxUlps(1)
        {
            for (size_t i = 0; i < NInputs; ++i)
            {
                snew_aligned(inputs[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(InputType));
            }
        }

        ~SimdFunctionTest()
        {
            for (size_t i = 0; i < NInputs; ++i)
            {
                sfree(inputs[i]);
            }
        }

        /* TODO Don't think it is useful to declare a generateInput()
           templated on InputKind, because the definition gets
           rejected unless the class template parameters (NInputs,
           NOutputs) are also specialized, and that code is more
           repetitive than this code here. */

        /*! Generate random reals in scale * (shift + [0,1)) */
        real generateInput(real, real scale, real shift)
        {
            return scale*(shift+gmx_rng_uniform_real(rng));
        }

        /*! Set a random bit in a 32-bit integer. Used for testing
            exclusion mask loads. */
        unsigned int generateInput(TestWithSingleBitsSet, real, unsigned int)
        {
            return 1 << (gmx_rng_uniform_uint32(rng) & 31);
        }

        /*! Get a random SIMD-boolean in a real variable. Different
            hardware encodes true and false in SIMD vectors in
            different ways. */
        real generateInput(booleanReal, real, real)
        {
            /* Generate random booleans according to the SIMD
             * hardware's view of true and false. */
            BitManipulater manipulater;
            /* When QPX is supported, or GMX_X86 exists, the following
             * lines will need to be versioned. */
            UnsignedIntWithSizeOfReal simdTrue  = ~0;
            UnsignedIntWithSizeOfReal simdFalse = 0;
            // Pretty wasteful of random bits, but no big deal.
            manipulater.i      = (gmx_rng_uniform_uint32(rng) <
                                  gmx_rng_uniform_uint32(rng)) ?
                simdTrue : simdFalse;
            return manipulater.r;
        }

        /*! This function gets defined differently in different
            specializations of the template class, to reflect the
            different use cases. */
        template<class SimdFunctionSet,
                 typename FunctionType>
        void callFunction(SimdFunctionSet &simdFunctionSet,
                          FunctionType     function,
                          real            *_result);

        /*! When a test fails, provide information about the inputs to
            help judge whether the output vectors are significantly wrong. */
        std::string const generateTraceString(size_t outputVectorIndex)
        {
            std::stringstream ss;
            ss << std::scientific;
            ss.precision(17);
            ss << "The above failure occurred while testing "
            << ::testing::UnitTest::GetInstance()->current_test_info()->name()
            << " with inputs \n";
            for (size_t i = 0; i < NInputs; ++i)
            {
                ss << (i+1) << ": ";
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    ss << inputs[i][j] << " ";
                }
                ss << "\n";
            }
            ss << "Error occured in the above element of output vector "
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
                 typename TestFunctionType,
                 typename InputKind>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType      testFunction,
                    InputKind = 0)
        {
            for (int k = 0; k < numberOfRepeats; ++k)
            {
                for (size_t i = 0; i < NInputs; ++i)
                {
                    // Make random values
                    for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                    {
                        // Declare a variable to force calling the
                        // correct form of generateInput
                        InputKind temp = 0;
                        inputs[i][j] = generateInput(temp, scale[i], shift[i]);
                    }
                }

                callFunction(referenceSimd, referenceFunction, referenceResult);
                callFunction(hardwareSimd, testFunction, testResult);

                for (size_t i = 0; i < NOutputs; ++i)
                {
                    ::testing::AssertionResult outputRealsAreEqual
                        = RealArraysAreEqual(referenceResult + i * GMX_SIMD_WIDTH_HERE,
                                             testResult + i * GMX_SIMD_WIDTH_HERE,
                                             GMX_SIMD_WIDTH_HERE,
                                             maxUlps);
                    EXPECT_TRUE(outputRealsAreEqual) << generateTraceString(i);
                }
            }
        }

        // Help the compiler find symbols in templated parent class
        using Parent::referenceSimd;
        using Parent::hardwareSimd;
        using Parent::rng;
        using Parent::referenceResult;
        using Parent::testResult;
        using Parent::numberOfRepeats;

        /*! A static constant seems to be found more regularly than a
            template parameter by the the compiler when compiling
            specializations of callFunctions(). */
        static const size_t numInputs = NInputs;

        // TODO reimplement all the data handling using
        // aligned-memory containers if/when we get such things.

        //! Pre-allocated input vectors to SIMD functions.
        InputType             *inputs[NInputs];
        // Variables to control the range of the random reals used for inputs.
        std::vector<InputType> shift;
        std::vector<real>      scale;
        // Maximum acceptable difference in "units in last place"
        // (ulps) for equality when comparing reals. Default is 1, and
        // produces standard Google Test behaviour. Larger values are
        // necessary for (e.g.) reduced-precision SIMD followed by
        // Newton-Raphson iterates, particularly in double precision.
        real maxUlps;
};

} // namespace
