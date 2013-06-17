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
 * \brief Defines general test fixture and supporting
 * classes for testing SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */

#ifndef GMX_SIMD_TESTS_GENERAL_H
#define GMX_SIMD_TESTS_GENERAL_H

#include "base.h"
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

namespace SIMDTests
{

/*! We need some types to use so that we can distinguish the type of
 * the test input data and the kind of data it contains. This lets us
 * call the correct function to generate the right kind of data. */

/*! \brief Describes unsigned integers with only a single bit set. */
class TestWithSingleBitsSet
{
    public:
        //! Constructor
        TestWithSingleBitsSet(unsigned int) {};
};

/*! \brief Describes 32-bit unsigned integers. */
class TestWith32BitInt
{
    public:
        //! Constructor
        TestWith32BitInt(unsigned int) {};
};

/*! \brief Describes reals that contain booleans that the current SIMD
 * hardware will interpret as either true or false. */
class booleanReal
{
    public:
        //! Constructor
        booleanReal(unsigned int) {};
};

/*! \brief Test fixture for testing SIMD functions
 *
 * \tparam InputType The type of the input vectors (real or unsigned int)
 *
 * The methods support the ability to loops over the number of repeats
 * specified on the command line with -numrepeats, generating random
 * input values, calling the reference and hardware SIMD functions on
 * the same inputs, and verifing that the output vector(s) are within
 * the specified tolerance of each other.
 *
 * It would be fractionally more elegant to define a class that does
 * testing of SIMD functions that has a proper constructor, rather
 * than use the prepare() method below. The prepare() function is
 * needed because a text fixture class must be
 * default-constructible. This would mean there is no need to derive a
 * text fixture from \::testing::Test at all - a typedef would do.
 */
template <typename InputType>
class SimdFunctionTest :
    public SimdFunctionTestBase
{
    public:
        //! Convenience typedef
        typedef SimdFunctionTestBase Parent;
        //! Convenience typedef
        typedef AlignedSimdAllocater<InputType> Allocater;
        using Parent::OutputAllocater;

        SimdFunctionTest() :
            Parent(), numInputs_(0), inputs_(),
            shift_(), scale_(),
            maxUlps_(1), bDonePrepare_(false)
        {
        }

        /*! Prepare the input and output vectors, and the values that
            will control the generation of the output. Ideally, this
            would take place in the constructor or SetUp(), but
            neither of those can take any arguments. One could also
            template the whole test fixture class to avoid needing
            arguments to a constructor or function call at all, but
            you'd still need a separate typedef for most of the TEST_F
            fixtures to make things work. */
        void prepare(size_t numInputs, size_t numOutputs, InputType shift = 0, real scale = 1.0)
        {
            numInputs_  = numInputs;
            numOutputs_ = numOutputs;

            OutputAllocater outputAllocater(numOutputs_);
            outputAllocater(referenceResult_);
            outputAllocater(testResult_);

            inputs_.resize(numInputs_);
            std::for_each(inputs_.begin(), inputs_.end(), Allocater(numInputs_));

            shift_.resize(numInputs_, shift);
            scale_.resize(numInputs_, scale);

            bDonePrepare_ = true;
        }

        ~SimdFunctionTest()
        {
            for (size_t i = 0; i < numInputs_; ++i)
            {
                sfree(inputs_[i]);
            }
        }

        /* TODO Don't think it is useful to declare a generateInput()
           templated on InputKind, because the definition gets
           rejected unless the class template parameters (numInputs_,
           NOutputs) are also specialized, and that code is more
           repetitive than this code here. */

        /*! Generate random reals in scale * (shift + [0,1)) */
        real generateInput(real, real shift, real scale)
        {
            return scale*(shift+gmx_rng_uniform_real(rng_));
        }

        /*! Generate a random bit in a 32-bit integer. Used for
            testing exclusion mask loads. */
        unsigned int generateInput(TestWithSingleBitsSet, unsigned int, real)
        {
            return 1 << (gmx_rng_uniform_uint32(rng_) & 31);
        }

        /*! Generate random 32-bit integer. Used for testing exclusion
            mask loads. */
        unsigned int generateInput(TestWith32BitInt, unsigned int, real)
        {
            unsigned int max32BitInt = 0xFFFF;
            return (gmx_rng_uniform_uint32(rng_) & max32BitInt);
        }

        /*! Get a random SIMD-boolean in a real variable. Different
            hardware encodes true and false in SIMD vectors in
            different ways. */
        real generateInput(booleanReal, real, real)
        {
            return hardwareSimd_.getSimdBool(gmx_rng_uniform_uint32(rng_) <
                                             gmx_rng_uniform_uint32(rng_));
        }

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
            for (size_t i = 0; i < numInputs_; ++i)
            {
                ss << (i+1) << ": ";
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    ss << inputs_[i][j] << " ";
                }
                ss << "\n";
            }
            ss << "Error occured in the above element of output vector "
            << (outputVectorIndex+1);
            return ss.str();
        }

        template <typename InputKind>
        void generateTheInputs(InputKind inputKind)
        {
            for (size_t i = 0; i < numInputs_; ++i)
            {
                // Make random values
                for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                {
                    inputs_[i][j] =
                        generateInput(inputKind,
                                      shift_[i],
                                      scale_[i]);
                }
            }
        }

        void testTheOutputs(int output_simd_width)
        {
            for (size_t i = 0; i < numOutputs_; ++i)
            {
                ::testing::AssertionResult outputRealsAreEqual
                    = RealArraysAreEqual(referenceResult_ + i * GMX_SIMD_WIDTH_HERE,
                                         testResult_ + i * GMX_SIMD_WIDTH_HERE,
                                         output_simd_width,
                                         maxUlps_);
                EXPECT_TRUE(outputRealsAreEqual) << generateTraceString(i);
            }
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
                 typename TestFunctionType,
                 typename InputKind>
        void RunTest(ReferenceFunctionType referenceFunction,
                     TestFunctionType      testFunction,
                     InputKind             inputKind = 0)
        {
            GMX_ASSERT(bDonePrepare_, "Must call prepare() before RunTest()");
            for (int k = 0; k < g_numberOfRepeats; ++k)
            {
                generateTheInputs(inputKind);
                callFunction(referenceSimd_, referenceFunction, inputs_, referenceResult_);
                callFunction(hardwareSimd_, testFunction, inputs_, testResult_);
                testTheOutputs(GMX_SIMD_WIDTH_HERE);
            }
        }

        //! Number of input vectors from the SIMD function under test
        size_t numInputs_;

        // \todo reimplement the data handling using aligned-memory
        // containers if/when we get such things.

        //! Pre-allocated input vectors to SIMD functions.
        std::vector<InputType *> inputs_;
        /*! Help control the range of the random reals used for
            inputs. If individual input vectors need different
            sub-ranges (e.g. the first input vector can have any sign,
            but the second must be non-negative), then the
            corresponding vector elements should be set. This should
            be done *after* the call to prepare(). */
        std::vector<InputType>   shift_;
        /*! Help control the range of the random reals used for
            inputs. */
        std::vector<real>        scale_;
        //! Maximum acceptable difference in "units in last place"
        // (ulps) for equality when comparing reals. Default is 1, and
        // produces standard Google Test behaviour. Larger values are
        // necessary for (e.g.) reduced-precision SIMD followed by
        // Newton-Raphson iterates, particularly in double precision.
        real maxUlps_;

    protected:
        //! Support preventing incorrect use of RunTest() without prepare()
        bool bDonePrepare_;
};

}      // namespace

#endif
