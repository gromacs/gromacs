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
 * text fixture from ::testing::Test at all - a typedef would do.
 */

template <typename InputType>
class SimdFunctionTestNew :
    public SimdFunctionTestBaseNew
{
    public:
        typedef SimdFunctionTestBaseNew Parent;
        // TODO fix this later
        //#ifdef GMX_X86
        typedef TestFunctions::x86SimdFunctionSet HardwareSimdFunctionSet;
        //#endif
        typedef AlignedSimdAllocater<InputType> Allocater;
        using Parent::OutputAllocater;

        SimdFunctionTestNew() :
            Parent(), numInputs_(0), inputs_(),
            shift_(), scale_(),
            maxUlps(1), bDonePrepare(false)
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
            numInputs_ = numInputs;
            numOutputs_ = numOutputs;

            OutputAllocater outputAllocater(numOutputs_);
            outputAllocater(referenceResult_);
            outputAllocater(testResult_);

            inputs_.resize(numInputs_);
            std::for_each(inputs_.begin(), inputs_.end(), Allocater(numInputs_));

            shift_.resize(numInputs_, shift);
            scale_.resize(numInputs_, scale);

            bDonePrepare = true;
        }

        ~SimdFunctionTestNew()
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

        /*! Set a random bit in a 32-bit integer. Used for testing
            exclusion mask loads. */
        unsigned int generateInput(TestWithSingleBitsSet, unsigned int, real)
        {
            return 1 << (gmx_rng_uniform_uint32(rng_) & 31);
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
            manipulater.i      = (gmx_rng_uniform_uint32(rng_) <
                                  gmx_rng_uniform_uint32(rng_)) ?
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
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType      testFunction,
                    InputKind             inputKind = 0)
        {
            GMX_ASSERT(bDonePrepare, "Must call prepare() before Tester()");
            for (int k = 0; k < numberOfRepeats_; ++k)
            {
                for (size_t i = 0; i < numInputs_; ++i)
                {
                    // Make random values
                    for (unsigned int j = 0; j != GMX_SIMD_WIDTH_HERE; ++j)
                    {
                        inputs_[i][j] = generateInput(inputKind, shift_[i], scale_[i]);
                    }
                }

                callFunction(referenceSimd_, referenceFunction, referenceResult_);
                callFunction(hardwareSimd_, testFunction, testResult_);

                for (size_t i = 0; i < numOutputs_; ++i)
                {
                    ::testing::AssertionResult outputRealsAreEqual
                        = RealArraysAreEqual(referenceResult_ + i * GMX_SIMD_WIDTH_HERE,
                                             testResult_ + i * GMX_SIMD_WIDTH_HERE,
                                             GMX_SIMD_WIDTH_HERE,
                                             maxUlps);
                    EXPECT_TRUE(outputRealsAreEqual) << generateTraceString(i);
                }
            }
        }

        //! Number of input vectors from the SIMD function under test
        size_t numInputs_;

        // TODO reimplement the data handling using aligned-memory
        // containers if/when we get such things.

        //! Pre-allocated input vectors to SIMD functions.
        std::vector<InputType *> inputs_;
        /*! Variables to control the range of the random reals used
            for inputs. If individual input vectors need different
            sub-ranges (e.g. the first input vector can have any sign,
            but the second must be non-negative), then the
            corresponding vector elements should be set. This should
            be done *after* the call to prepare(). */
        std::vector<InputType>   shift_;
        std::vector<real>        scale_;
        //! Maximum acceptable difference in "units in last place"
        // (ulps) for equality when comparing reals. Default is 1, and
        // produces standard Google Test behaviour. Larger values are
        // necessary for (e.g.) reduced-precision SIMD followed by
        // Newton-Raphson iterates, particularly in double precision.
        real maxUlps;

    private:
        //! Support preventing incorrect use of Tester() without prepare()
        bool bDonePrepare;
};

} // namespace
