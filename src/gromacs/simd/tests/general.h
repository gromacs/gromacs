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

/*! Need some types to use in templates so we can distinguish the type
 * of the test input data and the kind of data it contains. */

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

template <size_t NInputs,
          size_t NOutputs,
          typename InputType>
class SimdFunctionTest :
            public SimdFunctionTestBase<NOutputs>
{
    public:
        typedef SimdFunctionTestBase<NOutputs> Parent;

        SimdFunctionTest() : Parent(), shift(numInputs), scale(numInputs),
                                       maxUlps(1)
        {
            for (size_t i = 0; i < numInputs; ++i)
            {
                snew_aligned(inputs[i],
                             GMX_SIMD_WIDTH_HERE,
                             GMX_SIMD_WIDTH_HERE * sizeof(InputType));
            }
        }

        ~SimdFunctionTest()
        {
            for (size_t i = 0; i < numInputs; ++i)
            {
                sfree(inputs[i]);
            }
        }

        //! TODO Don't think it is useful to declare a generateInput()
        //! templated on InputKind, because the definition gets
        //! rejected unless the class template parameters (NInputs,
        //! NOutputs) are also specialized, and that code is more
        //! repetitive than this code here.
        /*! Generate random floating-point numbers in scale * (shift +
         * [0,1)) */
        real generateInput(real, real scale, real shift)
        {
            return scale*(shift+gmx_rng_uniform_real(rng));
        }

        /*! Generate random floating-point numbers in scale * (shift +
         * [0,MAX_INT]) */
        /* Unused
        unsigned int generateInput(unsigned int, real scale, unsigned int shift)
        {
            return (unsigned int) scale*(shift+gmx_rng_uniform_uint32(rng));
        }
        */

        /*! Set a random bit in a 32-bit integer. */
        //! TODO Find a way to get a guaranteed 32-bit integer?
        unsigned int generateInput(TestWithSingleBitsSet, real, unsigned int)
        {
            return 1 << (gmx_rng_uniform_uint32(rng) & 31);
        }

        /*! Get a random SIMD-boolean in a real variable. */
        real generateInput(booleanReal, real, real)
        {
            /* Generate random booleans according to the SIMD
             * hardware's view of true and false. */
            BitManipulater manipulater;
            /* When QPX is supported, or GMX_X86 exists, the following
             * lines will need to be versioned. */
            UnsignedIntWithSizeOfReal simdTrue = ~0;
            UnsignedIntWithSizeOfReal simdFalse = 0;
            // Pretty wasteful of random bits, but no big deal.
            manipulater.i      = (gmx_rng_uniform_uint32(rng) <
                                  gmx_rng_uniform_uint32(rng)) ?
                simdTrue : simdFalse;
            return manipulater.r;
        }

        //! This function gets defined differently in different
        // specializations.
        template<class SimdFunctionSet,
                 typename FunctionType>
        void callFunction(SimdFunctionSet &simdFunctionSet,
                          FunctionType function,
                          real *_result);

        template<typename ReferenceFunctionType,
                 typename TestFunctionType,
                 typename InputKind>
        void Tester(ReferenceFunctionType referenceFunction,
                    TestFunctionType testFunction,
                    InputKind = 0)
        {
            for (int k = 0; k < SIMD_TEST_NUM_REAL_VARIABLES; ++k)
            {
                for (size_t i = 0; i < numInputs; ++i)
                {
                    // make random values
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

                EXPECT_TRUE(RealArraysAreEqual(referenceResult,
                                               testResult,
                                               GMX_SIMD_WIDTH_HERE * NOutputs,
                                               maxUlps));
            }
        }

        using Parent::referenceSimd;
        using Parent::hardwareSimd;
        using Parent::rng;
        using Parent::referenceResult;
        using Parent::testResult;
        
        static const size_t numInputs = NInputs;

        //! TODO reimplement all the data handling using
        // aligned-memory containers if/when we get such things.
        InputType *inputs[NInputs];
        // Variables to control the range of the random reals used for inputs
        std::vector<InputType> shift;
        std::vector<real> scale;
        // Maximum acceptable difference in "units in last place"
        // (ulps) for equality when comparing reals. Default is 1, and
        // produces standard Google Test behaviour. Larger values are
        // necessary for (e.g.) reduced-precision SIMD followed by
        // Newton-Raphson iterates, particularly in double precision.
        real maxUlps;
};

} // namespace
