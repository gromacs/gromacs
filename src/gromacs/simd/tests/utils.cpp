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
 * Utility class and function for testing of SIMD functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */
#include "utils.h"
#include <gtest/internal/gtest-internal.h>
#include <sstream>

namespace SIMDTests
{

/*! \brief Class derived from GoogleTest to provide a way to tune
 * whether two floats should compare as equal. This needs to vary to
 * suit the requirements of the test.
 *
 * For example, "SIMD" functions that use table lookups,
 * rational-function approximations, or Newton-Raphson iterations will
 * not be binary exact with the reference versions, because generally
 * we don't need the accuracy in GROMACS.
 *
 * Also caters to the x86 quirk that SIMD-true is a NaN, and the IEEE
 * rules (which are otherwise followed by GoogleTest) specify that any
 * comparison with any NaN must be false (even with another NaN that
 * is bitwise identical).
 */
template <typename RawType>
class FloatingPoint : public ::testing::internal::FloatingPoint<RawType>
{
    public:
        //! Constructor
        FloatingPoint(const RawType &x) : ::testing::internal::FloatingPoint<RawType>(x)
        {
        }

        //! Returns true iff this number is at most kMaxUlps * scale
        // ULP's away from rhs.  In particular, this function:
        //
        //   - returns false if either number is (or both are) NAN.
        //   - treats really large numbers as almost equal to infinity.
        //   - thinks +0.0 and -0.0 are 0 DLP's apart.
        //
        // Based on Google Tests's FloatingPoint<RawType>::AlmostEquals
        bool ScaleableAlmostEquals(const ::testing::internal::FloatingPoint<RawType> &rhs, real scaleMaxUlps) const
        {
            typedef ::testing::internal::FloatingPoint<real>::Bits UnsignedIntWithSizeOfReal;

            // TODO fix this later
//#ifdef GMX_X86
            /* The IEEE standard says that any comparison operation
               involving a NaN must return false. But x86 SIMD uses
               0xffffffff as true, and that is a NaN. So we need the
               ability to test for that, first. */
            UnsignedIntWithSizeOfReal simdTrue = ~0;
            if ((simdTrue == this->bits()) &&
                (simdTrue == rhs.bits())) { return true; }
//#endif
            if (this->is_nan() || rhs.is_nan()) {return false; }

            return ::testing::internal::FloatingPoint<RawType>::DistanceBetweenSignAndMagnitudeNumbers(this->bits(), rhs.bits())
                   <= this->kMaxUlps * scaleMaxUlps;
        }
};

::testing::AssertionResult
RealArraysAreEqual(const real         *expected,
                   const real         *actual,
                   const unsigned long length,
                   const real          scaleMaxUlps)
{
    bool              bAllWereEqual(true);
    std::stringstream errorString;
    typedef ::testing::internal::FloatingPoint<real>::Bits UnsignedIntWithSizeOfReal;

    for (unsigned long i = 0; i < length; ++i)
    {
        const FloatingPoint<real> lhs(expected[i]), rhs(actual[i]);
        if (!lhs.ScaleableAlmostEquals(rhs, scaleMaxUlps))
        {
            if (!bAllWereEqual)
            {
                /* Print one error per line when there is more than
                 * one error */
                errorString << ",\n";
            }
            bAllWereEqual = false;
            BitManipulater manipulater_expected, manipulater_actual;
            manipulater_expected.r = expected[i];
            manipulater_actual.r   = actual[i];
            real difference = expected[i] - actual[i];
            errorString << std::scientific;
            errorString.precision(17);
            errorString << "element [" << std::dec << (i+1) << "] was not equal: expected "
            << expected[i] << " (0x" << std::hex << manipulater_expected.i
            << ") but was " << std::dec << actual[i] << " (0x"
            << std::hex << manipulater_actual.i
            << ") and the difference was " << std::dec << difference;
        }
    }

    if (bAllWereEqual)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure() << errorString.str();
    }
}

} // namespace
