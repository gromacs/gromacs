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
#include "utils.h"
#include <gtest/internal/gtest-internal.h>
#include <sstream>

namespace SIMDTests
{

template <typename RawType>
class FloatingPoint : public ::testing::internal::FloatingPoint<RawType>
{
    public:
        FloatingPoint(const RawType &x) : ::testing::internal::FloatingPoint<RawType>(x)
        {
        }

        // Returns true iff this number is at most kMaxUlps * scale
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

            // The IEEE standard says that any comparison operation
            // involving a NAN must return false. But x86 SIMD uses
            // 0xffffffff as true, and that is a NAN. So we need the
            // ability to test for that, first.
            if (((UnsignedIntWithSizeOfReal) ~0 == this->bits()) &&
                ((UnsignedIntWithSizeOfReal) ~0 == rhs.bits())) { return true; }
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
            manipulater_actual.r = actual[i];
            errorString << "element [" << i << "] was not equal: expected "
                        << expected[i] << " (0x" << std::hex << manipulater_expected.i
                        << ") but was "
                        << actual[i] << " (0x" << std::hex << manipulater_actual.i << ")";
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
