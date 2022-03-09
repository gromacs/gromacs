/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * This implements basic nblib test systems
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_TESTHELPERS_H
#define NBLIB_TESTHELPERS_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/conftest.h"
#include "testutils/refdata.h"

#include "nblib/box.h"
#include "nblib/vector.h"

namespace nblib
{

namespace test
{

/*! \internal \brief
 *  Simple test harness for checking 3D vectors like coordinates, velocities,
 *  forces against reference data
 *
 */
class RefDataChecker
{
public:
    RefDataChecker() : checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(
                gmx::test::FloatingPointTolerance(1e-8, 1.0e-12, 1e-8, 1.0e-12, 200, 100, true));
        checker_.setDefaultTolerance(tolerance);
    }

    RefDataChecker(real relativeFloatingPointTolerance) : checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(gmx::test::FloatingPointTolerance(
                1e-6, 1.0e-9, relativeFloatingPointTolerance, relativeFloatingPointTolerance, 200, 100, true));
        checker_.setDefaultTolerance(tolerance);
    }

    //! Compare a given input array of cartesians, reals, integers, etc with the reference data
    template<class T>
    void testArrays(gmx::ArrayRef<T> tArray, const std::string& testString)
    {
        checker_.checkSequence(tArray.begin(), tArray.end(), testString.c_str());
    }

    void testReal(real value, const std::string& testName)
    {
        checker_.checkReal(value, testName.c_str());
    }

private:
    gmx::test::TestReferenceData    refData_;
    gmx::test::TestReferenceChecker checker_;
};

//! Macros to compare floats and doubles with a specified tolerance
/// \cond DO_NOT_DOCUMENT
#if GMX_DOUBLE
#    define EXPECT_FLOAT_DOUBLE_EQ_TOL(value, refFloat, refDouble, tolerance) \
        EXPECT_DOUBLE_EQ_TOL(value, refDouble, tolerance)
#    define ASSERT_FLOAT_DOUBLE_EQ_TOL(value, refFloat, refDouble, tolerance) \
        ASSERT_DOUBLE_EQ_TOL(value, refDouble, tolerance)
#else
#    define EXPECT_FLOAT_DOUBLE_EQ_TOL(value, refFloat, refDouble, tolerance) \
        EXPECT_FLOAT_EQ_TOL(value, refFloat, tolerance)
#    define ASSERT_FLOAT_DOUBLE_EQ_TOL(value, refFloat, refDouble, tolerance) \
        ASSERT_FLOAT_EQ_TOL(value, refFloat, tolerance)
#endif
/// \endcond

} // namespace test
} // namespace nblib
#endif // NBLIB_TESTHELPERS_H
