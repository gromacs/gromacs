/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
/*! \libinternal \file
 * \brief
 * This implements basic nblib test systems
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_TESTS_TESTHELPERS_H
#define NBLIB_TESTS_TESTHELPERS_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "nblib/box.h"
#include "nblib/vector.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace nblib
{

namespace test
{

//! Compare between two instances of the Box object
bool operator==(const Box& a, const Box& b);

/*! \internal \brief
 *  Simple test harness for checking 3D vectors like coordinates, velocities,
 *  forces against reference data
 *
 */
class Vector3DTest
{
public:
    Vector3DTest() : checker_(refData_.rootChecker())
    {
        gmx::test::FloatingPointTolerance tolerance(
                gmx::test::FloatingPointTolerance(1e-8, 1.0e-12, 1e-8, 1.0e-12, 200, 100, true));
        checker_.setDefaultTolerance(tolerance);
    }

    explicit Vector3DTest(const gmx::test::FloatingPointTolerance& tolerance) :
        checker_(refData_.rootChecker())
    {
        checker_.setDefaultTolerance(tolerance);
    }

    //! Compare a given input vector of cartesians with the reference data
    void testVectors(gmx::ArrayRef<Vec3> forces, const std::string& testString)
    {
        checker_.checkSequence(forces.begin(), forces.end(), testString.c_str());
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
#endif // NBLIB_TESTS_TESTHELPERS_H
