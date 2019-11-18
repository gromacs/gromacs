/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Extra GoogleMock matchers for unit tests.
 *
 * This file provides the usual kind of GoogleMock matchers that
 * extend the usefulness of GoogleMock EXPECT_THAT constructs to the
 * kinds of containers of reals commonly used. This means that test
 * code can write one-liners rather than loops over whole containers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTMATCHERS_H
#define GMX_TESTUTILS_TESTMATCHERS_H

#include <memory>

#include <gmock/gmock.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/*! \brief Make matcher for floats for use with GoogleMock that compare
 * equal when \c tolerance is satisifed.
 *
 * Used like
 *
 *   EXPECT_THAT(testFloats, Pointwise(FloatEq(tolerance), referenceFloats));
 */
testing::Matcher<std::tuple<float, float>> FloatEq(const FloatingPointTolerance& tolerance);

/*! \brief Make matcher for doubles for use with GoogleMock that compare
 * equal when \c tolerance is satisifed.
 *
 * Used like
 *
 *   EXPECT_THAT(testDoubles, Pointwise(DoubleEq(tolerance), referenceDoubles));
 */
testing::Matcher<std::tuple<double, double>> DoubleEq(const FloatingPointTolerance& tolerance);

/*! \brief Make matcher for reals for use with GoogleMock that compare
 * equal when \c tolerance is satisifed.
 *
 * Used like
 *
 *   EXPECT_THAT(testReals, Pointwise(RealEq(tolerance), referenceReals));
 */
testing::Matcher<std::tuple<real, real>> RealEq(const FloatingPointTolerance& tolerance);

/*! \brief Make matcher for RVecs for use with GoogleMock that compare
 * equal when \c tolerance is satisifed.
 *
 * Used like
 *
 *   EXPECT_THAT(testRVecs, Pointwise(RVecEq(tolerance), referenceRVecs));
 */
testing::Matcher<std::tuple<RVec, RVec>> RVecEq(const FloatingPointTolerance& tolerance);

} // namespace test
} // namespace gmx

#endif
