/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include "gmxapi/version.h"

#include <climits>

#include <gtest/gtest.h>

#include "testingconfiguration.h"

namespace gmxapi
{

namespace testing
{

namespace
{

using gmxapi::Version;

/* Copy header version info with intentionally sloppy type-ing to try to catch
 * worst-case scenarios and unexpected behavior. Otherwise the isAtLeast function
 * uses major(), minor(), and patch() so testing them might be superfluous.
 */
//! \cond
const int current_major = gmxapi::c_majorVersion;
const int current_minor = gmxapi::c_minorVersion;
const int current_patch = gmxapi::c_patchVersion;
//! \endcond

/*!
 * \brief Check basic Version interface functionality.
 */
TEST_F(GmxApiTest, SaneVersionComparisons)
{
    EXPECT_TRUE(Version::isAtLeast(0, 0, 0));
    EXPECT_FALSE(Version::isAtLeast(SHRT_MAX, SHRT_MAX, SHRT_MAX));
    EXPECT_TRUE(Version::isAtLeast(current_major, current_minor, current_patch));
    EXPECT_FALSE(Version::isAtLeast(current_major + 1, current_minor, current_patch));
    EXPECT_FALSE(Version::isAtLeast(current_major, current_minor + 1, current_patch));
    EXPECT_FALSE(Version::isAtLeast(current_major, current_minor, current_patch + 1));
}

/*!
 * \brief Check whether gmxapi correctly advertises or refutes feature availability.
 *
 * Check for correct responses from the Version API for features or
 * functionality not (yet) guaranteed by the current API version.
 * If a feature is available, it is expected to conform to the API specification
 * for the library Version::release(). As we discover features that break
 * forward-compatibility of the API, we will have to provide developer documentation
 * or sample code for build-time CMake feature checks.
 *
 * This is the test for pre-0.1 features leading up to that specification.
 *
 * \internal
 * Designed but unimplemented features should be tested with ``EXPECT_FALSE``
 * until they are implemented, then toggled to ``EXPECT_TRUE`` as implemented as
 * extensions of the current API spec. (There aren't any yet.)
 */
TEST_F(GmxApiTest, VersionNamed0_1_Features)
{
    EXPECT_FALSE(Version::hasFeature(""));
    EXPECT_FALSE(Version::hasFeature("nonexistent feature"));
}

} // end anonymous namespace

} // namespace testing

} // namespace gmxapi
