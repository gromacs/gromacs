/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include <climits>

#include "gmxapi/version.h"
#include <gtest/gtest.h>

namespace
{

using gmxapi::Version;

// Copy header version info with intentionally sloppy type-ing to try to catch
// worst-case scenarios and unexpected behavior. Otherwise the isAtLeast function
// uses major(), minor(), and patch() so testing them might be superfluous.

const int current_major = gmxapi::GMXAPI_MAJOR;
const int current_minor = gmxapi::GMXAPI_MINOR;
const int current_patch = gmxapi::GMXAPI_PATCH;

/*!
 * \brief Check basic Version interface functionality.
 */
TEST(VersionTest, SaneComparisons)
{
    ASSERT_TRUE(Version::isAtLeast(0,
                                   0,
                                   0));
    // Should negative version values cause a sanity-check error?
    ASSERT_TRUE(Version::isAtLeast(-1,
                                   -1,
                                   -1));
    ASSERT_FALSE(Version::isAtLeast(SHRT_MAX,
                                    SHRT_MAX,
                                    SHRT_MAX));
    ASSERT_TRUE(Version::isAtLeast(current_major,
                                   current_minor,
                                   current_patch));
    ASSERT_FALSE(Version::isAtLeast(current_major + 1,
                                    current_minor,
                                    current_patch));
    ASSERT_FALSE(Version::isAtLeast(current_major,
                                    current_minor + 1,
                                    current_patch));
    ASSERT_FALSE(Version::isAtLeast(current_major,
                                    current_minor,
                                    current_patch + 1));
}

/*!
 * \brief Check whether gmxapi correctly advertises or refutes feature availability.
 *
 * A few unimplemented features are tests with ``ASSERT_FALSE`` just for sanity
 * checking and to give an idea of near-term targeted named features. If the
 * feature is available, it is expected to conform to the API specification
 * for the library Version::release(). As we discover features that break
 * forward-compatibility of the API, we will have to provide developer documentation
 * or sample code for build-time CMake feature checks.
 */
TEST(VersionTest, Named0_1_Features)
{
    ASSERT_FALSE(Version::hasFeature(""));
    ASSERT_FALSE(Version::hasFeature("MD_plugin_restraint_force"));
    ASSERT_FALSE(Version::hasFeature("MD_plugin_restraint_callback"));
    ASSERT_FALSE(Version::hasFeature("MD_plugin_mpi_domain_decomposition"));
    ASSERT_FALSE(Version::hasFeature("MD_stop_signal"));
    ASSERT_FALSE(Version::hasFeature("MD_set_final_trajectory_step"));
    ASSERT_FALSE(Version::hasFeature("gmxapi_communicator_from_client"));
    ASSERT_FALSE(Version::hasFeature("gmxapi_simulation_from_tpr"));
}

} // end anonymous namespace
