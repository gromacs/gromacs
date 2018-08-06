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
#include "gmxapi/version.h"
#include <gtest/gtest.h>


using namespace gmxapi;

const int current_major = GMXAPI_MAJOR;
const int current_minor = GMXAPI_MINOR;
const int current_patch = GMXAPI_PATCH;

namespace
{

// The is_at_least functin uses major(), minor(), and patch()
// so testing them is maybe superfluous

TEST(VersionTest, SaneComparisons)
{
    ASSERT_TRUE(Version::is_at_least(0, 0, 0));
    ASSERT_FALSE(Version::is_at_least(-1, -1, -1));
    ASSERT_TRUE(Version::is_at_least(current_major, current_minor, current_patch));
    ASSERT_FALSE(Version::is_at_least(current_major + 1, current_minor, current_patch));
    ASSERT_FALSE(Version::is_at_least(current_major, current_minor + 1, current_patch));
    ASSERT_FALSE(Version::is_at_least(current_major, current_minor, current_patch + 1));
}

} // end anonymous namespace
