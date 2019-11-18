/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Tests for gmx::DefaultInitializationAllocator used in std::vector
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/defaultinitializationallocator.h"

#include <gtest/gtest.h>

#include "gromacs/utility/gmxassert.h"


namespace gmx
{
namespace test
{
namespace
{

TEST(DefaultInitializationAllocator, PerformsValueInitialization)
{
    std::vector<int, DefaultInitializationAllocator<int>> v;

    v.resize(1, 2);
    EXPECT_EQ(v[0], 2);
}

TEST(DefaultInitializationAllocator, PerformsNoInitialization)
{
    std::vector<int, DefaultInitializationAllocator<int>> v{ 1, 2, 3 };

    const int* oldData = v.data();
    v.resize(0);
    v.resize(3);
    GMX_RELEASE_ASSERT(v.data() == oldData,
                       "According to the C++ standard std::vector will not reallocate when the "
                       "capacity is sufficient");
    // The allocation is the same, so the default initialization should
    // not have changed the contents
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
}

} // namespace
} // namespace test
} // namespace gmx
