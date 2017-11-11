/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2017, by the GROMACS development team, led by
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
 * \brief Tests for gmx::AlignedAllocator and gmx::PageAlignedAllocator.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/alignedallocator.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"

namespace gmx
{
namespace test
{

//! Declare allocator types to test.
using AllocatorTypesToTest = ::testing::Types<AlignedAllocator<real>,
                                              PageAlignedAllocator<real>,
                                              AlignedAllocator<int>,
                                              PageAlignedAllocator<int>,
                                              AlignedAllocator<RVec>,
                                              PageAlignedAllocator<RVec>
                                              >;

TYPED_TEST_CASE(AllocatorTest, AllocatorTypesToTest);

} // namespace
} // namespace

// Includes tests common to all allocation policies.
#include "gromacs/utility/tests/alignedallocator-impl.h"

namespace gmx
{
namespace test
{

TYPED_TEST(AllocatorTest, StatelessAllocatorUsesNoMemory)
{
    using value_type = typename TypeParam::value_type;
    EXPECT_EQ(sizeof(std::vector<value_type>),
              sizeof(std::vector<value_type, TypeParam>));
}

} // namespace
} // namespace
