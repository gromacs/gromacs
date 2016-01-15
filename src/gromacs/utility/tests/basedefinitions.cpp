/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief Tests for base definitions (only alignment attributes for now)
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/basedefinitions.h"

#include <cstdint>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

namespace gmx
{

TEST(BasedefinitionsTest, GmxAlignedDeclaresAlignedVariable)
{
    GMX_ALIGNED(real, 2)  r1;
    GMX_ALIGNED(real, 4)  r2;
    GMX_ALIGNED(real, 8)  r3;

    std::uint64_t addr1 = reinterpret_cast<std::uint64_t>(&r1);
    std::uint64_t addr2 = reinterpret_cast<std::uint64_t>(&r2);
    std::uint64_t addr3 = reinterpret_cast<std::uint64_t>(&r3);

    EXPECT_EQ(0, addr1 % 2);
    EXPECT_EQ(0, addr2 % 4);
    EXPECT_EQ(0, addr3 % 8);

    GMX_ALIGNED(int, 2)   i1;
    GMX_ALIGNED(int, 4)   i2;
    GMX_ALIGNED(int, 8)   i3;

    addr1 = reinterpret_cast<std::uint64_t>(&i1);
    addr2 = reinterpret_cast<std::uint64_t>(&i2);
    addr3 = reinterpret_cast<std::uint64_t>(&i3);

    EXPECT_EQ(0, addr1 % 2);
    EXPECT_EQ(0, addr2 % 4);
    EXPECT_EQ(0, addr3 % 8);
}

}
