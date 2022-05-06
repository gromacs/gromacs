/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Tests for atoms datastructures
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/topology/atoms.h"

#include <optional>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

namespace
{

TEST(PdbAtomEntryTest, CanCreateBasicEntry)
{
    PdbAtomEntry testEntry(PdbRecordType::Atom, 1, ' ', "CA");
    EXPECT_EQ(testEntry.type(), PdbRecordType::Atom);
    EXPECT_EQ(testEntry.atomSerialNumber(), 1);
    EXPECT_EQ(testEntry.altloc(), ' ');
    EXPECT_STREQ(testEntry.atomName().c_str(), "CA");
    EXPECT_FALSE(testEntry.occupancy().has_value());
    EXPECT_FALSE(testEntry.bFactor().has_value());
    EXPECT_FALSE(testEntry.anisotropy().has_value());
}

TEST(PdbAtomEntryTest, CanCreateEntryWithOccupAndBfac)
{
    PdbAtomEntry testEntry(PdbRecordType::Atom, 1, ' ', "CA", 0.53, 42.1);
    EXPECT_EQ(testEntry.type(), PdbRecordType::Atom);
    EXPECT_EQ(testEntry.atomSerialNumber(), 1);
    EXPECT_EQ(testEntry.altloc(), ' ');
    EXPECT_STREQ(testEntry.atomName().c_str(), "CA");
    ASSERT_TRUE(testEntry.occupancy().has_value());
    EXPECT_FLOAT_EQ(*testEntry.occupancy(), 0.53);
    ASSERT_TRUE(testEntry.bFactor().has_value());
    EXPECT_FLOAT_EQ(*testEntry.bFactor(), 42.1);
    EXPECT_FALSE(testEntry.anisotropy().has_value());
}

TEST(PdbAtomEntryTest, CanCreateFullEntry)
{
    std::array<real, 6> uij = { 1, 2, 3, 4, 5, 6 };
    PdbAtomEntry        testEntry(PdbRecordType::Atom, 1, ' ', "CA", 0.53, 42.1, uij);
    EXPECT_EQ(testEntry.type(), PdbRecordType::Atom);
    EXPECT_EQ(testEntry.atomSerialNumber(), 1);
    EXPECT_EQ(testEntry.altloc(), ' ');
    EXPECT_STREQ(testEntry.atomName().c_str(), "CA");
    ASSERT_TRUE(testEntry.occupancy().has_value());
    EXPECT_FLOAT_EQ(*testEntry.occupancy(), 0.53);
    ASSERT_TRUE(testEntry.bFactor().has_value());
    EXPECT_FLOAT_EQ(*testEntry.bFactor(), 42.1);
    ASSERT_TRUE(testEntry.anisotropy().has_value());
    auto accessedMatrix = *testEntry.anisotropy();
    EXPECT_THAT(accessedMatrix, ::testing::Pointwise(::testing::Eq(), uij));
}

} // namespace

} // namespace test

} // namespace gmx
