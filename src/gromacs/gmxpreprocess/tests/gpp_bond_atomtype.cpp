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
 * Test routines that handle check handling of bond atom types during preprocessing.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/enumerationhelpers.h"

class PreprocessingBondAtomTypeTest : public ::testing::Test
{
public:
    PreprocessingBondAtomTypeTest() {}

    int addType(const char* name);

    ~PreprocessingBondAtomTypeTest() override {}

protected:
    PreprocessingBondAtomType bat_;
};

int PreprocessingBondAtomTypeTest::addType(const char* name)
{
    return bat_.addBondAtomType(name);
}

TEST_F(PreprocessingBondAtomTypeTest, EmptyOnCreate)
{
    EXPECT_EQ(bat_.size(), 0);
}

TEST_F(PreprocessingBondAtomTypeTest, IndexOutOfRangeInvalid)
{
    EXPECT_FALSE(bat_.isSet(-1));
    EXPECT_FALSE(bat_.isSet(0));
}

TEST_F(PreprocessingBondAtomTypeTest, AddTypeWorks)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_TRUE(bat_.isSet(0));
    EXPECT_EQ(bat_.size(), 1);
}

TEST_F(PreprocessingBondAtomTypeTest, AddMultipleTypesWorks)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_TRUE(bat_.isSet(0));
    EXPECT_EQ(bat_.size(), 1);
    EXPECT_EQ(addType("Bar"), 1);
    EXPECT_TRUE(bat_.isSet(1));
    EXPECT_EQ(bat_.size(), 2);
}

TEST_F(PreprocessingBondAtomTypeTest, CannotAddDuplicateEntry)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_TRUE(bat_.isSet(0));
    EXPECT_EQ(bat_.size(), 1);
    EXPECT_EQ(addType("Bar"), 1);
    EXPECT_TRUE(bat_.isSet(1));
    EXPECT_EQ(bat_.size(), 2);
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_FALSE(bat_.isSet(3));
    EXPECT_EQ(bat_.size(), 2);
}

TEST_F(PreprocessingBondAtomTypeTest, ReturnsCorrectIndexOnDuplicateType)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_TRUE(bat_.isSet(0));
    EXPECT_EQ(bat_.size(), 1);
    EXPECT_EQ(addType("Bar"), 1);
    EXPECT_TRUE(bat_.isSet(1));
    EXPECT_EQ(bat_.size(), 2);
    EXPECT_EQ(addType("BAT"), 2);
    EXPECT_TRUE(bat_.isSet(2));
    EXPECT_EQ(bat_.size(), 3);
    EXPECT_EQ(addType("Bar"), 1);
    EXPECT_FALSE(bat_.isSet(4));
    EXPECT_EQ(bat_.size(), 3);
}

TEST_F(PreprocessingBondAtomTypeTest, CorrectNameFound)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_EQ(bat_.bondAtomTypeFromName("Foo"), 0);
}

TEST_F(PreprocessingBondAtomTypeTest, WrongNameNotFound)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_FALSE(bat_.bondAtomTypeFromName("Bar").has_value());
}

TEST_F(PreprocessingBondAtomTypeTest, CorrectNameFromTypeNumber)
{
    EXPECT_EQ(addType("Foo"), 0);
    EXPECT_EQ(addType("Bar"), 1);
    EXPECT_EQ(bat_.atomNameFromBondAtomType(0), "Foo");
    EXPECT_EQ(bat_.atomNameFromBondAtomType(1), "Bar");
}

TEST_F(PreprocessingBondAtomTypeTest, NoNameFromIncorrectTypeNumber)
{
    EXPECT_FALSE(bat_.atomNameFromBondAtomType(-1).has_value());
}
