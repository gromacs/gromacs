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
 * Test routines that handle check handling of atom types during preprocessing.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/gpp_atomtype.h"

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/enumerationhelpers.h"

namespace gmx
{
namespace test
{
namespace
{

class PreprocessingAtomTypesTest : public ::testing::Test
{
public:
    PreprocessingAtomTypesTest() : nb_({}, {}) {}

    int addType(const char* name, int bondAtomType, int atomNumber)
    {
        return atypes_.addType(atom_, name, nb_, bondAtomType, atomNumber);
    }

    ~PreprocessingAtomTypesTest() override {}

protected:
    PreprocessingAtomTypes atypes_;
    t_atom                 atom_;
    InteractionOfType      nb_;
};

TEST_F(PreprocessingAtomTypesTest, EmptyOnCreate)
{
    EXPECT_EQ(atypes_.size(), 0);
}

TEST_F(PreprocessingAtomTypesTest, IndexOutOfRangeInvalid)
{
    EXPECT_FALSE(atypes_.isSet(-1));
    EXPECT_FALSE(atypes_.isSet(0));
}

TEST_F(PreprocessingAtomTypesTest, AddTypeWorks)
{
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_TRUE(atypes_.isSet(0));
    EXPECT_EQ(atypes_.size(), 1);
}

TEST_F(PreprocessingAtomTypesTest, AddMultipleTypesWorks)
{
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_TRUE(atypes_.isSet(0));
    EXPECT_EQ(atypes_.size(), 1);
    EXPECT_EQ(addType("Bar", 3, 4), 1);
    EXPECT_TRUE(atypes_.isSet(1));
    EXPECT_EQ(atypes_.size(), 2);
}

TEST_F(PreprocessingAtomTypesTest, CannotAddDuplicateEntry)
{
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_TRUE(atypes_.isSet(0));
    EXPECT_EQ(atypes_.size(), 1);
    EXPECT_EQ(addType("Bar", 3, 4), 1);
    EXPECT_TRUE(atypes_.isSet(1));
    EXPECT_EQ(atypes_.size(), 2);
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_FALSE(atypes_.isSet(3));
    EXPECT_EQ(atypes_.size(), 2);
}

TEST_F(PreprocessingAtomTypesTest, CorrectNameFound)
{
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_EQ(atypes_.atomTypeFromName("Foo"), 0);
}

TEST_F(PreprocessingAtomTypesTest, WrongNameNotFound)
{
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_FALSE(atypes_.atomTypeFromName("Bar").has_value());
}

TEST_F(PreprocessingAtomTypesTest, CorrectNameFromTypeNumber)
{
    EXPECT_EQ(addType("Foo", 1, 2), 0);
    EXPECT_EQ(addType("Bar", 3, 4), 1);
    EXPECT_EQ(atypes_.atomNameFromAtomType(0), "Foo");
    EXPECT_EQ(atypes_.atomNameFromAtomType(1), "Bar");
}

TEST_F(PreprocessingAtomTypesTest, NoNameFromIncorrectTypeNumber)
{
    EXPECT_FALSE(atypes_.atomNameFromAtomType(-1).has_value());
}

} // namespace
} // namespace test
} // namespace gmx
