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
 * \brief
 * Test routines that handle check handling of bond atom types during preprocessing.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"

#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/enumerationhelpers.h"

class PreprocessingBondAtomTypeTest : public ::testing::Test
{
    public:
        PreprocessingBondAtomTypeTest()
        {
            open_symtab(&symtab);
            bat = init_bond_atomtype();
        }

        void addType(const char *name);

        ~PreprocessingBondAtomTypeTest() override
        {
            done_symtab(&symtab);
            done_bond_atomtype(&bat);
        }
    protected:
        gpp_bond_atomtype *bat;
        t_symtab           symtab;
};

void PreprocessingBondAtomTypeTest::addType(const char *name)
{
    add_bond_atomtype(bat,
                      &symtab,
                      const_cast<char *>(name));
}

// TODO activate those once we have a proper type for gpp_bond_atomtype.
//TEST_F(PreprocessingBondAtomTypeTest, EmptyOnCreate)
//{
//    EXPECT_EQ(bat->nr, 0);
//}
//
//TEST_F(PreprocessingBondAtomTypeTest, AddEntryWorks)
//{
//    addType("Foo");
//    EXPECT_EQ(bat->nr,1);
//}

TEST_F(PreprocessingBondAtomTypeTest, CorrectNameFound)
{
    addType("Foo");
    EXPECT_EQ(get_bond_atomtype_type("Foo", bat), 0);
}

TEST_F(PreprocessingBondAtomTypeTest, WrongNameNotFound)
{
    addType("Foo");
    EXPECT_EQ(get_bond_atomtype_type("Bar", bat), NOTSET);
}

TEST_F(PreprocessingBondAtomTypeTest, CorrectNameFromTypeNumber)
{
    addType("Foo");
    addType("Bar");
    EXPECT_STREQ(get_bond_atomtype_name(0, bat), "Foo");
    EXPECT_STREQ(get_bond_atomtype_name(1, bat), "Bar");
}

/* This tests the current behaviour that duplicate types can be added,
 * even though this is not intended and the calling code right now circumvents
 * this by having a check before a new type is added.
 */
TEST_F(PreprocessingBondAtomTypeTest, CanAddDuplicateType)
{
    addType("Foo");
    addType("Bar");
    addType("Foo");
    EXPECT_STREQ(get_bond_atomtype_name(0, bat), "Foo");
    EXPECT_STREQ(get_bond_atomtype_name(1, bat), "Bar");
    EXPECT_STREQ(get_bond_atomtype_name(2, bat), "Foo");
}

TEST_F(PreprocessingBondAtomTypeTest, NoNameFromIncorrectTypeNumber)
{
    EXPECT_EQ(get_bond_atomtype_name(-1, bat), nullptr);
}
