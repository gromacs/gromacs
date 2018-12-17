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
/*! \internal \file
 * \brief
 * Tests for symbol table
 *
 * \author
 */
#include "gmxpre.h"

#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"

#include <cstdio>
#include <cstring>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{

namespace test
{

namespace
{

TEST(SymtabTest, EmptyOnOpen)
{
    SymbolTable symtab;

    EXPECT_EQ(0, symtab.size());
    open_symtab(&symtab);
    EXPECT_EQ(0, symtab.size());
}

TEST(SymtabTest, EmptyInT_Topology)
{
    t_topology top;
    EXPECT_EQ(0, top.symtab.size());
    open_symtab(&top.symtab);
    EXPECT_EQ(0, top.symtab.size());
}

TEST(SymtabTest, EmptyInT_TopologyPointer)
{
    t_topology *top = new t_topology;
    EXPECT_EQ(0, top->symtab.size());
    open_symtab(&top->symtab);
    EXPECT_EQ(0, top->symtab.size());
}

TEST(SymtabTest, SymbolPtrInvalidIfNotAssigned)
{
    SymbolPtr symbol;
    EXPECT_FALSE(symbol.isSet());
}

TEST(SymtabTest, CanAddSingleEntry)
{
    SymbolTable symtab;
    SymbolPtr   symbol = put_symtab(&symtab, "Test Entry");
    ASSERT_TRUE(symbol.isSet());
    EXPECT_EQ(1, symtab.size());
    EXPECT_EQ(symbol.number(), 1);
    EXPECT_STREQ(symbol->c_str(), "Test Entry");
}

TEST(SymtabTest, CanAddMultipleEntries)
{
    SymbolTable            symtab;
    // Whitespace gets removed!
    std::vector<SymbolPtr> symbol;
    symbol.push_back(put_symtab(&symtab, "Test Entry"));
    EXPECT_EQ(1, symtab.size());
    ASSERT_TRUE(symbol[0].isSet());
    EXPECT_EQ(symbol[0].number(), 1);
    EXPECT_STREQ(symbol[0]->c_str(), "Test Entry");

    symbol.push_back(put_symtab(&symtab, "Test Entry 2"));
    EXPECT_EQ(2, symtab.size());
    ASSERT_TRUE(symbol[1].isSet());
    EXPECT_EQ(symbol[1].number(), 2);
    EXPECT_STREQ(symbol[1]->c_str(), "Test Entry 2");
}

TEST(SymtabTest, NoDuplicateEntries)
{
    SymbolTable            symtab;

    std::vector<SymbolPtr> symbol;
    symbol.push_back(put_symtab(&symtab, "Test Entry"));
    EXPECT_EQ(1, symtab.size());
    ASSERT_TRUE(symbol[0].isSet());
    EXPECT_EQ(symbol[0].number(), 1);
    EXPECT_STREQ((*symbol[0]).c_str(), "Test Entry");

    symbol.push_back(put_symtab(&symtab, "Test Entry 2"));
    EXPECT_EQ(2, symtab.size());
    ASSERT_TRUE(symbol[1].isSet());
    EXPECT_EQ(symbol[1].number(), 2);
    EXPECT_STREQ((*symbol[1]).c_str(), "Test Entry 2");

    symbol.push_back(put_symtab(&symtab, "Test Entry"));
    EXPECT_EQ(2, symtab.size());
    ASSERT_TRUE(symbol[2].isSet());
    EXPECT_EQ(symbol[2].number(), 1);
    EXPECT_STREQ((*symbol[2]).c_str(), "Test Entry");
}

TEST(SymtabTest, GetEntryByNumber)
{
    SymbolTable            symtab;

    std::vector<SymbolPtr> symbol;
    symbol.push_back(put_symtab(&symtab, "Test Entry"));
    EXPECT_EQ(1, symtab.size());
    ASSERT_TRUE(symbol[0].isSet());
    EXPECT_EQ(symbol[0].number(), 1);
    EXPECT_STREQ((*symbol[0]).c_str(), "Test Entry");

    symbol.push_back(put_symtab(&symtab, "Test Entry 2"));
    EXPECT_EQ(2, symtab.size());
    ASSERT_TRUE(symbol[1].isSet());
    EXPECT_EQ(symbol[1].number(), 2);
    EXPECT_STREQ((*symbol[1]).c_str(), "Test Entry 2");

    symbol.push_back(get_symtab_handle(&symtab, 1));
    ASSERT_TRUE(symbol[2].isSet());
    EXPECT_STREQ("Test Entry", symbol[2]->c_str());
}

TEST(SymtabTest, GetEntryByName)
{
    SymbolTable            symtab;

    std::vector<SymbolPtr> symbol;
    symbol.push_back(put_symtab(&symtab, "Test Entry"));
    EXPECT_EQ(1, symtab.size());
    ASSERT_TRUE(symbol[0].isSet());
    EXPECT_EQ(symbol[0].number(), 1);
    EXPECT_STREQ((*symbol[0]).c_str(), "Test Entry");

    symbol.push_back(put_symtab(&symtab, "Test Entry 2"));
    EXPECT_EQ(2, symtab.size());
    ASSERT_TRUE(symbol[1].isSet());
    EXPECT_EQ(symbol[1].number(), 2);
    EXPECT_STREQ((*symbol[1]).c_str(), "Test Entry 2");

    std::string entry = "Test Entry";
    symbol.push_back(lookup_symtab(symtab, entry));
    ASSERT_TRUE(symbol[2].isSet());
    EXPECT_STREQ("Test Entry", symbol[2]->c_str());
}


} // namespace

} // namespace test

} // namespace gmx
