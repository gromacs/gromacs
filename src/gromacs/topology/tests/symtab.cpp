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
 * Tests for legacy symbol table
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/topology/symtab.h"

#include <memory>

#include <gtest/gtest.h>

#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"

namespace gmx
{

namespace test
{

namespace
{

class SymtabTest : public ::testing::Test
{
public:
    SymtabTest() { open_symtab(&symtab_); }
    ~SymtabTest() override
    {
        done_symtab(&symtab_);
        EXPECT_EQ(symtab_.nr, 0);
        EXPECT_EQ(symtab_.symbuf, nullptr);
    }

    //! Get handle to symbol table.
    t_symtab* symtab() { return &symtab_; }
    //! Dump symtab. Similar to pr_symtab function.
    void dumpSymtab();

private:
    //! Get reference checker using lazy initialization
    TestReferenceChecker* checker()
    {
        if (!checker_)
        {
            checker_ = std::make_unique<TestReferenceChecker>(data_.rootChecker());
        }
        return checker_.get();
    }
    //! The symbol table being tested.
    t_symtab symtab_;
    //! Handler for reference data.
    TestReferenceData data_;
    //! Handler for checking reference data.
    std::unique_ptr<TestReferenceChecker> checker_;
};

void SymtabTest::dumpSymtab()
{
    int                      nr     = symtab_.nr;
    t_symbuf*                symbuf = symtab_.symbuf;
    std::vector<std::string> symtabDump;
    int                      pos = 0;
    while (symbuf != nullptr)
    {
        int i;
        for (i = 0; (i < symbuf->bufsize) && (i < nr); i++)
        {
            symtabDump.emplace_back(formatString("Symtab[%d]=\"%s\"", pos++, symbuf->buf[i]));
        }
        nr -= i;
        symbuf = symbuf->next;
    }
    checker()->checkSequence(symtabDump.begin(), symtabDump.end(), "Complete dump of SymbolTable");
}

/*! \brief
 * Helper that compares an input to a handle obtained from symtab lookup.
 *
 * \param[in] symtab Symbol table that contains the entries.
 * \param[in] symbol The entry obtained from placing a string in the symbol table.
 * \param[in] index  Index into symtab corresponding to an entry.
 * \returns Whether to \p symbol and the entry returned by \p index are the same pointer.
 */
bool entriesAreEqual(t_symtab* symtab, char** symbol, int index)
{
    return symbol == get_symtab_handle(symtab, index);
}

/*! \brief
 * Helper function to check internal consistency of symtab lookup.
 *
 * Checks that placing an entry resulted in valid symbol table, and that
 * the index obtained from a call to lookup_symtab returns the correct entry.
 *
 * \param[in] symtab Symbol table that contains the entries.
 * \param[in] symbol The entry obtained from placing a string in the symbol table.
 */
void compareSymtabLookupAndHandle(t_symtab* symtab, char** symbol)
{
    ASSERT_NE(symtab->symbuf, nullptr);
    auto index = lookup_symtab(symtab, symbol);
    EXPECT_TRUE(entriesAreEqual(symtab, symbol, index));
}
/*! \brief
 *  Check that symbols obtained from symtab compare correctly.
 *
 *  Helper function to find out if two entries obtained by a symtab lookup
 *  are equivalent or not, according to testing criteria.
 *
 *  \param[in] symtab Symbol table that contains the entries.
 *  \param[in] firstSymbol Handle into symtab obtained from placing string in symtab.
 *  \param[in] otherSymbol Other handle from obtained from separate string deposit.
 *  \param[in] expectedOutcome If the handles should result in equal entries or not.
 */
void compareDifferentHandles(t_symtab* symtab, char** firstSymbol, char** otherSymbol, bool expectedOutcome)
{
    ASSERT_NE(symtab->symbuf, nullptr);
    auto firstIndex = lookup_symtab(symtab, firstSymbol);
    auto otherIndex = lookup_symtab(symtab, otherSymbol);
    EXPECT_EQ(expectedOutcome, entriesAreEqual(symtab, firstSymbol, otherIndex));
    EXPECT_EQ(expectedOutcome, entriesAreEqual(symtab, otherSymbol, firstIndex));
}

TEST_F(SymtabTest, EmptyOnOpen)
{
    ASSERT_EQ(0, symtab()->nr);
    ASSERT_EQ(nullptr, symtab()->symbuf);
}

TEST_F(SymtabTest, AddSingleEntry)
{
    auto fooSymbol = put_symtab(symtab(), "Foo");
    ASSERT_EQ(1, symtab()->nr);
    compareSymtabLookupAndHandle(symtab(), fooSymbol);
    EXPECT_STREQ("Foo", *fooSymbol);
}

TEST_F(SymtabTest, AddTwoDistinctEntries)
{
    auto fooSymbol = put_symtab(symtab(), "Foo");
    auto barSymbol = put_symtab(symtab(), "Bar");
    ASSERT_EQ(2, symtab()->nr);

    compareSymtabLookupAndHandle(symtab(), fooSymbol);
    compareSymtabLookupAndHandle(symtab(), barSymbol);

    EXPECT_NE(fooSymbol, barSymbol);
    compareDifferentHandles(symtab(), fooSymbol, barSymbol, false);
    EXPECT_STREQ("Foo", *fooSymbol);
    EXPECT_STREQ("Bar", *barSymbol);
}

TEST_F(SymtabTest, TryToAddDuplicates)
{
    auto fooSymbol = put_symtab(symtab(), "Foo");
    auto barSymbol = put_symtab(symtab(), "Bar");
    ASSERT_EQ(2, symtab()->nr);

    compareSymtabLookupAndHandle(symtab(), fooSymbol);
    compareSymtabLookupAndHandle(symtab(), barSymbol);

    EXPECT_NE(fooSymbol, barSymbol);
    compareDifferentHandles(symtab(), fooSymbol, barSymbol, false);
    EXPECT_STREQ("Foo", *fooSymbol);
    EXPECT_STREQ("Bar", *barSymbol);

    // Insert a duplicate element
    auto anotherFooSymbol = put_symtab(symtab(), "Foo");
    ASSERT_EQ(2, symtab()->nr);

    // Check for correct post-conditions
    EXPECT_EQ(fooSymbol, anotherFooSymbol);
    EXPECT_STREQ("Foo", *anotherFooSymbol);
    EXPECT_STREQ("Foo", *fooSymbol);
    EXPECT_STREQ("Bar", *barSymbol);

    // Check for correct behaviours with new and old symbols
    compareDifferentHandles(symtab(), fooSymbol, anotherFooSymbol, true);
    compareDifferentHandles(symtab(), barSymbol, anotherFooSymbol, false);
    compareDifferentHandles(symtab(), fooSymbol, barSymbol, false);
}

TEST_F(SymtabTest, AddLargeNumberOfEntries)
{
    int                 numStringsToAdd = 7; // Larger than c_maxBufSize limit for size of symbuf.
    std::vector<char**> symbolsAdded;
    symbolsAdded.reserve(numStringsToAdd);
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        symbolsAdded.push_back(put_symtab(symtab(), toString(i).c_str()));
    }
    ASSERT_EQ(numStringsToAdd, symtab()->nr);
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        EXPECT_STREQ(toString(i).c_str(), *symbolsAdded[i]);
        compareSymtabLookupAndHandle(symtab(), symbolsAdded[i]);
    }
    // Add something unrelated and check that indices still work afterward.
    auto foobarSymbol = put_symtab(symtab(), "foobar");
    ASSERT_EQ(numStringsToAdd + 1, symtab()->nr);
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        EXPECT_STREQ(toString(i).c_str(), *symbolsAdded[i]);
        compareSymtabLookupAndHandle(symtab(), symbolsAdded[i]);
    }
    compareSymtabLookupAndHandle(symtab(), foobarSymbol);

    // Now dump the symtab to see that we can reproduce it if needed.
    dumpSymtab();
}

TEST_F(SymtabTest, NoDuplicatesInLargeTable)
{
    int halfOfStringsToAdd   = 7; // Larger than c_maxBufSize limit for size of symbuf.
    int totalNumStringsToAdd = 2 * halfOfStringsToAdd;
    std::vector<char**> symbolsAdded;
    symbolsAdded.reserve(halfOfStringsToAdd);
    for (int i = 0; i < halfOfStringsToAdd; ++i)
    {
        symbolsAdded.push_back(put_symtab(symtab(), toString(i).c_str()));
    }
    ASSERT_EQ(halfOfStringsToAdd, symtab()->nr);

    // We now try to mess around in the symtab.
    auto bazSymbol = put_symtab(symtab(), "baz");
    ASSERT_EQ(halfOfStringsToAdd + 1, symtab()->nr);
    compareSymtabLookupAndHandle(symtab(), bazSymbol);

    // Now try to add more symbols, also including those that are already there.
    for (int i = 0; i < totalNumStringsToAdd; i++)
    {
        symbolsAdded.push_back(put_symtab(symtab(), toString(i).c_str()));
    }
    ASSERT_EQ(totalNumStringsToAdd + 1, symtab()->nr);

    //! Check that entries that should be equal are, and new ones are not.
    for (int i = 0; i < halfOfStringsToAdd; i++)
    {
        compareDifferentHandles(symtab(), symbolsAdded[i], symbolsAdded[halfOfStringsToAdd + i], true);
        compareDifferentHandles(symtab(), symbolsAdded[i], symbolsAdded[2 * halfOfStringsToAdd + i], false);
        compareDifferentHandles(symtab(), symbolsAdded[i], bazSymbol, false);
    }
    compareSymtabLookupAndHandle(symtab(), bazSymbol);
    EXPECT_STREQ("baz", *bazSymbol);
    dumpSymtab();
}

} // namespace

} // namespace test

} // namespace gmx
