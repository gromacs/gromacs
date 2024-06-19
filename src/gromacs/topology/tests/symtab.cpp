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
 * Tests for legacy symbol table and replacement.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/topology/symtab.h"

#include <cstdio>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

class StringTableTest : public ::testing::Test
{
public:
    StringTableTest() {}
    //! Get handle to symbol table.
    StringTableBuilder& builder() { return stringTableBuilder_; }
    /* \brief
     * Check human readable format of the table.
     * \todo change when table writing to users is done also with serializer.
     * \parm[in] table The string table to check the serialization.
     */
    void checkTable(const StringTable& table);

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
    //! Symbol table for testing purposes.
    StringTableBuilder stringTableBuilder_;
    //! Handler for reference data.
    TestReferenceData data_;
    //! Handler for checking reference data.
    std::unique_ptr<TestReferenceChecker> checker_;
};

void StringTableTest::checkTable(const StringTable& table)
{
    TestFileManager files;
    std::string     filename(files.getTemporaryFilePath("table.txt").string());
    FILE*           fp = fopen(filename.c_str(), "w");
    table.printStringTableStorageToFile(fp, 4, "Test title");
    fclose(fp);
    const std::string text = TextReader::readFileToString(filename);
    checker()->checkTextBlock(text, "Output");
}

/*! \brief
 *  Check that symbols obtained from symtab compare correctly.
 *
 *  Helper function to find out if two entries obtained by a symtab lookup
 *  are equivalent or not, according to testing criteria.
 *  Checks that the indices match before finalizing storage.
 *
 *  \param[in] builder     StringTableBuilder table that contains the entries to validate.
 *  \param[in] firstSymbol Handle into \p builder obtained from placing string in \p builder.
 *  \param[in] otherSymbol Other handle from obtained from separate string deposit.
 *  \param[in] expectedOutcome If the handles should result in equal entries or not.
 */
static void compareDifferentIndices(const StringTableBuilder& builder,
                                    const StringTableEntry&   firstSymbol,
                                    const StringTableEntry&   otherSymbol,
                                    bool                      expectedOutcome)
{
    EXPECT_EQ(expectedOutcome, (firstSymbol == otherSymbol));
    auto firstIndex = builder.findEntryByName(*firstSymbol);
    auto otherIndex = builder.findEntryByName(*otherSymbol);
    EXPECT_EQ(expectedOutcome, (firstIndex == otherIndex))
            << "Expected was " << expectedOutcome << " firstIndex is " << firstIndex
            << " otherIndex is " << otherIndex;
}

/*! \brief
 * Helper to obtain the integer index from an entry.
 *
 * As the index is only used during (de-) serialization, use this machinery
 * to obtain it.
 *
 * \param[in] symbol Single StringTableEntry to obtain the index of.
 * \returns Integer index to be used to obtain value from StringTable.
 */
static int readIndexFromSerializer(const StringTableEntry& symbol)
{
    gmx::InMemorySerializer writer;
    symbol.serialize(&writer);
    auto                      buffer = writer.finishAndGetBuffer();
    gmx::InMemoryDeserializer reader(buffer, false);
    int                       index = 0;
    reader.doInt(&index);
    return index;
}

/*! \brief
 * Helper function to check that a string in matches when looked up in non finalized table.
 *
 * Checks that a string looked up by using the index in the symbol table matches
 * the string stored in the wrapper object obtained by entering a string.
 *
 * \param[in] symtab Symbol table that contains the entries.
 * \param[in] symbol The entry obtained from placing a string in the symbol table.
 * \param[in] string The string the entry should match.
 */
static void stringMatches(const StringTable& symtab, const StringTableEntry& symbol, const char* string)
{
    int  index          = readIndexFromSerializer(symbol);
    auto entryFromIndex = symtab.at(index);

    EXPECT_EQ(*entryFromIndex, string)
            << "Index is " << index << " Entry from index is " << entryFromIndex->c_str();
}


TEST_F(StringTableTest, AddSingleEntry)
{
    builder().addString("foo");
    StringTable table = builder().build();
    checkTable(table);
}

TEST_F(StringTableTest, CanAccessWithAt)
{
    builder().addString("foo");
    StringTable table = builder().build();
    EXPECT_NO_THROW(table.at(0));
    checkTable(table);
}

TEST_F(StringTableTest, CanAccessWithBracket)
{
    builder().addString("foo");
    StringTable table = builder().build();
    checkTable(table);
    auto entry = table[0];
    EXPECT_EQ(*entry, "foo");
}

TEST_F(StringTableTest, ThrowsOutOfRange)
{
    builder().addString("foo");
    StringTable table = builder().build();
    EXPECT_THROW(table.at(1), InternalError);
    checkTable(table);
}

TEST_F(StringTableTest, StringCompareIsCorrect)
{
    auto        fooSymbol = builder().addString("foo");
    StringTable table     = builder().build();
    stringMatches(table, fooSymbol, "foo");
    checkTable(table);
}

TEST_F(StringTableTest, AddTwoDistinctEntries)
{
    auto fooSymbol = builder().addString("foo");
    auto barSymbol = builder().addString("Bar");

    EXPECT_FALSE(fooSymbol == barSymbol);
    compareDifferentIndices(builder(), fooSymbol, barSymbol, false);
    EXPECT_TRUE("foo" == *fooSymbol);
    EXPECT_TRUE("Bar" == *barSymbol);
    auto table = builder().build();
    stringMatches(table, fooSymbol, "foo");
    stringMatches(table, barSymbol, "Bar");
    checkTable(table);
}

TEST_F(StringTableTest, TryToAddDuplicates)
{
    auto fooSymbol = builder().addString("foo");
    auto barSymbol = builder().addString("Bar");

    EXPECT_FALSE(fooSymbol == barSymbol);
    EXPECT_FALSE(fooSymbol->empty());
    compareDifferentIndices(builder(), fooSymbol, barSymbol, false);
    EXPECT_TRUE("foo" == *fooSymbol);
    EXPECT_TRUE("Bar" == *barSymbol);

    // Insert a duplicate element
    auto anotherFooSymbol = builder().addString("foo");
    // Insert element with different case
    auto capitalFooSymbol = builder().addString("Foo");

    // Check that no duplicate is made
    EXPECT_TRUE(fooSymbol == anotherFooSymbol);
    // Check case sensitivity
    EXPECT_FALSE(fooSymbol == capitalFooSymbol);

    // Check that underlying representation is same
    EXPECT_TRUE("foo" == *anotherFooSymbol);
    EXPECT_TRUE("foo" == *fooSymbol);
    EXPECT_FALSE(*fooSymbol == *capitalFooSymbol);
    EXPECT_TRUE("Bar" == *barSymbol);

    // Check for correct behaviours with new and old symbols
    compareDifferentIndices(builder(), fooSymbol, anotherFooSymbol, true);
    compareDifferentIndices(builder(), barSymbol, anotherFooSymbol, false);
    compareDifferentIndices(builder(), fooSymbol, barSymbol, false);
    compareDifferentIndices(builder(), fooSymbol, capitalFooSymbol, false);
    auto table = builder().build();
    checkTable(table);
}

TEST_F(StringTableTest, AddLargeNumberOfEntries)
{
    int                           numStringsToAdd = 7; // Random number of strings.
    std::vector<StringTableEntry> symbolsAdded;
    symbolsAdded.reserve(numStringsToAdd);
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        symbolsAdded.push_back(builder().addString(toString(i)));
    }
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        EXPECT_TRUE(toString(i) == *symbolsAdded[i]) << "index is " << i;
    }
    // Add something unrelated and check that indices still work afterward.
    builder().addString("foobar");
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        EXPECT_TRUE(toString(i) == *symbolsAdded[i]) << "index is " << i;
    }
    auto table = builder().build();
    for (int i = 0; i < numStringsToAdd; ++i)
    {
        stringMatches(table, symbolsAdded[i], toString(i).c_str());
    }
    checkTable(table);
}

TEST_F(StringTableTest, NoDuplicatesInLargeTable)
{
    int                           halfOfStringsToAdd   = 7; // Random number of strings.
    int                           totalNumStringsToAdd = 2 * halfOfStringsToAdd;
    std::vector<StringTableEntry> symbolsAdded;
    symbolsAdded.reserve(halfOfStringsToAdd);
    for (int i = 0; i < halfOfStringsToAdd; ++i)
    {
        symbolsAdded.push_back(builder().addString(toString(i)));
    }

    // We now try to mess around in the symtab.
    auto bazSymbol = builder().addString("baz");

    // Now try to add more symbols, also including those that are already there.
    for (int i = 0; i < totalNumStringsToAdd; i++)
    {
        symbolsAdded.push_back(builder().addString(toString(i)));
    }

    //! Check that entries that should be equal are, and new ones are not.
    for (int i = 0; i < halfOfStringsToAdd; i++)
    {
        compareDifferentIndices(builder(), symbolsAdded[i], symbolsAdded[halfOfStringsToAdd + i], true);
        compareDifferentIndices(builder(), symbolsAdded[i], symbolsAdded[2 * halfOfStringsToAdd + i], false);
        compareDifferentIndices(builder(), symbolsAdded[i], bazSymbol, false);
    }
    EXPECT_TRUE("baz" == *bazSymbol);
    symbolsAdded.emplace_back(bazSymbol);
    auto table = builder().build();
    checkTable(table);
}


TEST_F(StringTableTest, CanWriteToBuffer)
{
    builder().addString("foo");
    builder().addString("bar");
    builder().addString("baz");
    auto               finalTable = builder().build();
    InMemorySerializer writer;
    finalTable.serializeStringTable(&writer);

    auto buffer = writer.finishAndGetBuffer();
    EXPECT_EQ(buffer.size(), 37); // 4 (size) + 3*(8 (string size) + 3*1 (char size) )
}

TEST_F(StringTableTest, Roundtrip)
{
    // First generate a buffer from a string table
    builder().addString("foo");
    builder().addString("bar");
    builder().addString("baz");
    auto               finalTable = builder().build();
    InMemorySerializer writer;
    finalTable.serializeStringTable(&writer);

    auto buffer = writer.finishAndGetBuffer();
    EXPECT_EQ(buffer.size(), 37); // 4 (size) + 3*(8 (string size) + 3*1 (char size) )

    // Now try to make a new table from it.
    InMemoryDeserializer reader(buffer, false);
    StringTable          readInTable(&reader);
    EXPECT_EQ(*(finalTable.at(0)), *(readInTable.at(0)));
    EXPECT_EQ(*(finalTable.at(1)), *(readInTable.at(1)));
    EXPECT_EQ(*(finalTable.at(2)), *(readInTable.at(2)));
}

TEST_F(StringTableTest, RoundtripWithCorrectStringIndices)
{
    std::vector<StringTableEntry> testEntries;
    // First generate a buffer from a string table
    testEntries.emplace_back(builder().addString("foo"));
    testEntries.emplace_back(builder().addString("bar"));
    testEntries.emplace_back(builder().addString("baz"));
    auto               finalTable = builder().build();
    InMemorySerializer writer;
    finalTable.serializeStringTable(&writer);
    for (const auto& stringEntry : testEntries)
    {
        stringEntry.serialize(&writer);
    }

    auto buffer = writer.finishAndGetBuffer();
    EXPECT_EQ(buffer.size(), 49); // 4 (size) + 3*(8 (string size) + 3*1 (char size) + 3*4 (int size))

    // Now try to make a new table from it.
    InMemoryDeserializer          reader(buffer, false);
    StringTable                   readInTable(&reader);
    std::vector<StringTableEntry> deserializedEntries;
    for (Index gmx_unused i = 0; i < gmx::ssize(testEntries); i++)
    {
        deserializedEntries.emplace_back(readStringTableEntry(&reader, readInTable));
    }
    EXPECT_EQ(*(finalTable.at(0)), *(deserializedEntries[0]));
    EXPECT_EQ(*(finalTable.at(1)), *(deserializedEntries[1]));
    EXPECT_EQ(*(finalTable.at(2)), *(deserializedEntries[2]));
}

TEST_F(StringTableTest, CanCopyToLegacyTable)
{
    auto fooSymbol = builder().addString("foo");
    auto barSymbol = builder().addString("Bar");

    StringTable finalTable = builder().build();

    t_symtab legacySymtab;
    open_symtab(&legacySymtab);
    finalTable.copyToLegacySymtab(&legacySymtab);
    int fooEntryIndex = readIndexFromSerializer(fooSymbol);
    int barEntryIndex = readIndexFromSerializer(barSymbol);
    EXPECT_STREQ(finalTable.at(fooEntryIndex)->c_str(), *get_symtab_handle(&legacySymtab, fooEntryIndex));
    EXPECT_STREQ(finalTable.at(barEntryIndex)->c_str(), *get_symtab_handle(&legacySymtab, barEntryIndex));
    done_symtab(&legacySymtab);
}

namespace
{

class LegacySymtabTest : public ::testing::Test
{
public:
    LegacySymtabTest() { open_symtab(&symtab_); }
    ~LegacySymtabTest() override
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

void LegacySymtabTest::dumpSymtab()
{
    int                      nr     = symtab_.nr;
    t_symbuf*                symbuf = symtab_.symbuf;
    std::vector<std::string> symtabDump;
    int                      pos = 0;
    while (symbuf != nullptr)
    {
        int i = 0;
        for (; (i < symbuf->bufsize) && (i < nr); i++)
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

TEST_F(LegacySymtabTest, EmptyOnOpen)
{
    ASSERT_EQ(0, symtab()->nr);
    ASSERT_EQ(nullptr, symtab()->symbuf);
}

TEST_F(LegacySymtabTest, AddSingleEntry)
{
    auto* fooSymbol = put_symtab(symtab(), "Foo");
    ASSERT_EQ(1, symtab()->nr);
    compareSymtabLookupAndHandle(symtab(), fooSymbol);
    EXPECT_STREQ("Foo", *fooSymbol);
}

TEST_F(LegacySymtabTest, AddTwoDistinctEntries)
{
    auto* fooSymbol = put_symtab(symtab(), "Foo");
    auto* barSymbol = put_symtab(symtab(), "Bar");
    ASSERT_EQ(2, symtab()->nr);

    compareSymtabLookupAndHandle(symtab(), fooSymbol);
    compareSymtabLookupAndHandle(symtab(), barSymbol);

    EXPECT_NE(fooSymbol, barSymbol);
    compareDifferentHandles(symtab(), fooSymbol, barSymbol, false);
    EXPECT_STREQ("Foo", *fooSymbol);
    EXPECT_STREQ("Bar", *barSymbol);
}

TEST_F(LegacySymtabTest, TryToAddDuplicates)
{
    auto* fooSymbol = put_symtab(symtab(), "Foo");
    auto* barSymbol = put_symtab(symtab(), "Bar");
    ASSERT_EQ(2, symtab()->nr);

    compareSymtabLookupAndHandle(symtab(), fooSymbol);
    compareSymtabLookupAndHandle(symtab(), barSymbol);

    EXPECT_NE(fooSymbol, barSymbol);
    compareDifferentHandles(symtab(), fooSymbol, barSymbol, false);
    EXPECT_STREQ("Foo", *fooSymbol);
    EXPECT_STREQ("Bar", *barSymbol);

    // Insert a duplicate element
    auto* anotherFooSymbol = put_symtab(symtab(), "Foo");
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

TEST_F(LegacySymtabTest, AddLargeNumberOfEntries)
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
    auto* foobarSymbol = put_symtab(symtab(), "foobar");
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

TEST_F(LegacySymtabTest, NoDuplicatesInLargeTable)
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
    auto* bazSymbol = put_symtab(symtab(), "baz");
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
