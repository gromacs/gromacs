/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/*! \file
 * \brief
 * Declares modern and legacy symbol table used to store strings of characters.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 * \ingroup module_topology
 * \inlibraryapi
 */
#ifndef GMX_TOPOLOGY_SYMTAB_H
#define GMX_TOPOLOGY_SYMTAB_H

#include <cstdio>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

struct t_commrec;
struct t_fileio;

namespace gmx
{
class ISerializer;
namespace test
{
class StringTableTest;
} // namespace test
} // namespace gmx

//! Convenience typedef for pair stored in map.
using StringTablePair = std::pair<std::string, int>;
//! Convenience typedef for string reference wrapper.
using StringReference = std::reference_wrapper<const std::string>;

class StringTableBuilder;
class StringTableEntry;
/*! \brief
 * A class to store strings for lookup.
 *
 * We store the strings in a dedicated object to avoid
 * wrong usage of the flat string vector, and forcing people
 * to use an object that can only be constructed from the transitional
 * StringTableBuilder or filled during file IO.
 *
 * Note that strings are stripped of trailing and leading whitespace.
 */
class StringTable
{
public:
    //! Constructor used to generate table object from file reading.
    StringTable(gmx::ISerializer* serializer);
    //! Can move construct.
    StringTable(StringTable&&) = default;
    //! Can move assign.
    StringTable& operator=(StringTable&&) = default;
    //! No copy constructor.
    StringTable(const StringTable&) = delete;
    //! No copy assign.
    StringTable& operator=(const StringTable&) = delete;
    /*! \brief
     *  Access string at \p index.
     *
     *  \returns Entry type that constains both the string and the index,
     *           with the index needed during (de-)serialization.
     *  \throws  On index being out of range.
     */
    StringTableEntry at(gmx::Index index) const;
    //! Bracket operator.
    StringTableEntry operator[](gmx::Index index) const;
    //! Handle file IO.
    void serializeStringTable(gmx::ISerializer* serializer);

    //! Print human readable format of storage.
    void printStringTableStorageToFile(FILE* fp, int indent, const char* title) const;

    /*! \brief
     * Copy data in new datastructure to legacy version.
     *
     * The legacy datastructures need to be already initialized.
     *
     * \param[in] symtab Legacy symbol table to add entries to.
     */
    void copyToLegacySymtab(struct t_symtab* symtab) const;

    friend class StringTableBuilder;

private:
    /*! \brief
     * Private constructor so that only builder can create the final table.
     *
     * \param[in] table A vector of strings to be stored in the table.
     */
    StringTable(const std::vector<std::string>& table) : table_(table) {}

    //! The table is stored as a vector of strings.
    std::vector<std::string> table_;
};

/*! \brief
 * Helper class to access members in StringTable.
 *
 * This class is a wrapper around a string reference to access
 * the actual entry in the table, as well as an index used for
 * serializing the datastructure.
 *
 * This also provides efficient comparison calls between different entries.
 */
class StringTableEntry
{
public:
    //! Copy construct.
    StringTableEntry(const StringTableEntry&) = default;
    //! Move construct.
    StringTableEntry(StringTableEntry&&) noexcept = default;
    //! Copy assign.
    StringTableEntry& operator=(const StringTableEntry&) = default;
    //! Move assign.
    StringTableEntry& operator=(StringTableEntry&&) = default;

    //! Compare entries by indices. Same string should always have same index.
    bool operator==(const StringTableEntry& o) const { return tableIndex_ == o.tableIndex_; }
    //! Unequal comparison.
    bool operator!=(const StringTableEntry& o) const { return !(*this == o); }
    //! Access to underlying view.
    const std::string& operator*() const { return entry_; }
    //! Access to underlying view.
    const std::string* operator->() const { return &entry_.get(); }
    //! Serialize index.
    void serialize(gmx::ISerializer* serializer) const;

    // We only allow construction from the places that are known to create
    // valid objects for us.
    friend StringTableEntry readStringTableEntry(gmx::ISerializer* serializer, const StringTable& table);
    friend class StringTableBuilder;
    friend class StringTable;

private:
    //! Only allow construct with all information present.
    StringTableEntry(StringReference entry, int tableIndex) : entry_(entry), tableIndex_(tableIndex)
    {
    }
    //! The actual string reference that is stored.
    StringReference entry_;
    //! The index into the table.
    int tableIndex_ = -1;
};

/*! \brief
 * De-serialize StringTableEntry using the index into the \p table.
 *
 * \param[in] serializer  The object containing the serialized index.
 * \param[in] table       The storage object holding all strings.
 * \returns The entry into the Table as StringTableEntry.
 */
StringTableEntry readStringTableEntry(gmx::ISerializer* serializer, const StringTable& table);

/*! \libinternal \brief
 * Builds a memory efficient storage for strings of characters.
 *
 * Allows storing strings of characters with unique entries.
 */
class StringTableBuilder
{
public:
    /*! \brief
     * Place new unique string in storage object.
     *
     * Enters new string into the underlying storage or recovers existing entry.
     * \param[in] theString New string to enter.
     * \returns New entry object with reference to string and index into storage.
     *          The reference is only valid while the builder is in use, and becomes
     *          invalidated when generating the StringTable.
     */
    StringTableEntry addString(const std::string& theString);
    //! Find matching entry in storage by name as string.
    int findEntryByName(const std::string& name) const;
    /*! \brief
     * Build the StringTable from the internal map of strings.
     *
     * The unique indices returned from addString() can be used
     * to index into the returned StringTable. Clears the
     * temporary storage so that the StringTableBuilder can be re-used to
     * build a distinct StringTable.
     */
    StringTable build();

private:
    //! Storage object for entries.
    std::unordered_map<std::string, int> map_;
};

// Below this is the legacy code for the old symbol table, only used in
// deprecated datastructures.
/*! \libinternal \brief
 * Legacy symbol table entry as linked list.
 */
struct t_symbuf
{
    //! Number of entries in this item
    int bufsize;
    //! Storage for strings in this item.
    char** buf;
    //! Next item in linked list.
    struct t_symbuf* next;
};

/* \libinternal \brief
 * Legacy symbol table.
 */
struct t_symtab
{
    //! Total number of entries stored.
    int nr;
    //! First item in linked list of storage elements.
    t_symbuf* symbuf;
};

/*
 * This module handles symbol table manipulation. All text strings
 * needed by an application are allocated only once. All references
 * to these text strings use handles returned from the put_symtab()
 * routine. These handles can easily be converted to address independent
 * values by invoking lookup_symtab(). So when writing structures to
 * a file which contains text strings, this value can be written in stead
 * of the text string or its address. This value can easily be converted
 * back to a text string handle by get_symtab_handle().
 */

//! Initialises the symbol table symtab.
void open_symtab(t_symtab* symtab);

/*! \brief
 * Undoes the effect of open_symtab()
 *
 * After invoking this function, no value can be added to the
 * symbol table, only values can be retrieved using get_symtab_handle().
 *
 * Note that this does no work.
 * \param[inout] symtab Symbol table to close.
 */
void close_symtab(t_symtab* symtab);

/*! \brief Returns a deep copy of \c symtab. */
t_symtab* duplicateSymtab(const t_symtab* symtab);

//! Frees the space allocated by the symbol table itself.
void free_symtab(t_symtab* symtab);

//! Frees the space allocated by the symbol table, including all entries in it.
void done_symtab(t_symtab* symtab);

/*! \brief
 * Enters a string into the symbol table.
 *
 * If the string \p name was not present before, a reference to a copy is returned,
 * else a reference to the earlier entered value is returned. Strings are trimmed of spaces.
 *
 * \param[inout] symtab Symbol table to add string to.
 * \param[in] name String to add.
 * \returns Pointer to entry of string in symtab.
 */
char** put_symtab(t_symtab* symtab, const char* name);

/*! \brief
 * Returns unique handle for \p name.
 *
 * Looks up the string pointer \p name in the symbol table and returns the
 * index in it to the matching entry. Gives fatal error if \p name is
 * not found. \p name has to be entered first using put_symtab().
 *
 * \param[in] symtab Symbol table to search.
 * \param[in] name String pointer into \p symtab.
 * \returns Unique index to position in symbol table.
 */
int lookup_symtab(t_symtab* symtab, char** name);

/*! \brief
 * Returns text string corresponding to \p index.
 *
 * \p index needs to be value obtained from call to lookup_symtab().
 * get_symtab_handle() and lookup_symtab() are inverse functions.
 *
 * \param[in] symtab Symbol table to search.
 * \param[in] index  Entry to find in table.
 * \returns String pointer into \p symtab corresponding to the entry.
 */
char** get_symtab_handle(t_symtab* symtab, int index);

/*! \brief
 * Prints human readable form of \p symtab.
 *
 * \param[in] fp File to print to.
 * \param[in] indent Number of spaces to use for indentation.
 * \param[in] title Name for header text.
 * \param[in] symtab Symbol table to print out.
 */
void pr_symtab(FILE* fp, int indent, const char* title, t_symtab* symtab);

#endif
