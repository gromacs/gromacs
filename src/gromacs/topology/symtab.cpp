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
/*! \internal \file
 * \brief
 * Implements new and legacy symbol table routines.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include "gromacs/topology/symtab.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/txtdump.h"

StringTableEntry StringTableBuilder::addString(const std::string& theString)
{
    int         size     = map_.size();
    std::string stripped = gmx::stripString(theString);

    const auto foundEntry = map_.insert(StringTablePair(stripped, size));
    return StringTableEntry(foundEntry.first->first, foundEntry.first->second);
}

int StringTableBuilder::findEntryByName(const std::string& name) const
{
    auto foundEntry = map_.find(name);
    if (foundEntry != map_.end())
    {
        return foundEntry->second;
    }
    else
    {
        GMX_THROW(gmx::InternalError(
                gmx::formatString("Could not find string \"%s\" in SymbolTable", name.c_str())));
    }
}

StringTable StringTableBuilder::build()
{
    std::vector<std::string> table(map_.size());
    for (const auto& entry : map_)
    {
        table[entry.second] = entry.first;
    }
    map_.clear();
    return StringTable(table);
}

void StringTable::printStringTableStorageToFile(FILE* fp, int indent, const char* title) const
{
    indent = pr_title_n(fp, indent, title, table_.size());
    int i  = 0;
    for (const auto& entry : table_)
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s[%d]=\"%s\"\n", title, i++, entry.c_str());
    }
}

StringTable::StringTable(gmx::ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(serializer->reading(), "Can not use writing serializer to read string table");
    int nr = 0;
    serializer->doInt(&nr);
    table_.resize(nr);
    for (auto& entry : table_)
    {
        serializer->doString(&entry);
    }
}

void StringTable::serializeStringTable(gmx::ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(!serializer->reading(),
                       "Can not use reading serializer to write string table");
    int nr = table_.size();
    serializer->doInt(&nr);
    for (auto& entry : table_)
    {
        serializer->doString(&entry);
    }
}

StringTableEntry StringTable::at(gmx::Index index) const
{
    if (index >= gmx::ssize(table_))
    {
        GMX_THROW(gmx::InternalError("Can't read beyond last entry"));
    }
    return StringTableEntry(table_[index], index);
}

StringTableEntry StringTable::operator[](gmx::Index index) const
{
    GMX_ASSERT(index < gmx::ssize(table_), "Can't read beyond last entry");
    return StringTableEntry(table_[index], index);
}

void StringTableEntry::serialize(gmx::ISerializer* serializer) const
{
    GMX_RELEASE_ASSERT(!serializer->reading(),
                       "Can not use reading serializer to write string index");
    int entry = tableIndex_;
    serializer->doInt(&entry);
}

StringTableEntry readStringTableEntry(gmx::ISerializer* serializer, const StringTable& table)
{
    GMX_RELEASE_ASSERT(serializer->reading(), "Can not use writing serializer to read string index");
    int entry = 0;
    serializer->doInt(&entry);
    return table.at(entry);
}

void StringTable::copyToLegacySymtab(struct t_symtab* symtab) const
{
    for (const auto& entry : table_)
    {
        put_symtab(symtab, entry.c_str());
    }
}

// Old code for legacy data structure starts below.
//! Maximum size of character string in table.
constexpr int c_trimSize = 1024;
//! Maximum number of entries in each element of the linked list.
constexpr int c_maxBufSize = 5;

/*! \brief
 * Remove leading and trailing whitespace from string and enforce maximum length.
 *
 * \param[in]    s      String to trim.
 * \param[inout] out    String to return.
 * \param[in]    maxlen Maximum string length to use.
 * \returns New pruned string.
 */
static char* trim_string(const char* s, char* out, int maxlen)
{
    int len = 0, i = 0;

    if (strlen(s) > static_cast<size_t>(maxlen - 1))
    {
        gmx_fatal(FARGS, "String '%s' (%zu) is longer than buffer (%d).\n", s, strlen(s), maxlen - 1);
    }

    for (; (*s) == ' '; s++) {}
    for (len = strlen(s); (len > 0); len--)
    {
        if (s[len - 1] != ' ')
        {
            break;
        }
    }
    if (len >= c_trimSize)
    {
        len = c_trimSize - 1;
    }
    for (i = 0; i < len; i++)
    {
        out[i] = *(s++);
    }
    out[i] = 0;
    return out;
}

int lookup_symtab(t_symtab* symtab, char** name)
{
    int       base   = 0;
    t_symbuf* symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        const int index = name - symbuf->buf;
        if ((index >= 0) && (index < symbuf->bufsize))
        {
            return index + base;
        }
        else
        {
            base += symbuf->bufsize;
            symbuf = symbuf->next;
        }
    }
    gmx_fatal(FARGS, "symtab lookup \"%s\" not found", *name);
}

char** get_symtab_handle(t_symtab* symtab, int name)
{
    t_symbuf* symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        if (name < symbuf->bufsize)
        {
            return &(symbuf->buf[name]);
        }
        else
        {
            name -= symbuf->bufsize;
            symbuf = symbuf->next;
        }
    }
    gmx_fatal(FARGS, "symtab get_symtab_handle %d not found", name);
}

//! Returns a new initialized entry into the symtab linked list.
static t_symbuf* new_symbuf()
{
    t_symbuf* symbuf = nullptr;

    snew(symbuf, 1);
    symbuf->bufsize = c_maxBufSize;
    snew(symbuf->buf, symbuf->bufsize);
    symbuf->next = nullptr;

    return symbuf;
}

/*! \brief
 * Low level function to enter new string into legacy symtab.
 *
 * \param[inout] symtab Symbol table to add entry to.
 * \param[in]    name   New string to add to symtab.
 * \returns Pointer to new entry in the legacy symbol table, or to existing entry if it already existed.
 */
static char** enter_buf(t_symtab* symtab, char* name)
{
    bool bCont = false;

    if (symtab->symbuf == nullptr)
    {
        symtab->symbuf = new_symbuf();
    }

    t_symbuf* symbuf = symtab->symbuf;
    do
    {
        for (int i = 0; (i < symbuf->bufsize); i++)
        {
            if (symbuf->buf[i] == nullptr)
            {
                symtab->nr++;
                symbuf->buf[i] = gmx_strdup(name);
                return &(symbuf->buf[i]);
            }
            else if (strcmp(symbuf->buf[i], name) == 0)
            {
                return &(symbuf->buf[i]);
            }
        }
        if (symbuf->next != nullptr)
        {
            symbuf = symbuf->next;
            bCont  = TRUE;
        }
        else
        {
            bCont = FALSE;
        }
    } while (bCont);

    symbuf->next = new_symbuf();
    symbuf       = symbuf->next;

    symtab->nr++;
    symbuf->buf[0] = gmx_strdup(name);
    return &(symbuf->buf[0]);
}

char** put_symtab(t_symtab* symtab, const char* name)
{
    char buf[1024];

    return enter_buf(symtab, trim_string(name, buf, 1023));
}

void open_symtab(t_symtab* symtab)
{
    symtab->nr     = 0;
    symtab->symbuf = nullptr;
}

void close_symtab(t_symtab gmx_unused* symtab) {}

// TODO this will go away when we use a
// std::list<std::vector<std::string>>> for t_symtab.
t_symtab* duplicateSymtab(const t_symtab* symtab)
{
    t_symtab* copySymtab = nullptr;
    snew(copySymtab, 1);
    open_symtab(copySymtab);
    t_symbuf* symbuf = symtab->symbuf;
    if (symbuf != nullptr)
    {
        snew(copySymtab->symbuf, 1);
    }
    t_symbuf* copySymbuf = copySymtab->symbuf;
    while (symbuf != nullptr)
    {
        snew(copySymbuf->buf, symbuf->bufsize);
        copySymbuf->bufsize = symbuf->bufsize;
        for (int i = 0; (i < symbuf->bufsize) && (i < symtab->nr); i++)
        {
            if (symbuf->buf[i])
            {
                copySymbuf->buf[i] = gmx_strdup(symbuf->buf[i]);
            }
        }
        symbuf = symbuf->next;
        if (symbuf != nullptr)
        {
            snew(copySymbuf->next, 1);
            copySymbuf = copySymbuf->next;
        }
    }
    copySymtab->nr = symtab->nr;
    return copySymtab;
}

void done_symtab(t_symtab* symtab)
{
    close_symtab(symtab);
    t_symbuf* symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        int i = 0;
        for (; (i < symbuf->bufsize) && (i < symtab->nr); i++)
        {
            sfree(symbuf->buf[i]);
        }
        symtab->nr -= i;
        sfree(symbuf->buf);
        t_symbuf* freeptr = symbuf;
        symbuf            = symbuf->next;
        sfree(freeptr);
    }
    symtab->symbuf = nullptr;
    if (symtab->nr != 0)
    {
        gmx_incons("Freeing symbol table (symtab) structure");
    }
}

void free_symtab(t_symtab* symtab)
{
    close_symtab(symtab);
    t_symbuf* symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        symtab->nr -= std::min(symbuf->bufsize, symtab->nr);
        t_symbuf* freeptr = symbuf;
        symbuf            = symbuf->next;
        sfree(freeptr);
    }
    symtab->symbuf = nullptr;
    if (symtab->nr != 0)
    {
        gmx_incons("Freeing symbol table (symtab) structure");
    }
}

void pr_symtab(FILE* fp, int indent, const char* title, t_symtab* symtab)
{
    if (available(fp, symtab, indent, title))
    {
        indent           = pr_title_n(fp, indent, title, symtab->nr);
        int       i      = 0;
        int       nr     = symtab->nr;
        t_symbuf* symbuf = symtab->symbuf;
        while (symbuf != nullptr)
        {
            int j = 0;
            for (; (j < symbuf->bufsize) && (j < nr); j++)
            {
                pr_indent(fp, indent);
                (void)fprintf(fp, "%s[%d]=\"%s\"\n", title, i++, symbuf->buf[j]);
            }
            nr -= j;
            symbuf = symbuf->next;
        }
        if (nr != 0)
        {
            gmx_incons("Printing symbol table (symtab) structure");
        }
    }
}
