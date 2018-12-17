/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "symtab.h"

#include <string>
#include <algorithm>

#include <cstdio>
#include <cstring>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"
#include "gromacs/utility/exceptions.h"

/*! \brief
 * Returns a string without leading or trailing  whitespaces truncated to BUFSIZE positions.
 *
 * \param[in] s Reference to input string.
 * \returns Processed output string.
 */
static std::string trim_string(const std::string &s)
{
    const std::string &chars     = "\t\n\v\f\r ";
    std::string        newString = s;
    newString.erase(0, newString.find_first_not_of(chars));
    newString.erase(newString.find_last_not_of(chars) + 1);

    return newString;
}
const
SymbolPtr lookup_symtab(const SymbolTable &symtab, const std::string &name)
{

    const SymbolTable::const_iterator search = symtab.find(name);
    if (search != symtab.end())
    {
        return search;
    }
    else
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("symtab lookup \"%s\" not found", name.c_str())));
    }
}
const
SymbolPtr get_symtab_handle(SymbolTable *symtab, const int &name)
{
    for (SymbolTable::iterator it = symtab->begin(); it != symtab->end(); ++it)
    {
        if (it->second == name)
        {
            return it;
        }
    }

    GMX_THROW(gmx::InternalError(gmx::formatString("symtab geSymbolTable_handle %d not found", name)));
}

static SymbolPtr enter_buf(SymbolTable *symtab, const std::string &name)
{
    int         lastValue = symtab->size();

    const auto &entry = symtab->insert(std::pair<std::string, int>(name, lastValue + 1));

    return entry.first;
}
const
SymbolPtr put_symtab(SymbolTable *symtab, const std::string &name)
{
    return enter_buf(symtab, trim_string(name));
}

void open_symtab(SymbolTable *symtab)
{
    symtab->clear();
}

void close_symtab(SymbolTable gmx_unused *symtab)
{
}

void done_symtab(SymbolTable *symtab)
{
    open_symtab(symtab);
}

void free_symtab(SymbolTable *symtab)
{
    open_symtab(symtab);
}

void pr_symtab(FILE *fp, int indent, const char *title, const SymbolTable &symtab)
{
    indent = pr_title_n(fp, indent, title, symtab.size());
    int i      = 0;
    for (const auto &entry : symtab)
    {
        pr_indent(fp, indent);
        (void) fprintf(fp, "%s[%d]=\"%s\"\n", title, i++, entry.first.c_str());
    }
}
