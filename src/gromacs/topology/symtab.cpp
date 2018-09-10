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

#include <cstdio>
#include <cstring>

#include <algorithm>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

#define BUFSIZE         1024
#define TABLESIZE       5

static char *trim_string(const char *s, char *out, int maxlen)
/*
 * Returns a pointer to a static area which contains a copy
 * of s without leading or trailing spaces. Strings are
 * truncated to BUFSIZE positions.
 *
 * TODO This partially duplicates code in trim(), but perhaps
 * replacing symtab with a std::map is a better fix.
 */
{
    int len, i;

    if (strlen(s) > static_cast<size_t>(maxlen-1))
    {
        gmx_fatal(FARGS, "String '%s' (%zu) is longer than buffer (%d).\n",
                  s, strlen(s), maxlen-1);
    }

    for (; (*s) == ' '; s++)
    {
        ;
    }
    for (len = strlen(s); (len > 0); len--)
    {
        if (s[len-1] != ' ')
        {
            break;
        }
    }
    if (len >= BUFSIZE)
    {
        len = BUFSIZE-1;
    }
    for (i = 0; i < len; i++)
    {
        out[i] = *(s++);
    }
    out[i] = 0;
    return out;
}

int lookup_symtab(t_symtab *symtab, char **name)
{
    int       base;
    t_symbuf *symbuf;

    base   = 0;
    symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        const int index = name-symbuf->buf;
        if ( ( index >= 0 ) && ( index < symbuf->bufsize ) )
        {
            return index+base;
        }
        else
        {
            base  += symbuf->bufsize;
            symbuf = symbuf->next;
        }
    }
    gmx_fatal(FARGS, "symtab lookup \"%s\" not found", *name);
}

char **get_symtab_handle(t_symtab *symtab, int name)
{
    t_symbuf *symbuf;

    symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        if (name < symbuf->bufsize)
        {
            return &(symbuf->buf[name]);
        }
        else
        {
            name  -= symbuf->bufsize;
            symbuf = symbuf->next;
        }
    }
    gmx_fatal(FARGS, "symtab get_symtab_handle %d not found", name);
}

static t_symbuf *new_symbuf()
{
    t_symbuf *symbuf;

    snew(symbuf, 1);
    symbuf->bufsize = TABLESIZE;
    snew(symbuf->buf, symbuf->bufsize);
    symbuf->next = nullptr;

    return symbuf;
}

static char **enter_buf(t_symtab *symtab, char *name)
{
    int          i;
    t_symbuf    *symbuf;
    gmx_bool     bCont;

    if (symtab->symbuf == nullptr)
    {
        symtab->symbuf = new_symbuf();
    }

    symbuf = symtab->symbuf;
    do
    {
        for (i = 0; (i < symbuf->bufsize); i++)
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
    }
    while (bCont);

    symbuf->next = new_symbuf();
    symbuf       = symbuf->next;

    symtab->nr++;
    symbuf->buf[0] = gmx_strdup(name);
    return &(symbuf->buf[0]);
}

char **put_symtab(t_symtab *symtab, const char *name)
{
    char buf[1024];

    return enter_buf(symtab, trim_string(name, buf, 1023));
}

void open_symtab(t_symtab *symtab)
{
    symtab->nr     = 0;
    symtab->symbuf = nullptr;
}

void close_symtab(t_symtab gmx_unused *symtab)
{
}

// TODO this will go away when we use a
// std::list<std::vector<std::string>>> for t_symtab.
t_symtab *duplicateSymtab(const t_symtab *symtab)
{
    t_symtab *copySymtab;
    snew(copySymtab, 1);
    open_symtab(copySymtab);
    t_symbuf *symbuf = symtab->symbuf;
    if (symbuf != nullptr)
    {
        snew(copySymtab->symbuf, 1);
    }
    t_symbuf *copySymbuf = copySymtab->symbuf;
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

void done_symtab(t_symtab *symtab)
{
    int       i;
    t_symbuf *symbuf, *freeptr;

    close_symtab(symtab);
    symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        for (i = 0; (i < symbuf->bufsize) && (i < symtab->nr); i++)
        {
            sfree(symbuf->buf[i]);
        }
        symtab->nr -= i;
        sfree(symbuf->buf);
        freeptr = symbuf;
        symbuf  = symbuf->next;
        sfree(freeptr);
    }
    symtab->symbuf = nullptr;
    if (symtab->nr != 0)
    {
        gmx_incons("Freeing symbol table (symtab) structure");
    }
}

void free_symtab(t_symtab *symtab)
{
    t_symbuf *symbuf, *freeptr;

    close_symtab(symtab);
    symbuf = symtab->symbuf;
    while (symbuf != nullptr)
    {
        symtab->nr -= std::min(symbuf->bufsize, symtab->nr);
        freeptr     = symbuf;
        symbuf      = symbuf->next;
        sfree(freeptr);
    }
    symtab->symbuf = nullptr;
    if (symtab->nr != 0)
    {
        gmx_incons("Freeing symbol table (symtab) structure");
    }
}

void pr_symtab(FILE *fp, int indent, const char *title, t_symtab *symtab)
{
    int       i, j, nr;
    t_symbuf *symbuf;

    if (available(fp, symtab, indent, title))
    {
        indent = pr_title_n(fp, indent, title, symtab->nr);
        i      = 0;
        nr     = symtab->nr;
        symbuf = symtab->symbuf;
        while (symbuf != nullptr)
        {
            for (j = 0; (j < symbuf->bufsize) && (j < nr); j++)
            {
                pr_indent(fp, indent);
                (void) fprintf(fp, "%s[%d]=\"%s\"\n", title, i++, symbuf->buf[j]);
            }
            nr    -= j;
            symbuf = symbuf->next;
        }
        if (nr != 0)
        {
            gmx_incons("Printing symbol table (symtab) structure");
        }
    }
}
