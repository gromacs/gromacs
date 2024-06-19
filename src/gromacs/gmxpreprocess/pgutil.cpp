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
/* This file is completely threadsafe - keep it that way! */

#include "gmxpre.h"

#include "pgutil.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>

#include "gromacs/topology/atoms.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/snprintf.h"

#define BUFSIZE 1024
static void atom_not_found(int         fatal_errno,
                           const char* file,
                           int         line,
                           const char* atomname,
                           int         resind,
                           const char* resname,
                           const char* bondtype,
                           bool        bAllowMissing)
{
    char message_buffer[BUFSIZE];
    if (strcmp(bondtype, "check") != 0)
    {
        if (0 != strcmp(bondtype, "atom"))
        {
            snprintf(message_buffer,
                     1024,
                     "Residue %d named %s of a molecule in the input file was mapped\n"
                     "to an entry in the topology database, but the atom %s used in\n"
                     "an interaction of type %s in that entry is not found in the\n"
                     "input file. Perhaps your atom and/or residue naming needs to be\n"
                     "fixed.\n",
                     resind + 1,
                     resname,
                     atomname,
                     bondtype);
        }
        else
        {
            snprintf(message_buffer,
                     1024,
                     "Residue %d named %s of a molecule in the input file was mapped\n"
                     "to an entry in the topology database, but the atom %s used in\n"
                     "that entry is not found in the input file. Perhaps your atom\n"
                     "and/or residue naming needs to be fixed.\n",
                     resind + 1,
                     resname,
                     atomname);
        }
        if (bAllowMissing)
        {
            gmx_warning("WARNING: %s", message_buffer);
        }
        else
        {
            gmx_fatal(fatal_errno, file, line, "%s", message_buffer);
        }
    }
}

int search_atom(const char*              type,
                int                      start,
                const t_atoms*           atoms,
                const char*              bondtype,
                bool                     bAllowMissing,
                gmx::ArrayRef<const int> cyclicBondsIndex)
{
    int                                i, resind = -1;
    bool                               bPrevious, bNext, bOverring;
    int                                natoms = atoms->nr;
    t_atom*                            at     = atoms->atom;
    char** const*                      anm    = atoms->atomname;
    gmx::ArrayRef<const int>::iterator cyclicBondsIterator;

    bPrevious = (strchr(type, '-') != nullptr);
    bNext     = (strchr(type, '+') != nullptr);

    if (!bPrevious)
    {
        resind = at[start].resind;
        if (bNext)
        {
            /* The next residue */
            type++;
            bOverring = !cyclicBondsIndex.empty()
                        && (cyclicBondsIterator =
                                    std::find(cyclicBondsIndex.begin(), cyclicBondsIndex.end(), resind))
                                   != cyclicBondsIndex.end();
            if (bOverring && ((cyclicBondsIterator - cyclicBondsIndex.begin()) & 1))
            {
                resind = *(--cyclicBondsIterator);
                return search_res_atom(type, resind, atoms, bondtype, false);
            }
            else
            {
                while ((start < natoms) && (at[start].resind == resind))
                {
                    start++;
                }
                if (start < natoms)
                {
                    resind = at[start].resind;
                }
            }
        }

        for (i = start; (i < natoms) && (bNext || (at[i].resind == resind)); i++)
        {
            if (anm[i] && gmx_strcasecmp(type, *(anm[i])) == 0)
            {
                return i;
            }
        }
        if (!(bNext && at[start].resind == at[natoms - 1].resind))
        {
            atom_not_found(FARGS, type, at[start].resind, *atoms->resinfo[resind].name, bondtype, bAllowMissing);
        }
    }
    else
    {
        /* The previous residue */
        type++;
        resind    = at[start].resind;
        bOverring = !cyclicBondsIndex.empty()
                    && (cyclicBondsIterator =
                                std::find(cyclicBondsIndex.begin(), cyclicBondsIndex.end(), resind))
                               != cyclicBondsIndex.end();

        if (bOverring && !((cyclicBondsIterator - cyclicBondsIndex.begin()) & 1))
        {
            resind = *(++cyclicBondsIterator);
            return search_res_atom(type, resind, atoms, bondtype, false);
        }
        else
        {
            while ((start >= 0) && (at[start].resind == resind))
            {
                start--;
            }
            if (start >= 0)
            {
                resind = at[start].resind;
                start++;
            }
        }
        for (i = start - 1; (i >= 0) && (at[i].resind == resind); i--)
        {
            if (gmx_strcasecmp(type, *(anm[i])) == 0)
            {
                return i;
            }
        }
        if (start > 0)
        {
            atom_not_found(FARGS, type, at[start].resind, *atoms->resinfo[resind].name, bondtype, bAllowMissing);
        }
    }
    return -1;
}

int search_res_atom(const char* type, int resind, const t_atoms* atoms, const char* bondtype, bool bAllowMissing)
{
    int i;

    for (i = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].resind == resind)
        {
            return search_atom(type, i, atoms, bondtype, bAllowMissing, gmx::ArrayRef<const int>());
        }
    }

    return -1;
}
