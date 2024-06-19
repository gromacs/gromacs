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

#include "add_par.h"

#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "hackblock.h"

void add_param(InteractionsOfType* ps, int ai, int aj, gmx::ArrayRef<const real> c, const char* s)
{
    if ((ai < 0) || (aj < 0))
    {
        gmx_fatal(FARGS, "Trying to add impossible atoms: ai=%d, aj=%d", ai, aj);
    }
    std::vector<int>  atoms = { ai, aj };
    std::vector<real> forceParm(c.begin(), c.end());

    ps->interactionTypes.emplace_back(atoms, forceParm, s ? s : "");
}

void add_cmap_param(InteractionsOfType* ps, int ai, int aj, int ak, int al, int am, const char* s)
{
    std::vector<int> atoms = { ai, aj, ak, al, am };
    ps->interactionTypes.emplace_back(atoms, gmx::ArrayRef<const real>{}, s ? s : "");
}

void add_vsite2_param(InteractionsOfType* ps, int ai, int aj, int ak, real c0)
{
    std::vector<int>  atoms     = { ai, aj, ak };
    std::vector<real> forceParm = { c0 };
    ps->interactionTypes.emplace_back(atoms, forceParm);
}

void add_vsite3_param(InteractionsOfType* ps, int ai, int aj, int ak, int al, real c0, real c1)
{
    std::vector<int>  atoms     = { ai, aj, ak, al };
    std::vector<real> forceParm = { c0, c1 };
    ps->interactionTypes.emplace_back(atoms, forceParm);
}

void add_vsite3_atoms(InteractionsOfType* ps, int ai, int aj, int ak, int al, bool bSwapParity)
{
    std::vector<int> atoms = { ai, aj, ak, al };
    ps->interactionTypes.emplace_back(atoms, gmx::ArrayRef<const real>{});

    if (bSwapParity)
    {
        ps->interactionTypes.back().setForceParameter(1, -1);
    }
}

void add_vsite4_atoms(InteractionsOfType* ps, int ai, int aj, int ak, int al, int am)
{
    std::vector<int> atoms = { ai, aj, ak, al, am };
    ps->interactionTypes.emplace_back(atoms, gmx::ArrayRef<const real>{});
}

int search_jtype(const PreprocessResidue& localPpResidue, const char* name, bool bNterm)
{
    int    niter, jmax;
    size_t k, kmax, minstrlen;
    char * rtpname, searchname[12];

    strcpy(searchname, name);

    /* Do a best match comparison */
    /* for protein N-terminus, allow renaming of H1, H2 and H3 to H */
    if (bNterm && (strlen(searchname) == 2) && (searchname[0] == 'H')
        && ((searchname[1] == '1') || (searchname[1] == '2') || (searchname[1] == '3')))
    {
        niter = 2;
    }
    else
    {
        niter = 1;
    }
    kmax = 0;
    jmax = -1;
    for (int iter = 0; (iter < niter && jmax == -1); iter++)
    {
        if (iter == 1)
        {
            /* Try without the hydrogen number in the N-terminus */
            searchname[1] = '\0';
        }
        for (int j = 0; (j < localPpResidue.natom()); j++)
        {
            rtpname = *(localPpResidue.atomname[j]);
            if (gmx_strcasecmp(searchname, rtpname) == 0)
            {
                jmax = j;
                kmax = strlen(searchname);
                break;
            }
            if (iter == niter - 1)
            {
                minstrlen = std::min(strlen(searchname), strlen(rtpname));
                for (k = 0; k < minstrlen; k++)
                {
                    if (searchname[k] != rtpname[k])
                    {
                        break;
                    }
                }
                if (k > kmax)
                {
                    kmax = k;
                    jmax = j;
                }
            }
        }
    }
    if (jmax == -1)
    {
        gmx_fatal(FARGS,
                  "Atom %s not found in rtp database in residue %s",
                  searchname,
                  localPpResidue.resname.c_str());
    }
    if (kmax != strlen(searchname))
    {
        gmx_fatal(FARGS,
                  "Atom %s not found in rtp database in residue %s, "
                  "it looks a bit like %s",
                  searchname,
                  localPpResidue.resname.c_str(),
                  *(localPpResidue.atomname[jmax]));
    }
    return jmax;
}
