/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "add_par.h"

#include <cstring>

#include <algorithm>

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void add_param(t_params *ps, int ai, int aj, const real *c, const char *s)
{
    int i;

    if ((ai < 0) || (aj < 0))
    {
        gmx_fatal(FARGS, "Trying to add impossible atoms: ai=%d, aj=%d", ai, aj);
    }
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    if (c != nullptr)
    {
        for (i = 0; (i < MAXFORCEPARAM); i++)
        {
            ps->param.back().c[i] = c[i];
        }
    }
    set_p_string(&(ps->param.back()), s);
}

void add_imp_param(t_params *ps, int ai, int aj, int ak, int al, real c0, real c1,
                   const char *s)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().al() = al;
    ps->param.back().c0() = c0;
    ps->param.back().c1() = c1;
    set_p_string(&(ps->param.back()), s);
}

void add_dih_param(t_params *ps, int ai, int aj, int ak, int al, real c0, real c1,
                   real c2, const char *s)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().al() = al;
    ps->param.back().c0() = c0;
    ps->param.back().c1() = c1;
    ps->param.back().c2() = c2;
    set_p_string(&(ps->param.back()), s);
}

void add_cmap_param(t_params *ps, int ai, int aj, int ak, int al, int am, const char *s)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().al() = al;
    ps->param.back().am() = am;
    set_p_string(&(ps->param.back()), s);
}

void add_vsite2_atoms(t_params *ps, int ai, int aj, int ak)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    set_p_string(&(ps->param.back()), "");
}

void add_vsite2_param(t_params *ps, int ai, int aj, int ak, real c0)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().c0() = c0;
    set_p_string(&(ps->param.back()), "");
}

void add_vsite3_param(t_params *ps, int ai, int aj, int ak, int al,
                      real c0, real c1)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().al() = al;
    ps->param.back().c0() = c0;
    ps->param.back().c1() = c1;
    set_p_string(&(ps->param.back()), "");
}

void add_vsite3_atoms(t_params *ps, int ai, int aj, int ak, int al, bool bSwapParity)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().al() = al;
    if (bSwapParity)
    {
        ps->param.back().c1() = -1;
    }
    set_p_string(&(ps->param.back()), "");
}

void add_vsite4_atoms(t_params *ps, int ai, int aj, int ak, int al, int am)
{
    ps->param.push_back(t_param());
    ps->param.back().ai() = ai;
    ps->param.back().aj() = aj;
    ps->param.back().ak() = ak;
    ps->param.back().al() = al;
    ps->param.back().am() = am;
    set_p_string(&(ps->param.back()), "");
}

int search_jtype(const t_restp &rtp, const char *name, bool bNterm)
{
    int    niter, iter, jmax;
    size_t k, kmax, minstrlen;
    char  *rtpname, searchname[12];

    strcpy(searchname, name);

    /* Do a best match comparison */
    /* for protein N-terminus, allow renaming of H1, H2 and H3 to H */
    if (bNterm && (strlen(searchname) == 2) && (searchname[0] == 'H') &&
        ( (searchname[1] == '1') || (searchname[1] == '2') ||
          (searchname[1] == '3') ) )
    {
        niter = 2;
    }
    else
    {
        niter = 1;
    }
    kmax = 0;
    jmax = -1;
    for (iter = 0; (iter < niter && jmax == -1); iter++)
    {
        if (iter == 1)
        {
            /* Try without the hydrogen number in the N-terminus */
            searchname[1] = '\0';
        }
        for (int j = 0; (j < rtp.natom()); j++)
        {
            rtpname = *(rtp.atomname[j]);
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
        gmx_fatal(FARGS, "Atom %s not found in rtp database in residue %s",
                  searchname, rtp.resname.c_str());
    }
    if (kmax != strlen(searchname))
    {
        gmx_fatal(FARGS, "Atom %s not found in rtp database in residue %s, "
                  "it looks a bit like %s",
                  searchname, rtp.resname.c_str(), *(rtp.atomname[jmax]));
    }
    return jmax;
}
