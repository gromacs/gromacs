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
#include "gmxpre.h"

#include "nrama.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static const char* const pp_pat[] = { "C", "N", "CA", "C", "N" };
#define NPP (sizeof(pp_pat) / sizeof(pp_pat[0]))

static bool d_comp(const t_dih& a, const t_dih& b)
{
    if (a.ai[1] < b.ai[1])
    {
        return true;
    }
    else if (a.ai[1] == b.ai[1])
    {
        return a.ai[2] < b.ai[2];
    }
    else
    {
        return false;
    }
}


static void calc_dihs(t_xrama* xr)
{
    int         i, t1, t2, t3;
    rvec        r_ij, r_kj, r_kl, m, n;
    t_dih*      dd;
    gmx_rmpbc_t gpbc = nullptr;

    gpbc = gmx_rmpbc_init(xr->idef, xr->pbcType, xr->natoms);
    gmx_rmpbc_apply(gpbc, xr->natoms, xr->box, xr->x);
    gmx_rmpbc_done(gpbc);

    for (i = 0; (i < xr->ndih); i++)
    {
        dd      = &(xr->dih[i]);
        dd->ang = dih_angle(xr->x[dd->ai[0]],
                            xr->x[dd->ai[1]],
                            xr->x[dd->ai[2]],
                            xr->x[dd->ai[3]],
                            nullptr,
                            r_ij,
                            r_kj,
                            r_kl,
                            m,
                            n,
                            &t1,
                            &t2,
                            &t3);
    }
}

gmx_bool new_data(t_xrama* xr)
{
    if (!read_next_x(xr->oenv, xr->traj, &xr->t, xr->x, xr->box))
    {
        return FALSE;
    }

    calc_dihs(xr);

    return TRUE;
}

static int find_atom(const char* find, char*** names, int start, int nr)
{
    int i;

    for (i = start; (i < nr); i++)
    {
        if (std::strcmp(find, *names[i]) == 0)
        {
            return i;
        }
    }
    return -1;
}

static void add_xr(t_xrama* xr, const int ff[5], const t_atoms* atoms)
{
    char buf[12];
    int  i;

    srenew(xr->dih, xr->ndih + 2LL);
    for (i = 0; (i < 4); i++)
    {
        xr->dih[xr->ndih].ai[i] = ff[i];
    }
    for (i = 0; (i < 4); i++)
    {
        xr->dih[xr->ndih + 1].ai[i] = ff[i + 1];
    }
    xr->ndih += 2;

    srenew(xr->pp, xr->npp + 1LL);
    xr->pp[xr->npp].iphi  = xr->ndih - 2;
    xr->pp[xr->npp].ipsi  = xr->ndih - 1;
    xr->pp[xr->npp].bShow = FALSE;
    sprintf(buf,
            "%s-%d",
            *atoms->resinfo[atoms->atom[ff[1]].resind].name,
            atoms->resinfo[atoms->atom[ff[1]].resind].nr);
    xr->pp[xr->npp].label = gmx_strdup(buf);
    xr->npp++;
}

static void get_dih(t_xrama* xr, const t_atoms* atoms)
{
    int    found, ff[NPP];
    int    i;
    size_t j;

    for (i = 0; (i < atoms->nr);)
    {
        found = i;
        for (j = 0; (j < NPP); j++)
        {
            if ((ff[j] = find_atom(pp_pat[j], atoms->atomname, found, atoms->nr)) == -1)
            {
                break;
            }
            found = ff[j] + 1;
        }
        if (j != NPP)
        {
            break;
        }
        add_xr(xr, ff, atoms);
        i = ff[0] + 1;
    }
    fprintf(stderr, "Found %d phi-psi combinations\n", xr->npp);
}

static void min_max(t_xrama* xr)
{
    int ai, i, j;

    xr->amin = xr->natoms;
    xr->amax = 0;
    for (i = 0; (i < xr->ndih); i++)
    {
        for (j = 0; (j < 4); j++)
        {
            MSVC_DIAGNOSTIC_IGNORE(28182) // false positive in 2019 (16.5.4)
            ai = xr->dih[i].ai[j];
            MSVC_DIAGNOSTIC_RESET
            if (ai < xr->amin)
            {
                xr->amin = ai;
            }
            else if (ai > xr->amax)
            {
                xr->amax = ai;
            }
        }
    }
}

static void get_dih_props(t_xrama* xr, const t_idef* idef, int mult)
{
    int      i, ft, ftype, nra;
    t_iatom* ia;
    t_dih *  dd, key;

    ia = idef->il[F_PDIHS].iatoms;
    for (i = 0; (i < idef->il[F_PDIHS].nr);)
    {
        ft    = ia[0];
        ftype = idef->functype[ft];
        nra   = interaction_function[ftype].nratoms;

        if (ftype != F_PDIHS)
        {
            gmx_incons("ftype is not a dihedral");
        }

        key.ai[1] = ia[2];
        key.ai[2] = ia[3];
        dd        = std::lower_bound(xr->dih, xr->dih + xr->ndih, key, d_comp);
        if (dd < xr->dih + xr->ndih && !d_comp(key, *dd))
        {
            dd->mult = idef->iparams[ft].pdihs.mult;
            dd->phi0 = idef->iparams[ft].pdihs.phiA;
        }

        i += nra + 1;
        ia += nra + 1LL;
    }
    /* Fill in defaults for values not in the topology */
    for (i = 0; (i < xr->ndih); i++)
    {
        if (xr->dih[i].mult == 0)
        {
            fprintf(stderr,
                    "Dihedral around %d,%d not found in topology. Using mult=%d\n",
                    xr->dih[i].ai[1],
                    xr->dih[i].ai[2],
                    mult);
            xr->dih[i].mult = mult;
            xr->dih[i].phi0 = 180;
        }
    }
}


t_topology* init_rama(gmx_output_env_t* oenv, const char* infile, const char* topfile, t_xrama* xr, int mult)
{
    t_topology* top;
    real        t;

    top = read_top(topfile, &xr->pbcType);

    /*get_dih2(xr,top->idef.functype,&(top->idef.bondeds),&(top->atoms));*/
    get_dih(xr, &(top->atoms));
    get_dih_props(xr, &(top->idef), mult);
    xr->natoms = read_first_x(oenv, &xr->traj, infile, &t, &(xr->x), xr->box);
    xr->idef   = &(top->idef);
    xr->oenv   = oenv;

    min_max(xr);
    calc_dihs(xr);

    return top;
}
