/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "nrama.h"
#include <math.h>
#include "sysstuff.h"
#include "gromacs/utility/smalloc.h"
#include "typedefs.h"
#include "bondf.h"
#include "gromacs/fileio/futil.h"
#include "gmx_fatal.h"
#include "rmpbc.h"

static const char *pp_pat[] = { "C", "N", "CA", "C", "N" };
#define NPP (sizeof(pp_pat)/sizeof(pp_pat[0]))

static int d_comp(const void *a, const void *b)
{
    t_dih *da, *db;

    da = (t_dih *)a;
    db = (t_dih *)b;

    if (da->ai[1] < db->ai[1])
    {
        return -1;
    }
    else if (da->ai[1] == db->ai[1])
    {
        return (da->ai[2] - db->ai[2]);
    }
    else
    {
        return 1;
    }
}


static void calc_dihs(t_xrama *xr)
{
    int          i, t1, t2, t3;
    rvec         r_ij, r_kj, r_kl, m, n;
    real         sign;
    t_dih       *dd;
    gmx_rmpbc_t  gpbc = NULL;

    gpbc = gmx_rmpbc_init(xr->idef, xr->ePBC, xr->natoms);
    gmx_rmpbc(gpbc, xr->natoms, xr->box, xr->x);
    gmx_rmpbc_done(gpbc);

    for (i = 0; (i < xr->ndih); i++)
    {
        dd      = &(xr->dih[i]);
        dd->ang = dih_angle(xr->x[dd->ai[0]], xr->x[dd->ai[1]],
                            xr->x[dd->ai[2]], xr->x[dd->ai[3]],
                            NULL,
                            r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);
    }
}

gmx_bool new_data(t_xrama *xr)
{
    if (!read_next_x(xr->oenv, xr->traj, &xr->t, xr->x, xr->box))
    {
        return FALSE;
    }

    calc_dihs(xr);

    return TRUE;
}

static int find_atom(const char *find, char ***names, int start, int nr)
{
    int i;

    for (i = start; (i < nr); i++)
    {
        if (strcmp(find, *names[i]) == 0)
        {
            return i;
        }
    }
    return -1;
}

static void add_xr(t_xrama *xr, int ff[5], t_atoms *atoms)
{
    char buf[12];
    int  i;

    srenew(xr->dih, xr->ndih+2);
    for (i = 0; (i < 4); i++)
    {
        xr->dih[xr->ndih].ai[i] = ff[i];
    }
    for (i = 0; (i < 4); i++)
    {
        xr->dih[xr->ndih+1].ai[i] = ff[i+1];
    }
    xr->ndih += 2;

    srenew(xr->pp, xr->npp+1);
    xr->pp[xr->npp].iphi  = xr->ndih-2;
    xr->pp[xr->npp].ipsi  = xr->ndih-1;
    xr->pp[xr->npp].bShow = FALSE;
    sprintf(buf, "%s-%d", *atoms->resinfo[atoms->atom[ff[1]].resind].name,
            atoms->resinfo[atoms->atom[ff[1]].resind].nr);
    xr->pp[xr->npp].label = strdup(buf);
    xr->npp++;
}

static void get_dih(t_xrama *xr, t_atoms *atoms)
{
    int    found, ff[NPP];
    int    i;
    size_t j;

    for (i = 0; (i < atoms->nr); )
    {
        found = i;
        for (j = 0; (j < NPP); j++)
        {
            if ((ff[j] = find_atom(pp_pat[j], atoms->atomname, found, atoms->nr)) == -1)
            {
                break;
            }
            found = ff[j]+1;
        }
        if (j != NPP)
        {
            break;
        }
        add_xr(xr, ff, atoms);
        i = ff[0]+1;
    }
    fprintf(stderr, "Found %d phi-psi combinations\n", xr->npp);
}

static int search_ff(int thisff[NPP], int ndih, int **ff)
{
    int      j, k;
    gmx_bool bFound = FALSE;

    for (j = 0; (j < ndih); j++)
    {
        bFound = TRUE;
        for (k = 1; (k <= 3); k++)
        {
            bFound = bFound && (thisff[k] == ff[j][k]);
        }
        if (bFound)
        {
            if (thisff[0] == -1)
            {
                ff[j][4] = thisff[4];
            }
            else
            {
                ff[j][0] = thisff[0];
            }
            break;
        }
    }
    if (!bFound)
    {
        for (k = 0; (k < 5); k++)
        {
            ff[ndih][k] = thisff[k];
        }
        ndih++;
    }
    return ndih;
}

static void min_max(t_xrama *xr)
{
    int ai, i, j;

    xr->amin = xr->natoms;
    xr->amax = 0;
    for (i = 0; (i < xr->ndih); i++)
    {
        for (j = 0; (j < 4); j++)
        {
            ai = xr->dih[i].ai[j];
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

static void get_dih_props(t_xrama *xr, t_idef *idef, int mult)
{
    int      i, ft, ftype, nra;
    t_iatom *ia;
    t_dih   *dd, key;

    ia = idef->il[F_PDIHS].iatoms;
    for (i = 0; (i < idef->il[F_PDIHS].nr); )
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
        if ((dd = (t_dih *)bsearch(&key, xr->dih, xr->ndih, (size_t)sizeof(key), d_comp))
            != NULL)
        {
            dd->mult = idef->iparams[ft].pdihs.mult;
            dd->phi0 = idef->iparams[ft].pdihs.phiA;
        }

        i  += nra+1;
        ia += nra+1;
    }
    /* Fill in defaults for values not in the topology */
    for (i = 0; (i < xr->ndih); i++)
    {
        if (xr->dih[i].mult == 0)
        {
            fprintf(stderr,
                    "Dihedral around %d,%d not found in topology. Using mult=%d\n",
                    xr->dih[i].ai[1], xr->dih[i].ai[2], mult);
            xr->dih[i].mult = mult;
            xr->dih[i].phi0 = 180;
        }
    }
}



t_topology *init_rama(const output_env_t oenv, const char *infile,
                      const char *topfile, t_xrama *xr, int mult)
{
    t_topology *top;
    int         ePBC;
    real        t;

    top = read_top(topfile, &xr->ePBC);

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
