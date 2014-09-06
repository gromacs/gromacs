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
#include "gmxpre.h"

#include "rmpbc.h"

#include "gromacs/fileio/trx.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    int      natoms;
    t_graph *gr;
} rmpbc_graph_t;

struct gmx_rmpbc {
    t_idef        *idef;
    int            natoms_init;
    int            ePBC;
    int            ngraph;
    rmpbc_graph_t *graph;
};

static t_graph *gmx_rmpbc_get_graph(gmx_rmpbc_t gpbc, int ePBC, int natoms)
{
    int            i;
    rmpbc_graph_t *gr;

    if (ePBC == epbcNONE
        || NULL == gpbc
        || NULL == gpbc->idef
        || gpbc->idef->ntypes <= 0)
    {
        return NULL;
    }

    gr = NULL;
    for (i = 0; i < gpbc->ngraph; i++)
    {
        if (natoms == gpbc->graph[i].natoms)
        {
            gr = &gpbc->graph[i];
        }
    }
    if (gr == NULL)
    {
        /* We'd like to check with the number of atoms in the topology,
         * but we don't have that available.
         * So we check against the number of atoms that gmx_rmpbc_init
         * was called with.
         */
        if (natoms > gpbc->natoms_init)
        {
            gmx_fatal(FARGS, "Structure or trajectory file has more atoms (%d) than the topology (%d)", natoms, gpbc->natoms_init);
        }
        gpbc->ngraph++;
        srenew(gpbc->graph, gpbc->ngraph);
        gr         = &gpbc->graph[gpbc->ngraph-1];
        gr->natoms = natoms;
        gr->gr     = mk_graph(NULL, gpbc->idef, 0, natoms, FALSE, FALSE);
    }

    return gr->gr;
}

gmx_rmpbc_t gmx_rmpbc_init(t_idef *idef, int ePBC, int natoms)
{
    gmx_rmpbc_t gpbc;

    snew(gpbc, 1);

    gpbc->natoms_init = natoms;

    /* This sets pbc when we now it,
     * otherwise we guess it from the instantaneous box in the trajectory.
     */
    gpbc->ePBC = ePBC;

    gpbc->idef = idef;
    if (gpbc->idef->ntypes <= 0)
    {
        fprintf(stderr,
                "\n"
                "WARNING: If there are molecules in the input trajectory file\n"
                "         that are broken across periodic boundaries, they\n"
                "         cannot be made whole (or treated as whole) without\n"
                "         you providing a run input file.\n\n");
    }

    return gpbc;
}

void gmx_rmpbc_done(gmx_rmpbc_t gpbc)
{
    int i;

    if (NULL != gpbc)
    {
        for (i = 0; i < gpbc->ngraph; i++)
        {
            done_graph(gpbc->graph[i].gr);
            sfree(gpbc->graph[i].gr);
        }
        if (gpbc->graph != NULL)
        {
            sfree(gpbc->graph);
        }
        sfree(gpbc);
    }
}

static int gmx_rmpbc_ePBC(gmx_rmpbc_t gpbc, matrix box)
{
    if (NULL != gpbc && gpbc->ePBC >= 0)
    {
        return gpbc->ePBC;
    }
    else
    {
        return guess_ePBC(box);
    }
}

void gmx_rmpbc(gmx_rmpbc_t gpbc, int natoms, matrix box, rvec x[])
{
    int      ePBC;
    t_graph *gr;

    ePBC = gmx_rmpbc_ePBC(gpbc, box);
    gr   = gmx_rmpbc_get_graph(gpbc, ePBC, natoms);
    if (gr != NULL)
    {
        mk_mshift(stdout, gr, ePBC, box, x);
        shift_self(gr, box, x);
    }
}

void gmx_rmpbc_copy(gmx_rmpbc_t gpbc, int natoms, matrix box, rvec x[], rvec x_s[])
{
    int      ePBC;
    t_graph *gr;
    int      i;

    ePBC = gmx_rmpbc_ePBC(gpbc, box);
    gr   = gmx_rmpbc_get_graph(gpbc, ePBC, natoms);
    if (gr != NULL)
    {
        mk_mshift(stdout, gr, ePBC, box, x);
        shift_x(gr, box, x, x_s);
    }
    else
    {
        for (i = 0; i < natoms; i++)
        {
            copy_rvec(x[i], x_s[i]);
        }
    }
}

void gmx_rmpbc_trxfr(gmx_rmpbc_t gpbc, t_trxframe *fr)
{
    int      ePBC;
    t_graph *gr;

    if (fr->bX && fr->bBox)
    {
        ePBC = gmx_rmpbc_ePBC(gpbc, fr->box);
        gr   = gmx_rmpbc_get_graph(gpbc, ePBC, fr->natoms);
        if (gr != NULL)
        {
            mk_mshift(stdout, gr, ePBC, fr->box, fr->x);
            shift_self(gr, fr->box, fr->x);
        }
    }
}

void rm_gropbc(t_atoms *atoms, rvec x[], matrix box)
{
    real dist;
    int  n, m, d;

    /* check periodic boundary */
    for (n = 1; (n < atoms->nr); n++)
    {
        for (m = DIM-1; m >= 0; m--)
        {
            dist = x[n][m]-x[n-1][m];
            if (fabs(dist) > 0.9*box[m][m])
            {
                if (dist >  0)
                {
                    for (d = 0; d <= m; d++)
                    {
                        x[n][d] -= box[m][d];
                    }
                }
                else
                {
                    for (d = 0; d <= m; d++)
                    {
                        x[n][d] += box[m][d];
                    }
                }
            }
        }
    }
}
