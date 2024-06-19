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

#include "gromacs/pbcutil/rmpbc.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <filesystem>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

typedef struct
{
    int      natoms;
    t_graph* gr;
} rmpbc_graph_t;

struct gmx_rmpbc
{
    const InteractionDefinitions* interactionDefinitions;
    const t_idef*                 idef;
    int                           natoms_init;
    PbcType                       pbcType;
    int                           ngraph;
    rmpbc_graph_t*                graph;
};

static t_graph* gmx_rmpbc_get_graph(gmx_rmpbc_t gpbc, PbcType pbcType, int natoms)
{
    int            i;
    rmpbc_graph_t* gr;

    if (pbcType == PbcType::No || nullptr == gpbc
        || (nullptr == gpbc->interactionDefinitions && (nullptr == gpbc->idef || gpbc->idef->ntypes <= 0)))
    {
        return nullptr;
    }

    gr = nullptr;
    for (i = 0; i < gpbc->ngraph; i++)
    {
        if (natoms == gpbc->graph[i].natoms)
        {
            gr = &gpbc->graph[i];
        }
    }
    if (gr == nullptr)
    {
        /* We'd like to check with the number of atoms in the topology,
         * but we don't have that available.
         * So we check against the number of atoms that gmx_rmpbc_init
         * was called with.
         */
        if (natoms > gpbc->natoms_init)
        {
            gmx_fatal(FARGS,
                      "Structure or trajectory file has more atoms (%d) than the topology (%d)",
                      natoms,
                      gpbc->natoms_init);
        }
        gpbc->ngraph++;
        srenew(gpbc->graph, gpbc->ngraph);
        gr         = &gpbc->graph[gpbc->ngraph - 1];
        gr->natoms = natoms;
        if (gpbc->interactionDefinitions)
        {
            gr->gr = mk_graph(nullptr, *gpbc->interactionDefinitions, natoms, FALSE, FALSE);
        }
        else
        {
            gr->gr = mk_graph(nullptr, gpbc->idef, natoms, FALSE, FALSE);
        }
    }

    return gr->gr;
}

gmx_rmpbc_t gmx_rmpbc_init(const InteractionDefinitions& idef, PbcType pbcType, int natoms)
{
    gmx_rmpbc_t gpbc;

    snew(gpbc, 1);

    gpbc->natoms_init = natoms;

    /* This sets pbc when we now it,
     * otherwise we guess it from the instantaneous box in the trajectory.
     */
    gpbc->pbcType = pbcType;

    gpbc->interactionDefinitions = &idef;

    return gpbc;
}

gmx_rmpbc_t gmx_rmpbc_init(const t_idef* idef, PbcType pbcType, int natoms)
{
    gmx_rmpbc_t gpbc;

    snew(gpbc, 1);

    gpbc->natoms_init = natoms;

    /* This sets pbc when we now it,
     * otherwise we guess it from the instantaneous box in the trajectory.
     */
    gpbc->pbcType = pbcType;

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

    if (nullptr != gpbc)
    {
        for (i = 0; i < gpbc->ngraph; i++)
        {
            delete gpbc->graph[i].gr;
        }
        if (gpbc->graph != nullptr)
        {
            sfree(gpbc->graph);
        }
        sfree(gpbc);
    }
}

static PbcType gmx_rmpbc_ePBC(gmx_rmpbc_t gpbc, const matrix box)
{
    if (nullptr != gpbc && gpbc->pbcType != PbcType::Unset)
    {
        return gpbc->pbcType;
    }
    else
    {
        return guessPbcType(box);
    }
}

void gmx_rmpbc_apply(gmx_rmpbc_t gpbc, int natoms, const matrix box, rvec x[])
{
    PbcType  pbcType;
    t_graph* gr;

    pbcType = gmx_rmpbc_ePBC(gpbc, box);
    gr      = gmx_rmpbc_get_graph(gpbc, pbcType, natoms);
    if (gr != nullptr)
    {
        mk_mshift(stdout, gr, pbcType, box, x);
        shift_self(gr, box, x);
    }
}

void gmx_rmpbc_copy(gmx_rmpbc_t gpbc, int natoms, const matrix box, rvec x[], rvec x_s[])
{
    PbcType  pbcType;
    t_graph* gr;
    int      i;

    pbcType = gmx_rmpbc_ePBC(gpbc, box);
    gr      = gmx_rmpbc_get_graph(gpbc, pbcType, natoms);
    if (gr != nullptr)
    {
        mk_mshift(stdout, gr, pbcType, box, x);
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

void gmx_rmpbc_trxfr(gmx_rmpbc_t gpbc, t_trxframe* fr)
{
    PbcType  pbcType;
    t_graph* gr;

    if (fr->bX && fr->bBox)
    {
        pbcType = gmx_rmpbc_ePBC(gpbc, fr->box);
        gr      = gmx_rmpbc_get_graph(gpbc, pbcType, fr->natoms);
        if (gr != nullptr)
        {
            mk_mshift(stdout, gr, pbcType, fr->box, fr->x);
            shift_self(gr, fr->box, fr->x);
        }
    }
}

void rm_gropbc(const t_atoms* atoms, rvec x[], const matrix box)
{
    real dist;
    int  n, m, d;

    /* check periodic boundary */
    for (n = 1; (n < atoms->nr); n++)
    {
        for (m = DIM - 1; m >= 0; m--)
        {
            dist = x[n][m] - x[n - 1][m];
            if (std::abs(dist) > 0.9 * box[m][m])
            {
                if (dist > 0)
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
