/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file definees functions for DD to write PDB files
 * e.g. when reporting problems.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "dump.h"

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/utility/cstringutil.h"

#include "domdec_internal.h"

void write_dd_grid_pdb(const char *fn, int64_t step,
                       gmx_domdec_t *dd, matrix box, gmx_ddbox_t *ddbox)
{
    rvec   grid_s[2], *grid_r = nullptr, cx, r;
    char   fname[STRLEN], buf[22];
    FILE  *out;
    int    a, i, d, z, y, x;
    matrix tric;
    real   vol;

    copy_rvec(dd->comm->cell_x0, grid_s[0]);
    copy_rvec(dd->comm->cell_x1, grid_s[1]);

    if (DDMASTER(dd))
    {
        snew(grid_r, 2*dd->nnodes);
    }

    dd_gather(dd, 2*sizeof(rvec), grid_s, DDMASTER(dd) ? grid_r : nullptr);

    if (DDMASTER(dd))
    {
        for (d = 0; d < DIM; d++)
        {
            for (i = 0; i < DIM; i++)
            {
                if (d == i)
                {
                    tric[d][i] = 1;
                }
                else
                {
                    if (d < ddbox->npbcdim && dd->nc[d] > 1)
                    {
                        tric[d][i] = box[i][d]/box[i][i];
                    }
                    else
                    {
                        tric[d][i] = 0;
                    }
                }
            }
        }
        sprintf(fname, "%s_%s.pdb", fn, gmx_step_str(step, buf));
        out = gmx_fio_fopen(fname, "w");
        gmx_write_pdb_box(out, dd->bScrewPBC ? epbcSCREW : epbcXYZ, box);
        a = 1;
        for (i = 0; i < dd->nnodes; i++)
        {
            vol = dd->nnodes/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
            for (d = 0; d < DIM; d++)
            {
                vol *= grid_r[i*2+1][d] - grid_r[i*2][d];
            }
            for (z = 0; z < 2; z++)
            {
                for (y = 0; y < 2; y++)
                {
                    for (x = 0; x < 2; x++)
                    {
                        cx[XX] = grid_r[i*2+x][XX];
                        cx[YY] = grid_r[i*2+y][YY];
                        cx[ZZ] = grid_r[i*2+z][ZZ];
                        mvmul(tric, cx, r);
                        gmx_fprintf_pdb_atomline(out, epdbATOM, a++, "CA", ' ', "GLY", ' ', i+1, ' ',
                                                 10*r[XX], 10*r[YY], 10*r[ZZ], 1.0, vol, "");
                    }
                }
            }
            for (d = 0; d < DIM; d++)
            {
                for (x = 0; x < 4; x++)
                {
                    switch (d)
                    {
                        case 0: y = 1 + i*8 + 2*x; break;
                        case 1: y = 1 + i*8 + 2*x - (x % 2); break;
                        case 2: y = 1 + i*8 + x; break;
                    }
                    fprintf(out, "%6s%5d%5d\n", "CONECT", y, y+(1<<d));
                }
            }
        }
        gmx_fio_fclose(out);
        sfree(grid_r);
    }
}

void write_dd_pdb(const char *fn, int64_t step, const char *title,
                  const gmx_mtop_t *mtop, const t_commrec *cr,
                  int natoms, const rvec x[], const matrix box)
{
    char          fname[STRLEN], buf[22];
    FILE         *out;
    int           resnr;
    const char   *atomname, *resname;
    gmx_domdec_t *dd;

    dd = cr->dd;
    if (natoms == -1)
    {
        natoms = dd->comm->atomRanges.end(DDAtomRanges::Type::Vsites);
    }

    sprintf(fname, "%s_%s_n%d.pdb", fn, gmx_step_str(step, buf), cr->sim_nodeid);

    out = gmx_fio_fopen(fname, "w");

    fprintf(out, "TITLE     %s\n", title);
    gmx_write_pdb_box(out, dd->bScrewPBC ? epbcSCREW : epbcXYZ, box);
    int molb = 0;
    for (int i = 0; i < natoms; i++)
    {
        int  ii = dd->globalAtomIndices[i];
        mtopGetAtomAndResidueName(mtop, ii, &molb, &atomname, &resnr, &resname, nullptr);
        int  c;
        real b;
        if (i < dd->comm->atomRanges.end(DDAtomRanges::Type::Zones))
        {
            c = 0;
            while (i >= dd->atomGrouping().subRange(0, dd->comm->zones.cg_range[c + 1]).end())
            {
                c++;
            }
            b = c;
        }
        else if (i < dd->comm->atomRanges.end(DDAtomRanges::Type::Vsites))
        {
            b = dd->comm->zones.n;
        }
        else
        {
            b = dd->comm->zones.n + 1;
        }
        gmx_fprintf_pdb_atomline(out, epdbATOM, ii+1, atomname, ' ', resname, ' ', resnr, ' ',
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ], 1.0, b, "");
    }
    fprintf(out, "TER\n");

    gmx_fio_fclose(out);
}
