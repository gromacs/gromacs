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
/* \internal \file
 *
 * \brief Implements DD cell-size related functions.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Roland Schulz <roland.schulz@intel.com>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "cellsizes.h"

#include "config.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

#include "atomdistribution.h"
#include "domdec_internal.h"
#include "utility.h"

static void set_pme_maxshift(gmx_domdec_t *dd, gmx_ddpme_t *ddpme,
                             gmx_bool bUniform, const gmx_ddbox_t *ddbox,
                             const real *cellFrac)
{
    gmx_domdec_comm_t *comm;
    int                nc, ns, s;
    int               *xmin, *xmax;
    real               range, pme_boundary;
    int                sh;

    comm = dd->comm;
    nc   = dd->nc[ddpme->dim];
    ns   = ddpme->nslab;

    if (!ddpme->dim_match)
    {
        /* PP decomposition is not along dim: the worst situation */
        sh = ns/2;
    }
    else if (ns <= 3 || (bUniform && ns == nc))
    {
        /* The optimal situation */
        sh = 1;
    }
    else
    {
        /* We need to check for all pme nodes which nodes they
         * could possibly need to communicate with.
         */
        xmin = ddpme->pp_min;
        xmax = ddpme->pp_max;
        /* Allow for atoms to be maximally 2/3 times the cut-off
         * out of their DD cell. This is a reasonable balance between
         * between performance and support for most charge-group/cut-off
         * combinations.
         */
        range  = 2.0/3.0*comm->cutoff/ddbox->box_size[ddpme->dim];
        /* Avoid extra communication when we are exactly at a boundary */
        range *= 0.999;

        sh = 1;
        for (s = 0; s < ns; s++)
        {
            /* PME slab s spreads atoms between box frac. s/ns and (s+1)/ns */
            pme_boundary = static_cast<real>(s)/ns;
            while (sh + 1 < ns &&
                   ((s - (sh + 1) >= 0 &&
                     cellFrac[xmax[s - (sh + 1)     ] + 1]     + range > pme_boundary) ||
                    (s - (sh + 1) <  0 &&
                     cellFrac[xmax[s - (sh + 1) + ns] + 1] - 1 + range > pme_boundary)))
            {
                sh++;
            }
            pme_boundary = static_cast<real>(s+1)/ns;
            while (sh + 1 < ns &&
                   ((s + (sh + 1) <  ns &&
                     cellFrac[xmin[s + (sh + 1)    ]  ]     - range < pme_boundary) ||
                    (s + (sh + 1) >= ns &&
                     cellFrac[xmin[s +(sh + 1) - ns]  ] + 1 - range < pme_boundary)))
            {
                sh++;
            }
        }
    }

    ddpme->maxshift = sh;

    if (debug)
    {
        fprintf(debug, "PME slab communication range for dim %d is %d\n",
                ddpme->dim, ddpme->maxshift);
    }
}

static void check_box_size(const gmx_domdec_t *dd,
                           const gmx_ddbox_t  *ddbox)
{
    int d, dim;

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];
        if (dim < ddbox->nboundeddim &&
            ddbox->box_size[dim]*ddbox->skew_fac[dim] <
            dd->nc[dim]*dd->comm->cellsize_limit*DD_CELL_MARGIN)
        {
            gmx_fatal(FARGS, "The %c-size of the box (%f) times the triclinic skew factor (%f) is smaller than the number of DD cells (%d) times the smallest allowed cell size (%f)\n",
                      dim2char(dim), ddbox->box_size[dim], ddbox->skew_fac[dim],
                      dd->nc[dim], dd->comm->cellsize_limit);
        }
    }
}

real grid_jump_limit(const gmx_domdec_comm_t *comm,
                     real                     cutoff,
                     int                      dim_ind)
{
    real grid_jump_limit;

    /* The distance between the boundaries of cells at distance
     * x+-1,y+-1 or y+-1,z+-1 is limited by the cut-off restrictions
     * and by the fact that cells should not be shifted by more than
     * half their size, such that cg's only shift by one cell
     * at redecomposition.
     */
    grid_jump_limit = comm->cellsize_limit;
    if (!comm->bVacDLBNoLimit)
    {
        if (comm->bPMELoadBalDLBLimits)
        {
            cutoff = std::max(cutoff, comm->PMELoadBal_max_cutoff);
        }
        grid_jump_limit = std::max(grid_jump_limit,
                                   cutoff/comm->cd[dim_ind].numPulses());
    }

    return grid_jump_limit;
}

/* This function should be used for moving the domain boudaries during DLB,
 * for obtaining the minimum cell size. It checks the initially set limit
 * comm->cellsize_min, for bonded and initial non-bonded cut-offs,
 * and, possibly, a longer cut-off limit set for PME load balancing.
 */
static real cellsize_min_dlb(gmx_domdec_comm_t *comm, int dim_ind, int dim)
{
    real cellsize_min;

    cellsize_min = comm->cellsize_min[dim];

    if (!comm->bVacDLBNoLimit)
    {
        /* The cut-off might have changed, e.g. by PME load balacning,
         * from the value used to set comm->cellsize_min, so check it.
         */
        cellsize_min = std::max(cellsize_min, comm->cutoff/comm->cd[dim_ind].np_dlb);

        if (comm->bPMELoadBalDLBLimits)
        {
            /* Check for the cut-off limit set by the PME load balancing */
            cellsize_min = std::max(cellsize_min, comm->PMELoadBal_max_cutoff/comm->cd[dim_ind].np_dlb);
        }
    }

    return cellsize_min;
}

/* Set the domain boundaries. Use for static (or no) load balancing,
 * and also for the starting state for dynamic load balancing.
 * setmode determine if and where the boundaries are stored, use enum above.
 * Returns the number communication pulses in npulse.
 */
gmx::ArrayRef < const std::vector < real>>
set_dd_cell_sizes_slb(gmx_domdec_t *dd, const gmx_ddbox_t *ddbox,
                      int setmode, ivec npulse)
{
    gmx_domdec_comm_t *comm = dd->comm;

    gmx::ArrayRef < std::vector < real>> cell_x_master;
    if (setmode == setcellsizeslbMASTER)
    {
        cell_x_master = dd->ma->cellSizesBuffer;
    }

    rvec cellsize_min;
    for (int d = 0; d < DIM; d++)
    {
        cellsize_min[d] = ddbox->box_size[d]*ddbox->skew_fac[d];
        npulse[d]       = 1;
        if (dd->nc[d] == 1 || comm->slb_frac[d] == nullptr)
        {
            /* Uniform grid */
            real cell_dx = ddbox->box_size[d]/dd->nc[d];
            switch (setmode)
            {
                case setcellsizeslbMASTER:
                    for (int j = 0; j < dd->nc[d]+1; j++)
                    {
                        cell_x_master[d][j] = ddbox->box0[d] + j*cell_dx;
                    }
                    break;
                case setcellsizeslbLOCAL:
                    comm->cell_x0[d] = ddbox->box0[d] + (dd->ci[d]  )*cell_dx;
                    comm->cell_x1[d] = ddbox->box0[d] + (dd->ci[d]+1)*cell_dx;
                    break;
                default:
                    break;
            }
            real cellsize = cell_dx*ddbox->skew_fac[d];
            while (cellsize*npulse[d] < comm->cutoff)
            {
                npulse[d]++;
            }
            cellsize_min[d] = cellsize;
        }
        else
        {
            /* Statically load balanced grid */
            /* Also when we are not doing a master distribution we determine
             * all cell borders in a loop to obtain identical values
             * to the master distribution case and to determine npulse.
             */
            gmx::ArrayRef<real> cell_x;
            std::vector<real>   cell_x_buffer;
            if (setmode == setcellsizeslbMASTER)
            {
                cell_x = cell_x_master[d];
            }
            else
            {
                cell_x_buffer.resize(dd->nc[d] + 1);
                cell_x = cell_x_buffer;
            }
            cell_x[0] = ddbox->box0[d];
            for (int j = 0; j < dd->nc[d]; j++)
            {
                real cell_dx  = ddbox->box_size[d]*comm->slb_frac[d][j];
                cell_x[j+1]   = cell_x[j] + cell_dx;
                real cellsize = cell_dx*ddbox->skew_fac[d];
                while (cellsize*npulse[d] < comm->cutoff &&
                       npulse[d] < dd->nc[d]-1)
                {
                    npulse[d]++;
                }
                cellsize_min[d] = std::min(cellsize_min[d], cellsize);
            }
            if (setmode == setcellsizeslbLOCAL)
            {
                comm->cell_x0[d] = cell_x[dd->ci[d]];
                comm->cell_x1[d] = cell_x[dd->ci[d]+1];
            }
        }
        /* The following limitation is to avoid that a cell would receive
         * some of its own home charge groups back over the periodic boundary.
         * Double charge groups cause trouble with the global indices.
         */
        if (d < ddbox->npbcdim &&
            dd->nc[d] > 1 && npulse[d] >= dd->nc[d])
        {
            char error_string[STRLEN];

            sprintf(error_string,
                    "The box size in direction %c (%f) times the triclinic skew factor (%f) is too small for a cut-off of %f with %d domain decomposition cells, use 1 or more than %d %s or increase the box size in this direction",
                    dim2char(d), ddbox->box_size[d], ddbox->skew_fac[d],
                    comm->cutoff,
                    dd->nc[d], dd->nc[d],
                    dd->nnodes > dd->nc[d] ? "cells" : "ranks");

            if (setmode == setcellsizeslbLOCAL)
            {
                gmx_fatal_collective(FARGS, dd->mpi_comm_all, DDMASTER(dd),
                                     "%s", error_string);
            }
            else
            {
                gmx_fatal(FARGS, "%s", error_string);
            }
        }
    }

    if (!isDlbOn(comm))
    {
        copy_rvec(cellsize_min, comm->cellsize_min);
    }

    for (int d = 0; d < comm->npmedecompdim; d++)
    {
        set_pme_maxshift(dd, &comm->ddpme[d],
                         comm->slb_frac[dd->dim[d]] == nullptr, ddbox,
                         comm->ddpme[d].slb_dim_f);
    }

    return cell_x_master;
}


static void dd_cell_sizes_dlb_root_enforce_limits(gmx_domdec_t *dd,
                                                  int d, int dim, RowMaster *rowMaster,
                                                  const gmx_ddbox_t *ddbox,
                                                  gmx_bool bUniform, int64_t step, real cellsize_limit_f, int range[])
{
    gmx_domdec_comm_t *comm;
    real               halfway, cellsize_limit_f_i, region_size;
    gmx_bool           bLastHi  = FALSE;
    int                nrange[] = {range[0], range[1]};

    region_size = rowMaster->cellFrac[range[1]] - rowMaster->cellFrac[range[0]];

    comm = dd->comm;

    const int           ncd       = dd->nc[dim];

    const bool          dimHasPbc = (dim < ddbox->npbcdim);

    gmx::ArrayRef<real> cell_size = rowMaster->buf_ncd;

    if (debug)
    {
        fprintf(debug, "enforce_limits: %d %d\n", range[0], range[1]);
    }

    /* First we need to check if the scaling does not make cells
     * smaller than the smallest allowed size.
     * We need to do this iteratively, since if a cell is too small,
     * it needs to be enlarged, which makes all the other cells smaller,
     * which could in turn make another cell smaller than allowed.
     */
    for (int i = range[0]; i < range[1]; i++)
    {
        rowMaster->isCellMin[i] = false;
    }
    int nmin = 0;
    int nmin_old;
    do
    {
        nmin_old = nmin;
        /* We need the total for normalization */
        real fac = 0;
        for (int i = range[0]; i < range[1]; i++)
        {
            if (!rowMaster->isCellMin[i])
            {
                fac += cell_size[i];
            }
        }
        fac = ( region_size - nmin*cellsize_limit_f)/fac; /* substracting cells already set to cellsize_limit_f */
        /* Determine the cell boundaries */
        for (int i = range[0]; i < range[1]; i++)
        {
            if (!rowMaster->isCellMin[i])
            {
                cell_size[i] *= fac;
                if (!dimHasPbc && (i == 0 || i == dd->nc[dim] - 1))
                {
                    cellsize_limit_f_i = 0;
                }
                else
                {
                    cellsize_limit_f_i = cellsize_limit_f;
                }
                if (cell_size[i] < cellsize_limit_f_i)
                {
                    rowMaster->isCellMin[i] = true;
                    cell_size[i]            = cellsize_limit_f_i;
                    nmin++;
                }
            }
            rowMaster->cellFrac[i + 1] = rowMaster->cellFrac[i] + cell_size[i];
        }
    }
    while (nmin > nmin_old);

    const int i  = range[1] - 1;
    cell_size[i] = rowMaster->cellFrac[i + 1] - rowMaster->cellFrac[i];
    /* For this check we should not use DD_CELL_MARGIN,
     * but a slightly smaller factor,
     * since rounding could get use below the limit.
     */
    if (dimHasPbc &&
        cell_size[i] < cellsize_limit_f*DD_CELL_MARGIN2/DD_CELL_MARGIN)
    {
        char buf[22];
        gmx_fatal(FARGS, "step %s: the dynamic load balancing could not balance dimension %c: box size %f, triclinic skew factor %f, #cells %d, minimum cell size %f\n",
                  gmx_step_str(step, buf),
                  dim2char(dim), ddbox->box_size[dim], ddbox->skew_fac[dim],
                  ncd, comm->cellsize_min[dim]);
    }

    rowMaster->dlbIsLimited = (nmin > 0) || (range[0] > 0) || (range[1] < ncd);

    if (!bUniform)
    {
        /* Check if the boundary did not displace more than halfway
         * each of the cells it bounds, as this could cause problems,
         * especially when the differences between cell sizes are large.
         * If changes are applied, they will not make cells smaller
         * than the cut-off, as we check all the boundaries which
         * might be affected by a change and if the old state was ok,
         * the cells will at most be shrunk back to their old size.
         */
        for (int i = range[0]+1; i < range[1]; i++)
        {
            halfway = 0.5*(rowMaster->oldCellFrac[i] + rowMaster->oldCellFrac[i - 1]);
            if (rowMaster->cellFrac[i] < halfway)
            {
                rowMaster->cellFrac[i] = halfway;
                /* Check if the change also causes shifts of the next boundaries */
                for (int j = i + 1; j < range[1]; j++)
                {
                    if (rowMaster->cellFrac[j] < rowMaster->cellFrac[j - 1] + cellsize_limit_f)
                    {
                        rowMaster->cellFrac[j] =  rowMaster->cellFrac[j - 1] + cellsize_limit_f;
                    }
                }
            }
            halfway = 0.5*(rowMaster->oldCellFrac[i] + rowMaster->oldCellFrac[i + 1]);
            if (rowMaster->cellFrac[i] > halfway)
            {
                rowMaster->cellFrac[i] = halfway;
                /* Check if the change also causes shifts of the next boundaries */
                for (int j = i - 1; j >= range[0] + 1; j--)
                {
                    if (rowMaster->cellFrac[j] > rowMaster->cellFrac[j + 1] - cellsize_limit_f)
                    {
                        rowMaster->cellFrac[j] = rowMaster->cellFrac[j + 1] - cellsize_limit_f;
                    }
                }
            }
        }
    }

    /* nrange is defined as [lower, upper) range for new call to enforce_limits */
    /* find highest violation of LimLo (a) and the following violation of LimHi (thus the lowest following) (b)
     * then call enforce_limits for (oldb,a), (a,b). In the next step: (b,nexta). oldb and nexta can be the boundaries.
     * for a and b nrange is used */
    if (d > 0)
    {
        /* Take care of the staggering of the cell boundaries */
        if (bUniform)
        {
            for (int i = range[0]; i < range[1]; i++)
            {
                rowMaster->bounds[i].cellFracLowerMax = rowMaster->cellFrac[i];
                rowMaster->bounds[i].cellFracUpperMin = rowMaster->cellFrac[i + 1];
            }
        }
        else
        {
            for (int i = range[0] + 1; i < range[1]; i++)
            {
                const RowMaster::Bounds &bounds = rowMaster->bounds[i];

                bool                     bLimLo = (rowMaster->cellFrac[i] < bounds.boundMin);
                bool                     bLimHi = (rowMaster->cellFrac[i] > bounds.boundMax);
                if (bLimLo && bLimHi)
                {
                    /* Both limits violated, try the best we can */
                    /* For this case we split the original range (range) in two parts and care about the other limitiations in the next iteration. */
                    rowMaster->cellFrac[i] = 0.5*(bounds.boundMin + bounds.boundMax);
                    nrange[0]              = range[0];
                    nrange[1]              = i;
                    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);

                    nrange[0] = i;
                    nrange[1] = range[1];
                    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);

                    return;
                }
                else if (bLimLo)
                {
                    /* rowMaster->cellFrac[i] = rowMaster->boundMin[i]; */
                    nrange[1] = i;  /* only store violation location. There could be a LimLo violation following with an higher index */
                    bLastHi   = FALSE;
                }
                else if (bLimHi && !bLastHi)
                {
                    bLastHi = TRUE;
                    if (nrange[1] < range[1])   /* found a LimLo before */
                    {
                        rowMaster->cellFrac[nrange[1]] = rowMaster->bounds[1].boundMin;
                        dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);
                        nrange[0] = nrange[1];
                    }
                    rowMaster->cellFrac[i] = rowMaster->bounds[i].boundMax;
                    nrange[1]              = i;
                    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);
                    nrange[0] = i;
                    nrange[1] = range[1];
                }
            }
            if (nrange[1] < range[1])   /* found last a LimLo */
            {
                rowMaster->cellFrac[nrange[1]] = rowMaster->bounds[nrange[1]].boundMin;
                dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);
                nrange[0] = nrange[1];
                nrange[1] = range[1];
                dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);
            }
            else if (nrange[0] > range[0]) /* found at least one LimHi */
            {
                dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, nrange);
            }
        }
    }
}


static void set_dd_cell_sizes_dlb_root(gmx_domdec_t *dd,
                                       int d, int dim, RowMaster *rowMaster,
                                       const gmx_ddbox_t *ddbox,
                                       gmx_bool bDynamicBox,
                                       gmx_bool bUniform, int64_t step)
{
    gmx_domdec_comm_t *comm    = dd->comm;
    constexpr real     c_relax = 0.5;
    int                range[] = { 0, 0 };

    /* Convert the maximum change from the input percentage to a fraction */
    const real          change_limit = comm->dlb_scale_lim*0.01;

    const int           ncd          = dd->nc[dim];

    const bool          bPBC = (dim < ddbox->npbcdim);

    gmx::ArrayRef<real> cell_size = rowMaster->buf_ncd;

    /* Store the original boundaries */
    for (int i = 0; i < ncd + 1; i++)
    {
        rowMaster->oldCellFrac[i] = rowMaster->cellFrac[i];
    }
    if (bUniform)
    {
        for (int i = 0; i < ncd; i++)
        {
            cell_size[i] = 1.0/ncd;
        }
    }
    else if (dd_load_count(comm) > 0)
    {
        real load_aver  = comm->load[d].sum_m/ncd;
        real change_max = 0;
        real load_i;
        real change;
        for (int i = 0; i < ncd; i++)
        {
            /* Determine the relative imbalance of cell i */
            load_i         = comm->load[d].load[i*comm->load[d].nload+2];
            real imbalance = (load_i - load_aver)/(load_aver > 0 ? load_aver : 1);
            /* Determine the change of the cell size using underrelaxation */
            change     = -c_relax*imbalance;
            change_max = std::max(change_max, std::max(change, -change));
        }
        /* Limit the amount of scaling.
         * We need to use the same rescaling for all cells in one row,
         * otherwise the load balancing might not converge.
         */
        real sc = c_relax;
        if (change_max > change_limit)
        {
            sc *= change_limit/change_max;
        }
        for (int i = 0; i < ncd; i++)
        {
            /* Determine the relative imbalance of cell i */
            load_i         = comm->load[d].load[i*comm->load[d].nload+2];
            real imbalance = (load_i - load_aver)/(load_aver > 0 ? load_aver : 1);
            /* Determine the change of the cell size using underrelaxation */
            change       = -sc*imbalance;
            cell_size[i] = (rowMaster->cellFrac[i + 1] - rowMaster->cellFrac[i])*(1 + change);
        }
    }

    real cellsize_limit_f  = cellsize_min_dlb(comm, d, dim)/ddbox->box_size[dim];
    cellsize_limit_f      *= DD_CELL_MARGIN;
    real dist_min_f_hard   = grid_jump_limit(comm, comm->cutoff, d)/ddbox->box_size[dim];
    real dist_min_f        = dist_min_f_hard * DD_CELL_MARGIN;
    if (ddbox->tric_dir[dim])
    {
        cellsize_limit_f /= ddbox->skew_fac[dim];
        dist_min_f       /= ddbox->skew_fac[dim];
    }
    if (bDynamicBox && d > 0)
    {
        dist_min_f *= DD_PRES_SCALE_MARGIN;
    }
    if (d > 0 && !bUniform)
    {
        /* Make sure that the grid is not shifted too much */
        for (int i = 1; i < ncd; i++)
        {
            const RowMaster::Bounds &boundsNeighbor = rowMaster->bounds[i - 1];
            RowMaster::Bounds       &bounds         = rowMaster->bounds[i];

            if (bounds.cellFracUpperMin - boundsNeighbor.cellFracLowerMax < 2 * dist_min_f_hard)
            {
                gmx_incons("Inconsistent DD boundary staggering limits!");
            }
            bounds.boundMin = boundsNeighbor.cellFracLowerMax + dist_min_f;
            real space       = rowMaster->cellFrac[i] - (boundsNeighbor.cellFracLowerMax + dist_min_f);
            if (space > 0)
            {
                bounds.boundMin += 0.5*space;
            }
            bounds.boundMax  = bounds.cellFracUpperMin - dist_min_f;
            space            = rowMaster->cellFrac[i] - (bounds.cellFracUpperMin - dist_min_f);
            if (space < 0)
            {
                rowMaster->bounds[i].boundMax += 0.5*space;
            }
            if (debug)
            {
                fprintf(debug,
                        "dim %d boundary %d %.3f < %.3f < %.3f < %.3f < %.3f\n",
                        d, i,
                        boundsNeighbor.cellFracLowerMax + dist_min_f,
                        bounds.boundMin, rowMaster->cellFrac[i], bounds.boundMax,
                        bounds.cellFracUpperMin - dist_min_f);
            }
        }
    }
    range[1]                 = ncd;
    rowMaster->cellFrac[0]   = 0;
    rowMaster->cellFrac[ncd] = 1;
    dd_cell_sizes_dlb_root_enforce_limits(dd, d, dim, rowMaster, ddbox, bUniform, step, cellsize_limit_f, range);


    /* After the checks above, the cells should obey the cut-off
     * restrictions, but it does not hurt to check.
     */
    for (int i = 0; i < ncd; i++)
    {
        if (debug)
        {
            fprintf(debug, "Relative bounds dim %d  cell %d: %f %f\n",
                    dim, i, rowMaster->cellFrac[i], rowMaster->cellFrac[i + 1]);
        }

        if ((bPBC || (i != 0 && i != dd->nc[dim]-1)) &&
            rowMaster->cellFrac[i + 1] - rowMaster->cellFrac[i] <
            cellsize_limit_f/DD_CELL_MARGIN)
        {
            char buf[22];
            fprintf(stderr,
                    "\nWARNING step %s: direction %c, cell %d too small: %f\n",
                    gmx_step_str(step, buf), dim2char(dim), i,
                    (rowMaster->cellFrac[i + 1] - rowMaster->cellFrac[i])
                    *ddbox->box_size[dim]*ddbox->skew_fac[dim]);
        }
    }

    int pos = ncd + 1;
    /* Store the cell boundaries of the lower dimensions at the end */
    for (int d1 = 0; d1 < d; d1++)
    {
        rowMaster->cellFrac[pos++] = comm->cellsizesWithDlb[d1].fracLower;
        rowMaster->cellFrac[pos++] = comm->cellsizesWithDlb[d1].fracUpper;
    }

    if (d < comm->npmedecompdim)
    {
        /* The master determines the maximum shift for
         * the coordinate communication between separate PME nodes.
         */
        set_pme_maxshift(dd, &comm->ddpme[d], bUniform, ddbox, rowMaster->cellFrac.data());
    }
    rowMaster->cellFrac[pos++] = comm->ddpme[0].maxshift;
    if (d >= 1)
    {
        rowMaster->cellFrac[pos++] = comm->ddpme[1].maxshift;
    }
}

static void relative_to_absolute_cell_bounds(gmx_domdec_t      *dd,
                                             const gmx_ddbox_t *ddbox,
                                             int                dimind)
{
    gmx_domdec_comm_t        *comm      = dd->comm;
    const DDCellsizesWithDlb &cellsizes = comm->cellsizesWithDlb[dimind];

    /* Set the cell dimensions */
    int dim            = dd->dim[dimind];
    comm->cell_x0[dim] = cellsizes.fracLower*ddbox->box_size[dim];
    comm->cell_x1[dim] = cellsizes.fracUpper*ddbox->box_size[dim];
    if (dim >= ddbox->nboundeddim)
    {
        comm->cell_x0[dim] += ddbox->box0[dim];
        comm->cell_x1[dim] += ddbox->box0[dim];
    }
}

static void distribute_dd_cell_sizes_dlb(gmx_domdec_t *dd,
                                         int d, int dim,
                                         gmx::ArrayRef<real> cellFracRow,
                                         const gmx_ddbox_t *ddbox)
{
    gmx_domdec_comm_t &comm = *dd->comm;

#if GMX_MPI
    /* Each node would only need to know two fractions,
     * but it is probably cheaper to broadcast the whole array.
     */
    MPI_Bcast(cellFracRow.data(), ddCellFractionBufferSize(dd, d)*sizeof(real), MPI_BYTE,
              0, comm.mpi_comm_load[d]);
#endif
    /* Copy the fractions for this dimension from the buffer */
    comm.cellsizesWithDlb[d].fracLower = cellFracRow[dd->ci[dim]    ];
    comm.cellsizesWithDlb[d].fracUpper = cellFracRow[dd->ci[dim] + 1];
    /* The whole array was communicated, so set the buffer position */
    int pos = dd->nc[dim] + 1;
    for (int d1 = 0; d1 <= d; d1++)
    {
        if (d1 < d)
        {
            /* Copy the cell fractions of the lower dimensions */
            comm.cellsizesWithDlb[d1].fracLower = cellFracRow[pos++];
            comm.cellsizesWithDlb[d1].fracUpper = cellFracRow[pos++];
        }
        relative_to_absolute_cell_bounds(dd, ddbox, d1);
    }
    /* Convert the communicated shift from float to int */
    comm.ddpme[0].maxshift = gmx::roundToInt(cellFracRow[pos++]);
    if (d >= 1)
    {
        comm.ddpme[1].maxshift = gmx::roundToInt(cellFracRow[pos++]);
    }
}

static void set_dd_cell_sizes_dlb_change(gmx_domdec_t *dd,
                                         const gmx_ddbox_t *ddbox,
                                         gmx_bool bDynamicBox,
                                         gmx_bool bUniform, int64_t step)
{
    for (int d = 0; d < dd->ndim; d++)
    {
        const int dim         = dd->dim[d];
        bool      isRowMember = true;
        bool      isRowRoot   = true;
        for (int d1 = d; d1 < dd->ndim; d1++)
        {
            if (dd->ci[dd->dim[d1]] > 0)
            {
                if (d1 != d)
                {
                    isRowMember = false;
                }
                isRowRoot = false;
            }
        }
        if (isRowMember)
        {
            DDCellsizesWithDlb  &cellsizes = dd->comm->cellsizesWithDlb[d];
            gmx::ArrayRef<real>  cellFracRow;

            if (isRowRoot)
            {
                RowMaster *rowMaster = cellsizes.rowMaster.get();
                set_dd_cell_sizes_dlb_root(dd, d, dim, rowMaster,
                                           ddbox, bDynamicBox, bUniform, step);
                cellFracRow = rowMaster->cellFrac;
            }
            else
            {
                cellFracRow = cellsizes.fracRow;
            }
            distribute_dd_cell_sizes_dlb(dd, d, dim, cellFracRow, ddbox);
        }
    }
}

static void set_dd_cell_sizes_dlb_nochange(gmx_domdec_t      *dd,
                                           const gmx_ddbox_t *ddbox)
{
    int d;

    /* This function assumes the box is static and should therefore
     * not be called when the box has changed since the last
     * call to dd_partition_system.
     */
    for (d = 0; d < dd->ndim; d++)
    {
        relative_to_absolute_cell_bounds(dd, ddbox, d);
    }
}



static void set_dd_cell_sizes_dlb(gmx_domdec_t *dd,
                                  const gmx_ddbox_t *ddbox, gmx_bool bDynamicBox,
                                  gmx_bool bUniform, gmx_bool bDoDLB, int64_t step,
                                  gmx_wallcycle_t wcycle)
{
    gmx_domdec_comm_t *comm;
    int                dim;

    comm = dd->comm;

    if (bDoDLB)
    {
        wallcycle_start(wcycle, ewcDDCOMMBOUND);
        set_dd_cell_sizes_dlb_change(dd, ddbox, bDynamicBox, bUniform, step);
        wallcycle_stop(wcycle, ewcDDCOMMBOUND);
    }
    else if (bDynamicBox)
    {
        set_dd_cell_sizes_dlb_nochange(dd, ddbox);
    }

    /* Set the dimensions for which no DD is used */
    for (dim = 0; dim < DIM; dim++)
    {
        if (dd->nc[dim] == 1)
        {
            comm->cell_x0[dim] = 0;
            comm->cell_x1[dim] = ddbox->box_size[dim];
            if (dim >= ddbox->nboundeddim)
            {
                comm->cell_x0[dim] += ddbox->box0[dim];
                comm->cell_x1[dim] += ddbox->box0[dim];
            }
        }
    }
}

void set_dd_cell_sizes(gmx_domdec_t *dd,
                       const gmx_ddbox_t *ddbox, gmx_bool bDynamicBox,
                       gmx_bool bUniform, gmx_bool bDoDLB, int64_t step,
                       gmx_wallcycle_t wcycle)
{
    gmx_domdec_comm_t *comm = dd->comm;

    /* Copy the old cell boundaries for the cg displacement check */
    copy_rvec(comm->cell_x0, comm->old_cell_x0);
    copy_rvec(comm->cell_x1, comm->old_cell_x1);

    if (isDlbOn(comm))
    {
        if (DDMASTER(dd))
        {
            check_box_size(dd, ddbox);
        }
        set_dd_cell_sizes_dlb(dd, ddbox, bDynamicBox, bUniform, bDoDLB, step, wcycle);
    }
    else
    {
        ivec numPulses;
        set_dd_cell_sizes_slb(dd, ddbox, setcellsizeslbLOCAL, numPulses);

        /* Check if the change in cell size requires a different number
         * of communication pulses and if so change the number.
         */
        for (int d = 0; d < dd->ndim; d++)
        {
            gmx_domdec_comm_dim_t &cd           = comm->cd[d];
            int                    numPulsesDim = numPulses[dd->dim[d]];
            if (cd.numPulses() != numPulsesDim)
            {
                if (debug)
                {
                    fprintf(debug, "Changing the number of halo communication pulses along dim %c from %d to %d\n",
                            dim2char(dd->dim[d]), cd.numPulses(), numPulsesDim);
                }
                cd.ind.resize(numPulsesDim);
            }
        }
    }

    if (debug)
    {
        for (int d = 0; d < DIM; d++)
        {
            fprintf(debug, "cell_x[%d] %f - %f skew_fac %f\n",
                    d, comm->cell_x0[d], comm->cell_x1[d], ddbox->skew_fac[d]);
        }
    }
}
