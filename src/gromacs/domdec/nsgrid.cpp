/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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

#include "nsgrid.h"

#include <cmath>
#include <cstdio>

#include "gromacs/domdec/dlb.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"

/*! \brief The extent of the neighborsearch grid is a bit larger than sqrt(3)
 * to account for less dense regions at the edges of the system.
 */
constexpr real c_stdDevFactor = 2.0;

/***********************************
 *         Grid Routines
 ***********************************/

static void calc_x_av_stddev(int n, rvec* x, rvec av, rvec stddev)
{
    dvec s1, s2;
    int  i, d;

    clear_dvec(s1);
    clear_dvec(s2);

    for (i = 0; i < n; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            s1[d] += x[i][d];
            s2[d] += x[i][d] * x[i][d];
        }
    }

    dsvmul(1.0 / n, s1, s1);
    dsvmul(1.0 / n, s2, s2);

    for (d = 0; d < DIM; d++)
    {
        av[d]     = s1[d];
        stddev[d] = std::sqrt(s2[d] - s1[d] * s1[d]);
    }
}

static void get_nsgrid_boundaries_vac(real av, real stddev, real* bound0, real* bound1, real* bdens0, real* bdens1)
{
    /* Set the grid to 2 times the standard deviation of
     * the charge group centers in both directions.
     * For a uniform bounded distribution the width is sqrt(3)*stddev,
     * so all charge groups fall within the width.
     * For a sphere stddev is r/sqrt(5): 99.2% falls within the width.
     * For a Gaussian distribution 98% fall within the width.
     */
    *bound0 = av - c_stdDevFactor * stddev;
    *bound1 = av + c_stdDevFactor * stddev;

    *bdens0 = av - c_gridStdDevFactor * stddev;
    *bdens1 = av + c_gridStdDevFactor * stddev;
}

static void dd_box_bounds_to_ns_bounds(real box0, real box_size, real* gr0, real* gr1)
{
    real av, stddev;

    /* Redetermine av and stddev from the DD box boundaries */
    av     = box0 + 0.5 * box_size;
    stddev = 0.5 * box_size / c_gridStdDevFactor;

    *gr0 = av - c_stdDevFactor * stddev;
    *gr1 = av + c_stdDevFactor * stddev;
}

void get_nsgrid_boundaries(int           nboundeddim,
                           matrix        box,
                           gmx_domdec_t* dd,
                           gmx_ddbox_t*  ddbox,
                           gmx::RVec*    gr0,
                           gmx::RVec*    gr1,
                           int           ncg,
                           rvec*         cgcm,
                           rvec          grid_x0,
                           rvec          grid_x1)
{
    rvec av, stddev;
    real bdens0, bdens1;
    int  d;

    if (nboundeddim < DIM)
    {
        calc_x_av_stddev(ncg, cgcm, av, stddev);
    }

    for (d = 0; d < DIM; d++)
    {
        if (d < nboundeddim)
        {
            grid_x0[d] = (gr0 != nullptr ? (*gr0)[d] : 0);
            grid_x1[d] = (gr1 != nullptr ? (*gr1)[d] : box[d][d]);
        }
        else
        {
            if (ddbox == nullptr)
            {
                get_nsgrid_boundaries_vac(av[d], stddev[d], &grid_x0[d], &grid_x1[d], &bdens0, &bdens1);
            }
            else
            {
                /* Temporary fix which uses global ddbox boundaries
                 * for unbounded dimensions.
                 * Should be replaced by local boundaries, which makes
                 * the ns grid smaller and does not require global comm.
                 */
                dd_box_bounds_to_ns_bounds(ddbox->box0[d], ddbox->box_size[d], &grid_x0[d], &grid_x1[d]);
                bdens0 = grid_x0[d];
                bdens1 = grid_x1[d];
            }
            /* Check for a DD cell not at a lower edge */
            if (dd != nullptr && gr0 != nullptr && dd->ci[d] > 0)
            {
                grid_x0[d] = (*gr0)[d];
                bdens0     = (*gr0)[d];
            }
            /* Check for a DD cell not at a higher edge */
            if (dd != nullptr && gr1 != nullptr && dd->ci[d] < dd->numCells[d] - 1)
            {
                grid_x1[d] = (*gr1)[d];
                bdens1     = (*gr1)[d];
            }
        }

        if (debug)
        {
            fprintf(debug, "Set grid boundaries dim %d: %f %f\n", d, grid_x0[d], grid_x1[d]);
        }
    }
}
