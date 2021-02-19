/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2019,2021, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_NSGRID_H
#define GMX_MDLIB_NSGRID_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct gmx_domdec_zones_t;
struct gmx_ddbox_t;
struct t_forcerec;

/*! \brief Used when estimating the interaction density.
 *
 * GRID_STDDEV_FAC * stddev estimates the interaction density. The
 * value sqrt(3) == 1.73205080757 gives a uniform load for a
 * rectangular 3D block of charge groups. For a sphere, it is not a
 * bad approximation for 4x1x1 up to 4x2x2.
 *
 * \todo It would be nicer to use sqrt(3) here, when all code that
 * includes this file is in C++, which will let us cope with the
 * std::sqrt<T> on Windows. */
static const real GRID_STDDEV_FAC = 1.73205080757;

/*! \brief The extent of the neighborsearch grid is a bit larger than sqrt(3)
 * to account for less dense regions at the edges of the system.
 */
static const real NSGRID_STDDEV_FAC = 2.0;

#define NSGRID_SIGNAL_MOVED_FAC 4
/* A cell index of NSGRID_SIGNAL_MOVED_FAC*ncells signals
 * that a charge group moved to another DD domain.
 */

void get_nsgrid_boundaries(int                  nboundeddim,
                           matrix               box,
                           struct gmx_domdec_t* dd,
                           gmx_ddbox_t*         ddbox,
                           rvec*                gr0,
                           rvec*                gr1,
                           int                  ncg,
                           rvec*                cgcm,
                           rvec                 grid_x0,
                           rvec                 grid_x1,
                           real*                grid_density);
/* Return the ns grid boundaries grid_x0 and grid_x1
 * and the estimate for the grid density.
 * For non-bounded dimensions the boundaries are determined
 * from the average and std.dev. of cgcm.
 * The are determined from box, unless gr0!=NULL or gr1!=NULL,
 * then they are taken from gr0 or gr1.
 * With dd and unbounded dimensions, the proper grid borders for cells
 * on the edges are determined from cgcm.
 */

#endif
