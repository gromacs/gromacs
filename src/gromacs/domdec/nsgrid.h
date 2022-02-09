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
#ifndef GMX_MDLIB_NSGRID_H
#define GMX_MDLIB_NSGRID_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct gmx_ddbox_t;

/*! \brief Used when estimating the interaction density.
 *
 * c_gridStdDevFactor * stddev estimates the interaction density. The
 * value sqrt(3) == 1.73205080757 gives a uniform load for a
 * rectangular 3D block of charge groups. For a sphere, it is not a
 * bad approximation for 4x1x1 up to 4x2x2.
 *
 * \todo It would be nicer to use sqrt(3) here, when all code that
 * includes this file is in C++, which will let us cope with the
 * std::sqrt<T> on Windows. */
constexpr real c_gridStdDevFactor = 1.73205080757;

void get_nsgrid_boundaries(int                  nboundeddim,
                           matrix               box,
                           struct gmx_domdec_t* dd,
                           gmx_ddbox_t*         ddbox,
                           gmx::RVec*           gr0,
                           gmx::RVec*           gr1,
                           int                  ncg,
                           rvec*                cgcm,
                           rvec                 grid_x0,
                           rvec                 grid_x1);
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
