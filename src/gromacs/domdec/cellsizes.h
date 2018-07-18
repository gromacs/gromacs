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
 * \brief Declares DD cell-size related functions.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_CELLSIZES_H
#define GMX_DOMDEC_DOMDEC_CELLSIZES_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"

struct gmx_ddbox_t;
struct gmx_domdec_comm_t;
struct gmx_domdec_t;

/*! \brief Options for setting up a regular, possibly static load balanced, cell grid geometry */
enum
{
    setcellsizeslbLOCAL,      //!< Set cell sizes locally on each rank
    setcellsizeslbMASTER,     //!< Set cell sizes on master rank only
    setcellsizeslbPULSE_ONLY  //!< Only set the communication pulses, not the cell sizes
};

/*! \brief Returns the minimum allowed distance between lower and upper bounds of zones along dimension dim_ind */
real grid_jump_limit(const gmx_domdec_comm_t *comm,
                     real                     cutoff,
                     int                      dim_ind);

/*! \brief Sets up an initial, non-staggered grid geometry, possibly using static load balancing
 *
 * The number of communication pulses per dimension is returned in numPulses.
 * When setmode==setcellsizeslbMASTER, the cell boundaries per dimension are
 * returned, otherwise an empty arrayref is returned.
 */
gmx::ArrayRef < const std::vector < real>>
set_dd_cell_sizes_slb(gmx_domdec_t      *dd,
                      const gmx_ddbox_t *ddbox,
                      int                setmode,
                      ivec               numPulses);

/*! \brief General cell size adjustment, possibly applying dynamic load balancing */
void set_dd_cell_sizes(gmx_domdec_t *dd,
                       const gmx_ddbox_t *ddbox, gmx_bool bDynamicBox,
                       gmx_bool bUniform, gmx_bool bDoDLB, int64_t step,
                       gmx_wallcycle_t wcycle);

#endif
