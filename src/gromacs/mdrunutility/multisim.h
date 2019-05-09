/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief Declares the multi-simulation support routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_MULTISIM_H
#define GMX_MDRUNUTILITY_MULTISIM_H

#include <string>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxmpi.h"

struct gmx_multisim_t;

/*! \brief Initializes multi-simulations.
 *
 * Splits the communication into multidirs.size() separate
 * simulations, if >1, and creates a communication structure between
 * the master these simulations. */
gmx_multisim_t *init_multisystem(MPI_Comm                         comm,
                                 gmx::ArrayRef<const std::string> multidirs);

//! Cleans up multi-system handler.
void done_multisim(gmx_multisim_t *ms);

//! Are we doing multiple independent simulations?
static bool inline isMultiSim(const gmx_multisim_t *ms)
{
    return ms != nullptr;
}

//! Are we the master simulation of a possible multi-simulation?
bool isMasterSim(const gmx_multisim_t *ms);

/*! \brief Are we the master rank (of the master simulation, for a multi-sim).
 *
 * This rank prints the remaining run time etc. */
bool isMasterSimMasterRank(const gmx_multisim_t *ms,
                           bool                  isMaster);

#endif
