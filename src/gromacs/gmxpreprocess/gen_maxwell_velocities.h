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

#ifndef GMX_MAXWELL_VELOCITIES
#define GMX_MAXWELL_VELOCITIES

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;

namespace gmx
{
class MDLogger;
}

/*! \brief
 * Generate Maxwellian velocities.
 *
 * \param[in] tempi Temperature to generate around
 * \param[in] seed  Random number generator seed. A new one is
 *                  generated if this is -1.
 * \param[in] mtop  Molecular Topology
 * \param[out] v    Velocities
 * \param[in] logger Handle to logging interface.
 */
void maxwell_speed(real tempi, int seed, gmx_mtop_t* mtop, rvec v[], const gmx::MDLogger& logger);

/*! \brief
 * Remove the center of mass motion in a set of coordinates.
 *
 * \param[in]  logger Handle to logging interface.
 * \param[in]  natoms Number of atoms
 * \param[in]  mass   Atomic masses
 * \param[in]  x      Coordinates
 * \param[out] v      Velocities
 */
void stop_cm(const gmx::MDLogger& logger, int natoms, real mass[], rvec x[], rvec v[]);

#endif
