/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
#ifndef GMX_EWALD_PME_GATHER_H
#define GMX_EWALD_PME_GATHER_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class PmeAtomComm;
struct gmx_pme_t;
struct splinedata_t;

namespace gmx
{
template<typename T>
class ArrayRef;
}

/*! Gather the forces from the grid
 *
 * \param[in] pme          General PME settings
 * \param[in] grid         The grid with potential values to gather the forces from
 * \param[in] clearForces  Whether the forces in \p atc should be cleared (set instead of accumulate)
 * \param[in,out] atc      Contains atom indices, charges and force output buffer
 * \param[in] spline       The spline coefficients for the atoms
 * \param[in] scaleFactor  Factor to scale the forces with
 */
void gather_f_bsplines(const gmx_pme_t&          pme,
                       gmx::ArrayRef<const real> grid,
                       bool                      clearForces,
                       PmeAtomComm*              atc,
                       const splinedata_t&       spline,
                       real                      scaleFactor);

/*! Computes and returns the potential energy of the PME mesh part
 *
 * Computes and return the potential energy of the PME mesh part using the potential values
 * stored in \p pme and the charges \p atc. This is only useful when the charges in \p atc
 * are different from those used to compute the potential (otherwise is simpler and
 * computationally cheaper to get the potential from the solve function).
 *
 * \note: Does not support MPI or OpenMP parallelization (is asserted)
 *
 * \param[in] pme   General PME settings
 * \param[in] grid  The grid with potential values to gather the forces from
 * \param[in] atc   Contains atom indices, charges and force output buffer
 */
real gather_energy_bsplines(const gmx_pme_t& pme, gmx::ArrayRef<const real> grid, const PmeAtomComm& atc);

#endif
