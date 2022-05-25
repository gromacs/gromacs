/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

/*! \libinternal \file
 * \brief Defines a struct useful for transferring the PME output
 * values
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_OUTPUT_H
#define GMX_EWALD_PME_OUTPUT_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

// TODO There's little value in computing the Coulomb and LJ virial
// separately, so we should simplify that.
// TODO The matrices might be best as a view, but not currently
// possible. Use mdspan?
struct PmeOutput
{
    //!< Host staging area for PME forces
    gmx::ArrayRef<gmx::RVec> forces_;
    //!< True if forces have been staged other false (when forces are reduced on the GPU).
    bool haveForceOutput_ = false;
    //!< Host staging area for PME coulomb energy
    real coulombEnergy_ = 0;
    //!< Host staging area for PME coulomb virial contributions
    matrix coulombVirial_ = { { 0 } };
    //!< Host staging area for PME coulomb dVdl.
    real coulombDvdl_ = 0;
    //!< Host staging area for PME LJ dVdl.
    real lennardJonesDvdl_ = 0;
    //!< Host staging area for PME LJ energy
    real lennardJonesEnergy_ = 0;
    //!< Host staging area for PME LJ virial contributions
    matrix lennardJonesVirial_ = { { 0 } };
};

#endif
