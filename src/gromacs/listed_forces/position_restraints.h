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
/*! \libinternal \file
 *
 * \brief This file contains declarations necessary for low-level
 * functions for computing energies and forces for position
 * restraints.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_listed_forces
 */

#ifndef GMX_LISTED_FORCES_POSITION_RESTRAINTS_H
#define GMX_LISTED_FORCES_POSITION_RESTRAINTS_H

#include <cstdio>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

struct gmx_enerdata_t;
struct gmx_wallcycle;
struct t_forcerec;
class InteractionDefinitions;
struct t_pbc;

namespace gmx
{
class ForceWithVirial;
}

/*! \brief Helper function that wraps calls to posres
 *
 * \returns the potential
 */
real posres_wrapper(gmx::ArrayRef<const int>                  iatoms,
                    gmx::ArrayRef<const t_iparams>            iparamsPosres,
                    const t_pbc&                              pbc,
                    const rvec*                               x,
                    gmx::ArrayRef<const real>                 lambda,
                    const t_forcerec*                         fr,
                    const gmx::ArrayRef<const unsigned short> refScaleComIndices,
                    gmx::ArrayRef<gmx::RVec>                  centersOfMassScaledBuffer,
                    gmx::ArrayRef<gmx::RVec>                  centersOfMassBScaledBuffer,
                    gmx::ArrayRef<rvec4>                      forces,
                    gmx::RVec*                                virial,
                    real*                                     dvdl);

/*! \brief Helper function that wraps calls to posres for free-energy
    pertubation */
void posres_wrapper_lambda(struct gmx_wallcycle*                     wcycle,
                           const InteractionDefinitions&             idef,
                           const t_pbc&                              pbc,
                           const rvec                                x[],
                           gmx_enerdata_t*                           enerd,
                           gmx::ArrayRef<const real>                 lambda,
                           const t_forcerec*                         fr,
                           const gmx::ArrayRef<const unsigned short> refScaleComIndices,
                           gmx::ArrayRef<gmx::RVec>                  centersOfMassScaledBuffer,
                           gmx::ArrayRef<gmx::RVec>                  centersOfMassBScaledBuffer);

/*! \brief Helper function that wraps calls to fbposres for free-energy perturbation
 *
 * \returns the potential
 */
real fbposres_wrapper(gmx::ArrayRef<const int>                  iatoms,
                      gmx::ArrayRef<const t_iparams>            iparamsFBPosres,
                      const t_pbc&                              pbc,
                      const rvec*                               x,
                      const t_forcerec*                         fr,
                      const gmx::ArrayRef<const unsigned short> refScaleComIndices,
                      gmx::ArrayRef<gmx::RVec>                  centersOfMassScaledBuffer,
                      gmx::ArrayRef<rvec4>                      forces,
                      gmx::RVec*                                virial);

#endif
