/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief
 * Declares the nbnxn pair interaction kernel dispatcher.
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef _nbnxn_kernel_cpu_h
#define _nbnxn_kernel_cpu_h

#include "gromacs/fda/FDA.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct interaction_const_t;
struct nbnxn_atomdata_t;
struct nonbonded_verlet_group_t;

/*! \brief Dispatches the non-bonded N versus M atom cluster CPU kernels.
 *
 * OpenMP parallelization is performed within this function.
 * Energy reduction, but not force and shift force reduction, is performed
 * within this function.
 *
 * \param[in,out] nbvg          The group (local/non-local) to compute interaction for
 * \param[in]     nbat          The atomdata for the interactions
 * \param[in]     ic            Non-bonded interaction constants
 * \param[in]     shiftVectors  The PBC shift vectors
 * \param[in]     forceFlags    Flags that tell what to compute
 * \param[in]     clearF        Enum that tells if to clear the force output buffer
 * \param[out]    fshift        Shift force output buffer
 * \param[out]    vCoulomb      Output buffer for Coulomb energies
 * \param[out]    vVdw          Output buffer for Van der Waals energies
 */
void
nbnxn_kernel_cpu(nonbonded_verlet_group_t  *nbvg,
                 const nbnxn_atomdata_t    *nbat,
                 const interaction_const_t *ic,
                 rvec                      *shiftVectors,
                 int                        forceFlags,
                 int                        clearF,
                 real                      *fshift,
                 real                      *vCoulomb,
                 real                      *vVdw
#ifdef BUILD_WITH_FDA
                 ,
                 FDA                       *fda,
                 int                       *cellInv
#endif
                );

#endif
