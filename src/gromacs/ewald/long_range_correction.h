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
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Erik Lindahl <erik@kth.se>
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_LONG_RANGE_CORRECTION_H
#define GMX_EWALD_LONG_RANGE_CORRECTION_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct t_commrec;
struct t_forcerec;
struct t_inputrec;
enum class EwaldGeometry : int;

namespace gmx
{
template<typename>
class ArrayRef;
}

/*! \brief Calculate long-range Ewald correction terms.
 *
 * Calculate correction for electrostatic surface dipole terms.
 */
void ewald_LRcorrection(int                            numAtomsLocal,
                        const t_commrec*               commrec,
                        int                            numThreads,
                        int                            thread,
                        real                           epsilonR,
                        gmx::ArrayRef<const double>    qsum,
                        EwaldGeometry                  ewaldGeometry,
                        real                           epsilonSurface,
                        bool                           havePbcXY2Walls,
                        real                           wallEwaldZfac,
                        gmx::ArrayRef<const real>      chargeA,
                        gmx::ArrayRef<const real>      chargeB,
                        bool                           bHaveChargePerturbed,
                        gmx::ArrayRef<const gmx::RVec> coords,
                        const matrix                   box,
                        gmx::ArrayRef<const gmx::RVec> mu_tot,
                        gmx::ArrayRef<gmx::RVec>       forces,
                        real*                          Vcorr_q,
                        real                           lambda_q,
                        real*                          dvdlambda_q);

#endif
