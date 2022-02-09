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
/*! \libinternal \defgroup module_ewald Ewald-family treatments of long-ranged forces
 * \ingroup group_mdrun
 *
 * \brief Computes energies and forces for long-ranged interactions
 * using the Ewald decomposition. Includes plain Ewald, PME, P3M for
 * Coulomb, PME for Lennard-Jones, load-balancing for PME, and
 * supporting code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Roland Schulz <roland@rschulz.eu>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Christan Wennberg <cwennberg@kth.se>
 */

/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the plain-Ewald long-ranged part,
 * and the correction for overall system charge for all Ewald-family
 * methods.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_EWALD_H
#define GMX_EWALD_EWALD_H

#include <stdio.h>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_complex;
enum class FreeEnergyPerturbationType : int;

namespace gmx
{
template<typename>
class ArrayRef;
}

struct gmx_ewald_tab_t
{
    gmx_ewald_tab_t(const t_inputrec& ir, FILE* fp);

    ~gmx_ewald_tab_t();

    int nx;
    int ny;
    int nz;
    int kmax;

    std::vector<t_complex> tab_xy;
    std::vector<t_complex> tab_qxyz;
};

/*! \brief Do the long-ranged part of an Ewald calculation */
real do_ewald(bool                           havePbcXY2Walls,
              real                           wallEwaldZfac,
              real                           epsilonR,
              FreeEnergyPerturbationType     freeEnergyPerturbationType,
              gmx::ArrayRef<const gmx::RVec> coords,
              gmx::ArrayRef<gmx::RVec>       forces,
              gmx::ArrayRef<const real>      chargeA,
              gmx::ArrayRef<const real>      chargeB,
              const matrix                   box,
              const t_commrec*               commrec,
              int                            natoms,
              matrix                         lrvir,
              real                           ewaldcoeff,
              real                           lambda,
              real*                          dvdlambda,
              gmx_ewald_tab_t*               et);

/*! \brief Calculate the correction to the Ewald sum, due to a net system
 * charge.
 *
 * Should only be called on one thread. */
real ewald_charge_correction(const t_commrec*            commrec,
                             real                        epsilonR,
                             real                        ewaldcoeffQ,
                             gmx::ArrayRef<const double> qsum,
                             real                        lambda,
                             const matrix                box,
                             real*                       dvdlambda,
                             tensor                      vir);

#endif
