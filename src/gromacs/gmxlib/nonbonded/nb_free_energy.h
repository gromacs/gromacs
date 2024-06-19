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

#ifndef GMX_GMXLIB_NONBONDED_NB_FREE_ENERGY_H
#define GMX_GMXLIB_NONBONDED_NB_FREE_ENERGY_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct t_forcerec;
struct t_nrnb;
struct t_nblist;
struct interaction_const_t;
namespace gmx
{
template<typename>
class ArrayRef;
template<typename>
class ArrayRefWithPadding;
} // namespace gmx

/*! \brief The non-bonded free-energy kernel
 *
 * Note that this uses a regular atom pair, not cluster pair, list.
 *
 * \throws InvalidInputError when an excluded pair is beyond the rcoulomb with reaction-field.
 */
void gmx_nb_free_energy_kernel(const t_nblist&                                  nlist,
                               const gmx::ArrayRefWithPadding<const gmx::RVec>& coords,
                               bool                                             useSimd,
                               int                                              ntype,
                               const interaction_const_t&                       ic,
                               gmx::ArrayRef<const gmx::RVec>                   shiftvec,
                               gmx::ArrayRef<const real>                        nbfp,
                               gmx::ArrayRef<const real>                        nbfp_grid,
                               gmx::ArrayRef<const real>                        chargeA,
                               gmx::ArrayRef<const real>                        chargeB,
                               gmx::ArrayRef<const int>                         typeA,
                               gmx::ArrayRef<const int>                         typeB,
                               int                                              flags,
                               gmx::ArrayRef<const real>                        lambda,
                               t_nrnb* gmx_restrict                             nrnb,
                               gmx::ArrayRefWithPadding<gmx::RVec>              threadForceBuffer,
                               rvec*               threadForceShiftBuffer,
                               gmx::ArrayRef<real> threadVc,
                               gmx::ArrayRef<real> threadVv,
                               gmx::ArrayRef<real> threadDvdl);

#endif
