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

/*! \file
 * \internal
 *
 * \brief Declares and defines the AtomPairlist class.
 *
 * This class in currently only used for perturbed interactions.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_FREEENERGYKERNEL_H
#define GMX_NBNXM_FREEENERGYKERNEL_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct t_forcerec;
struct t_nrnb;
struct interaction_const_t;
namespace gmx
{
template<typename>
class ArrayRef;
template<typename>
class ArrayRefWithPadding;
class AtomPairlist;
class StepWorkload;

/*! \brief The non-bonded free-energy kernel
 *
 * \param[in] nlist      A plain atom pair list containing only perturbed interactions
 * \param[in] coords     The list of local and non-local coordinates
 * \param[in] useSimd    Whether to use the SIMD intrinsics kernel
 * \param[in] ntype      The number of non-bonded atom types
 * \param[in] ic         The non-bonded interactions constants
 * \param[in] shiftvec   The periodic shift vectors
 * \param[in] nbfp       The matrix of LJ parameters between atom types
 * \param[in] nbfp_grid  The matrix of LJ parameters for PME grid corrections
 * \param[in] chargeA    List of charges per atom for state A
 * \param[in] chargeB    List of charges per atom for state B
 * \param[in] typeA      List of LJ atom types per atom for state A
 * \param[in] typeB      List of LJ atom types per atom for state B
 * \param[in] computeForeignLambda  When true, only compute energies and dV/dlambda, igore \p stepWork
 * \param[in] stepWork   Tells what to compute, ignored and can be nullptr when
 *                       \p computeForeignLambda=true
 * \param[in] lambda     List of free-energy lambdas for different components
 * \param[in,out] nrnb   Flop counters to add to
 * \param[in,out] threadForceBuffer  Thread-local force buffer to accumulate into
 * \param[in,out] threadForceShiftBuffer  Thread-local shift force buffer to accumulate into
 * \param[in,out] threadVc  Thread-local Coulomb energy group pair buffer to accumulate into
 * \param[in,out] threadVv  Thread-local VdW energy group pair buffer to accumulate into
 * \param[in,out] threadDvdl  Thread-local dV/dlambda component buffer to accumulate into
 *
 * \throws InvalidInputError when an excluded pair is beyond the rcoulomb with reaction-field.
 */
void gmx_nb_free_energy_kernel(const AtomPairlist&                              nlist,
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
                               bool                                computeForeignLambda,
                               const StepWorkload*                 stepWork,
                               gmx::ArrayRef<const real>           lambda,
                               t_nrnb* gmx_restrict                nrnb,
                               gmx::ArrayRefWithPadding<gmx::RVec> threadForceBuffer,
                               rvec*                               threadForceShiftBuffer,
                               gmx::ArrayRef<real>                 threadVc,
                               gmx::ArrayRef<real>                 threadVv,
                               gmx::ArrayRef<real>                 threadDvdl);

} // namespace gmx

#endif
