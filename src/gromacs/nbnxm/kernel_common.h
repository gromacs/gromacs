/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

/*! \internal \file
 *
 * \brief
 * Declares the nbnxm pair interaction kernel function types and kind counts, also declares utility functions used in nbnxm_kernel.cpp.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBXNM_KERNEL_COMMON_H
#define GMX_NBXNM_KERNEL_COMMON_H

#include "gromacs/math/vectypes.h"
/* nbnxn_atomdata_t and nbnxn_pairlist_t could be forward declared, but that requires modifications in all SIMD kernel files */
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/utility/real.h"

#include "pairlist.h"

struct interaction_const_t;
enum class CoulombInteractionType : int;
enum class VanDerWaalsType : int;
enum class InteractionModifiers : int;
enum class LongRangeVdW : int;

namespace Nbnxm
{
enum class EwaldExclusionType : int;
}

//! \brief Kinds of electrostatic treatments in SIMD Verlet kernels
enum class CoulombKernelType : int
{
    ReactionField,
    Table,
    TableTwin,
    Ewald,
    EwaldTwin,
    Count
};

//! \brief Whether have a separate cut-off check for VDW interactions
enum class VdwCutoffCheck : int
{
    No,
    Yes
};

//! \brief Kind of Lennard-Jones Ewald treatments in NBNxM SIMD kernels
enum class LJEwald : int
{
    None,
    CombGeometric
};

enum class EnergyOutput : int
{
    None,
    System,
    GroupPairs
};

/*! \brief Pair-interaction kernel type that also calculates energies.
 */
typedef void(NbnxmKernelFunc)(const NbnxnPairlistCpu*    nbl,
                              const nbnxn_atomdata_t*    nbat,
                              const interaction_const_t* ic,
                              const rvec*                shift_vec,
                              nbnxn_atomdata_output_t*   out);

//! \brief Lookup function for Coulomb kernel type
CoulombKernelType getCoulombKernelType(Nbnxm::EwaldExclusionType ewaldExclusionType,
                                       CoulombInteractionType    coulombInteractionType,
                                       bool                      haveEqualCoulombVwdRadii);

/*! \brief Kinds of Van der Waals treatments in NBNxM SIMD kernels
 *
 * The \p LJCUT_COMB refers to the LJ combination rule for the short range.
 * The \p EWALDCOMB refers to the combination rule for the grid part.
 * \p vdwktNR is the number of VdW treatments for the SIMD kernels.
 * \p vdwktNR_ref is the number of VdW treatments for the C reference kernels.
 * These two numbers differ, because currently only the reference kernels
 * support LB combination rules for the LJ-Ewald grid part.
 */
enum
{
    vdwktLJCUT_COMBGEOM,
    vdwktLJCUT_COMBLB,
    vdwktLJCUT_COMBNONE,
    vdwktLJFORCESWITCH,
    vdwktLJPOTSWITCH,
    vdwktLJEWALDCOMBGEOM,
    vdwktLJEWALDCOMBLB,
    vdwktNR = vdwktLJEWALDCOMBLB,
    vdwktNR_ref
};

//! \brief Lookup function for Vdw kernel type
int getVdwKernelType(Nbnxm::KernelType    kernelType,
                     LJCombinationRule    ljCombinationRule,
                     VanDerWaalsType      vanDerWaalsType,
                     InteractionModifiers interactionModifiers,
                     LongRangeVdW         longRangeVdW);

/*! \brief Clears the force buffer.
 *
 * Either the whole buffer is cleared or only the parts used
 * by thread/task \p outputIndex when nbat->bUseBufferFlags is set.
 *
 * \param[in,out] nbat         The Nbnxm atom data
 * \param[in]     outputIndex  The index of the output object to clear
 */
void clearForceBuffer(nbnxn_atomdata_t* nbat, int outputIndex);

/*! \brief Clears the shift forces.
 */
void clear_fshift(real* fshift);

/*! \brief Reduces the collected energy terms over the pair-lists/threads.
 */
void reduce_energies_over_lists(const nbnxn_atomdata_t* nbat, int nlist, real* Vvdw, real* Vc);

#endif
