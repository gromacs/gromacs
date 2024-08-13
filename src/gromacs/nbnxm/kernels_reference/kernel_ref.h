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
 * Declares CPU reference kernels
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_KERNELS_REFERENCE_KERNEL_REF_H
#define GMX_NBNXM_KERNELS_REFERENCE_KERNEL_REF_H

#include "gromacs/nbnxm/kernel_common.h"

namespace gmx
{

//! All the different CPU reference kernel functions.
//! \{
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJ_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJFsw_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJPsw_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJ_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJFsw_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJPsw_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_ref;

NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJ_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJFsw_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJPsw_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJFsw_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJPsw_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_ref;

NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJFsw_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJPsw_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecRF_VdwLJEwCombLB_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJFsw_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJPsw_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_ref;
NbnxmKernelFunc nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_ref;
//! \}

#ifdef INCLUDE_KERNELFUNCTION_TABLES

/*! \brief Declare and define the kernel function pointer lookup tables.
 *
 * The minor index of the array goes over both the LJ combination rules,
 * which is only supported by plain cut-off, and the LJ switch/PME functions.
 * For the C reference kernels, unlike the SIMD kernels, there is not much
 * advantage in using combination rules, so we (re-)use the same kernel.
 */
//! \{
static NbnxmKernelFunc* const nbnxn_kernel_noener_ref[static_cast<int>(CoulombKernelType::Count)][vdwktNR_ref] = {
    { nbnxn_kernel_ElecRF_VdwLJ_F_ref,
      nbnxn_kernel_ElecRF_VdwLJ_F_ref,
      nbnxn_kernel_ElecRF_VdwLJ_F_ref,
      nbnxn_kernel_ElecRF_VdwLJFsw_F_ref,
      nbnxn_kernel_ElecRF_VdwLJPsw_F_ref,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_ref,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_ref },
    { nbnxn_kernel_ElecQSTab_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTab_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTab_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTab_VdwLJFsw_F_ref,
      nbnxn_kernel_ElecQSTab_VdwLJPsw_F_ref,
      nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_F_ref,
      nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_F_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_ref }
};

static NbnxmKernelFunc* const nbnxn_kernel_ener_ref[static_cast<int>(CoulombKernelType::Count)][vdwktNR_ref] = {
    { nbnxn_kernel_ElecRF_VdwLJ_VF_ref,
      nbnxn_kernel_ElecRF_VdwLJ_VF_ref,
      nbnxn_kernel_ElecRF_VdwLJ_VF_ref,
      nbnxn_kernel_ElecRF_VdwLJFsw_VF_ref,
      nbnxn_kernel_ElecRF_VdwLJPsw_VF_ref,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_ref,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_ref },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJFsw_VF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJPsw_VF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VF_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_ref }
};

static NbnxmKernelFunc* const nbnxn_kernel_energrp_ref[static_cast<int>(CoulombKernelType::Count)][vdwktNR_ref] = {
    { nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecRF_VdwLJFsw_VgrpF_ref,
      nbnxn_kernel_ElecRF_VdwLJPsw_VgrpF_ref,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VgrpF_ref,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_VgrpF_ref },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJFsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJPsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VgrpF_ref,
      nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VgrpF_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_ref },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_ref,
      nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_ref }
};
//! \}

#endif /* INCLUDE_KERNELFUNCTION_TABLES */

} // namespace gmx

#endif
