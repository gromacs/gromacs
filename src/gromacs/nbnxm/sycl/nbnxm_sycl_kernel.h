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
/*! \internal \file
 * \brief
 * Declares nbnxn sycl helper functions
 *
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SYCL_NBNXM_SYCL_KERNEL_H
#define GMX_NBNXM_SYCL_NBNXM_SYCL_KERNEL_H

// Forward declarations
namespace gmx
{
enum class InteractionLocality;
class StepWorkload;
struct NbnxmGpu;

// Ensure any changes are in sync with device_management_sycl.cpp
#define SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_8 (GMX_GPU_NB_CLUSTER_SIZE == 4)
#define SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_32 (GMX_GPU_NB_CLUSTER_SIZE == 8)
#define SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_64 \
    (GMX_GPU_NB_CLUSTER_SIZE == 8 && !(GMX_SYCL_HIPSYCL && !GMX_HIPSYCL_HAVE_HIP_TARGET))

/*! \brief Launch SYCL NBNXM kernel.
 *
 * \param nb Non-bonded parameters.
 * \param stepWork Workload flags for the current step.
 * \param iloc Interaction locality.
 * \param doPrune Whether to do neighborlist pruning.
 */
void launchNbnxmKernel(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc, bool doPrune);

} // namespace gmx

#endif // GMX_NBNXM_SYCL_NBNXM_SYCL_KERNEL_H
