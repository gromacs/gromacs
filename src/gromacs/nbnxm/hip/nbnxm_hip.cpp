/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Define HIP implementation (stubs) for GPU execution for NBNXM module
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

void gpu_copy_xq_to_gpu(NbnxmGpu* /* nb */, const nbnxn_atomdata_t* /* nbdata */, gmx::AtomLocality /* aloc */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_launch_kernel(NbnxmGpu* /*    nb */,
                       const gmx::StepWorkload& /* stepWork */,
                       gmx::InteractionLocality /* iloc */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_launch_kernel_pruneonly(NbnxmGpu* /* nb */, gmx::InteractionLocality /* iloc */, int /* numParts */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_launch_cpyback(NbnxmGpu* /* nb */,
                        nbnxn_atomdata_t* /* nbatom */,
                        const gmx::StepWorkload& /* stepWork */,
                        gmx::AtomLocality /* aloc */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

bool gpu_try_finish_task(NbnxmGpu* /* nb */,
                         const gmx::StepWorkload& /* stepWork */,
                         gmx::AtomLocality /* aloc */,
                         real* /* e_lj */,
                         real* /* e_el */,
                         gmx::ArrayRef<gmx::RVec> /* shiftForces */,
                         GpuTaskCompletion /* completionKind */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return false;
}

float gpu_wait_finish_task(NbnxmGpu* /* nb */,
                           const gmx::StepWorkload& /* stepWork */,
                           gmx::AtomLocality /* aloc */,
                           real* /* e_lj */,
                           real* /* e_el */,
                           gmx::ArrayRef<gmx::RVec> /* shiftForces */,
                           gmx_wallcycle* /* wcycle */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return 0.0;
}

void nbnxn_gpu_init_x_to_nbat_x(const GridSet& /* gridSet */, NbnxmGpu* /* gpu_nbv */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void nbnxn_gpu_x_to_nbat_x(const Grid& /* grid */,
                           NbnxmGpu* /* gpu_nbv */,
                           DeviceBuffer<gmx::RVec> /* d_x */,
                           GpuEventSynchronizer* /* xReadyOnDevice */,
                           gmx::AtomLocality /* locality */,
                           int /* gridId */,
                           int /* numColumnsMax */,
                           bool /* mustInsertNonLocalDependency */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void nbnxnInsertNonlocalGpuDependency(NbnxmGpu* /* nb */,
                                      gmx::InteractionLocality /* interactionLocality */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void setupGpuShortRangeWorkLow(NbnxmGpu* /* nb */,
                               const gmx::ListedForcesGpu* /* listedForcesGpu */,
                               gmx::InteractionLocality /* iLocality */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

bool haveGpuShortRangeWork(const NbnxmGpu* /* nb */, gmx::InteractionLocality /* interactionLocality */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return false;
}
} // namespace gmx
