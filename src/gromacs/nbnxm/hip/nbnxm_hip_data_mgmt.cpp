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
 *  \brief Define HIP implementation (stubs) for GPU data transfer for NBNXM module
 *
 *  \author Paul bauer <paul.bauer.q@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

NbnxmGpu* gpu_init(const gmx::DeviceStreamManager& /* deviceStreamManager */,
                   const interaction_const_t* /* ic */,
                   const PairlistParams& /* listParams */,
                   const nbnxn_atomdata_t* /* nbat */,
                   /* true if both local and non-local are done on GPU */
                   bool /* bLocalAndNonlocal */)
{
    return nullptr;
}

void gpu_init_pairlist(NbnxmGpu* /* nb */, const NbnxnPairlistGpu* /* h_nblist */, gmx::InteractionLocality /* iloc */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_init_atomdata(NbnxmGpu* /* nb */, const nbnxn_atomdata_t* /* nbat */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_pme_loadbal_update_param(nonbonded_verlet_t* /* nbv */, const interaction_const_t& /* ic */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_upload_shiftvec(NbnxmGpu* /* nb */, const nbnxn_atomdata_t* /* nbatom */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_clear_outputs(NbnxmGpu* /* nb */, bool /* computeVirial */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

void gpu_free(NbnxmGpu* /* nb */) {}

struct gmx_wallclock_gpu_nbnxn_t* gpu_get_timings(NbnxmGpu* /* nb */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return nullptr;
}

void gpu_reset_timings(struct nonbonded_verlet_t* /* nbv */)
{
    GMX_ASSERT(false, "Not implemented yet");
}

int gpu_min_ci_balanced(NbnxmGpu* /* nb */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return -1;
}

bool gpu_is_kernel_ewald_analytical(const NbnxmGpu* /* nb */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return false;
}

NBAtomDataGpu* gpuGetNBAtomData(NbnxmGpu* /* nb */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return nullptr;
}

DeviceBuffer<gmx::RVec> gpu_get_f(NbnxmGpu* /* nb */)
{
    GMX_ASSERT(false, "Not implemented yet");
    return {};
}

} // namespace gmx
