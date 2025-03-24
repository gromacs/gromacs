/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief This file contains internal implementation stubs
 * for performing the PME calculations on GPU.
 *
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "config.h"

#include "gromacs/utility/gmxassert.h"

#include "pme_gpu_internal.h"

int pme_gpu_get_atoms_per_warp(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
    return 0;
}

void pme_gpu_synchronize(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_spread(PmeGpu* /* pmeGpu */,
                    GpuEventSynchronizer* /* xReadyOnDevice */,
                    gmx::ArrayRef<PmeAndFftGrids> /* h_grids */,
                    bool /* computeSplines */,
                    bool /* spreadCharges */,
                    real /* lambda */,
                    bool /* useGpuDirectComm */,
                    gmx::PmeCoordinateReceiverGpu* /* pmeCoordinateReceiverGpu */,
                    bool /* useMdGpuGraph */,
                    gmx_wallcycle* /* wcycle */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_solve(PmeGpu* /* pmeGpu */,
                   int /* gridIndex */,
                   t_complex* /* h_grid */,
                   GridOrdering /* gridOrdering */,
                   bool /* computeEnergyAndVirial */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_gather(PmeGpu* /* pmeGpu */,
                    gmx::ArrayRef<PmeAndFftGrids> /* h_grids */,
                    float /* lambda */,
                    gmx_wallcycle* /* wcycle */,
                    bool /* computeVirial */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_set_kernelparam_coordinates(const PmeGpu* /* pmeGpu */, DeviceBuffer<gmx::RVec> /* d_x */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

DeviceBuffer<gmx::RVec> pme_gpu_get_kernelparam_forces(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
    return {};
}

void pme_gpu_set_kernelparam_useNvshmem(const PmeGpu* /* pmeGpu */, bool /* useNvshmem */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

GpuEventSynchronizer* pme_gpu_get_forces_ready_synchronizer(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
    return nullptr;
}

void pme_gpu_getEnergyAndVirial(const gmx_pme_t& /* pme */, float /* lambda */, PmeOutput* /* output */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

PmeOutput pme_gpu_getOutput(gmx_pme_t* /* pme */, bool /* computeEnergyAndVirial */, real /* lambdaQ */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
    return {};
}

void pme_gpu_update_input_box(PmeGpu* /* pmeGpu */, const matrix /* box */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_get_real_grid_sizes(const PmeGpu* /* pmeGpu */, gmx::IVec* /* gridSize */, gmx::IVec* /* paddedGridSize */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_reinit(gmx_pme_t* /* pme */,
                    const DeviceContext* /* deviceContext */,
                    const DeviceStream* /* deviceStream */,
                    const PmeGpuProgram* /* pmeGpuProgram */,
                    bool /* useMdGpuGraph */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_destroy(PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_reinit_atoms(PmeGpu* /* pmeGpu */,
                          int /* nAtoms */,
                          const real* /* chargesA */,
                          const real* /* chargesB = nullptr */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_sync_spread_grid(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_3dfft(const PmeGpu* /* pmeGpu */, enum gmx_fft_direction /* direction */, int /* gridIndex = 0 */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_clear_grids(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_clear_energy_virial(const PmeGpu* /* pmeGpu */, bool /* gpuGraphWithSeparatePmeRank */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

int pme_gpu_get_atom_data_block_size()
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
    return -1;
}

void pme_gpu_update_timings(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_reinit_timings(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_reset_timings(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

void pme_gpu_get_timings(const PmeGpu* /* pmeGpu */, gmx_wallclock_gpu_pme_t* /* timings */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
}

bool pme_gpu_stream_query(const PmeGpu* /* pmeGpu */)
{
    GMX_RELEASE_ASSERT(false, "HIP PME not implemented yet");
    return false;
}
