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
 * \brief Implements high-level PME GPU functions which do not require GPU framework-specific code.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "config.h"

#include <list>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_coordinate_receiver_gpu.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "pme_gpu_internal.h"
#include "pme_gpu_settings.h"
#include "pme_gpu_timings.h"
#include "pme_gpu_types_host.h"
#include "pme_grid.h"
#include "pme_internal.h"
#include "pme_solve.h"

/*! \brief
 * Finds out if PME is currently running on GPU.
 *
 * \todo The GPU module should not be constructed (or at least called)
 * when it is not active, so there should be no need to check whether
 * it is active. An assertion that this is true makes sense.
 *
 * \param[in] pme  The PME structure.
 * \returns        True if PME runs on GPU currently, false otherwise.
 */
static inline bool pme_gpu_active(const gmx_pme_t* pme)
{
    return (pme != nullptr) && (pme->runMode != PmeRunMode::CPU);
}

void pme_gpu_reset_timings(const gmx_pme_t* pme)
{
    if (pme_gpu_active(pme))
    {
        pme_gpu_reset_timings(pme->gpu);
    }
}

void pme_gpu_get_timings(const gmx_pme_t* pme, gmx_wallclock_gpu_pme_t* timings)
{
    if (pme_gpu_active(pme))
    {
        pme_gpu_get_timings(pme->gpu, timings);
    }
}

int pme_gpu_get_block_size(const gmx_pme_t* pme)
{

    if (!pme || !pme_gpu_active(pme))
    {
        return 0;
    }
    else
    {
        return pme_gpu_get_atom_data_block_size();
    }
}

/*! \brief
 * A convenience wrapper for launching either the GPU or CPU FFT.
 *
 * \param[in] pme            The PME structure.
 * \param[in] gridIndex      The grid index - should currently always be 0.
 * \param[in] dir            The FFT direction enum.
 * \param[in] wcycle         The wallclock counter.
 */
void inline parallel_3dfft_execute_gpu_wrapper(gmx_pme_t*             pme,
                                               const int              gridIndex,
                                               enum gmx_fft_direction dir,
                                               gmx_wallcycle*         wcycle)
{
    if (pme_gpu_settings(pme->gpu).performGPUFFT)
    {
        // use a separate sub-counter for GPU FFT launch
        // this is specially important for PME decomposition where distributed FFT
        // implementations are used
        wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);
        wallcycle_sub_start(wcycle, WallCycleSubCounter::LaunchGpuPmeFft);
        pme_gpu_3dfft(pme->gpu, dir, gridIndex);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuPmeFft);
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
    }
    else
    {
        wallcycle_start(wcycle, WallCycleCounter::PmeFftMixedMode);
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
        for (int thread = 0; thread < pme->nthread; thread++)
        {
            gmx_parallel_3dfft_execute(pme->gridsCoulomb[gridIndex].pfft_setup.get(), dir, thread, wcycle);
        }
        wallcycle_stop(wcycle, WallCycleCounter::PmeFftMixedMode);
    }
}

/* The PME computation code split into a few separate functions. */

void pme_gpu_prepare_computation(gmx_pme_t*               pme,
                                 const matrix             box,
                                 gmx_wallcycle*           wcycle,
                                 const gmx::StepWorkload& stepWork)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");
    GMX_ASSERT(pme->nnodes > 0, "");
    GMX_ASSERT(pme->nnodes == 1 || pme->ndecompdim > 0, "");

    PmeGpu* pmeGpu = pme->gpu;
    // TODO these flags are only here to honor the CPU PME code, and probably should be removed
    pmeGpu->settings.useGpuForceReduction = stepWork.useGpuPmeFReduction;

    bool shouldUpdateBox = false;
    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            shouldUpdateBox |= (pmeGpu->common->previousBox[i][j] != box[i][j]);
            pmeGpu->common->previousBox[i][j] = box[i][j];
        }
    }

    if (stepWork.haveDynamicBox || shouldUpdateBox) // || is to make the first computation always update
    {
        wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);
        pme_gpu_update_input_box(pmeGpu, box);
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);

        if (!pme_gpu_settings(pmeGpu).performGPUSolve)
        {
            // TODO remove code duplication and add test coverage
            matrix scaledBox;
            pmeGpu->common->boxScaler->scaleBox(box, scaledBox);
            gmx::invertBoxMatrix(scaledBox, pme->recipbox);
            pme->boxVolume = scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ];
        }
    }
}

void pme_gpu_launch_spread(gmx_pme_t*                     pme,
                           GpuEventSynchronizer*          xReadyOnDevice,
                           gmx_wallcycle*                 wcycle,
                           const real                     lambdaQ,
                           const bool                     useGpuDirectComm,
                           gmx::PmeCoordinateReceiverGpu* pmeCoordinateReceiverGpu,
                           const bool                     useMdGpuGraph)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");
    GMX_ASSERT(xReadyOnDevice || !pme->bPPnode, "Need a valid xReadyOnDevice on PP+PME ranks.");
    GMX_ASSERT(pme->doCoulomb, "Only Coulomb PME can be run on GPU.");

    PmeGpu* pmeGpu = pme->gpu;

    GMX_ASSERT(pmeGpu->common->ngrids == 1 || (pmeGpu->common->ngrids == 2 && pme->bFEP_q),
               "If not decoupling Coulomb interactions there should only be one FEP grid. If "
               "decoupling Coulomb interactions there should be two grids.");

    /* PME on GPU can currently manage two grids:
     * grid_index=0: Coulomb PME with charges in the normal state or from FEP state A.
     * grid_index=1: Coulomb PME with charges from FEP state B.
     */
    /* Spread the coefficients on a grid */
    const bool computeSplines = true;
    const bool spreadCharges  = true;
    pme_gpu_spread(pmeGpu,
                   xReadyOnDevice,
                   pme->gridsCoulomb,
                   computeSplines,
                   spreadCharges,
                   lambdaQ,
                   useGpuDirectComm,
                   pmeCoordinateReceiverGpu,
                   useMdGpuGraph,
                   wcycle);
}

void pme_gpu_launch_complex_transforms(gmx_pme_t* pme, gmx_wallcycle* wcycle, const gmx::StepWorkload& stepWork)
{
    PmeGpu*     pmeGpu   = pme->gpu;
    const auto& settings = pmeGpu->settings;
    // There's no support for computing energy without virial, or vice versa
    const bool computeEnergyAndVirial = stepWork.computeEnergy || stepWork.computeVirial;

    if (!settings.performGPUFFT)
    {
        wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmeGridD2hCopy);
        pme_gpu_sync_spread_grid(pme->gpu);
        wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmeGridD2hCopy);
    }

    try
    {
        /* The 3dffts and the solve are done in a loop to simplify things, even if this means that
         * there will be two kernel launches for solve. */
        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            /* do R2C 3D-FFT */
            t_complex* cfftgrid = pme->gridsCoulomb[gridIndex].cfftgrid;
            parallel_3dfft_execute_gpu_wrapper(pme, gridIndex, GMX_FFT_REAL_TO_COMPLEX, wcycle);

            /* solve in k-space for our local cells */
            if (settings.performGPUSolve)
            {
                const auto gridOrdering = GridOrdering::XYZ;
                wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);
                pme_gpu_solve(pmeGpu, gridIndex, cfftgrid, gridOrdering, computeEnergyAndVirial);
                wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
            }
            else
            {
                wallcycle_start(wcycle, WallCycleCounter::PmeSolveMixedMode);
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
                for (int thread = 0; thread < pme->nthread; thread++)
                {
                    pme->pmeSolve->solveCoulombYZX(
                            *pme, cfftgrid, pme->boxVolume, computeEnergyAndVirial, thread);
                }
                wallcycle_stop(wcycle, WallCycleCounter::PmeSolveMixedMode);
            }

            parallel_3dfft_execute_gpu_wrapper(pme, gridIndex, GMX_FFT_COMPLEX_TO_REAL, wcycle);
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void pme_gpu_launch_gather(gmx_pme_t* pme, gmx_wallcycle gmx_unused* wcycle, const real lambdaQ, const bool computeVirial)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    if (!pme_gpu_settings(pme->gpu).performGPUGather)
    {
        return;
    }

    pme_gpu_gather(pme->gpu, pme->gridsCoulomb, lambdaQ, wcycle, computeVirial);
}

//! Accumulate the \c forcesToAdd to \c f, using the available threads.
static void sum_forces(gmx::ArrayRef<gmx::RVec> f, gmx::ArrayRef<const gmx::RVec> forceToAdd)
{
    const int end = forceToAdd.size();

    int gmx_unused nt = gmx_omp_nthreads_get(ModuleMultiThread::Pme);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (int i = 0; i < end; i++)
    {
        f[i] += forceToAdd[i];
    }
}

//! Reduce quantities from \c output to \c forceWithVirial and \c enerd.
static void pme_gpu_reduce_outputs(const bool            computeEnergyAndVirial,
                                   const PmeOutput&      output,
                                   gmx_wallcycle*        wcycle,
                                   gmx::ForceWithVirial* forceWithVirial,
                                   gmx_enerdata_t*       enerd)
{
    wallcycle_start(wcycle, WallCycleCounter::PmeGpuFReduction);
    GMX_ASSERT(forceWithVirial, "Invalid force pointer");

    if (computeEnergyAndVirial)
    {
        GMX_ASSERT(enerd, "Invalid energy output manager");
        forceWithVirial->addVirialContribution(output.coulombVirial_);
        enerd->term[F_COUL_RECIP] += output.coulombEnergy_;
        enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Coul] += output.coulombDvdl_;
    }
    if (output.haveForceOutput_)
    {
        sum_forces(forceWithVirial->force_, output.forces_);
    }
    wallcycle_stop(wcycle, WallCycleCounter::PmeGpuFReduction);
}

bool pme_gpu_try_finish_task(gmx_pme_t*               pme,
                             const gmx::StepWorkload& stepWork,
                             gmx_wallcycle*           wcycle,
                             gmx::ForceWithVirial*    forceWithVirial,
                             gmx_enerdata_t*          enerd,
                             const real               lambdaQ,
                             GpuTaskCompletion        completionKind)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");
    GMX_ASSERT(!pme->gpu->settings.useGpuForceReduction,
               "GPU force reduction should not be active on the pme_gpu_try_finish_task() path");

    // First, if possible, check whether all tasks on the stream have
    // completed, and return fast if not. Accumulate to wcycle the
    // time needed for that checking, but do not yet record that the
    // gather has occurred.
    bool           needToSynchronize      = true;
    constexpr bool c_streamQuerySupported = GMX_GPU_CUDA;

    // TODO: implement c_streamQuerySupported with an additional GpuEventSynchronizer per stream (#2521)
    if ((completionKind == GpuTaskCompletion::Check) && c_streamQuerySupported)
    {
        wallcycle_start_nocount(wcycle, WallCycleCounter::WaitGpuPmeGather);
        // Query the PME stream for completion of all tasks enqueued and
        // if we're not done, stop the timer before early return.
        const bool pmeGpuDone = pme_gpu_stream_query(pme->gpu);
        wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmeGather);

        if (!pmeGpuDone)
        {
            return false;
        }
        needToSynchronize = false;
    }

    wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmeGather);
    // If the above check passed, then there is no need to make an
    // explicit synchronization call.
    if (needToSynchronize)
    {
        // Synchronize the whole PME stream at once, including D2H result transfers.
        pme_gpu_synchronize(pme->gpu);
    }
    pme_gpu_update_timings(pme->gpu);
    // There's no support for computing energy without virial, or vice versa
    const bool computeEnergyAndVirial = stepWork.computeEnergy || stepWork.computeVirial;
    PmeOutput  output                 = pme_gpu_getOutput(
            pme, computeEnergyAndVirial, pme->gpu->common->ngrids > 1 ? lambdaQ : 1.0);
    wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmeGather);

    GMX_ASSERT(pme->gpu->settings.useGpuForceReduction == !output.haveForceOutput_,
               "When forces are reduced on the CPU, there needs to be force output");
    pme_gpu_reduce_outputs(computeEnergyAndVirial, output, wcycle, forceWithVirial, enerd);

    return true;
}

// This is used by PME-only ranks
PmeOutput pme_gpu_wait_finish_task(gmx_pme_t*     pme,
                                   const bool     computeEnergyAndVirial,
                                   const real     lambdaQ,
                                   gmx_wallcycle* wcycle)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmeGather);

    // Synchronize the whole PME stream at once, including D2H result transfers
    // if there are outputs we need to wait for at this step; we still call getOutputs
    // for uniformity and because it sets the PmeOutput.haveForceOutput_.
    if (!pme->gpu->settings.useGpuForceReduction || computeEnergyAndVirial)
    {
        pme_gpu_synchronize(pme->gpu);
    }

    PmeOutput output = pme_gpu_getOutput(
            pme, computeEnergyAndVirial, pme->gpu->common->ngrids > 1 ? lambdaQ : 1.0);
    wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmeGather);
    return output;
}

// This is used when not using the alternate-waiting reduction
void pme_gpu_wait_and_reduce(gmx_pme_t*               pme,
                             const gmx::StepWorkload& stepWork,
                             gmx_wallcycle*           wcycle,
                             gmx::ForceWithVirial*    forceWithVirial,
                             gmx_enerdata_t*          enerd,
                             const real               lambdaQ)
{
    // There's no support for computing energy without virial, or vice versa
    const bool computeEnergyAndVirial = stepWork.computeEnergy || stepWork.computeVirial;
    PmeOutput  output                 = pme_gpu_wait_finish_task(
            pme, computeEnergyAndVirial, pme->gpu->common->ngrids > 1 ? lambdaQ : 1.0, wcycle);
    GMX_ASSERT(pme->gpu->settings.useGpuForceReduction == !output.haveForceOutput_,
               "When forces are reduced on the CPU, there needs to be force output");
    pme_gpu_reduce_outputs(computeEnergyAndVirial, output, wcycle, forceWithVirial, enerd);
}

void pme_gpu_reinit_computation(const gmx_pme_t* pme, const bool gpuGraphWithSeparatePmeRank, gmx_wallcycle* wcycle)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

    pme_gpu_update_timings(pme->gpu);

    pme_gpu_clear_grids(pme->gpu);
    pme_gpu_clear_energy_virial(pme->gpu, gpuGraphWithSeparatePmeRank);

    wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
}

DeviceBuffer<gmx::RVec> pme_gpu_get_device_f(const gmx_pme_t* pme)
{
    if (!pme || !pme_gpu_active(pme))
    {
        return DeviceBuffer<gmx::RVec>{};
    }
    return pme_gpu_get_kernelparam_forces(pme->gpu);
}

void pme_gpu_set_device_x(const gmx_pme_t* pme, DeviceBuffer<gmx::RVec> d_x)
{
    GMX_ASSERT(pme != nullptr, "Null pointer is passed as a PME to the set coordinates function.");
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    pme_gpu_set_kernelparam_coordinates(pme->gpu, d_x);
}

GpuEventSynchronizer* pme_gpu_get_f_ready_synchronizer(const gmx_pme_t* pme)
{
    if (!pme || !pme_gpu_active(pme))
    {
        return nullptr;
    }

    return pme_gpu_get_forces_ready_synchronizer(pme->gpu);
}

void pme_gpu_use_nvshmem(PmeGpu* pmeGpu, bool useNvshmem)
{
    pmeGpu->useNvshmem = useNvshmem;
    pme_gpu_set_kernelparam_useNvshmem(pmeGpu, useNvshmem);
    if (useNvshmem)
    {
        pmeGpu->nvshmemParams = std::make_unique<PmeNvshmemHost>();
    }
}
