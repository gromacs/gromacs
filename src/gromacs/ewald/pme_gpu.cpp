/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
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
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "pme_gpu_internal.h"
#include "pme_grid.h"
#include "pme_internal.h"
#include "pme_solve.h"

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

int pme_gpu_get_padding_size(const gmx_pme_t* pme)
{

    if (!pme || !pme_gpu_active(pme))
    {
        return 0;
    }
    else
    {
        return pme_gpu_get_atom_data_alignment(pme->gpu);
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
                                               gmx_wallcycle_t        wcycle)
{
    GMX_ASSERT(gridIndex == 0, "Only single grid supported");
    if (pme_gpu_performs_FFT(pme->gpu))
    {
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_PME);
        pme_gpu_3dfft(pme->gpu, dir, gridIndex);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);
    }
    else
    {
        wallcycle_start(wcycle, ewcPME_FFT_MIXED_MODE);
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
        for (int thread = 0; thread < pme->nthread; thread++)
        {
            gmx_parallel_3dfft_execute(pme->pfft_setup[gridIndex], dir, thread, wcycle);
        }
        wallcycle_stop(wcycle, ewcPME_FFT_MIXED_MODE);
    }
}

/* The PME computation code split into a few separate functions. */

void pme_gpu_prepare_computation(gmx_pme_t*     pme,
                                 bool           needToUpdateBox,
                                 const matrix   box,
                                 gmx_wallcycle* wcycle,
                                 int            flags,
                                 bool           useGpuForceReduction)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");
    GMX_ASSERT(pme->nnodes > 0, "");
    GMX_ASSERT(pme->nnodes == 1 || pme->ndecompdim > 0, "");

    PmeGpu* pmeGpu                = pme->gpu;
    pmeGpu->settings.currentFlags = flags;
    // TODO these flags are only here to honor the CPU PME code, and probably should be removed
    pmeGpu->settings.useGpuForceReduction = useGpuForceReduction;

    bool shouldUpdateBox = false;
    for (int i = 0; i < DIM; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            shouldUpdateBox |= (pmeGpu->common->previousBox[i][j] != box[i][j]);
            pmeGpu->common->previousBox[i][j] = box[i][j];
        }
    }

    if (needToUpdateBox || shouldUpdateBox) // || is to make the first computation always update
    {
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_PME);
        pme_gpu_update_input_box(pmeGpu, box);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);

        if (!pme_gpu_performs_solve(pmeGpu))
        {
            // TODO remove code duplication and add test coverage
            matrix scaledBox;
            pmeGpu->common->boxScaler->scaleBox(box, scaledBox);
            gmx::invertBoxMatrix(scaledBox, pme->recipbox);
            pme->boxVolume = scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ];
        }
    }
}

void pme_gpu_launch_spread(gmx_pme_t* pme, GpuEventSynchronizer* xReadyOnDevice, gmx_wallcycle* wcycle)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");
    GMX_ASSERT(xReadyOnDevice || !pme->bPPnode || (GMX_GPU != GMX_GPU_CUDA),
               "Need a valid xReadyOnDevice on PP+PME ranks with CUDA.");

    PmeGpu* pmeGpu = pme->gpu;

    const unsigned int gridIndex = 0;
    real*              fftgrid   = pme->fftgrid[gridIndex];
    if (pmeGpu->settings.currentFlags & GMX_PME_SPREAD)
    {
        /* Spread the coefficients on a grid */
        const bool computeSplines = true;
        const bool spreadCharges  = true;
        wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
        wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_PME);
        pme_gpu_spread(pmeGpu, xReadyOnDevice, gridIndex, fftgrid, computeSplines, spreadCharges);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME);
        wallcycle_stop(wcycle, ewcLAUNCH_GPU);
    }
}

void pme_gpu_launch_complex_transforms(gmx_pme_t* pme, gmx_wallcycle* wcycle)
{
    PmeGpu*    pmeGpu                 = pme->gpu;
    const bool computeEnergyAndVirial = (pmeGpu->settings.currentFlags & GMX_PME_CALC_ENER_VIR) != 0;
    const bool performBackFFT = (pmeGpu->settings.currentFlags & (GMX_PME_CALC_F | GMX_PME_CALC_POT)) != 0;
    const unsigned int gridIndex = 0;
    t_complex*         cfftgrid  = pme->cfftgrid[gridIndex];

    if (pmeGpu->settings.currentFlags & GMX_PME_SPREAD)
    {
        if (!pme_gpu_performs_FFT(pmeGpu))
        {
            wallcycle_start(wcycle, ewcWAIT_GPU_PME_SPREAD);
            pme_gpu_sync_spread_grid(pme->gpu);
            wallcycle_stop(wcycle, ewcWAIT_GPU_PME_SPREAD);
        }
    }

    try
    {
        if (pmeGpu->settings.currentFlags & GMX_PME_SOLVE)
        {
            /* do R2C 3D-FFT */
            parallel_3dfft_execute_gpu_wrapper(pme, gridIndex, GMX_FFT_REAL_TO_COMPLEX, wcycle);

            /* solve in k-space for our local cells */
            if (pme_gpu_performs_solve(pmeGpu))
            {
                const auto gridOrdering = pme_gpu_uses_dd(pmeGpu) ? GridOrdering::YZX : GridOrdering::XYZ;
                wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
                wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_PME);
                pme_gpu_solve(pmeGpu, cfftgrid, gridOrdering, computeEnergyAndVirial);
                wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME);
                wallcycle_stop(wcycle, ewcLAUNCH_GPU);
            }
            else
            {
                wallcycle_start(wcycle, ewcPME_SOLVE_MIXED_MODE);
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
                for (int thread = 0; thread < pme->nthread; thread++)
                {
                    solve_pme_yzx(pme, cfftgrid, pme->boxVolume, computeEnergyAndVirial,
                                  pme->nthread, thread);
                }
                wallcycle_stop(wcycle, ewcPME_SOLVE_MIXED_MODE);
            }
        }

        if (performBackFFT)
        {
            parallel_3dfft_execute_gpu_wrapper(pme, gridIndex, GMX_FFT_COMPLEX_TO_REAL, wcycle);
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void pme_gpu_launch_gather(const gmx_pme_t* pme, gmx_wallcycle gmx_unused* wcycle, PmeForceOutputHandling forceTreatment)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    if (!pme_gpu_performs_gather(pme->gpu))
    {
        return;
    }

    wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
    wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_PME);
    const unsigned int gridIndex = 0;
    real*              fftgrid   = pme->fftgrid[gridIndex];
    pme_gpu_gather(pme->gpu, forceTreatment, reinterpret_cast<float*>(fftgrid));
    wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME);
    wallcycle_stop(wcycle, ewcLAUNCH_GPU);
}

//! Accumulate the \c forcesToAdd to \c f, using the available threads.
static void sum_forces(gmx::ArrayRef<gmx::RVec> f, gmx::ArrayRef<const gmx::RVec> forceToAdd)
{
    const int end = forceToAdd.size();

    int gmx_unused nt = gmx_omp_nthreads_get(emntPME);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (int i = 0; i < end; i++)
    {
        f[i] += forceToAdd[i];
    }
}

//! Reduce quantities from \c output to \c forceWithVirial and \c enerd.
static void pme_gpu_reduce_outputs(const int             flags,
                                   const PmeOutput&      output,
                                   gmx_wallcycle*        wcycle,
                                   gmx::ForceWithVirial* forceWithVirial,
                                   gmx_enerdata_t*       enerd)
{
    wallcycle_start(wcycle, ewcPME_GPU_F_REDUCTION);
    GMX_ASSERT(forceWithVirial, "Invalid force pointer");

    const bool haveComputedEnergyAndVirial = (flags & GMX_PME_CALC_ENER_VIR) != 0;
    if (haveComputedEnergyAndVirial)
    {
        GMX_ASSERT(enerd, "Invalid energy output manager");
        forceWithVirial->addVirialContribution(output.coulombVirial_);
        enerd->term[F_COUL_RECIP] += output.coulombEnergy_;
    }
    if (output.haveForceOutput_)
    {
        sum_forces(forceWithVirial->force_, output.forces_);
    }
    wallcycle_stop(wcycle, ewcPME_GPU_F_REDUCTION);
}

bool pme_gpu_try_finish_task(gmx_pme_t*            pme,
                             const int             flags,
                             gmx_wallcycle*        wcycle,
                             gmx::ForceWithVirial* forceWithVirial,
                             gmx_enerdata_t*       enerd,
                             GpuTaskCompletion     completionKind)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");
    GMX_ASSERT(!pme->gpu->settings.useGpuForceReduction,
               "GPU force reduction should not be active on the pme_gpu_try_finish_task() path");

    // First, if possible, check whether all tasks on the stream have
    // completed, and return fast if not. Accumulate to wcycle the
    // time needed for that checking, but do not yet record that the
    // gather has occured.
    bool           needToSynchronize      = true;
    constexpr bool c_streamQuerySupported = (GMX_GPU == GMX_GPU_CUDA);
    // TODO: implement c_streamQuerySupported with an additional GpuEventSynchronizer per stream (#2521)
    if ((completionKind == GpuTaskCompletion::Check) && c_streamQuerySupported)
    {
        wallcycle_start_nocount(wcycle, ewcWAIT_GPU_PME_GATHER);
        // Query the PME stream for completion of all tasks enqueued and
        // if we're not done, stop the timer before early return.
        const bool pmeGpuDone = pme_gpu_stream_query(pme->gpu);
        wallcycle_stop(wcycle, ewcWAIT_GPU_PME_GATHER);

        if (!pmeGpuDone)
        {
            return false;
        }
        needToSynchronize = false;
    }

    wallcycle_start(wcycle, ewcWAIT_GPU_PME_GATHER);
    // If the above check passed, then there is no need to make an
    // explicit synchronization call.
    if (needToSynchronize)
    {
        // Synchronize the whole PME stream at once, including D2H result transfers.
        pme_gpu_synchronize(pme->gpu);
    }
    pme_gpu_update_timings(pme->gpu);
    PmeOutput output = pme_gpu_getOutput(*pme, flags);
    wallcycle_stop(wcycle, ewcWAIT_GPU_PME_GATHER);

    GMX_ASSERT(pme->gpu->settings.useGpuForceReduction == !output.haveForceOutput_,
               "When forces are reduced on the CPU, there needs to be force output");
    pme_gpu_reduce_outputs(flags, output, wcycle, forceWithVirial, enerd);

    return true;
}

// This is used by PME-only ranks
PmeOutput pme_gpu_wait_finish_task(gmx_pme_t* pme, const int flags, gmx_wallcycle* wcycle)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    wallcycle_start(wcycle, ewcWAIT_GPU_PME_GATHER);

    // Synchronize the whole PME stream at once, including D2H result transfers
    // if there are outputs we need to wait for at this step; we still call getOutputs
    // for uniformity and because it sets the PmeOutput.haveForceOutput_.
    const bool haveComputedEnergyAndVirial = (flags & GMX_PME_CALC_ENER_VIR) != 0;
    if (!pme->gpu->settings.useGpuForceReduction || haveComputedEnergyAndVirial)
    {
        pme_gpu_synchronize(pme->gpu);
    }

    PmeOutput output = pme_gpu_getOutput(*pme, flags);
    wallcycle_stop(wcycle, ewcWAIT_GPU_PME_GATHER);
    return output;
}

// This is used when not using the alternate-waiting reduction
void pme_gpu_wait_and_reduce(gmx_pme_t*            pme,
                             const int             flags,
                             gmx_wallcycle*        wcycle,
                             gmx::ForceWithVirial* forceWithVirial,
                             gmx_enerdata_t*       enerd)
{
    PmeOutput output = pme_gpu_wait_finish_task(pme, flags, wcycle);
    GMX_ASSERT(pme->gpu->settings.useGpuForceReduction == !output.haveForceOutput_,
               "When forces are reduced on the CPU, there needs to be force output");
    pme_gpu_reduce_outputs(flags, output, wcycle, forceWithVirial, enerd);
}

void pme_gpu_reinit_computation(const gmx_pme_t* pme, gmx_wallcycle* wcycle)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU);
    wallcycle_sub_start_nocount(wcycle, ewcsLAUNCH_GPU_PME);

    pme_gpu_update_timings(pme->gpu);

    pme_gpu_clear_grids(pme->gpu);
    pme_gpu_clear_energy_virial(pme->gpu);

    wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME);
    wallcycle_stop(wcycle, ewcLAUNCH_GPU);
}

DeviceBuffer<float> pme_gpu_get_device_x(const gmx_pme_t* pme)
{
    GMX_ASSERT((pme && pme_gpu_active(pme)),
               "PME GPU coordinates buffer was requested from uninitialized PME module");
    return pme_gpu_get_kernelparam_coordinates(pme->gpu);
}

void* pme_gpu_get_device_f(const gmx_pme_t* pme)
{
    if (!pme || !pme_gpu_active(pme))
    {
        return nullptr;
    }
    return pme_gpu_get_kernelparam_forces(pme->gpu);
}

void pme_gpu_set_device_x(const gmx_pme_t* pme, DeviceBuffer<float> d_x)
{
    GMX_ASSERT(pme != nullptr, "Null pointer is passed as a PME to the set coordinates function.");
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    pme_gpu_set_kernelparam_coordinates(pme->gpu, d_x);
}

void* pme_gpu_get_device_stream(const gmx_pme_t* pme)
{
    if (!pme || !pme_gpu_active(pme))
    {
        return nullptr;
    }
    return pme_gpu_get_stream(pme->gpu);
}

void* pme_gpu_get_device_context(const gmx_pme_t* pme)
{
    if (!pme || !pme_gpu_active(pme))
    {
        return nullptr;
    }
    return pme_gpu_get_context(pme->gpu);
}

GpuEventSynchronizer* pme_gpu_get_f_ready_synchronizer(const gmx_pme_t* pme)
{
    if (!pme || !pme_gpu_active(pme))
    {
        return nullptr;
    }

    return pme_gpu_get_forces_ready_synchronizer(pme->gpu);
}
