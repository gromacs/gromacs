/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

#include "gromacs/ewald/pme.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "pme-gpu-internal.h"
#include "pme-grid.h"
#include "pme-internal.h"
#include "pme-solve.h"

bool pme_gpu_task_enabled(const gmx_pme_t *pme)
{
    return pme && pme->useGPU;
}

void pme_gpu_reset_timings(const gmx_pme_t *pme)
{
    if (pme_gpu_active(pme))
    {
        pme_gpu_reset_timings(pme->gpu);
    }
}

void pme_gpu_get_timings(const gmx_pme_t *pme, gmx_wallclock_gpu_pme_t *timings)
{
    if (pme_gpu_active(pme))
    {
        pme_gpu_get_timings(pme->gpu, timings);
    }
}

/*! \brief \internal
 * A convenience wrapper for launching either the GPU or CPU FFT.
 *
 * \param[in] pme            The PME structure.
 * \param[in] grid_index     The grid index - should currently always be 0.
 * \param[in] dir            The FFT direction enum.
 * \param[in] wcycle         The wallclock counter.
 */
void parallel_3dfft_execute_gpu_wrapper(gmx_pme_t              *pme,
                                        const int               grid_index,
                                        enum gmx_fft_direction  dir,
                                        gmx_wallcycle_t         wcycle)
{
    GMX_ASSERT(grid_index == 0, "Only single grid supported");
    if (pme_gpu_performs_FFT(pme->gpu))
    {
        wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_PME_FFT);
        pme_gpu_3dfft(pme->gpu, dir, grid_index);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME_FFT);
    }
    else
    {
        wallcycle_start(wcycle, ewcPME_FFT);
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
        for (int thread = 0; thread < pme->nthread; thread++)
        {
            gmx_parallel_3dfft_execute(pme->pfft_setup[grid_index], dir, thread, wcycle);
        }
        wallcycle_stop(wcycle, ewcPME_FFT);
    }
}

/* Finally, the actual PME step code.
 * Together, they are a GPU counterpart to gmx_pme_do, albeit cut down due to unsupported features
 */
// TODO: add gmx_ prefix?

void pme_gpu_launch_everything_but_gather(gmx_pme_t            *pme,
                                          int                   nAtoms,
                                          const rvec           *x,
                                          const real           *charges,
                                          const matrix          box,
                                          gmx_wallcycle_t       wcycle,
                                          int                   flags)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    wallcycle_start(wcycle, ewcLAUNCH_GPU_PME);

    pme_gpu_t           *pmeGpu = pme->gpu;

    real                *grid        = nullptr;
    real                *fftgrid;
    t_complex           *cfftgrid;
    bool                 bFirst;
    const bool           bCalcEnerVir            = flags & GMX_PME_CALC_ENER_VIR;
    const bool           bBackFFT                = flags & (GMX_PME_CALC_F | GMX_PME_CALC_POT);

    GMX_ASSERT(pme->nnodes > 0, "");
    GMX_ASSERT(pme->nnodes == 1 || pme->ndecompdim > 0, "");

    bFirst = true;

    wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_PME_INIT);
    if (pme->gpu->settings.needToUpdateAtoms)
    {
        /* This only does a one-time atom data init at the first MD step.
         * Later, pme_gpu_reinit_atoms is called when needed after gmx_pme_recv_coeffs_coords.
         */
        pme_gpu_reinit_atoms(pmeGpu, nAtoms, charges);
        pme->gpu->settings.needToUpdateAtoms = false;
    }
    pme_gpu_start_step(pmeGpu, box, x);
    wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME_INIT);

    const unsigned int grid_index = 0;

    /* Unpack structure */
    const pmegrids_t *pmegrid     = &pme->pmegrid[grid_index];
    fftgrid     = pme->fftgrid[grid_index];
    cfftgrid    = pme->cfftgrid[grid_index];

    grid = pmegrid->grid.grid;

    // no bBackFFT, no bCalcF checks

    if (flags & GMX_PME_SPREAD)
    {
        /* Spread the coefficients on a grid */
        wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_PME_SPREAD);
        pme_gpu_spread(pmeGpu, grid_index, grid, bFirst, true);
        wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME_SPREAD);

        //if (!pme->bUseThreads)
        {
            if (!pme_gpu_performs_wrapping(pmeGpu))
            {
                wrap_periodic_pmegrid(pme, grid);
            }

            /* sum contributions to local grid from other nodes */
#if GMX_MPI
            if (pme->nnodes > 1)
            {
                gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_FORWARD);
                where();
            }
#endif
            if (!pme_gpu_performs_FFT(pmeGpu))
            {
                copy_pmegrid_to_fftgrid(pme, grid, fftgrid, grid_index);
            }
        }
    }

    try
    {
        if (flags & GMX_PME_SOLVE)
        {
            /* do 3d-fft */
            parallel_3dfft_execute_gpu_wrapper(pme, grid_index, GMX_FFT_REAL_TO_COMPLEX,
                                               wcycle);

            /* solve in k-space for our local cells */
            if (pme_gpu_performs_solve(pmeGpu))
            {
                wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_PME_SOLVE);
                pme_gpu_solve(pmeGpu, cfftgrid, bCalcEnerVir, pme_gpu_uses_dd(pmeGpu));
                wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME_SOLVE);
            }
            else
            {
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
                for (int thread = 0; thread < pme->nthread; thread++)
                {
                    solve_pme_yzx(pme, cfftgrid,
                                  box[XX][XX]*box[YY][YY]*box[ZZ][ZZ],
                                  bCalcEnerVir, pme->nthread, thread);
                }
            }
        }

        if (bBackFFT)
        {
            /* do 3d-invfft */
            parallel_3dfft_execute_gpu_wrapper(pme, grid_index, GMX_FFT_COMPLEX_TO_REAL, wcycle);

            if (!pme_gpu_performs_FFT(pmeGpu) || !pme_gpu_performs_gather(pmeGpu))
            {
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
                for (int thread = 0; thread < pme->nthread; thread++)
                {
                    copy_fftgrid_to_pmegrid(pme, fftgrid, grid, grid_index, pme->nthread, thread);
                }
            }
        }
    } GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if (bBackFFT)
    {
        /* distribute local grid to all nodes */
        if (!pme_gpu_performs_wrapping(pmeGpu))
        {
            unwrap_periodic_pmegrid(pme, grid);
        }
    }

    wallcycle_stop(wcycle, ewcLAUNCH_GPU_PME);
}

void pme_gpu_launch_gather(const gmx_pme_t                 *pme,
                           gmx_wallcycle_t gmx_unused       wcycle,
                           rvec                            *forces,
                           bool                             overwriteForces)
{
    if (!pme_gpu_performs_gather(pme->gpu))
    {
        return;
    }

    wallcycle_start_nocount(wcycle, ewcLAUNCH_GPU_PME);
    wallcycle_sub_start(wcycle, ewcsLAUNCH_GPU_PME_GATHER);
    pme_gpu_gather(pme->gpu, reinterpret_cast<float *>(forces), overwriteForces, nullptr);
    wallcycle_sub_stop(wcycle, ewcsLAUNCH_GPU_PME_GATHER);
    wallcycle_stop(wcycle, ewcLAUNCH_GPU_PME);
}

void pme_gpu_get_results(const gmx_pme_t *pme,
                         gmx_wallcycle_t  wcycle,
                         matrix           vir_q,
                         real            *energy_q,
                         int              flags)
{
    GMX_ASSERT(pme->useGPU, "pme_gpu_get_results should not be called on the CPU PME run.");

    const bool       bCalcEnerVir            = flags & GMX_PME_CALC_ENER_VIR;
    const bool       bCalcF                  = flags & GMX_PME_CALC_F;

    // FIXME: the wallcycle hierarchy is messy - we pause ewcFORCE just to exclude the PME GPU sync time
    wallcycle_stop(wcycle, ewcFORCE);

    wallcycle_start(wcycle, ewcWAIT_GPU_PME);
    wallcycle_sub_start(wcycle, ewcsWAIT_GPU_PME);
    pme_gpu_finish_step(pme->gpu, bCalcF, bCalcEnerVir);
    wallcycle_sub_stop(wcycle, ewcsWAIT_GPU_PME);
    wallcycle_stop(wcycle, ewcWAIT_GPU_PME);

    wallcycle_start_nocount(wcycle, ewcFORCE);

    if (bCalcEnerVir)
    {
        if (pme->doCoulomb)
        {
            pme_gpu_get_energy_virial(pme->gpu, energy_q, vir_q);
            if (debug)
            {
                fprintf(debug, "Electrostatic PME mesh energy [GPU]: %g\n", *energy_q);
            }
        }
        else
        {
            *energy_q = 0;
        }
    }
    /* No bCalcF code since currently forces are copied to the output host buffer with no transformation. */
}
