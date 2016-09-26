/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief This file contains internal function implementations
 * for performing the PME calculations on GPU.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "pme-internal.h"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <string>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "pme-gpu.h"

void pme_gpu_set_io_ranges(pme_gpu_t *pmeGPU, rvec *coordinates, rvec *forces)
{
    pmeGPU->io.h_forces      = reinterpret_cast<float *>(forces);
    pmeGPU->io.h_coordinates = reinterpret_cast<float *>(coordinates);
}

void pme_gpu_get_energy_virial(const pme_gpu_t *pmeGPU, real *energy, matrix virial)
{
    assert(energy);
    size_t j = 0;
    virial[XX][XX] = 0.25 * pmeGPU->io.h_virialAndEnergy[j++];
    virial[YY][YY] = 0.25 * pmeGPU->io.h_virialAndEnergy[j++];
    virial[ZZ][ZZ] = 0.25 * pmeGPU->io.h_virialAndEnergy[j++];
    virial[XX][YY] = virial[YY][XX] = 0.25 * pmeGPU->io.h_virialAndEnergy[j++];
    virial[XX][ZZ] = virial[ZZ][XX] = 0.25 * pmeGPU->io.h_virialAndEnergy[j++];
    virial[YY][ZZ] = virial[ZZ][YY] = 0.25 * pmeGPU->io.h_virialAndEnergy[j++];
    *energy        = 0.5 * pmeGPU->io.h_virialAndEnergy[j++];
}

void pme_gpu_start_step(pme_gpu_t *pmeGPU, const matrix box)
{
    pme_gpu_copy_coordinates(pmeGPU);

    const size_t   boxMemorySize        = sizeof(matrix);
    const gmx_bool haveToUpdateUnitCell = memcmp(pmeGPU->previousBox, box, boxMemorySize);
    /* There could be a pressure coupling check here, but this is more straightforward.
     * This is an exact comparison of float values though.
     */
    if (haveToUpdateUnitCell)
    {
        memcpy(pmeGPU->previousBox, box, boxMemorySize);

        pmeGPU->kernelParams.step.boxVolume = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
        assert(pmeGPU->kernelParams.step.boxVolume != 0.0f);

#if GMX_DOUBLE
        assert("PME is single-precision only on GPU. You shouldn't be seeing this message!");
#else
        matrix recipBox;
        gmx::invertBoxMatrix(box, recipBox);
        /* The GPU recipBox is transposed as compared to the CPU recipBox.
         * Spread uses matrix columns (while solve and gather use rows).
         * There is no particular reason for this; it might be further rethought/optimized for better access patterns.
         */
        const real newRecipBox[DIM][DIM] =
        {
            {recipBox[XX][XX], recipBox[YY][XX], recipBox[ZZ][XX]},
            {             0.0, recipBox[YY][YY], recipBox[ZZ][YY]},
            {             0.0,              0.0, recipBox[ZZ][ZZ]}
        };
        memcpy(pmeGPU->kernelParams.step.recipBox, newRecipBox, boxMemorySize);
#endif
    }
}

/*! \brief \libinternal
 * The PME GPU reinitialization function that is called both at the end of any MD step and on any load balancing step.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
void pme_gpu_reinit_step(const pme_gpu_t *pmeGPU)
{
    pme_gpu_clear_grids(pmeGPU);
    pme_gpu_clear_energy_virial(pmeGPU);
}

void pme_gpu_finish_step(const pme_gpu_t *pmeGPU, const gmx_bool bCalcF, const gmx_bool bCalcEnerVir)
{
    /* Needed for copy back as well as timing events */
    pme_gpu_synchronize(pmeGPU);

    if (bCalcF)
    {
        pme_gpu_sync_output_forces(pmeGPU);
    }
    if (bCalcEnerVir)
    {
        pme_gpu_sync_energy_virial(pmeGPU);
    }
    pme_gpu_update_timings(pmeGPU);
    pme_gpu_reinit_step(pmeGPU);
}

/*! \brief \libinternal
 * Copies the grid sizes for overlapping (used in the PME wrap/unwrap).
 *
 * \param[in] pmeGPU             The PME GPU structure.
 */
void pme_gpu_copy_wrap_zones(const pme_gpu_t *pmeGPU)
{
    const int nx      = pmeGPU->kernelParams.grid.localGridSize[XX];
    const int ny      = pmeGPU->kernelParams.grid.localGridSize[YY];
    const int nz      = pmeGPU->kernelParams.grid.localGridSize[ZZ];
    const int overlap = pmeGPU->common->pme_order - 1;

    /* Cell counts in the 7 overlapped grid parts */
    const int zoneSizes_h[PME_GPU_OVERLAP_ZONES_COUNT][DIM] =
    {
        {     nx,        ny,   overlap},
        {     nx,   overlap,        nz},
        {overlap,        ny,        nz},
        {     nx,   overlap,   overlap},
        {overlap,        ny,   overlap},
        {overlap,   overlap,        nz},
        {overlap,   overlap,   overlap}
    };
    /* The X is never used on the GPU, actually */
    int zoneSizesYZ_h[PME_GPU_OVERLAP_ZONES_COUNT * 2];
    for (size_t i = 0; i < PME_GPU_OVERLAP_ZONES_COUNT; i++)
    {
        zoneSizesYZ_h[2 * i    ] = zoneSizes_h[i][YY];
        zoneSizesYZ_h[2 * i + 1] = zoneSizes_h[i][ZZ];
    }
    int cellsAccumCount_h[PME_GPU_OVERLAP_ZONES_COUNT];
    for (int i = 0; i < PME_GPU_OVERLAP_ZONES_COUNT; i++)
    {
        cellsAccumCount_h[i] = zoneSizes_h[i][XX] * zoneSizes_h[i][YY] * zoneSizes_h[i][ZZ];
    }
    /* Accumulation */
    for (int i = 1; i < PME_GPU_OVERLAP_ZONES_COUNT; i++)
    {
        cellsAccumCount_h[i] = cellsAccumCount_h[i] + cellsAccumCount_h[i - 1];
    }
    memcpy((void *)pmeGPU->kernelParams.grid.overlapSizes, zoneSizesYZ_h, sizeof(zoneSizesYZ_h));
    memcpy((void *)pmeGPU->kernelParams.grid.overlapCellCounts, cellsAccumCount_h, sizeof(cellsAccumCount_h));
}

/*! \brief \libinternal
 * (Re-)initializes all the PME GPU data related to the grid size and cut-off.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
void pme_gpu_reinit_grids(pme_gpu_t *pmeGPU)
{
    pmeGPU->kernelParams.grid.ewaldFactor = (M_PI * M_PI) / (pmeGPU->common->ewaldcoeff_q * pmeGPU->common->ewaldcoeff_q);

    /* The grid size variants */
    for (int i = 0; i < DIM; i++)
    {
        pmeGPU->kernelParams.grid.localGridSize[i]       = pmeGPU->common->nk[i];
        pmeGPU->kernelParams.grid.localGridSizeFP[i]     = (float)pmeGPU->kernelParams.grid.localGridSize[i];
        pmeGPU->kernelParams.grid.localGridSizePadded[i] = pmeGPU->common->pmegrid_n[i];
    }

    pme_gpu_copy_wrap_zones(pmeGPU);
    pme_gpu_realloc_and_copy_fract_shifts(pmeGPU);
    pme_gpu_realloc_and_copy_bspline_values(pmeGPU);
    pme_gpu_realloc_grids(pmeGPU);
    pme_gpu_reinit_3dfft(pmeGPU); // FIXME: this is a memory leak - allocating new plans each time!
}

/* Several GPU functions that refer to the CPU PME data live here.
 * We would like to keep these away from the GPU-framework specific code for clarity,
 * as well as compilation issues with MPI.
 */

/*! \brief \libinternal
 * Copies everything useful from the PME CPU to the PME GPU structure.
 * The goal is to minimize interaction with the PME CPU structure in the GPU code.
 *
 * \param[in] pme         The PME structure.
 */
void pme_gpu_fetch_shared_data(const gmx_pme_t *pme)
{
    /* TODO: Consider refactoring the CPU PME code to use the same structure,
     * so that this function becomes 2 lines */
    pme_gpu_t *pmeGPU             = pme->gpu;
    pmeGPU->common->ngrids        = pme->ngrids;
    pmeGPU->common->epsilon_r     = pme->epsilon_r;
    pmeGPU->common->ewaldcoeff_q  = pme->ewaldcoeff_q;
    pmeGPU->common->nk[XX]        = pme->nkx;
    pmeGPU->common->nk[YY]        = pme->nky;
    pmeGPU->common->nk[ZZ]        = pme->nkz;
    pmeGPU->common->pmegrid_n[XX] = pme->pmegrid_nx;
    pmeGPU->common->pmegrid_n[YY] = pme->pmegrid_ny;
    pmeGPU->common->pmegrid_n[ZZ] = pme->pmegrid_nz;
    pmeGPU->common->pme_order     = pme->pme_order;
    for (int i = 0; i < DIM; i++)
    {
        pmeGPU->common->bsp_mod[i].assign(pme->bsp_mod[i], pme->bsp_mod[i] + pmeGPU->common->nk[i]);
    }
    const int magic = 5; /* Cell count (with neighboring cells) */
    pmeGPU->common->fsh[XX].assign(pme->fshx, pme->fshx + magic * pme->nkx);
    pmeGPU->common->fsh[YY].assign(pme->fshy, pme->fshy + magic * pme->nky);
    pmeGPU->common->fsh[ZZ].assign(pme->fshz, pme->fshz + magic * pme->nkz);
    pmeGPU->common->nn[XX].assign(pme->nnx, pme->nnx + magic * pme->nkx);
    pmeGPU->common->nn[YY].assign(pme->nny, pme->nny + magic * pme->nky);
    pmeGPU->common->nn[ZZ].assign(pme->nnz, pme->nnz + magic * pme->nkz);
}

/*! \brief \libinternal
 * Finds out if PME with given inputs is possible to run on GPU.
 *
 * \param[in]  pme          The PME structure.
 * \param[out] error        The error message if the input is not supported on GPU.
 * \returns                 TRUE if this PME input is possible to run on GPU, FALSE otherwise.
 */
gmx_bool pme_gpu_check_restrictions(const gmx_pme_t *pme,
                                    std::string     &error)
{
    if (pme->nnodes != 1)
    {
        error.append("PME is only implemented for a single rank on GPU. ");
    }
    if (pme->pme_order != 4)
    {
        error.append("PME is only implemented for the interpolation order of 4 on GPU. ");
    }
    if (pme->bFEP)
    {
        error.append("PME is only implemented for a single grid on GPU. ");
    }
    if (pme->doLJ)
    {
        error.append("PME LJ is not implemented on GPU. ");
    }
#if GMX_DOUBLE
    {
        error.append("PME is only implemented for single precision on GPU. ");
    }
#endif
#if GMX_GPU != GMX_GPU_CUDA
    {
        error.append("PME is only implemented for CUDA framework on GPU. ");
    }
#endif

    return error.empty();
}

void pme_gpu_reinit(gmx_pme_t *pme, const gmx_hw_info_t *hwinfo, const gmx_gpu_opt_t *gpu_opt)
{
    if (!gmx_pme_gpu_enabled(pme))
    {
        return;
    }

    pme_gpu_t     *pmeGPU    = pme->gpu;
    const gmx_bool firstInit = !pmeGPU;
    if (firstInit) /* One-time initialization */
    {
        std::string error;
        pme->bGPU = pme_gpu_check_restrictions(pme, error);
        if (!pme->bGPU)
        {
            gmx_fatal(FARGS, error.c_str());
        }

        snew(pmeGPU, 1);
        pme->gpu       = pmeGPU;
        pmeGPU->common = new pme_shared_t;

        /* Some permanent settings are set here */
        pmeGPU->settings.bGPUSingle = (pme->nnodes == 1);
        /* A convenience variable. */
        pmeGPU->settings.bGPUFFT = !pme_gpu_uses_dd(pmeGPU) && !getenv("GMX_PME_GPU_FFTW");
        /* GPU FFT will only used for a single rank. */
        pmeGPU->settings.bGPUSolve = TRUE;
        /* pmeGPU->settings.bGPUFFT - CPU solve with the CPU FFTW is definitely broken at the moment - 20160511 */
        pmeGPU->settings.bGPUGather = TRUE;
        /* CPU gather has got to be broken as well due to different theta/dtheta layout. */
        pmeGPU->settings.bNeedToUpdateAtoms = TRUE;
        /* For the delayed atom data init hack */

        pme_gpu_init_specific(pmeGPU, hwinfo, gpu_opt);
        pme_gpu_init_sync_events(pmeGPU);
        pme_gpu_init_timings(pmeGPU);
        pme_gpu_alloc_energy_virial(pmeGPU);
    }
    pme_gpu_fetch_shared_data(pme);
    /* After this call nothing in the GPU code should refer to the gmx_pme_t *pme - until the next pme_gpu_reinit */

    if (firstInit)
    {
        assert(pmeGPU->common->epsilon_r != 0.0f);
        pmeGPU->kernelParams.constants.elFactor = ONE_4PI_EPS0 / pmeGPU->common->epsilon_r;
        // assert(pmeGPU->common->ngrids == 1);
        // this assert will fail now because PME CPU is stupid and has 2 grids minimum
    }

    pme_gpu_reinit_grids(pmeGPU);
    pme_gpu_reinit_step(pmeGPU);
}

void pme_gpu_destroy(pme_gpu_t *pmeGPU)
{
    stopGpuProfiler();

    /* Free lots of data */
    pme_gpu_free_energy_virial(pmeGPU);
    pme_gpu_free_bspline_values(pmeGPU);
    pme_gpu_free_forces(pmeGPU);
    pme_gpu_free_coordinates(pmeGPU);
    pme_gpu_free_coefficients(pmeGPU);
    pme_gpu_free_spline_data(pmeGPU);
    pme_gpu_free_grid_indices(pmeGPU);
    pme_gpu_free_fract_shifts(pmeGPU);
    pme_gpu_free_grids(pmeGPU);

    pme_gpu_destroy_3dfft(pmeGPU); /* FIXME: Why are some things freed and some destroyed?
                                      No logic in naming */
    pme_gpu_destroy_sync_events(pmeGPU);
    pme_gpu_destroy_timings(pmeGPU);

    /* Free the GPU-framework specific data last */
    pme_gpu_destroy_specific(pmeGPU);

    delete pmeGPU->common;
    sfree(pmeGPU);
}

void pme_gpu_reinit_atoms(pme_gpu_t *pmeGPU, const int nAtoms, real *coefficients)
{
    pmeGPU->kernelParams.atoms.nAtoms = nAtoms;
    const int      alignment = 8; // FIXME: this is particlesPerBlock
    pmeGPU->nAtomsPadded = ((nAtoms + alignment - 1) / alignment) * alignment;
    int            nAtomsAlloc   = PME_GPU_USE_PADDING ? pmeGPU->nAtomsPadded : nAtoms;
    const gmx_bool haveToRealloc = (pmeGPU->nAtomsAlloc < nAtomsAlloc); /* This check might be redundant, but is logical */
    pmeGPU->nAtomsAlloc = nAtomsAlloc;

    pmeGPU->io.h_coefficients = reinterpret_cast<float *>(coefficients);
    pme_gpu_realloc_and_copy_coefficients(pmeGPU); /* Could also be checked for haveToRealloc, but the copy always needs to be performed */

    if (haveToRealloc)
    {
        pme_gpu_realloc_coordinates(pmeGPU);
        pme_gpu_realloc_forces(pmeGPU);
        pme_gpu_realloc_spline_data(pmeGPU);
        pme_gpu_realloc_grid_indices(pmeGPU);
    }
}
