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

/*! \internal \file
 * \brief This file contains internal CUDA function implementations
 * for performing the PME calculations on GPU.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include "config.h"

/* GPU initialization includes */
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/logger.h"

/* The rest */
#include "pme.h"

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/pmalloc_cuda.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"
#include "pme-3dfft.cuh"

void pme_gpu_synchronize(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaStreamSynchronize(pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "Failed to synchronize the PME GPU stream!");
}

void pme_gpu_alloc_energy_virial(const pme_gpu_t *pmeGPU)
{
    const size_t energyAndVirialSize = PME_GPU_VIRIAL_AND_ENERGY_COUNT * sizeof(float);
    cudaError_t  stat                = cudaMalloc((void **)&pmeGPU->kernelParams.constants.virialAndEnergy, energyAndVirialSize);
    CU_RET_ERR(stat, "cudaMalloc failed on PME energy and virial");
    pmalloc((void **)&pmeGPU->io.h_virialAndEnergy, energyAndVirialSize);
}

void pme_gpu_free_energy_virial(pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaFree(pmeGPU->kernelParams.constants.virialAndEnergy);
    CU_RET_ERR(stat, "cudaFree failed on PME energy and virial");
    pmeGPU->kernelParams.constants.virialAndEnergy = NULL;
    pfree(pmeGPU->io.h_virialAndEnergy);
    pmeGPU->io.h_virialAndEnergy = NULL;
}

void pme_gpu_clear_energy_virial(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams.constants.virialAndEnergy, 0,
                                       PME_GPU_VIRIAL_AND_ENERGY_COUNT * sizeof(float), pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME energies/virial cudaMemsetAsync error");
}

void pme_gpu_realloc_and_copy_bspline_values(const pme_gpu_t *pmeGPU)
{
    const int splineValuesOffset[DIM] = {
        0,
        pmeGPU->kernelParams.grid.localGridSize[XX],
        pmeGPU->kernelParams.grid.localGridSize[XX] + pmeGPU->kernelParams.grid.localGridSize[YY]
    };
    memcpy((void *)&pmeGPU->kernelParams.grid.splineValuesOffset, &splineValuesOffset, sizeof(splineValuesOffset));

    const int newSplineValuesSize = pmeGPU->kernelParams.grid.localGridSize[XX] +
        pmeGPU->kernelParams.grid.localGridSize[YY] +
        pmeGPU->kernelParams.grid.localGridSize[ZZ];
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.grid.splineValuesArray, NULL, sizeof(float),
                        &pmeGPU->archSpecific->splineValuesSize, &pmeGPU->archSpecific->splineValuesSizeAlloc, newSplineValuesSize, pmeGPU->archSpecific->pmeStream, true);

    for (int i = 0; i < DIM; i++)
    {
        /* Reallocate the host buffer */
        const size_t modSize = pmeGPU->kernelParams.grid.localGridSize[i] * sizeof(float);
        if ((pmeGPU->io.h_splineValues[i] == NULL) || (pmeGPU->io.splineValuesSizes[i] < modSize))
        {
            pfree(pmeGPU->io.h_splineValues[i]);
            pmalloc((void **)&pmeGPU->io.h_splineValues[i], modSize);
        }
        memcpy((void *)pmeGPU->io.h_splineValues[i], pmeGPU->common->bsp_mod[i].data(), modSize);
        /* TODO: use pinning here as well! */
        /* FIXME: no need for separate buffers */
        cu_copy_H2D_async(pmeGPU->kernelParams.grid.splineValuesArray + splineValuesOffset[i], pmeGPU->io.h_splineValues[i], modSize, pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_free_bspline_values(const pme_gpu_t *pmeGPU)
{
    for (int i = 0; i < DIM; i++)
    {
        pfree(pmeGPU->io.h_splineValues[i]);
    }
    cu_free_buffered(pmeGPU->kernelParams.grid.splineValuesArray, &pmeGPU->archSpecific->splineValuesSize, &pmeGPU->archSpecific->splineValuesSizeAlloc);
}

void pme_gpu_realloc_forces(const pme_gpu_t *pmeGPU)
{
    const size_t newForcesSize = pmeGPU->nAtomsAlloc * DIM;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.atoms.forces, NULL, sizeof(float),
                        &pmeGPU->archSpecific->forcesSize, &pmeGPU->archSpecific->forcesSizeAlloc, newForcesSize, pmeGPU->archSpecific->pmeStream, true);
}

void pme_gpu_free_forces(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams.atoms.forces, &pmeGPU->archSpecific->forcesSize, &pmeGPU->archSpecific->forcesSizeAlloc);
}

void pme_gpu_copy_input_forces(const pme_gpu_t *pmeGPU)
{
    GMX_ASSERT(pmeGPU->io.h_forces, "NULL host forces pointer in PME GPU");
    const size_t forcesSize = DIM * pmeGPU->kernelParams.atoms.nAtoms * sizeof(float);
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    cu_copy_H2D_async(pmeGPU->kernelParams.atoms.forces, pmeGPU->io.h_forces, forcesSize, pmeGPU->archSpecific->pmeStream);
}

void pme_gpu_sync_output_forces(const pme_gpu_t *pmeGPU)
{
    cudaStream_t s    = pmeGPU->archSpecific->pmeStream;
    cudaError_t  stat = cudaStreamWaitEvent(s, pmeGPU->archSpecific->syncForcesD2H, 0);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU forces");

    for (int i = 0; i < DIM * pmeGPU->kernelParams.atoms.nAtoms; i++)
    {
        GMX_ASSERT(!isnan(pmeGPU->io.h_forces[i]), "PME GPU - wrong forces produced.");
    }
}

void pme_gpu_realloc_coordinates(const pme_gpu_t *pmeGPU)
{
    const size_t newCoordinatesSize = pmeGPU->nAtomsAlloc * DIM;
    GMX_ASSERT(newCoordinatesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.atoms.coordinates, NULL, sizeof(float),
                        &pmeGPU->archSpecific->coordinatesSize, &pmeGPU->archSpecific->coordinatesSizeAlloc, newCoordinatesSize, pmeGPU->archSpecific->pmeStream, true);
#if PME_GPU_USE_PADDING
    const size_t paddingIndex = DIM * pmeGPU->kernelParams.atoms.nAtoms;
    const size_t paddingCount = DIM * pmeGPU->nAtomsAlloc - paddingIndex;
    if (paddingCount > 0)
    {
        cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams.atoms.coordinates + paddingIndex, 0, paddingCount * sizeof(float), pmeGPU->archSpecific->pmeStream);
        CU_RET_ERR(stat, "PME failed to clear the padded coordinates");
    }
#endif
}

void pme_gpu_copy_coordinates(const pme_gpu_t *pmeGPU)
{
    GMX_ASSERT(pmeGPU->io.h_coordinates, "Bad host-side coordinate buffer in PME GPU");
    cu_copy_H2D_async(pmeGPU->kernelParams.atoms.coordinates, pmeGPU->io.h_coordinates,
                      pmeGPU->kernelParams.atoms.nAtoms * DIM * sizeof(float), pmeGPU->archSpecific->pmeStream);
}

void pme_gpu_free_coordinates(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams.atoms.coordinates, &pmeGPU->archSpecific->coordinatesSize, &pmeGPU->archSpecific->coordinatesSizeAlloc);
}

void pme_gpu_realloc_and_copy_coefficients(const pme_gpu_t *pmeGPU)
{
    GMX_ASSERT(pmeGPU->io.h_coefficients, "Bad host-side charge buffer in PME GPU");
    const size_t newCoefficientsSize = pmeGPU->nAtomsAlloc;
    GMX_ASSERT(newCoefficientsSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.atoms.coefficients, NULL, sizeof(float),
                        &pmeGPU->archSpecific->coefficientsSize, &pmeGPU->archSpecific->coefficientsSizeAlloc, newCoefficientsSize, pmeGPU->archSpecific->pmeStream, true);
    cu_copy_H2D_async(pmeGPU->kernelParams.atoms.coefficients, pmeGPU->io.h_coefficients, pmeGPU->kernelParams.atoms.nAtoms * sizeof(float), pmeGPU->archSpecific->pmeStream);
#if PME_GPU_USE_PADDING
    const size_t paddingIndex = pmeGPU->kernelParams.atoms.nAtoms;
    const size_t paddingCount = pmeGPU->nAtomsAlloc - paddingIndex;
    if (paddingCount > 0)
    {
        cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams.atoms.coefficients + paddingIndex, 0, paddingCount * sizeof(float), pmeGPU->archSpecific->pmeStream);
        CU_RET_ERR(stat, "PME failed to clear the padded charges");
    }
#endif
}

void pme_gpu_free_coefficients(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams.atoms.coefficients, &pmeGPU->archSpecific->coefficientsSize, &pmeGPU->archSpecific->coefficientsSizeAlloc);
}

void pme_gpu_realloc_spline_data(const pme_gpu_t *pmeGPU)
{
    const int    order             = pmeGPU->common->pme_order;
    const int    alignment         = PME_SPREADGATHER_PARTICLES_PER_WARP;
    const size_t nAtomsPadded      = ((pmeGPU->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t newSplineDataSize = DIM * order * nAtomsPadded;
    GMX_ASSERT(newSplineDataSize > 0, "Bad number of atoms in PME GPU");
    /* Two arrays of the same size */
    int currentSizeTemp      = pmeGPU->archSpecific->splineDataSize;
    int currentSizeTempAlloc = pmeGPU->archSpecific->splineDataSizeAlloc;
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.atoms.theta, NULL, sizeof(float),
                        &currentSizeTemp, &currentSizeTempAlloc, newSplineDataSize, pmeGPU->archSpecific->pmeStream, true);
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.atoms.dtheta, NULL, sizeof(float),
                        &pmeGPU->archSpecific->splineDataSize, &pmeGPU->archSpecific->splineDataSizeAlloc, newSplineDataSize, pmeGPU->archSpecific->pmeStream, true);
}

void pme_gpu_free_spline_data(const pme_gpu_t *pmeGPU)
{
    /* Two arrays of the same size */
    cu_free_buffered(pmeGPU->kernelParams.atoms.theta);
    cu_free_buffered(pmeGPU->kernelParams.atoms.dtheta, &pmeGPU->archSpecific->splineDataSize, &pmeGPU->archSpecific->splineDataSizeAlloc);
}

void pme_gpu_realloc_grid_indices(const pme_gpu_t *pmeGPU)
{
    const size_t newIndicesSize = DIM * pmeGPU->nAtomsAlloc;
    GMX_ASSERT(newIndicesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.atoms.gridlineIndices, NULL, sizeof(int),
                        &pmeGPU->archSpecific->gridlineIndicesSize, &pmeGPU->archSpecific->gridlineIndicesSizeAlloc, newIndicesSize, pmeGPU->archSpecific->pmeStream, true);
}

void pme_gpu_free_grid_indices(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams.atoms.gridlineIndices, &pmeGPU->archSpecific->gridlineIndicesSize, &pmeGPU->archSpecific->gridlineIndicesSizeAlloc);
}

void pme_gpu_realloc_grids(pme_gpu_t *pmeGPU)
{
    // TODO: make tests to be assured this grid size is always suffcieint for copying the CPU grids
    // TODO: put the gridsize in the structure maybe?
    const int newGridSize = pmeGPU->kernelParams.grid.localGridSizePadded[XX] *
        pmeGPU->kernelParams.grid.localGridSizePadded[YY] *
        pmeGPU->kernelParams.grid.localGridSizePadded[ZZ];

    if (pmeGPU->archSpecific->performOutOfPlaceFFT)
    {
        /* Allocate a separate complex grid */
        int tempGridSize      = pmeGPU->archSpecific->gridSize;
        int tempGridSizeAlloc = pmeGPU->archSpecific->gridSizeAlloc;
        cu_realloc_buffered((void **)&pmeGPU->kernelParams.grid.fourierGrid, NULL, sizeof(float),
                            &tempGridSize, &tempGridSizeAlloc, newGridSize, pmeGPU->archSpecific->pmeStream, true);
    }
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.grid.realGrid, NULL, sizeof(float),
                        &pmeGPU->archSpecific->gridSize, &pmeGPU->archSpecific->gridSizeAlloc, newGridSize, pmeGPU->archSpecific->pmeStream, true);
    if (!pmeGPU->archSpecific->performOutOfPlaceFFT)
    {
        /* Using the same grid */
        pmeGPU->kernelParams.grid.fourierGrid = pmeGPU->kernelParams.grid.realGrid;
    }
}

void pme_gpu_free_grids(const pme_gpu_t *pmeGPU)
{
    if (pmeGPU->archSpecific->performOutOfPlaceFFT)
    {
        /* Free a separate complex grid of the same size */
        cu_free_buffered(pmeGPU->kernelParams.grid.fourierGrid);
    }
    cu_free_buffered(pmeGPU->kernelParams.grid.realGrid, &pmeGPU->archSpecific->gridSize, &pmeGPU->archSpecific->gridSizeAlloc);
}

void pme_gpu_clear_grids(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams.grid.realGrid, 0,
                                       pmeGPU->archSpecific->gridSize * sizeof(float), pmeGPU->archSpecific->pmeStream);
    /* Should the complex grid be cleared in some weird case? */
    CU_RET_ERR(stat, "cudaMemsetAsync on the PME grid error");
}

void pme_gpu_realloc_and_copy_fract_shifts(pme_gpu_t *pmeGPU)
{
    cudaStream_t s = pmeGPU->archSpecific->pmeStream;

    const int    nx = pmeGPU->kernelParams.grid.localGridSize[XX];
    const int    ny = pmeGPU->kernelParams.grid.localGridSize[YY];
    const int    nz = pmeGPU->kernelParams.grid.localGridSize[ZZ];

    const int    cellCount = 5;
    /* This is the number of neighbor cells that is also hardcoded in make_gridindex5_to_localindex and should be the same
     * TODO: eliminate stray fives!
     */

    const int fshOffset[DIM] = {0, cellCount * nx, cellCount * (nx + ny)};
    memcpy(pmeGPU->kernelParams.grid.fshOffset, &fshOffset, sizeof(fshOffset));

    const int    newFractShiftsSize  = cellCount * (nx + ny + nz);

    /* Two arrays, same size */
    int currentSizeTemp      = pmeGPU->archSpecific->fractShiftsSize;
    int currentSizeTempAlloc = pmeGPU->archSpecific->fractShiftsSizeAlloc;
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.grid.fshArray, NULL, sizeof(float),
                        &currentSizeTemp, &currentSizeTempAlloc,
                        newFractShiftsSize, pmeGPU->archSpecific->pmeStream, true);
    float *fshArray = pmeGPU->kernelParams.grid.fshArray;
    cu_realloc_buffered((void **)&pmeGPU->kernelParams.grid.nnArray, NULL, sizeof(int),
                        &pmeGPU->archSpecific->fractShiftsSize, &pmeGPU->archSpecific->fractShiftsSizeAlloc,
                        newFractShiftsSize, pmeGPU->archSpecific->pmeStream, true);
    int *nnArray = pmeGPU->kernelParams.grid.nnArray;

    /* TODO: pinning */

    for (int i = 0; i < DIM; i++)
    {
        pmeGPU->kernelParams.grid.fshOffset[i] = fshOffset[i];
        cu_copy_H2D_async(fshArray + fshOffset[i], pmeGPU->common->fsh[i].data(), cellCount * pmeGPU->kernelParams.grid.localGridSize[i] * sizeof(float), s);
        cu_copy_H2D_async(nnArray + fshOffset[i], pmeGPU->common->nn[i].data(), cellCount * pmeGPU->kernelParams.grid.localGridSize[i] * sizeof(int), s);
    }

    /* TODO: fix the textures code */
}

void pme_gpu_free_fract_shifts(const pme_gpu_t *pmeGPU)
{
    /* Two arrays, same size */
    cu_free_buffered(pmeGPU->kernelParams.grid.fshArray);
    cu_free_buffered(pmeGPU->kernelParams.grid.nnArray, &pmeGPU->archSpecific->fractShiftsSize, &pmeGPU->archSpecific->fractShiftsSizeAlloc);
}

void pme_gpu_sync_energy_virial(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaStreamWaitEvent(pmeGPU->archSpecific->pmeStream, pmeGPU->archSpecific->syncEnerVirD2H, 0);
    CU_RET_ERR(stat, "Error while waiting for PME solve");

    for (int j = 0; j < PME_GPU_VIRIAL_AND_ENERGY_COUNT; j++)
    {
        GMX_ASSERT(!isnan(pmeGPU->io.h_virialAndEnergy[j]), "PME GPU produces incorrect energy/virial.");
    }
}

void pme_gpu_sync_grid(const pme_gpu_t *pmeGPU, const gmx_fft_direction dir)
{
    /* FIXME: this function does not actually seem to be used when it should be, with CPU FFT? */
    bool syncGPUGrid = ((dir == GMX_FFT_REAL_TO_COMPLEX) ? true : pme_gpu_performs_solve(pmeGPU));
    if (syncGPUGrid)
    {
        cudaEvent_t syncEvent = (dir == GMX_FFT_REAL_TO_COMPLEX) ? pmeGPU->archSpecific->syncSpreadGridD2H : pmeGPU->archSpecific->syncSolveGridD2H;
        cudaError_t stat      = cudaStreamWaitEvent(pmeGPU->archSpecific->pmeStream, syncEvent, 0);
        CU_RET_ERR(stat, "Error while waiting for the PME GPU grid to be copied to CPU");
    }
}

void pme_gpu_init_specific(pme_gpu_t *pmeGPU, const gmx_hw_info_t *hwinfo, const gmx_gpu_opt_t *gpu_opt)
{
    /* FIXME: fix the GPU ID selection as well as initialization */
    int       gpuIndex = 0;
    char      gpu_err_str[STRLEN];
    GMX_ASSERT(hwinfo, "No hardware information");
    GMX_ASSERT(hwinfo->gpu_info.gpu_dev, "No GPU information");
    GMX_ASSERT(gpu_opt->dev_use, "No GPU information");
    /* Use GPU #0 for now since the code for GPU init has to be reworked anyway.
     * And don't forget to resurrect the external GMX_PME_GPU_ID env. variable.
     */
    pmeGPU->deviceInfo = &hwinfo->gpu_info.gpu_dev[gpu_opt->dev_use[gpuIndex]];
    const gmx::MDLogger temp;
    if (!init_gpu(temp, gpuIndex, gpu_err_str, &hwinfo->gpu_info, gpu_opt))
    {
        gmx_fatal(FARGS, "Could not select GPU %d for PME rank %d\n", pmeGPU->deviceInfo->id, gpuIndex);
    }

    /* Allocate the GPU-specific structure itself */
    pmeGPU->archSpecific = std::shared_ptr<pme_gpu_specific_t>(new pme_gpu_specific_t());

    pmeGPU->archSpecific->performOutOfPlaceFFT = true;
    /* This should give better performance, according to the cuFFT documentation.
     * The performance seems to be the same though.
     * Perhaps the limiting factor is using paddings/overlaps in the grid, which is also frowned upon.
     * PME could also try to pick up nice grid sizes (with factors of 2, 3, 5, 7).
     */

    pmeGPU->archSpecific->useTiming = (getenv("GMX_DISABLE_CUDA_TIMING") == NULL);
    /* This should also check for NB GPU being launched, and NB should check for PME GPU! */

    //pmeGPU->archSpecific->bUseTextureObjects = (pmeGPU->deviceInfo->prop.major >= 3);
    /* TODO: have to fix the GPU id selection */

    /* Creating a PME CUDA stream */
    cudaError_t stat;
#if GMX_CUDA_VERSION >= 5050
    int         highest_priority;
    int         lowest_priority;
    stat = cudaDeviceGetStreamPriorityRange(&lowest_priority, &highest_priority);
    CU_RET_ERR(stat, "PME cudaDeviceGetStreamPriorityRange failed");
    stat = cudaStreamCreateWithPriority(&pmeGPU->archSpecific->pmeStream,
                                        cudaStreamDefault, //cudaStreamNonBlocking,
                                        highest_priority);

    CU_RET_ERR(stat, "cudaStreamCreateWithPriority on the PME stream failed");
#else
    stat = cudaStreamCreate(&pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME cudaStreamCreate error");
#endif
}

void pme_gpu_destroy_specific(const pme_gpu_t *pmeGPU)
{
    /* Destroy the CUDA stream */
    cudaError_t stat = cudaStreamDestroy(pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME cudaStreamDestroy error");
}

void pme_gpu_init_sync_events(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat;
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncEnerVirD2H, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on syncEnerVirH2D failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncForcesD2H, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on syncForcesH2D failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncSpreadGridD2H, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on syncSpreadGridH2D failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncSolveGridD2H, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on syncSolveGridH2D failed");
}

void pme_gpu_destroy_sync_events(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat;
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncEnerVirD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncEnerVirH2D");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncForcesD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncForcesH2D");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncSpreadGridD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncpreadGridH2D");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncSolveGridD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncSolveGridH2D");
}

void pme_gpu_reinit_3dfft(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_performs_FFT(pmeGPU))
    {
        pmeGPU->archSpecific->pfft_setup_gpu.resize(0); // FIXME: reallocations
        for (int i = 0; i < pmeGPU->common->ngrids; i++)
        {
            pmeGPU->archSpecific->pfft_setup_gpu.push_back(std::unique_ptr<parallel_3dfft_gpu_t>(new parallel_3dfft_gpu_t(pmeGPU)));
        }
    }
}

void pme_gpu_destroy_3dfft(const pme_gpu_t *pmeGPU)
{
    pmeGPU->archSpecific->pfft_setup_gpu.resize(0);
}
