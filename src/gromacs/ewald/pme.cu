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
 *  \brief Implements PME GPU functions in CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

/* GPU initialization includes */
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/logger.h"

/* The rest */
#include "pme.h"

#include <assert.h>

#include "gromacs/gpu_utils/pmalloc_cuda.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"
#include "pme-3dfft.cuh"
#include "pme-grid.h"
#include "pme-solve.h"

/*! \brief \internal
 * Allocates the fixed size energy and virial buffer both on GPU and CPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_alloc_energy_virial(const gmx_pme_t *pme)
{
    const size_t energyAndVirialSize = PME_GPU_VIRIAL_AND_ENERGY_COUNT * sizeof(float);
    cudaError_t  stat                = cudaMalloc((void **)&pme->gpu->kernelParams.constants.virialAndEnergy, energyAndVirialSize);
    CU_RET_ERR(stat, "cudaMalloc failed on PME energy and virial");
    pmalloc((void **)&pme->gpu->virialAndEnergyHost, energyAndVirialSize);
}

/*! \brief \internal
 * Frees the energy and virial memory both on GPU and CPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_energy_virial(const gmx_pme_t *pme)
{
    cudaError_t stat = cudaFree(pme->gpu->kernelParams.constants.virialAndEnergy);
    CU_RET_ERR(stat, "cudaFree failed on PME energy and virial");
    pme->gpu->kernelParams.constants.virialAndEnergy = NULL;
    pfree(pme->gpu->virialAndEnergyHost);
    pme->gpu->virialAndEnergyHost = NULL;
}

/*! \brief \internal
 * Clears the energy and virial memory on GPU with 0. Should be called at the end of the energy/virial calculation step.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_clear_energy_virial(const gmx_pme_t *pme)
{
    cudaError_t stat = cudaMemsetAsync(pme->gpu->kernelParams.constants.virialAndEnergy, 0,
                                       PME_GPU_VIRIAL_AND_ENERGY_COUNT * sizeof(float), pme->gpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME energies/virial cudaMemsetAsync error");
}

/*! \brief \internal
 * Returns the output virial and energy of the PME solving.
 *
 * \param[in]  pme               The PME structure.
 * \param[out] energy            The output energy.
 * \param[out] virial            The output virial matrix.
 */
void pme_gpu_get_energy_virial(const gmx_pme_t *pme, real *energy, matrix virial)
{
    assert(energy);
    size_t j = 0;
    virial[XX][XX] = 0.25 * pme->gpu->virialAndEnergyHost[j++];
    virial[YY][YY] = 0.25 * pme->gpu->virialAndEnergyHost[j++];
    virial[ZZ][ZZ] = 0.25 * pme->gpu->virialAndEnergyHost[j++];
    virial[XX][YY] = virial[YY][XX] = 0.25 * pme->gpu->virialAndEnergyHost[j++];
    virial[XX][ZZ] = virial[ZZ][XX] = 0.25 * pme->gpu->virialAndEnergyHost[j++];
    virial[YY][ZZ] = virial[ZZ][YY] = 0.25 * pme->gpu->virialAndEnergyHost[j++];
    *energy        = 0.5 * pme->gpu->virialAndEnergyHost[j++];
}

/*! \brief \internal
 * Reallocates and copies the pre-computed B-spline values to the GPU.
 * FIXME: currently uses just a global memory, could be using texture memory/ldg.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_realloc_and_copy_bspline_values(const gmx_pme_t *pme)
{
    const int splineValuesOffset[DIM] = {0, pme->nkx, pme->nkx + pme->nky}; //?replace nkx
    memcpy(&pme->gpu->kernelParams.grid.splineValuesOffset, &splineValuesOffset, sizeof(splineValuesOffset));

    const int newSplineValuesSize = pme->nkx + pme->nky + pme->nkz;
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.grid.splineValuesArray, NULL, sizeof(float),
                        &pme->gpu->archSpecific->splineValuesSize, &pme->gpu->archSpecific->splineValuesSizeAlloc, newSplineValuesSize, pme->gpu->archSpecific->pmeStream, true);

    for (int i = 0; i < DIM; i++)
    {
        size_t       gridSize;
        switch (i)
        {
            case XX:
                gridSize = pme->nkx;
                break;

            case YY:
                gridSize = pme->nky;
                break;

            case ZZ:
                gridSize = pme->nkz;
                break;
        }
        size_t  modSize  = gridSize * sizeof(float);
        /* Reallocate the host buffer */
        if ((pme->gpu->splineValuesHost[i] == NULL) || (pme->gpu->splineValuesHostSizes[i] < modSize))
        {
            pfree(pme->gpu->splineValuesHost[i]);
            pmalloc((void **)&pme->gpu->splineValuesHost[i], modSize);
        }
        memcpy(pme->gpu->splineValuesHost[i], pme->bsp_mod[i], modSize);
        /* TODO: use pinning here as well! */
        cu_copy_H2D_async(pme->gpu->kernelParams.grid.splineValuesArray + splineValuesOffset[i], pme->gpu->splineValuesHost[i], modSize, pme->gpu->archSpecific->pmeStream);
    }
}

/*! \brief \internal
 * Frees the pre-computed B-spline values on the GPU (and the transfer CPU buffers).
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_bspline_values(const gmx_pme_t *pme)
{
    for (int i = 0; i < DIM; i++)
    {
        pfree(pme->gpu->splineValuesHost[i]);
    }
    cu_free_buffered(pme->gpu->kernelParams.grid.splineValuesArray, &pme->gpu->archSpecific->splineValuesSize, &pme->gpu->archSpecific->splineValuesSizeAlloc);
}

/*! \brief \internal
 * Copies the grid sizes for overlapping (used in the PME wrap/unwrap).
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_copy_wrap_zones(const gmx_pme_t *pme)
{
    const int nx      = pme->gpu->kernelParams.grid.localGridSize[XX];
    const int ny      = pme->gpu->kernelParams.grid.localGridSize[YY];
    const int nz      = pme->gpu->kernelParams.grid.localGridSize[ZZ];
    const int overlap = pme->pme_order - 1;

    /* Cell counts in the 7 overlapped grid parts */
    /* Is this correct? No Z alignment changes? */
    const int3 zoneSizes_h[PME_GPU_OVERLAP_ZONES_COUNT] =
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
    int2 zoneSizesYZ_h[PME_GPU_OVERLAP_ZONES_COUNT];
    for (int i = 0; i < PME_GPU_OVERLAP_ZONES_COUNT; i++)
    {
        zoneSizesYZ_h[i].x = zoneSizes_h[i].y;
        zoneSizesYZ_h[i].y = zoneSizes_h[i].z;
    }
    int cellsAccumCount_h[PME_GPU_OVERLAP_ZONES_COUNT];
    for (int i = 0; i < PME_GPU_OVERLAP_ZONES_COUNT; i++)
    {
        cellsAccumCount_h[i] = zoneSizes_h[i].x * zoneSizes_h[i].y * zoneSizes_h[i].z;
    }
    /* Accumulation */
    for (int i = 1; i < PME_GPU_OVERLAP_ZONES_COUNT; i++)
    {
        cellsAccumCount_h[i] = cellsAccumCount_h[i] + cellsAccumCount_h[i - 1];
    }
    memcpy(pme->gpu->kernelParams.grid.overlapSizes, zoneSizesYZ_h, sizeof(zoneSizesYZ_h));
    memcpy(pme->gpu->kernelParams.grid.overlapCellCounts, cellsAccumCount_h, sizeof(cellsAccumCount_h));
}

/*! \brief \internal
 * Reallocates the GPU buffer for the PME forces.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_realloc_forces(const gmx_pme_t *pme)
{
    const size_t newForcesSize = pme->gpu->archSpecific->nAtomsAlloc * DIM;
    assert(newForcesSize > 0);
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.atoms.forces, NULL, sizeof(float),
                        &pme->gpu->archSpecific->forcesSize, &pme->gpu->archSpecific->forcesSizeAlloc, newForcesSize, pme->gpu->archSpecific->pmeStream, true);
}

/*! \brief \internal
 * Frees the GPU buffer for the PME forces.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_forces(const gmx_pme_t *pme)
{
    cu_free_buffered(pme->gpu->kernelParams.atoms.forces, &pme->gpu->archSpecific->forcesSize, &pme->gpu->archSpecific->forcesSizeAlloc);
}

void pme_gpu_copy_input_forces(const gmx_pme_t *pme)
{
    GMX_ASSERT(pme->gpu->forcesHost, "NULL host forces pointer in PME GPU");
    const size_t forcesSize = DIM * pme->gpu->kernelParams.atoms.nAtoms * sizeof(float);
    assert(forcesSize > 0);
    cu_copy_H2D_async(pme->gpu->kernelParams.atoms.forces, pme->gpu->forcesHost, forcesSize, pme->gpu->archSpecific->pmeStream);
}

/*! \brief \internal
 * Waits for the PME GPU output forces copying to he CPU buffer (pme->gpu->forcesHost) to finish.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_sync_output_forces(const gmx_pme_t *pme)
{
    cudaStream_t s    = pme->gpu->archSpecific->pmeStream;
    cudaError_t  stat = cudaStreamWaitEvent(s, pme->gpu->archSpecific->syncForcesD2H, 0);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU forces");

    for (int i = 0; i < DIM * pme->gpu->kernelParams.atoms.nAtoms; i++)
    {
        GMX_ASSERT(!isnan(pme->gpu->forcesHost[i]), "PME GPU - wrong forces produced.");
    }
}

/*! \brief \internal
 * Reallocates the input coordinates buffer on the GPU (and clears the padded part if needed).
 *
 * \param[in] pme            The PME structure.
 *
 * Needs to be called on every DD step/in the beginning.
 */
void pme_gpu_realloc_coordinates(const gmx_pme_t *pme)
{
    const size_t newCoordinatesSize = pme->gpu->archSpecific->nAtomsAlloc * DIM;
    assert(newCoordinatesSize > 0);
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.atoms.coordinates, NULL, sizeof(float),
                        &pme->gpu->archSpecific->coordinatesSize, &pme->gpu->archSpecific->coordinatesSizeAlloc, newCoordinatesSize, pme->gpu->archSpecific->pmeStream, true);
#if PME_GPU_USE_PADDING
    const size_t paddingIndex = DIM * pme->gpu->kernelParams.atoms.nAtoms;
    const size_t paddingCount = DIM * pme->gpu->archSpecific->nAtomsAlloc - paddingIndex;
    if (paddingCount > 0)
    {
        cudaError_t stat = cudaMemsetAsync(pme->gpu->kernelParams.atoms.coordinates + paddingIndex, 0, paddingCount * sizeof(float), pme->gpu->archSpecific->pmeStream);
        CU_RET_ERR(stat, "PME failed to clear the padded coordinates");
    }
#endif
}

void pme_gpu_copy_coordinates(const gmx_pme_t *pme)
{
    assert(pme->gpu->coordinatesHost);
    cu_copy_H2D_async(pme->gpu->kernelParams.atoms.coordinates, pme->gpu->coordinatesHost, pme->gpu->kernelParams.atoms.nAtoms * DIM * sizeof(float), pme->gpu->archSpecific->pmeStream);
}

/*! \brief \internal
 * Frees the coordinates on the GPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_coordinates(const gmx_pme_t *pme)
{
    cu_free_buffered(pme->gpu->kernelParams.atoms.coordinates, &pme->gpu->archSpecific->coordinatesSize, &pme->gpu->archSpecific->coordinatesSizeAlloc);
}

/*! \brief \internal
 * Reallocates the buffer on the GPU and copies the charges/coefficients from the CPU buffer (pme->gpu->coefficientsHost). Clears the padded part if needed.
 *
 * \param[in] pme            The PME structure.
 *
 * Does not need to be done every MD step, only whenever the local charges change.
 * (So, in the beginning of the run, or on DD step).
 */
void pme_gpu_realloc_and_copy_coefficients(const gmx_pme_t *pme)
{
    assert(pme->gpu->coefficientsHost);
    const size_t newCoefficientsSize = pme->gpu->archSpecific->nAtomsAlloc;
    assert(newCoefficientsSize > 0);
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.atoms.coefficients, NULL, sizeof(float),
                        &pme->gpu->archSpecific->coefficientsSize, &pme->gpu->archSpecific->coefficientsSizeAlloc, newCoefficientsSize, pme->gpu->archSpecific->pmeStream, true);
    cu_copy_H2D_async(pme->gpu->kernelParams.atoms.coefficients, pme->gpu->coefficientsHost, pme->gpu->kernelParams.atoms.nAtoms * sizeof(float), pme->gpu->archSpecific->pmeStream);
#if PME_GPU_USE_PADDING
    const size_t paddingIndex = pme->gpu->kernelParams.atoms.nAtoms;
    const size_t paddingCount = pme->gpu->archSpecific->nAtomsAlloc - paddingIndex;
    if (paddingCount > 0)
    {
        cudaError_t stat = cudaMemsetAsync(pme->gpu->kernelParams.atoms.coefficients + paddingIndex, 0, paddingCount * sizeof(float), pme->gpu->archSpecific->pmeStream);
        CU_RET_ERR(stat, "PME failed to clear the padded charges");
    }
#endif
}

/*! \brief \internal
 * Frees the charges on the GPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_charges(const gmx_pme_t *pme)
{
    cu_free_buffered(pme->gpu->kernelParams.atoms.coefficients, &pme->gpu->archSpecific->coefficientsSize, &pme->gpu->archSpecific->coefficientsSizeAlloc);
}

/*! \brief \internal
 * Reallocates the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_realloc_spline_data(const gmx_pme_t *pme)
{
    const int    order             = pme->pme_order;
    const int    alignment         = PME_SPREADGATHER_PARTICLES_PER_WARP;
    const size_t nAtomsPadded      = ((pme->gpu->archSpecific->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t newSplineDataSize = DIM * order * nAtomsPadded;
    assert(newSplineDataSize > 0);
    /* Two arrays of the same size */
    int currentSizeTemp      = pme->gpu->archSpecific->splineDataSize;
    int currentSizeTempAlloc = pme->gpu->archSpecific->splineDataSizeAlloc;
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.atoms.theta, NULL, sizeof(float),
                        &currentSizeTemp, &currentSizeTempAlloc, newSplineDataSize, pme->gpu->archSpecific->pmeStream, true);
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.atoms.dtheta, NULL, sizeof(float),
                        &pme->gpu->archSpecific->splineDataSize, &pme->gpu->archSpecific->splineDataSizeAlloc, newSplineDataSize, pme->gpu->archSpecific->pmeStream, true);
}

/*! \brief \internal
 * Frees the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_spline_data(const gmx_pme_t *pme)
{
    /* Two arrays of the same size */
    cu_free_buffered(pme->gpu->kernelParams.atoms.theta);
    cu_free_buffered(pme->gpu->kernelParams.atoms.dtheta, &pme->gpu->archSpecific->splineDataSize, &pme->gpu->archSpecific->splineDataSizeAlloc);
}

/*! \brief \internal
 * Reallocates the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_realloc_grid_indices(const gmx_pme_t *pme)
{
    const size_t newIndicesSize = DIM * pme->gpu->archSpecific->nAtomsAlloc;
    assert(newIndicesSize > 0);
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.atoms.gridlineIndices, NULL, sizeof(int),
                        &pme->gpu->archSpecific->gridlineIndicesSize, &pme->gpu->archSpecific->gridlineIndicesSizeAlloc, newIndicesSize, pme->gpu->archSpecific->pmeStream, true);
}

/*! \brief \internal
 * Frees the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_grid_indices(const gmx_pme_t *pme)
{
    cu_free_buffered(pme->gpu->kernelParams.atoms.gridlineIndices, &pme->gpu->archSpecific->gridlineIndicesSize, &pme->gpu->archSpecific->gridlineIndicesSizeAlloc);
}

/*! \brief \internal
 * Reallocates the real space grid (and possibly the complex reciprocal grid) on the GPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_realloc_grids(const gmx_pme_t *pme)
{
    const int pnx         = pme->pmegrid_nx; //?
    const int pny         = pme->pmegrid_ny;
    const int pnz         = pme->pmegrid_nz;
    const int newGridSize = pnx * pny * pnz;

    if (pme->gpu->archSpecific->bOutOfPlaceFFT)
    {
        /* Allocate a separate complex grid */
        int tempGridSize      = pme->gpu->archSpecific->gridSize;
        int tempGridSizeAlloc = pme->gpu->archSpecific->gridSizeAlloc;
        cu_realloc_buffered((void **)&pme->gpu->kernelParams.grid.fourierGrid, NULL, sizeof(float),
                            &tempGridSize, &tempGridSizeAlloc, newGridSize, pme->gpu->archSpecific->pmeStream, true);
    }
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.grid.realGrid, NULL, sizeof(float),
                        &pme->gpu->archSpecific->gridSize, &pme->gpu->archSpecific->gridSizeAlloc, newGridSize, pme->gpu->archSpecific->pmeStream, true);
    if (!pme->gpu->archSpecific->bOutOfPlaceFFT)
    {
        /* Using the same grid */
        pme->gpu->kernelParams.grid.fourierGrid = pme->gpu->kernelParams.grid.realGrid;
    }
}

/*! \brief \internal
 * Frees the real space grid (and possibly the complex reciprocal grid) on the GPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_grids(const gmx_pme_t *pme)
{
    if (pme->gpu->archSpecific->bOutOfPlaceFFT)
    {
        /* Free a separate complex grid of the same size */
        cu_free_buffered(pme->gpu->kernelParams.grid.fourierGrid);
    }
    cu_free_buffered(pme->gpu->kernelParams.grid.realGrid, &pme->gpu->archSpecific->gridSize, &pme->gpu->archSpecific->gridSizeAlloc);
}

/*! \brief \internal
 * Clears the real space grid on the GPU. Should be called at the end of each MD step.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_clear_grids(const gmx_pme_t *pme)
{
    cudaError_t stat = cudaMemsetAsync(pme->gpu->kernelParams.grid.realGrid, 0, pme->gpu->archSpecific->gridSize * sizeof(float), pme->gpu->archSpecific->pmeStream);
    /* Should the complex grid be cleared in some weird case? */
    CU_RET_ERR(stat, "cudaMemsetAsync on the PME grid error");
}

/*! \brief \internal
 * Reallocates and copies the pre-computed fractional coordinates' shifts to the GPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_realloc_and_copy_fract_shifts(const gmx_pme_t *pme)
{
    cudaStream_t s = pme->gpu->archSpecific->pmeStream;

    const int    nx = pme->nkx; /* TODO: replace */
    const int    ny = pme->nky;
    const int    nz = pme->nkz;

    const int    cellCount = 5;
    /* This is the number of neighbor cells that is also hardcoded in make_gridindex5_to_localindex and should be the same */

    const int3 fshOffset = {0, cellCount * nx, cellCount * (nx + ny)};
    memcpy(&pme->gpu->kernelParams.grid.fshOffset, &fshOffset, sizeof(fshOffset));

    const int    newFractShiftsSize  = cellCount * (nx + ny + nz);

    /* Two arrays, same size */
    int currentSizeTemp      = pme->gpu->archSpecific->fractShiftsSize;
    int currentSizeTempAlloc = pme->gpu->archSpecific->fractShiftsSizeAlloc;
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.grid.fshArray, NULL, sizeof(float),
                        &currentSizeTemp, &currentSizeTempAlloc,
                        newFractShiftsSize, pme->gpu->archSpecific->pmeStream, true);
    float *fshArray = pme->gpu->kernelParams.grid.fshArray;
    cu_realloc_buffered((void **)&pme->gpu->kernelParams.grid.nnArray, NULL, sizeof(int),
                        &pme->gpu->archSpecific->fractShiftsSize, &pme->gpu->archSpecific->fractShiftsSizeAlloc,
                        newFractShiftsSize, pme->gpu->archSpecific->pmeStream, true);
    int *nnArray = pme->gpu->kernelParams.grid.nnArray;

    /* TODO: pinning */

    cu_copy_H2D_async(fshArray + fshOffset.x, pme->fshx, cellCount * nx * sizeof(float), s);
    cu_copy_H2D_async(fshArray + fshOffset.y, pme->fshy, cellCount * ny * sizeof(float), s);
    cu_copy_H2D_async(fshArray + fshOffset.z, pme->fshz, cellCount * nz * sizeof(float), s);

    cu_copy_H2D_async(nnArray + fshOffset.x, pme->nnx, cellCount * nx * sizeof(int), s);
    cu_copy_H2D_async(nnArray + fshOffset.y, pme->nny, cellCount * ny * sizeof(int), s);
    cu_copy_H2D_async(nnArray + fshOffset.z, pme->nnz, cellCount * nz * sizeof(int), s);

    /* TODO: fix the textures code */
}

/*! \brief \internal
 * Frees the pre-computed fractional coordinates' shifts on the GPU.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_free_fract_shifts(const gmx_pme_t *pme)
{
    /* Two arrays, same size */
    cu_free_buffered(pme->gpu->kernelParams.grid.fshArray);
    cu_free_buffered(pme->gpu->kernelParams.grid.nnArray, &pme->gpu->archSpecific->fractShiftsSize, &pme->gpu->archSpecific->fractShiftsSizeAlloc);
    /* TODO: unbind textures here! */
}

/*! \brief \internal
 * The PME GPU reinitialization function that is called both at the end of any MD step and on any load balancing step.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_reinit_step(const gmx_pme_t *pme)
{
    pme_gpu_clear_grids(pme);
    pme_gpu_clear_energy_virial(pme);
}

/*! \brief \internal
 * (Re-)initializes all the PME GPU data related to the grid size and cut-off.
 *
 * \param[in] pme            The PME structure.
 */
void pme_gpu_reinit_grids(const gmx_pme_t *pme)
{
    pme->gpu->kernelParams.grid.ewaldFactor = (M_PI * M_PI) / (pme->ewaldcoeff_q * pme->ewaldcoeff_q);

    /* The grid size variants */
    const int3   localGridSize = {pme->nkx, pme->nky, pme->nkz};
    memcpy(&pme->gpu->kernelParams.grid.localGridSize, &localGridSize, sizeof(localGridSize));
    const float3 localGridSizeFP = {(float)localGridSize.x, (float)localGridSize.y, (float)localGridSize.z};
    memcpy(&pme->gpu->kernelParams.grid.localGridSizeFP, &localGridSizeFP, sizeof(localGridSizeFP));
    const int3   localGridSizePadded = {pme->pmegrid_nx, pme->pmegrid_ny, pme->pmegrid_nz};
    memcpy(&pme->gpu->kernelParams.grid.localGridSizePadded, &localGridSizePadded, sizeof(localGridSizePadded));

    pme_gpu_copy_wrap_zones(pme);
    pme_gpu_realloc_and_copy_fract_shifts(pme);
    pme_gpu_realloc_and_copy_bspline_values(pme);
    pme_gpu_realloc_grids(pme);

    if (pme_gpu_performs_FFT(pme))
    {
        for (int i = 0; i < pme->ngrids; ++i)
        {
            pme_gpu_init_3dfft(&pme->gpu->archSpecific->pfft_setup_gpu[i], (int *)&localGridSize, pme);
        }
    }
}

/*! \brief \internal
 * Waits for the PME GPU output virial/energy copying to the intermediate CPU buffer to finish.
 *
 * \param[in] pme  The PME structure.
 */
void pme_gpu_sync_energy_virial(const gmx_pme_t *pme)
{
    cudaError_t stat = cudaStreamWaitEvent(pme->gpu->archSpecific->pmeStream, pme->gpu->archSpecific->syncEnerVirD2H, 0);
    CU_RET_ERR(stat, "Error while waiting for PME solve");

    for (int j = 0; j < PME_GPU_VIRIAL_AND_ENERGY_COUNT; j++)
    {
        GMX_ASSERT(!isnan(pme->gpu->virialAndEnergyHost[j]), "PME GPU produces incorrect energy/virial.");
    }
}

/*! \brief \internal
 * Waits for the PME GPU grid copying to the host-side buffer to finish.
 *
 * \param[in] pme  The PME structure.
 */
void pme_gpu_sync_grid(const gmx_pme_t *pme, const gmx_fft_direction dir)
{
    /* FIXME: this function does not actually seem to be used when it should be, with CPU FFT? */
    if (!pme_gpu_enabled(pme))
    {
        return;
    }

    gmx_bool syncGPUGrid = ((dir == GMX_FFT_REAL_TO_COMPLEX) ? true : pme_gpu_performs_solve(pme));
    if (syncGPUGrid)
    {
        cudaEvent_t syncEvent = (dir == GMX_FFT_REAL_TO_COMPLEX) ? pme->gpu->archSpecific->syncSpreadGridD2H : pme->gpu->archSpecific->syncSolveGridD2H;
        cudaError_t stat      = cudaStreamWaitEvent(pme->gpu->archSpecific->pmeStream, syncEvent, 0);
        CU_RET_ERR(stat, "Error while waiting for the PME GPU grid to be copied to CPU");
    }
}

/*! \brief \internal
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

    return error.empty();
}

/* The external PME GPU functions follow below */

void pme_gpu_reinit(gmx_pme_t *pme, const gmx_hw_info_t *hwinfo, const gmx_gpu_opt_t *gpu_opt)
{
    if (!pme_gpu_enabled(pme))
    {
        return;
    }

    const gmx_bool firstInit = !pme->gpu;
    if (firstInit)
    {
        std::string error;
        pme->bGPU = pme_gpu_check_restrictions(pme, error);
        if (!pme->bGPU)
        {
            gmx_fatal(FARGS, error.c_str());
        }

        snew(pme->gpu, 1);
        snew(pme->gpu->archSpecific, 1);

        cudaError_t stat;

        /* FIXME: fix the GPU ID selection as well as initialization */
        const int PMEGPURank = pme->nodeid;
        /* This is a node id within PME MPI communication group.
         * It doesn't really make much sense as a gpu index, right? */
        char      gpu_err_str[STRLEN];
        assert(hwinfo);
        assert(hwinfo->gpu_info.gpu_dev);
        assert(gpu_opt->dev_use);
        char *forcedGPUIdString = getenv("GMX_PME_GPU_ID");
        if (forcedGPUIdString)
        {
            int forcedGPUId = atoi(forcedGPUIdString); /* This is a CUDA GPU id */
            if (debug)
            {
                fprintf(debug, "PME GPU rank %d trying to use GPU %d\n", PMEGPURank, forcedGPUId);
            }
            stat = cudaSetDevice(forcedGPUId);
            CU_RET_ERR(stat, "PME failed to set the GPU device");
        }
        else
        {
            pme->gpu->deviceInfo = &hwinfo->gpu_info.gpu_dev[gpu_opt->dev_use[PMEGPURank]];
            const gmx::MDLogger temp;
            if (!init_gpu(temp, PMEGPURank, gpu_err_str, &hwinfo->gpu_info, gpu_opt))
            {
                gmx_fatal(FARGS, "Could not select GPU %d for PME rank %d\n", pme->gpu->deviceInfo->id, PMEGPURank);
            }
        }

        /* Some permanent settings are set here */

        pme->gpu->bGPUSingle = pme_gpu_enabled(pme) && (pme->nnodes == 1);
        /* A convenience variable. */

        pme->gpu->bGPUFFT = !pme_gpu_uses_dd(pme) && !getenv("GMX_PME_GPU_FFTW");
        /* cuFFT will only used for a single rank. */

        pme->gpu->bGPUSolve = TRUE;
        /* pme->gpu->archSpecific->bGPUFFT - CPU solve with the CPU FFTW is definitely broken at the moment - 20160511 */

        pme->gpu->bGPUGather = TRUE;
        /* CPU gather has got to be broken as well due to different theta/dtheta layout. */

        pme->gpu->bNeedToUpdateAtoms = TRUE;
        /* For the delayed atom data init */

        pme->gpu->archSpecific->bOutOfPlaceFFT = TRUE;
        /* This should give better performance, according to the cuFFT documentation.
         * The performance seems to be the same though.
         * Perhaps the limiting factor is using paddings/overlaps in the grid, which is also frowned upon.
         * PME could also try to pick up nice grid sizes (with factors of 2, 3, 5, 7)
         */

        pme->gpu->archSpecific->bTiming = (getenv("GMX_DISABLE_CUDA_TIMING") == NULL); /* This should also check for NB GPU being launched, and NB should check for PME GPU! */

        //pme->gpu->archSpecific->bUseTextureObjects = (pme->gpu->deviceInfo->prop.major >= 3);
        /* TODO: have to fix the GPU id selection, forced GPUIdHack?*/

        /* Creating a PME CUDA stream */
#if GMX_CUDA_VERSION >= 5050
        int highest_priority;
        int lowest_priority;
        stat = cudaDeviceGetStreamPriorityRange(&lowest_priority, &highest_priority);
        CU_RET_ERR(stat, "PME cudaDeviceGetStreamPriorityRange failed");
        stat = cudaStreamCreateWithPriority(&pme->gpu->archSpecific->pmeStream,
                                            cudaStreamDefault, //cudaStreamNonBlocking,
                                            highest_priority);

        CU_RET_ERR(stat, "cudaStreamCreateWithPriority on the PME stream failed");
#else
        stat = cudaStreamCreate(&pme->gpu->archSpecific->pmeStream);
        CU_RET_ERR(stat, "PME cudaStreamCreate error");
#endif

        /* Creating synchronization events */
        stat = cudaEventCreateWithFlags(&pme->gpu->archSpecific->syncEnerVirD2H, cudaEventDisableTiming);
        CU_RET_ERR(stat, "cudaEventCreate on syncEnerVirH2D failed");
        stat = cudaEventCreateWithFlags(&pme->gpu->archSpecific->syncForcesD2H, cudaEventDisableTiming);
        CU_RET_ERR(stat, "cudaEventCreate on syncForcesH2D failed");
        stat = cudaEventCreateWithFlags(&pme->gpu->archSpecific->syncSpreadGridD2H, cudaEventDisableTiming);
        CU_RET_ERR(stat, "cudaEventCreate on syncSpreadGridH2D failed");
        stat = cudaEventCreateWithFlags(&pme->gpu->archSpecific->syncSolveGridD2H, cudaEventDisableTiming);
        CU_RET_ERR(stat, "cudaEventCreate on syncSolveGridH2D failed");

        pme_gpu_init_timings(pme);

        pme_gpu_alloc_energy_virial(pme);

        assert(pme->epsilon_r != 0.0f);
        pme->gpu->kernelParams.constants.elFactor = ONE_4PI_EPS0 / pme->epsilon_r;

        snew(pme->gpu->archSpecific->pfft_setup_gpu, pme->ngrids);
    }

    pme_gpu_reinit_grids(pme);
    pme_gpu_reinit_step(pme);
}

void pme_gpu_destroy(gmx_pme_t *pme)
{
    if (!pme_gpu_enabled(pme))
    {
        return;
    }

    stopGpuProfiler();

    cudaError_t stat;

    /* Free lots of dynamic data */
    pme_gpu_free_energy_virial(pme);
    pme_gpu_free_bspline_values(pme);
    pme_gpu_free_forces(pme);
    pme_gpu_free_coordinates(pme);
    pme_gpu_free_charges(pme);
    pme_gpu_free_spline_data(pme);
    pme_gpu_free_grid_indices(pme);
    pme_gpu_free_fract_shifts(pme);
    pme_gpu_free_grids(pme);

    /* cuFFT cleanup */
    if (pme->gpu->archSpecific->pfft_setup_gpu)
    {
        for (int i = 0; i < pme->ngrids; i++)
        {
            pme_gpu_destroy_3dfft(pme->gpu->archSpecific->pfft_setup_gpu[i]);
        }
        sfree(pme->gpu->archSpecific->pfft_setup_gpu);
    }

    /* Free the synchronization events */
    stat = cudaEventDestroy(pme->gpu->archSpecific->syncEnerVirD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncEnerVirH2D");
    stat = cudaEventDestroy(pme->gpu->archSpecific->syncForcesD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncForcesH2D");
    stat = cudaEventDestroy(pme->gpu->archSpecific->syncSpreadGridD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncpreadGridH2D");
    stat = cudaEventDestroy(pme->gpu->archSpecific->syncSolveGridD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncSolveGridH2D");

    /* Free the timing events */
    pme_gpu_destroy_timings(pme);

    /* Destroy the CUDA stream */
    stat = cudaStreamDestroy(pme->gpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME cudaStreamDestroy error");

    /* Finally free the GPU structure itself */
    sfree(pme->gpu->archSpecific);
    sfree(pme->gpu);
    pme->gpu = NULL;
}

void pme_gpu_reinit_atoms(const gmx_pme_t *pme, const int nAtoms, float *coefficients)
{
    if (!pme_gpu_enabled(pme))
    {
        return;
    }

    pme->gpu->kernelParams.atoms.nAtoms = nAtoms;
    const int      alignment = 8; // FIXME: this is particlesPerBlock
    pme->gpu->archSpecific->nAtomsPadded = ((nAtoms + alignment - 1) / alignment) * alignment;
    int            nAtomsAlloc   = PME_GPU_USE_PADDING ? pme->gpu->archSpecific->nAtomsPadded : nAtoms;
    const gmx_bool haveToRealloc = (pme->gpu->archSpecific->nAtomsAlloc < nAtomsAlloc); /* This check might be redundant, but is logical */
    pme->gpu->archSpecific->nAtomsAlloc = nAtomsAlloc;

    pme->gpu->coefficientsHost = reinterpret_cast<float *>(coefficients);
    pme_gpu_realloc_and_copy_coefficients(pme); /* Could also be checked for haveToRealloc, but the copy always needs to be performed */

    if (haveToRealloc)
    {
        pme_gpu_realloc_coordinates(pme);
        pme_gpu_realloc_forces(pme);
        pme_gpu_realloc_spline_data(pme);
        pme_gpu_realloc_grid_indices(pme);
    }
}

void pme_gpu_finish_step(const gmx_pme_t *pme, const gmx_bool bCalcF, const gmx_bool bCalcEnerVir)
{
    if (!pme_gpu_enabled(pme))
    {
        return;
    }

    /* Needed for copy back as well as timing events */
    cudaError_t stat = cudaStreamSynchronize(pme->gpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "Failed to synchronize the PME GPU stream!");

    if (bCalcF)
    {
        pme_gpu_sync_output_forces(pme);
    }
    if (bCalcEnerVir)
    {
        pme_gpu_sync_energy_virial(pme);
    }
    pme_gpu_update_timings(pme);
    pme_gpu_reinit_step(pme);
}
