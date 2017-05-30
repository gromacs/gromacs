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
 * \brief This file contains internal CUDA function implementations
 * for performing the PME calculations on GPU.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cmath>

/* The rest */
#include "pme.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/pmalloc_cuda.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"
#include "pme-3dfft.cuh"
#include "pme-grid.h"

int pme_gpu_get_atom_data_alignment(const pme_gpu_t *pmeGPU)
{
    const int order = pmeGPU->common->pme_order;
    GMX_ASSERT(order > 0, "Invalid PME order");
    return PME_ATOM_DATA_ALIGNMENT;
}

int pme_gpu_get_atoms_per_warp(const pme_gpu_t *pmeGPU)
{
    const int order = pmeGPU->common->pme_order;
    GMX_ASSERT(order > 0, "Invalid PME order");
    return PME_SPREADGATHER_ATOMS_PER_WARP;
}

void pme_gpu_synchronize(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaStreamSynchronize(pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "Failed to synchronize the PME GPU stream!");
}

void pme_gpu_alloc_energy_virial(const pme_gpu_t *pmeGPU)
{
    const size_t energyAndVirialSize = c_virialAndEnergyCount * sizeof(float);
    cudaError_t  stat                = cudaMalloc((void **)&pmeGPU->kernelParams->constants.d_virialAndEnergy, energyAndVirialSize);
    CU_RET_ERR(stat, "cudaMalloc failed on PME energy and virial");
    pmalloc((void **)&pmeGPU->staging.h_virialAndEnergy, energyAndVirialSize);
}

void pme_gpu_free_energy_virial(pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaFree(pmeGPU->kernelParams->constants.d_virialAndEnergy);
    CU_RET_ERR(stat, "cudaFree failed on PME energy and virial");
    pmeGPU->kernelParams->constants.d_virialAndEnergy = nullptr;
    pfree(pmeGPU->staging.h_virialAndEnergy);
    pmeGPU->staging.h_virialAndEnergy = nullptr;
}

void pme_gpu_clear_energy_virial(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams->constants.d_virialAndEnergy, 0,
                                       c_virialAndEnergyCount * sizeof(float), pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME energy/virial cudaMemsetAsync error");
}

void pme_gpu_realloc_and_copy_bspline_values(const pme_gpu_t *pmeGPU)
{
    const int splineValuesOffset[DIM] = {
        0,
        pmeGPU->kernelParams->grid.realGridSize[XX],
        pmeGPU->kernelParams->grid.realGridSize[XX] + pmeGPU->kernelParams->grid.realGridSize[YY]
    };
    memcpy((void *)&pmeGPU->kernelParams->grid.splineValuesOffset, &splineValuesOffset, sizeof(splineValuesOffset));

    const int newSplineValuesSize = pmeGPU->kernelParams->grid.realGridSize[XX] +
        pmeGPU->kernelParams->grid.realGridSize[YY] +
        pmeGPU->kernelParams->grid.realGridSize[ZZ];
    const bool shouldRealloc = (newSplineValuesSize > pmeGPU->archSpecific->splineValuesSize);
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->grid.d_splineModuli, nullptr, sizeof(float),
                        &pmeGPU->archSpecific->splineValuesSize, &pmeGPU->archSpecific->splineValuesSizeAlloc, newSplineValuesSize, pmeGPU->archSpecific->pmeStream, true);
    if (shouldRealloc)
    {
        /* Reallocate the host buffer */
        pfree(pmeGPU->staging.h_splineModuli);
        pmalloc((void **)&pmeGPU->staging.h_splineModuli, newSplineValuesSize * sizeof(float));
    }
    for (int i = 0; i < DIM; i++)
    {
        memcpy(pmeGPU->staging.h_splineModuli + splineValuesOffset[i], pmeGPU->common->bsp_mod[i].data(), pmeGPU->common->bsp_mod[i].size() * sizeof(float));
    }
    /* TODO: pin original buffer instead! */
    cu_copy_H2D_async(pmeGPU->kernelParams->grid.d_splineModuli, pmeGPU->staging.h_splineModuli,
                      newSplineValuesSize * sizeof(float), pmeGPU->archSpecific->pmeStream);
}

void pme_gpu_free_bspline_values(const pme_gpu_t *pmeGPU)
{
    pfree(pmeGPU->staging.h_splineModuli);
    cu_free_buffered(pmeGPU->kernelParams->grid.d_splineModuli, &pmeGPU->archSpecific->splineValuesSize,
                     &pmeGPU->archSpecific->splineValuesSizeAlloc);
}

void pme_gpu_realloc_forces(const pme_gpu_t *pmeGPU)
{
    const size_t newForcesSize = pmeGPU->nAtomsAlloc * DIM;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->atoms.d_forces, nullptr, sizeof(float),
                        &pmeGPU->archSpecific->forcesSize, &pmeGPU->archSpecific->forcesSizeAlloc, newForcesSize, pmeGPU->archSpecific->pmeStream, true);
}

void pme_gpu_free_forces(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams->atoms.d_forces, &pmeGPU->archSpecific->forcesSize, &pmeGPU->archSpecific->forcesSizeAlloc);
}

void pme_gpu_copy_input_forces(const pme_gpu_t *pmeGPU, const float *h_forces)
{
    GMX_ASSERT(h_forces, "nullptr host forces pointer in PME GPU");
    const size_t forcesSize = DIM * pmeGPU->kernelParams->atoms.nAtoms * sizeof(float);
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    cu_copy_H2D_async(pmeGPU->kernelParams->atoms.d_forces, const_cast<float *>(h_forces), forcesSize, pmeGPU->archSpecific->pmeStream);
}

void pme_gpu_copy_output_forces(const pme_gpu_t *pmeGPU, float *h_forces)
{
    GMX_ASSERT(h_forces, "nullptr host forces pointer in PME GPU");
    const size_t forcesSize   = DIM * pmeGPU->kernelParams->atoms.nAtoms * sizeof(float);
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    cu_copy_D2H_async(h_forces, pmeGPU->kernelParams->atoms.d_forces, forcesSize, pmeGPU->archSpecific->pmeStream);
    cudaError_t stat = cudaEventRecord(pmeGPU->archSpecific->syncForcesD2H, pmeGPU->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME gather forces synchronization failure");
}

void pme_gpu_sync_output_forces(const pme_gpu_t *pmeGPU)
{
    cudaError_t  stat = cudaEventSynchronize(pmeGPU->archSpecific->syncForcesD2H);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU forces");
}

void pme_gpu_realloc_coordinates(const pme_gpu_t *pmeGPU)
{
    const size_t newCoordinatesSize = pmeGPU->nAtomsAlloc * DIM;
    GMX_ASSERT(newCoordinatesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->atoms.d_coordinates, nullptr, sizeof(float),
                        &pmeGPU->archSpecific->coordinatesSize, &pmeGPU->archSpecific->coordinatesSizeAlloc, newCoordinatesSize, pmeGPU->archSpecific->pmeStream, true);
    if (c_usePadding)
    {
        const size_t paddingIndex = DIM * pmeGPU->kernelParams->atoms.nAtoms;
        const size_t paddingCount = DIM * pmeGPU->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams->atoms.d_coordinates + paddingIndex, 0, paddingCount * sizeof(float), pmeGPU->archSpecific->pmeStream);
            CU_RET_ERR(stat, "PME failed to clear the padded coordinates");
        }
    }
}

void pme_gpu_copy_input_coordinates(const pme_gpu_t *pmeGPU, const rvec *h_coordinates)
{
    GMX_ASSERT(h_coordinates, "Bad host-side coordinate buffer in PME GPU");
#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision is supported");
    GMX_UNUSED_VALUE(h_coordinates);
#else
    cu_copy_H2D_async(pmeGPU->kernelParams->atoms.d_coordinates, const_cast<rvec *>(h_coordinates),
                      pmeGPU->kernelParams->atoms.nAtoms * sizeof(rvec), pmeGPU->archSpecific->pmeStream);
#endif
}

void pme_gpu_free_coordinates(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams->atoms.d_coordinates, &pmeGPU->archSpecific->coordinatesSize, &pmeGPU->archSpecific->coordinatesSizeAlloc);
}

void pme_gpu_realloc_and_copy_input_coefficients(const pme_gpu_t *pmeGPU, const float *h_coefficients)
{
    GMX_ASSERT(h_coefficients, "Bad host-side charge buffer in PME GPU");
    const size_t newCoefficientsSize = pmeGPU->nAtomsAlloc;
    GMX_ASSERT(newCoefficientsSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->atoms.d_coefficients, nullptr, sizeof(float),
                        &pmeGPU->archSpecific->coefficientsSize, &pmeGPU->archSpecific->coefficientsSizeAlloc,
                        newCoefficientsSize, pmeGPU->archSpecific->pmeStream, true);
    cu_copy_H2D_async(pmeGPU->kernelParams->atoms.d_coefficients, const_cast<float *>(h_coefficients),
                      pmeGPU->kernelParams->atoms.nAtoms * sizeof(float), pmeGPU->archSpecific->pmeStream);
    if (c_usePadding)
    {
        const size_t paddingIndex = pmeGPU->kernelParams->atoms.nAtoms;
        const size_t paddingCount = pmeGPU->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams->atoms.d_coefficients + paddingIndex, 0, paddingCount * sizeof(float), pmeGPU->archSpecific->pmeStream);
            CU_RET_ERR(stat, "PME failed to clear the padded charges");
        }
    }
}

void pme_gpu_free_coefficients(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams->atoms.d_coefficients, &pmeGPU->archSpecific->coefficientsSize, &pmeGPU->archSpecific->coefficientsSizeAlloc);
}

void pme_gpu_realloc_spline_data(const pme_gpu_t *pmeGPU)
{
    const int    order             = pmeGPU->common->pme_order;
    const int    alignment         = pme_gpu_get_atoms_per_warp(pmeGPU);
    const size_t nAtomsPadded      = ((pmeGPU->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const int    newSplineDataSize = DIM * order * nAtomsPadded;
    GMX_ASSERT(newSplineDataSize > 0, "Bad number of atoms in PME GPU");
    /* Two arrays of the same size */
    const bool shouldRealloc        = (newSplineDataSize > pmeGPU->archSpecific->splineDataSize);
    int        currentSizeTemp      = pmeGPU->archSpecific->splineDataSize;
    int        currentSizeTempAlloc = pmeGPU->archSpecific->splineDataSizeAlloc;
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->atoms.d_theta, nullptr, sizeof(float),
                        &currentSizeTemp, &currentSizeTempAlloc, newSplineDataSize, pmeGPU->archSpecific->pmeStream, true);
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->atoms.d_dtheta, nullptr, sizeof(float),
                        &pmeGPU->archSpecific->splineDataSize, &pmeGPU->archSpecific->splineDataSizeAlloc, newSplineDataSize, pmeGPU->archSpecific->pmeStream, true);
    // the host side reallocation
    if (shouldRealloc)
    {
        pfree(pmeGPU->staging.h_theta);
        pmalloc((void **)&pmeGPU->staging.h_theta, newSplineDataSize * sizeof(float));
        pfree(pmeGPU->staging.h_dtheta);
        pmalloc((void **)&pmeGPU->staging.h_dtheta, newSplineDataSize * sizeof(float));
    }
}

void pme_gpu_free_spline_data(const pme_gpu_t *pmeGPU)
{
    /* Two arrays of the same size */
    cu_free_buffered(pmeGPU->kernelParams->atoms.d_theta);
    cu_free_buffered(pmeGPU->kernelParams->atoms.d_dtheta, &pmeGPU->archSpecific->splineDataSize, &pmeGPU->archSpecific->splineDataSizeAlloc);
    pfree(pmeGPU->staging.h_theta);
    pfree(pmeGPU->staging.h_dtheta);
}

void pme_gpu_realloc_grid_indices(const pme_gpu_t *pmeGPU)
{
    const size_t newIndicesSize = DIM * pmeGPU->nAtomsAlloc;
    GMX_ASSERT(newIndicesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGPU->kernelParams->atoms.d_gridlineIndices, nullptr, sizeof(int),
                        &pmeGPU->archSpecific->gridlineIndicesSize, &pmeGPU->archSpecific->gridlineIndicesSizeAlloc, newIndicesSize, pmeGPU->archSpecific->pmeStream, true);
    pfree(pmeGPU->staging.h_gridlineIndices);
    pmalloc((void **)&pmeGPU->staging.h_gridlineIndices, newIndicesSize * sizeof(int));
}

void pme_gpu_free_grid_indices(const pme_gpu_t *pmeGPU)
{
    cu_free_buffered(pmeGPU->kernelParams->atoms.d_gridlineIndices, &pmeGPU->archSpecific->gridlineIndicesSize, &pmeGPU->archSpecific->gridlineIndicesSizeAlloc);
    pfree(pmeGPU->staging.h_gridlineIndices);
}

void pme_gpu_realloc_grids(pme_gpu_t *pmeGPU)
{
    auto     *kernelParamsPtr = pmeGPU->kernelParams.get();
    const int newRealGridSize = kernelParamsPtr->grid.realGridSizePadded[XX] *
        kernelParamsPtr->grid.realGridSizePadded[YY] *
        kernelParamsPtr->grid.realGridSizePadded[ZZ];
    const int newComplexGridSize = kernelParamsPtr->grid.complexGridSizePadded[XX] *
        kernelParamsPtr->grid.complexGridSizePadded[YY] *
        kernelParamsPtr->grid.complexGridSizePadded[ZZ] * 2;
    // Multiplied by 2 because we count complex grid size for complex numbers, but all allocations/pointers are float
    if (pmeGPU->archSpecific->performOutOfPlaceFFT)
    {
        /* 2 separate grids */
        cu_realloc_buffered((void **)&kernelParamsPtr->grid.d_fourierGrid, nullptr, sizeof(float),
                            &pmeGPU->archSpecific->complexGridSize, &pmeGPU->archSpecific->complexGridSizeAlloc,
                            newComplexGridSize, pmeGPU->archSpecific->pmeStream, true);
        cu_realloc_buffered((void **)&kernelParamsPtr->grid.d_realGrid, nullptr, sizeof(float),
                            &pmeGPU->archSpecific->realGridSize, &pmeGPU->archSpecific->realGridSizeAlloc,
                            newRealGridSize, pmeGPU->archSpecific->pmeStream, true);
    }
    else
    {
        /* A single buffer so that any grid will fit */
        const int newGridsSize = std::max(newRealGridSize, newComplexGridSize);
        cu_realloc_buffered((void **)&kernelParamsPtr->grid.d_realGrid, nullptr, sizeof(float),
                            &pmeGPU->archSpecific->realGridSize, &pmeGPU->archSpecific->realGridSizeAlloc,
                            newGridsSize, pmeGPU->archSpecific->pmeStream, true);
        kernelParamsPtr->grid.d_fourierGrid   = kernelParamsPtr->grid.d_realGrid;
        pmeGPU->archSpecific->complexGridSize = pmeGPU->archSpecific->realGridSize;
        // the size might get used later for copying the grid
    }
}

void pme_gpu_free_grids(const pme_gpu_t *pmeGPU)
{
    if (pmeGPU->archSpecific->performOutOfPlaceFFT)
    {
        cu_free_buffered(pmeGPU->kernelParams->grid.d_fourierGrid);
    }
    cu_free_buffered(pmeGPU->kernelParams->grid.d_realGrid,
                     &pmeGPU->archSpecific->realGridSize, &pmeGPU->archSpecific->realGridSizeAlloc);
}

void pme_gpu_clear_grids(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaMemsetAsync(pmeGPU->kernelParams->grid.d_realGrid, 0,
                                       pmeGPU->archSpecific->realGridSize * sizeof(float), pmeGPU->archSpecific->pmeStream);
    /* Should the complex grid be cleared in some weird case? */
    CU_RET_ERR(stat, "cudaMemsetAsync on the PME grid error");
}

void pme_gpu_realloc_and_copy_fract_shifts(pme_gpu_t *pmeGPU)
{
    pme_gpu_free_fract_shifts(pmeGPU);

    auto        *kernelParamsPtr = pmeGPU->kernelParams.get();

    const int    nx                  = kernelParamsPtr->grid.realGridSize[XX];
    const int    ny                  = kernelParamsPtr->grid.realGridSize[YY];
    const int    nz                  = kernelParamsPtr->grid.realGridSize[ZZ];
    const int    cellCount           = c_pmeNeighborUnitcellCount;
    const int    gridDataOffset[DIM] = {0, cellCount * nx, cellCount * (nx + ny)};

    memcpy(kernelParamsPtr->grid.tablesOffsets, &gridDataOffset, sizeof(gridDataOffset));

    const int    newFractShiftsSize  = cellCount * (nx + ny + nz);

    initParamLookupTable(kernelParamsPtr->grid.d_fractShiftsTable,
                         kernelParamsPtr->fractShiftsTableTexture,
                         &pme_gpu_get_fract_shifts_texref(),
                         pmeGPU->common->fsh.data(),
                         newFractShiftsSize,
                         pmeGPU->deviceInfo);

    initParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                         kernelParamsPtr->gridlineIndicesTableTexture,
                         &pme_gpu_get_gridline_texref(),
                         pmeGPU->common->nn.data(),
                         newFractShiftsSize,
                         pmeGPU->deviceInfo);
}

void pme_gpu_free_fract_shifts(const pme_gpu_t *pmeGPU)
{
    auto *kernelParamsPtr = pmeGPU->kernelParams.get();
    destroyParamLookupTable(kernelParamsPtr->grid.d_fractShiftsTable,
                            kernelParamsPtr->fractShiftsTableTexture,
                            &pme_gpu_get_fract_shifts_texref(),
                            pmeGPU->deviceInfo);
    destroyParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                            kernelParamsPtr->gridlineIndicesTableTexture,
                            &pme_gpu_get_gridline_texref(),
                            pmeGPU->deviceInfo);
}

void pme_gpu_sync_output_energy_virial(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaEventSynchronize(pmeGPU->archSpecific->syncEnerVirD2H);
    CU_RET_ERR(stat, "Error while waiting for PME solve output");

    for (int j = 0; j < c_virialAndEnergyCount; j++)
    {
        GMX_ASSERT(std::isfinite(pmeGPU->staging.h_virialAndEnergy[j]), "PME GPU produces incorrect energy/virial.");
    }
}

void pme_gpu_copy_input_gather_grid(const pme_gpu_t *pmeGpu, float *h_grid)
{
    const size_t gridSize = pmeGpu->archSpecific->realGridSize * sizeof(float);
    cu_copy_H2D_async(pmeGpu->kernelParams->grid.d_realGrid, h_grid, gridSize, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_output_spread_grid(const pme_gpu_t *pmeGpu, float *h_grid)
{
    const size_t gridSize = pmeGpu->archSpecific->realGridSize * sizeof(float);
    cu_copy_D2H_async(h_grid, pmeGpu->kernelParams->grid.d_realGrid, gridSize, pmeGpu->archSpecific->pmeStream);
    cudaError_t  stat = cudaEventRecord(pmeGpu->archSpecific->syncSpreadGridD2H, pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME spread grid sync event record failure");
}

void pme_gpu_copy_output_spread_atom_data(const pme_gpu_t *pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesSize     = DIM * nAtomsPadded * pmeGpu->common->pme_order * sizeof(float);
    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();
    cu_copy_D2H_async(pmeGpu->staging.h_dtheta, kernelParamsPtr->atoms.d_dtheta, splinesSize, pmeGpu->archSpecific->pmeStream);
    cu_copy_D2H_async(pmeGpu->staging.h_theta, kernelParamsPtr->atoms.d_theta, splinesSize, pmeGpu->archSpecific->pmeStream);
    cu_copy_D2H_async(pmeGpu->staging.h_gridlineIndices, kernelParamsPtr->atoms.d_gridlineIndices,
                      kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), pmeGpu->archSpecific->pmeStream);
    cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncSplineAtomDataD2H, pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME spread atom data sync event record failure");
}

void pme_gpu_copy_input_gather_atom_data(const pme_gpu_t *pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesSize     = DIM * nAtomsPadded * pmeGpu->common->pme_order * sizeof(float);
    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();
    if (c_usePadding)
    {
        const size_t gridlineIndicesSizePerAtom = DIM * sizeof(int);
        const size_t splineDataSizePerAtom      = pmeGpu->common->pme_order * DIM * sizeof(float);
        // TODO: could clear only the padding and not the whole thing, but this is a test-exclusive code anyway
        CU_RET_ERR(cudaMemsetAsync(kernelParamsPtr->atoms.d_gridlineIndices, 0, pmeGpu->nAtomsAlloc * gridlineIndicesSizePerAtom, pmeGpu->archSpecific->pmeStream),
                   "PME failed to clear the gridline indices");
        CU_RET_ERR(cudaMemsetAsync(kernelParamsPtr->atoms.d_dtheta, 0, pmeGpu->nAtomsAlloc * splineDataSizePerAtom, pmeGpu->archSpecific->pmeStream),
                   "PME failed to clear the spline derivatives");
        CU_RET_ERR(cudaMemsetAsync(kernelParamsPtr->atoms.d_theta, 0, pmeGpu->nAtomsAlloc * splineDataSizePerAtom, pmeGpu->archSpecific->pmeStream),
                   "PME failed to clear the spline values");
    }
    cu_copy_H2D_async(kernelParamsPtr->atoms.d_dtheta, pmeGpu->staging.h_dtheta, splinesSize, pmeGpu->archSpecific->pmeStream);
    cu_copy_H2D_async(kernelParamsPtr->atoms.d_theta, pmeGpu->staging.h_theta, splinesSize, pmeGpu->archSpecific->pmeStream);
    cu_copy_H2D_async(kernelParamsPtr->atoms.d_gridlineIndices, pmeGpu->staging.h_gridlineIndices,
                      kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_sync_spread_grid(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaEventSynchronize(pmeGPU->archSpecific->syncSpreadGridD2H);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU spread grid to be copied to the host");
}

void pme_gpu_sync_spline_atom_data(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaEventSynchronize(pmeGPU->archSpecific->syncSplineAtomDataD2H);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU atom data to be copied to the host");
}

void pme_gpu_sync_solve_grid(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat = cudaEventSynchronize(pmeGPU->archSpecific->syncSolveGridD2H);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU solve grid to be copied to the host");
    //should check for pme_gpu_performs_solve(pmeGPU)
}

void pme_gpu_init_internal(pme_gpu_t *pmeGPU)
{
    /* Allocate the target-specific structures */
    pmeGPU->archSpecific.reset(new pme_gpu_specific_t());
    pmeGPU->kernelParams.reset(new pme_gpu_kernel_params_t());

    pmeGPU->archSpecific->performOutOfPlaceFFT = true;
    /* This should give better performance, according to the cuFFT documentation.
     * The performance seems to be the same though.
     * TODO: PME could also try to pick up nice grid sizes (with factors of 2, 3, 5, 7).
     */

    pmeGPU->archSpecific->useTiming = (getenv("GMX_DISABLE_CUDA_TIMING") == nullptr) &&
        (getenv("GMX_DISABLE_GPU_TIMING") == nullptr);
    /* TODO: multiple CUDA streams on same GPU cause nonsense cudaEvent_t timings.
     * This should probably also check for gpuId exclusivity?
     */

    /* Creating a PME CUDA stream */
    cudaError_t stat;
    int         highest_priority, lowest_priority;
    stat = cudaDeviceGetStreamPriorityRange(&lowest_priority, &highest_priority);
    CU_RET_ERR(stat, "PME cudaDeviceGetStreamPriorityRange failed");
    stat = cudaStreamCreateWithPriority(&pmeGPU->archSpecific->pmeStream,
                                        cudaStreamDefault, //cudaStreamNonBlocking,
                                        highest_priority);
    CU_RET_ERR(stat, "cudaStreamCreateWithPriority on the PME stream failed");
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
    const auto  eventFlags = cudaEventDisableTiming;
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncEnerVirD2H, eventFlags);
    CU_RET_ERR(stat, "cudaEventCreate on syncEnerVirD2H failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncForcesD2H, eventFlags);
    CU_RET_ERR(stat, "cudaEventCreate on syncForcesD2H failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncSpreadGridD2H, eventFlags);
    CU_RET_ERR(stat, "cudaEventCreate on syncSpreadGridD2H failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncSplineAtomDataD2H, eventFlags);
    CU_RET_ERR(stat, "cudaEventCreate on syncSplineAtomDataD2H failed");
    stat = cudaEventCreateWithFlags(&pmeGPU->archSpecific->syncSolveGridD2H, eventFlags);
    CU_RET_ERR(stat, "cudaEventCreate on syncSolveGridD2H failed");
}

void pme_gpu_destroy_sync_events(const pme_gpu_t *pmeGPU)
{
    cudaError_t stat;
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncEnerVirD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncEnerVirD2H");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncForcesD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncForcesD2H");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncSpreadGridD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncSpreadGridD2H");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncSplineAtomDataD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncSplineAtomDataD2H");
    stat = cudaEventDestroy(pmeGPU->archSpecific->syncSolveGridD2H);
    CU_RET_ERR(stat, "cudaEventDestroy failed on syncSolveGridD2H");
}

void pme_gpu_reinit_3dfft(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_performs_FFT(pmeGPU))
    {
        pmeGPU->archSpecific->fftSetup.resize(0);
        for (int i = 0; i < pmeGPU->common->ngrids; i++)
        {
            pmeGPU->archSpecific->fftSetup.push_back(std::unique_ptr<GpuParallel3dFft>(new GpuParallel3dFft(pmeGPU)));
        }
    }
}

void pme_gpu_destroy_3dfft(const pme_gpu_t *pmeGPU)
{
    pmeGPU->archSpecific->fftSetup.resize(0);
}
