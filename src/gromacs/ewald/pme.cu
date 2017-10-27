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

int pme_gpu_get_atom_data_alignment(const PmeGpu *pmeGpu)
{
    const int order = pmeGpu->common->pme_order;
    GMX_ASSERT(order > 0, "Invalid PME order");
    return PME_ATOM_DATA_ALIGNMENT;
}

int pme_gpu_get_atoms_per_warp(const PmeGpu *pmeGpu)
{
    const int order = pmeGpu->common->pme_order;
    GMX_ASSERT(order > 0, "Invalid PME order");
    return PME_SPREADGATHER_ATOMS_PER_WARP;
}

void pme_gpu_synchronize(const PmeGpu *pmeGpu)
{
    cudaError_t stat = cudaStreamSynchronize(pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "Failed to synchronize the PME GPU stream!");
}

void pme_gpu_alloc_energy_virial(const PmeGpu *pmeGpu)
{
    const size_t energyAndVirialSize = c_virialAndEnergyCount * sizeof(float);
    cudaError_t  stat                = cudaMalloc((void **)&pmeGpu->kernelParams->constants.d_virialAndEnergy, energyAndVirialSize);
    CU_RET_ERR(stat, "cudaMalloc failed on PME energy and virial");
    pmalloc((void **)&pmeGpu->staging.h_virialAndEnergy, energyAndVirialSize);
}

void pme_gpu_free_energy_virial(PmeGpu *pmeGpu)
{
    cudaError_t stat = cudaFree(pmeGpu->kernelParams->constants.d_virialAndEnergy);
    CU_RET_ERR(stat, "cudaFree failed on PME energy and virial");
    pmeGpu->kernelParams->constants.d_virialAndEnergy = nullptr;
    pfree(pmeGpu->staging.h_virialAndEnergy);
    pmeGpu->staging.h_virialAndEnergy = nullptr;
}

void pme_gpu_clear_energy_virial(const PmeGpu *pmeGpu)
{
    cudaError_t stat = cudaMemsetAsync(pmeGpu->kernelParams->constants.d_virialAndEnergy, 0,
                                       c_virialAndEnergyCount * sizeof(float), pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME energy/virial cudaMemsetAsync error");
}

void pme_gpu_realloc_and_copy_bspline_values(const PmeGpu *pmeGpu)
{
    const int splineValuesOffset[DIM] = {
        0,
        pmeGpu->kernelParams->grid.realGridSize[XX],
        pmeGpu->kernelParams->grid.realGridSize[XX] + pmeGpu->kernelParams->grid.realGridSize[YY]
    };
    memcpy((void *)&pmeGpu->kernelParams->grid.splineValuesOffset, &splineValuesOffset, sizeof(splineValuesOffset));

    const int newSplineValuesSize = pmeGpu->kernelParams->grid.realGridSize[XX] +
        pmeGpu->kernelParams->grid.realGridSize[YY] +
        pmeGpu->kernelParams->grid.realGridSize[ZZ];
    const bool shouldRealloc = (newSplineValuesSize > pmeGpu->archSpecific->splineValuesSize);
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->grid.d_splineModuli, nullptr, sizeof(float),
                        &pmeGpu->archSpecific->splineValuesSize, &pmeGpu->archSpecific->splineValuesSizeAlloc, newSplineValuesSize, pmeGpu->archSpecific->pmeStream, true);
    if (shouldRealloc)
    {
        /* Reallocate the host buffer */
        pfree(pmeGpu->staging.h_splineModuli);
        pmalloc((void **)&pmeGpu->staging.h_splineModuli, newSplineValuesSize * sizeof(float));
    }
    for (int i = 0; i < DIM; i++)
    {
        memcpy(pmeGpu->staging.h_splineModuli + splineValuesOffset[i], pmeGpu->common->bsp_mod[i].data(), pmeGpu->common->bsp_mod[i].size() * sizeof(float));
    }
    /* TODO: pin original buffer instead! */
    cu_copy_H2D(pmeGpu->kernelParams->grid.d_splineModuli, pmeGpu->staging.h_splineModuli,
                newSplineValuesSize * sizeof(float), pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_free_bspline_values(const PmeGpu *pmeGpu)
{
    pfree(pmeGpu->staging.h_splineModuli);
    cu_free_buffered(pmeGpu->kernelParams->grid.d_splineModuli, &pmeGpu->archSpecific->splineValuesSize,
                     &pmeGpu->archSpecific->splineValuesSizeAlloc);
}

void pme_gpu_realloc_forces(PmeGpu *pmeGpu)
{
    const size_t newForcesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->atoms.d_forces, nullptr, sizeof(float),
                        &pmeGpu->archSpecific->forcesSize, &pmeGpu->archSpecific->forcesSizeAlloc, newForcesSize, pmeGpu->archSpecific->pmeStream, true);
    pmeGpu->staging.h_forces.reserve(pmeGpu->nAtomsAlloc);
    pmeGpu->staging.h_forces.resize(pmeGpu->kernelParams->atoms.nAtoms);
}

void pme_gpu_free_forces(const PmeGpu *pmeGpu)
{
    cu_free_buffered(pmeGpu->kernelParams->atoms.d_forces, &pmeGpu->archSpecific->forcesSize, &pmeGpu->archSpecific->forcesSizeAlloc);
}

void pme_gpu_copy_input_forces(PmeGpu *pmeGpu)
{
    const size_t forcesSize = DIM * pmeGpu->kernelParams->atoms.nAtoms * sizeof(float);
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    cu_copy_H2D(pmeGpu->kernelParams->atoms.d_forces, pmeGpu->staging.h_forces.data(), forcesSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_output_forces(PmeGpu *pmeGpu)
{
    const size_t forcesSize   = DIM * pmeGpu->kernelParams->atoms.nAtoms * sizeof(float);
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    cu_copy_D2H(pmeGpu->staging.h_forces.data(), pmeGpu->kernelParams->atoms.d_forces, forcesSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_realloc_coordinates(const PmeGpu *pmeGpu)
{
    const size_t newCoordinatesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newCoordinatesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->atoms.d_coordinates, nullptr, sizeof(float),
                        &pmeGpu->archSpecific->coordinatesSize, &pmeGpu->archSpecific->coordinatesSizeAlloc, newCoordinatesSize, pmeGpu->archSpecific->pmeStream, true);
    if (c_usePadding)
    {
        const size_t paddingIndex = DIM * pmeGpu->kernelParams->atoms.nAtoms;
        const size_t paddingCount = DIM * pmeGpu->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            cudaError_t stat = cudaMemsetAsync(pmeGpu->kernelParams->atoms.d_coordinates + paddingIndex, 0, paddingCount * sizeof(float), pmeGpu->archSpecific->pmeStream);
            CU_RET_ERR(stat, "PME failed to clear the padded coordinates");
        }
    }
}

void pme_gpu_copy_input_coordinates(const PmeGpu *pmeGpu, const rvec *h_coordinates)
{
    GMX_ASSERT(h_coordinates, "Bad host-side coordinate buffer in PME GPU");
#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision is supported");
    GMX_UNUSED_VALUE(h_coordinates);
#else
    cu_copy_H2D(pmeGpu->kernelParams->atoms.d_coordinates, const_cast<rvec *>(h_coordinates),
                pmeGpu->kernelParams->atoms.nAtoms * sizeof(rvec), pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
#endif
}

void pme_gpu_free_coordinates(const PmeGpu *pmeGpu)
{
    cu_free_buffered(pmeGpu->kernelParams->atoms.d_coordinates, &pmeGpu->archSpecific->coordinatesSize, &pmeGpu->archSpecific->coordinatesSizeAlloc);
}

void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu *pmeGpu, const float *h_coefficients)
{
    GMX_ASSERT(h_coefficients, "Bad host-side charge buffer in PME GPU");
    const size_t newCoefficientsSize = pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newCoefficientsSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->atoms.d_coefficients, nullptr, sizeof(float),
                        &pmeGpu->archSpecific->coefficientsSize, &pmeGpu->archSpecific->coefficientsSizeAlloc,
                        newCoefficientsSize, pmeGpu->archSpecific->pmeStream, true);
    cu_copy_H2D(pmeGpu->kernelParams->atoms.d_coefficients, const_cast<float *>(h_coefficients),
                pmeGpu->kernelParams->atoms.nAtoms * sizeof(float), pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
    if (c_usePadding)
    {
        const size_t paddingIndex = pmeGpu->kernelParams->atoms.nAtoms;
        const size_t paddingCount = pmeGpu->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            cudaError_t stat = cudaMemsetAsync(pmeGpu->kernelParams->atoms.d_coefficients + paddingIndex, 0, paddingCount * sizeof(float), pmeGpu->archSpecific->pmeStream);
            CU_RET_ERR(stat, "PME failed to clear the padded charges");
        }
    }
}

void pme_gpu_free_coefficients(const PmeGpu *pmeGpu)
{
    cu_free_buffered(pmeGpu->kernelParams->atoms.d_coefficients, &pmeGpu->archSpecific->coefficientsSize, &pmeGpu->archSpecific->coefficientsSizeAlloc);
}

void pme_gpu_realloc_spline_data(const PmeGpu *pmeGpu)
{
    const int    order             = pmeGpu->common->pme_order;
    const int    alignment         = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded      = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const int    newSplineDataSize = DIM * order * nAtomsPadded;
    GMX_ASSERT(newSplineDataSize > 0, "Bad number of atoms in PME GPU");
    /* Two arrays of the same size */
    const bool shouldRealloc        = (newSplineDataSize > pmeGpu->archSpecific->splineDataSize);
    int        currentSizeTemp      = pmeGpu->archSpecific->splineDataSize;
    int        currentSizeTempAlloc = pmeGpu->archSpecific->splineDataSizeAlloc;
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->atoms.d_theta, nullptr, sizeof(float),
                        &currentSizeTemp, &currentSizeTempAlloc, newSplineDataSize, pmeGpu->archSpecific->pmeStream, true);
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->atoms.d_dtheta, nullptr, sizeof(float),
                        &pmeGpu->archSpecific->splineDataSize, &pmeGpu->archSpecific->splineDataSizeAlloc, newSplineDataSize, pmeGpu->archSpecific->pmeStream, true);
    // the host side reallocation
    if (shouldRealloc)
    {
        pfree(pmeGpu->staging.h_theta);
        pmalloc((void **)&pmeGpu->staging.h_theta, newSplineDataSize * sizeof(float));
        pfree(pmeGpu->staging.h_dtheta);
        pmalloc((void **)&pmeGpu->staging.h_dtheta, newSplineDataSize * sizeof(float));
    }
}

void pme_gpu_free_spline_data(const PmeGpu *pmeGpu)
{
    /* Two arrays of the same size */
    cu_free_buffered(pmeGpu->kernelParams->atoms.d_theta);
    cu_free_buffered(pmeGpu->kernelParams->atoms.d_dtheta, &pmeGpu->archSpecific->splineDataSize, &pmeGpu->archSpecific->splineDataSizeAlloc);
    pfree(pmeGpu->staging.h_theta);
    pfree(pmeGpu->staging.h_dtheta);
}

void pme_gpu_realloc_grid_indices(const PmeGpu *pmeGpu)
{
    const size_t newIndicesSize = DIM * pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newIndicesSize > 0, "Bad number of atoms in PME GPU");
    cu_realloc_buffered((void **)&pmeGpu->kernelParams->atoms.d_gridlineIndices, nullptr, sizeof(int),
                        &pmeGpu->archSpecific->gridlineIndicesSize, &pmeGpu->archSpecific->gridlineIndicesSizeAlloc, newIndicesSize, pmeGpu->archSpecific->pmeStream, true);
    pfree(pmeGpu->staging.h_gridlineIndices);
    pmalloc((void **)&pmeGpu->staging.h_gridlineIndices, newIndicesSize * sizeof(int));
}

void pme_gpu_free_grid_indices(const PmeGpu *pmeGpu)
{
    cu_free_buffered(pmeGpu->kernelParams->atoms.d_gridlineIndices, &pmeGpu->archSpecific->gridlineIndicesSize, &pmeGpu->archSpecific->gridlineIndicesSizeAlloc);
    pfree(pmeGpu->staging.h_gridlineIndices);
}

void pme_gpu_realloc_grids(PmeGpu *pmeGpu)
{
    auto     *kernelParamsPtr = pmeGpu->kernelParams.get();
    const int newRealGridSize = kernelParamsPtr->grid.realGridSizePadded[XX] *
        kernelParamsPtr->grid.realGridSizePadded[YY] *
        kernelParamsPtr->grid.realGridSizePadded[ZZ];
    const int newComplexGridSize = kernelParamsPtr->grid.complexGridSizePadded[XX] *
        kernelParamsPtr->grid.complexGridSizePadded[YY] *
        kernelParamsPtr->grid.complexGridSizePadded[ZZ] * 2;
    // Multiplied by 2 because we count complex grid size for complex numbers, but all allocations/pointers are float
    if (pmeGpu->archSpecific->performOutOfPlaceFFT)
    {
        /* 2 separate grids */
        cu_realloc_buffered((void **)&kernelParamsPtr->grid.d_fourierGrid, nullptr, sizeof(float),
                            &pmeGpu->archSpecific->complexGridSize, &pmeGpu->archSpecific->complexGridSizeAlloc,
                            newComplexGridSize, pmeGpu->archSpecific->pmeStream, true);
        cu_realloc_buffered((void **)&kernelParamsPtr->grid.d_realGrid, nullptr, sizeof(float),
                            &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc,
                            newRealGridSize, pmeGpu->archSpecific->pmeStream, true);
    }
    else
    {
        /* A single buffer so that any grid will fit */
        const int newGridsSize = std::max(newRealGridSize, newComplexGridSize);
        cu_realloc_buffered((void **)&kernelParamsPtr->grid.d_realGrid, nullptr, sizeof(float),
                            &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc,
                            newGridsSize, pmeGpu->archSpecific->pmeStream, true);
        kernelParamsPtr->grid.d_fourierGrid   = kernelParamsPtr->grid.d_realGrid;
        pmeGpu->archSpecific->complexGridSize = pmeGpu->archSpecific->realGridSize;
        // the size might get used later for copying the grid
    }
}

void pme_gpu_free_grids(const PmeGpu *pmeGpu)
{
    if (pmeGpu->archSpecific->performOutOfPlaceFFT)
    {
        cu_free_buffered(pmeGpu->kernelParams->grid.d_fourierGrid);
    }
    cu_free_buffered(pmeGpu->kernelParams->grid.d_realGrid,
                     &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc);
}

void pme_gpu_clear_grids(const PmeGpu *pmeGpu)
{
    cudaError_t stat = cudaMemsetAsync(pmeGpu->kernelParams->grid.d_realGrid, 0,
                                       pmeGpu->archSpecific->realGridSize * sizeof(float), pmeGpu->archSpecific->pmeStream);
    /* Should the complex grid be cleared in some weird case? */
    CU_RET_ERR(stat, "cudaMemsetAsync on the PME grid error");
}

void pme_gpu_realloc_and_copy_fract_shifts(PmeGpu *pmeGpu)
{
    pme_gpu_free_fract_shifts(pmeGpu);

    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();

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
                         pmeGpu->common->fsh.data(),
                         newFractShiftsSize,
                         pmeGpu->deviceInfo);

    initParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                         kernelParamsPtr->gridlineIndicesTableTexture,
                         &pme_gpu_get_gridline_texref(),
                         pmeGpu->common->nn.data(),
                         newFractShiftsSize,
                         pmeGpu->deviceInfo);
}

void pme_gpu_free_fract_shifts(const PmeGpu *pmeGpu)
{
    auto *kernelParamsPtr = pmeGpu->kernelParams.get();
    destroyParamLookupTable(kernelParamsPtr->grid.d_fractShiftsTable,
                            kernelParamsPtr->fractShiftsTableTexture,
                            &pme_gpu_get_fract_shifts_texref(),
                            pmeGpu->deviceInfo);
    destroyParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                            kernelParamsPtr->gridlineIndicesTableTexture,
                            &pme_gpu_get_gridline_texref(),
                            pmeGpu->deviceInfo);
}

bool pme_gpu_stream_query(const PmeGpu *pmeGpu)
{
    return haveStreamTasksCompleted(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_input_gather_grid(const PmeGpu *pmeGpu, float *h_grid)
{
    const size_t gridSize = pmeGpu->archSpecific->realGridSize * sizeof(float);
    cu_copy_H2D(pmeGpu->kernelParams->grid.d_realGrid, h_grid, gridSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_output_spread_grid(const PmeGpu *pmeGpu, float *h_grid)
{
    const size_t gridSize = pmeGpu->archSpecific->realGridSize * sizeof(float);
    cu_copy_D2H(h_grid, pmeGpu->kernelParams->grid.d_realGrid, gridSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
    cudaError_t  stat = cudaEventRecord(pmeGpu->archSpecific->syncSpreadGridD2H, pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME spread grid sync event record failure");
}

void pme_gpu_copy_output_spread_atom_data(const PmeGpu *pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesSize     = DIM * nAtomsPadded * pmeGpu->common->pme_order * sizeof(float);
    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();
    cu_copy_D2H(pmeGpu->staging.h_dtheta, kernelParamsPtr->atoms.d_dtheta, splinesSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
    cu_copy_D2H(pmeGpu->staging.h_theta, kernelParamsPtr->atoms.d_theta, splinesSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
    cu_copy_D2H(pmeGpu->staging.h_gridlineIndices, kernelParamsPtr->atoms.d_gridlineIndices,
                kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_input_gather_atom_data(const PmeGpu *pmeGpu)
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
    cu_copy_H2D(kernelParamsPtr->atoms.d_dtheta, pmeGpu->staging.h_dtheta, splinesSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
    cu_copy_H2D(kernelParamsPtr->atoms.d_theta, pmeGpu->staging.h_theta, splinesSize, pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
    cu_copy_H2D(kernelParamsPtr->atoms.d_gridlineIndices, pmeGpu->staging.h_gridlineIndices,
                kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), pmeGpu->settings.transferKind, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_sync_spread_grid(const PmeGpu *pmeGpu)
{
    cudaError_t stat = cudaEventSynchronize(pmeGpu->archSpecific->syncSpreadGridD2H);
    CU_RET_ERR(stat, "Error while waiting for the PME GPU spread grid to be copied to the host");
}

void pme_gpu_init_internal(PmeGpu *pmeGpu)
{
    /* Allocate the target-specific structures */
    pmeGpu->archSpecific.reset(new PmeGpuSpecific());
    pmeGpu->kernelParams.reset(new PmeGpuKernelParams());

    pmeGpu->archSpecific->performOutOfPlaceFFT = true;
    /* This should give better performance, according to the cuFFT documentation.
     * The performance seems to be the same though.
     * TODO: PME could also try to pick up nice grid sizes (with factors of 2, 3, 5, 7).
     */

    /* WARNING: CUDA timings are incorrect with multiple streams.
     *          This is the main reason why they are disabled by default.
     */
    // TODO: Consider turning on by default when we can detect nr of streams.
    pmeGpu->archSpecific->useTiming = (getenv("GMX_ENABLE_GPU_TIMING") != nullptr);

    /* Creating a PME CUDA stream */
    cudaError_t stat;
    int         highest_priority, lowest_priority;
    stat = cudaDeviceGetStreamPriorityRange(&lowest_priority, &highest_priority);
    CU_RET_ERR(stat, "PME cudaDeviceGetStreamPriorityRange failed");
    stat = cudaStreamCreateWithPriority(&pmeGpu->archSpecific->pmeStream,
                                        cudaStreamDefault, //cudaStreamNonBlocking,
                                        highest_priority);
    CU_RET_ERR(stat, "cudaStreamCreateWithPriority on the PME stream failed");
}

void pme_gpu_destroy_specific(const PmeGpu *pmeGpu)
{
    /* Destroy the CUDA stream */
    cudaError_t stat = cudaStreamDestroy(pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME cudaStreamDestroy error");
}

void pme_gpu_init_sync_events(const PmeGpu *pmeGpu)
{
    const auto  eventFlags = cudaEventDisableTiming;
    CU_RET_ERR(cudaEventCreateWithFlags(&pmeGpu->archSpecific->syncSpreadGridD2H, eventFlags), "cudaEventCreate on syncSpreadGridD2H failed");
}

void pme_gpu_destroy_sync_events(const PmeGpu *pmeGpu)
{
    CU_RET_ERR(cudaEventDestroy(pmeGpu->archSpecific->syncSpreadGridD2H), "cudaEventDestroy failed on syncSpreadGridD2H");
}

void pme_gpu_reinit_3dfft(const PmeGpu *pmeGpu)
{
    if (pme_gpu_performs_FFT(pmeGpu))
    {
        pmeGpu->archSpecific->fftSetup.resize(0);
        for (int i = 0; i < pmeGpu->common->ngrids; i++)
        {
            pmeGpu->archSpecific->fftSetup.push_back(std::unique_ptr<GpuParallel3dFft>(new GpuParallel3dFft(pmeGpu)));
        }
    }
}

void pme_gpu_destroy_3dfft(const PmeGpu *pmeGpu)
{
    pmeGpu->archSpecific->fftSetup.resize(0);
}
