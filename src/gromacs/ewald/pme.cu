/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
#include "gromacs/gpu_utils/devicebuffer.h"
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
    clearDeviceBufferAsync(&pmeGpu->kernelParams->constants.d_virialAndEnergy, 0,
                           c_virialAndEnergyCount, pmeGpu->archSpecific->pmeStream);
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
    reallocateDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli, newSplineValuesSize,
                           &pmeGpu->archSpecific->splineValuesSize, &pmeGpu->archSpecific->splineValuesSizeAlloc, pmeGpu->archSpecific->pmeStream);
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
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli, pmeGpu->staging.h_splineModuli,
                       0, newSplineValuesSize,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_free_bspline_values(const PmeGpu *pmeGpu)
{
    pfree(pmeGpu->staging.h_splineModuli);
    freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli);
}

void pme_gpu_realloc_forces(PmeGpu *pmeGpu)
{
    const size_t newForcesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces, newForcesSize,
                           &pmeGpu->archSpecific->forcesSize, &pmeGpu->archSpecific->forcesSizeAlloc, pmeGpu->archSpecific->pmeStream);
    pmeGpu->staging.h_forces.reserve(pmeGpu->nAtomsAlloc);
    pmeGpu->staging.h_forces.resize(pmeGpu->kernelParams->atoms.nAtoms);
}

void pme_gpu_free_forces(const PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces);
}

void pme_gpu_copy_input_forces(PmeGpu *pmeGpu)
{
    GMX_ASSERT(pmeGpu->kernelParams->atoms.nAtoms > 0, "Bad number of atoms in PME GPU");
    float *h_forcesFloat = reinterpret_cast<float *>(pmeGpu->staging.h_forces.data());
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces, h_forcesFloat,
                       0, DIM * pmeGpu->kernelParams->atoms.nAtoms,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_output_forces(PmeGpu *pmeGpu)
{
    GMX_ASSERT(pmeGpu->kernelParams->atoms.nAtoms > 0, "Bad number of atoms in PME GPU");
    float *h_forcesFloat = reinterpret_cast<float *>(pmeGpu->staging.h_forces.data());
    copyFromDeviceBuffer(h_forcesFloat, &pmeGpu->kernelParams->atoms.d_forces,
                         0, DIM * pmeGpu->kernelParams->atoms.nAtoms,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_realloc_coordinates(const PmeGpu *pmeGpu)
{
    const size_t newCoordinatesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newCoordinatesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates, newCoordinatesSize,
                           &pmeGpu->archSpecific->coordinatesSize, &pmeGpu->archSpecific->coordinatesSizeAlloc, pmeGpu->archSpecific->pmeStream);
    if (c_usePadding)
    {
        const size_t paddingIndex = DIM * pmeGpu->kernelParams->atoms.nAtoms;
        const size_t paddingCount = DIM * pmeGpu->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            clearDeviceBufferAsync(&pmeGpu->kernelParams->atoms.d_coordinates, paddingIndex,
                                   paddingCount, pmeGpu->archSpecific->pmeStream);
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
    const float *h_coordinatesFloat = reinterpret_cast<const float *>(h_coordinates);
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates, h_coordinatesFloat,
                       0, pmeGpu->kernelParams->atoms.nAtoms * DIM,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
#endif
}

void pme_gpu_free_coordinates(const PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates);
}

void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu *pmeGpu, const float *h_coefficients)
{
    GMX_ASSERT(h_coefficients, "Bad host-side charge buffer in PME GPU");
    const size_t newCoefficientsSize = pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newCoefficientsSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients, newCoefficientsSize,
                           &pmeGpu->archSpecific->coefficientsSize, &pmeGpu->archSpecific->coefficientsSizeAlloc, pmeGpu->archSpecific->pmeStream);
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients, const_cast<float *>(h_coefficients),
                       0, pmeGpu->kernelParams->atoms.nAtoms,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    if (c_usePadding)
    {
        const size_t paddingIndex = pmeGpu->kernelParams->atoms.nAtoms;
        const size_t paddingCount = pmeGpu->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            clearDeviceBufferAsync(&pmeGpu->kernelParams->atoms.d_coefficients, paddingIndex,
                                   paddingCount, pmeGpu->archSpecific->pmeStream);
        }
    }
}

void pme_gpu_free_coefficients(const PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients);
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
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_theta, newSplineDataSize,
                           &currentSizeTemp, &currentSizeTempAlloc, pmeGpu->archSpecific->pmeStream);
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta, newSplineDataSize,
                           &pmeGpu->archSpecific->splineDataSize, &pmeGpu->archSpecific->splineDataSizeAlloc, pmeGpu->archSpecific->pmeStream);
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
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_theta);
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta);
    pfree(pmeGpu->staging.h_theta);
    pfree(pmeGpu->staging.h_dtheta);
}

void pme_gpu_realloc_grid_indices(const PmeGpu *pmeGpu)
{
    const size_t newIndicesSize = DIM * pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newIndicesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices, newIndicesSize,
                           &pmeGpu->archSpecific->gridlineIndicesSize, &pmeGpu->archSpecific->gridlineIndicesSizeAlloc, pmeGpu->archSpecific->pmeStream);
    pfree(pmeGpu->staging.h_gridlineIndices);
    pmalloc((void **)&pmeGpu->staging.h_gridlineIndices, newIndicesSize * sizeof(int));
}

void pme_gpu_free_grid_indices(const PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices);
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
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_fourierGrid, newComplexGridSize,
                               &pmeGpu->archSpecific->complexGridSize, &pmeGpu->archSpecific->complexGridSizeAlloc, pmeGpu->archSpecific->pmeStream);
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid, newRealGridSize,
                               &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc, pmeGpu->archSpecific->pmeStream);
    }
    else
    {
        /* A single buffer so that any grid will fit */
        const int newGridsSize = std::max(newRealGridSize, newComplexGridSize);
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid, newGridsSize,
                               &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc, pmeGpu->archSpecific->pmeStream);
        kernelParamsPtr->grid.d_fourierGrid   = kernelParamsPtr->grid.d_realGrid;
        pmeGpu->archSpecific->complexGridSize = pmeGpu->archSpecific->realGridSize;
        // the size might get used later for copying the grid
    }
}

void pme_gpu_free_grids(const PmeGpu *pmeGpu)
{
    if (pmeGpu->archSpecific->performOutOfPlaceFFT)
    {
        freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_fourierGrid);
    }
    freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid);
}

void pme_gpu_clear_grids(const PmeGpu *pmeGpu)
{
    clearDeviceBufferAsync(&pmeGpu->kernelParams->grid.d_realGrid, 0,
                           pmeGpu->archSpecific->realGridSize, pmeGpu->archSpecific->pmeStream);
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
                         pmeGpu->common->fsh.data(),
                         newFractShiftsSize,
                         pmeGpu->deviceInfo);

    initParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                         kernelParamsPtr->gridlineIndicesTableTexture,
                         pmeGpu->common->nn.data(),
                         newFractShiftsSize,
                         pmeGpu->deviceInfo);
}

void pme_gpu_free_fract_shifts(const PmeGpu *pmeGpu)
{
    auto *kernelParamsPtr = pmeGpu->kernelParams.get();
    destroyParamLookupTable(kernelParamsPtr->grid.d_fractShiftsTable,
                            kernelParamsPtr->fractShiftsTableTexture,
                            pmeGpu->deviceInfo);
    destroyParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                            kernelParamsPtr->gridlineIndicesTableTexture,
                            pmeGpu->deviceInfo);
}

bool pme_gpu_stream_query(const PmeGpu *pmeGpu)
{
    return haveStreamTasksCompleted(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_input_gather_grid(const PmeGpu *pmeGpu, float *h_grid)
{
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid, h_grid,
                       0, pmeGpu->archSpecific->realGridSize,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_output_spread_grid(const PmeGpu *pmeGpu, float *h_grid)
{
    copyFromDeviceBuffer(h_grid, &pmeGpu->kernelParams->grid.d_realGrid,
                         0, pmeGpu->archSpecific->realGridSize,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    cudaError_t  stat = cudaEventRecord(pmeGpu->archSpecific->syncSpreadGridD2H, pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME spread grid sync event record failure");
}

void pme_gpu_copy_output_spread_atom_data(const PmeGpu *pmeGpu)
{
    const int    alignment        = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded     = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesCount     = DIM * nAtomsPadded * pmeGpu->common->pme_order;
    auto        *kernelParamsPtr  = pmeGpu->kernelParams.get();
    copyFromDeviceBuffer(pmeGpu->staging.h_dtheta, &kernelParamsPtr->atoms.d_dtheta,
                         0, splinesCount,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_theta, &kernelParamsPtr->atoms.d_theta,
                         0, splinesCount,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_gridlineIndices, &kernelParamsPtr->atoms.d_gridlineIndices,
                         0, kernelParamsPtr->atoms.nAtoms * DIM,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_input_gather_atom_data(const PmeGpu *pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesCount    = DIM * nAtomsPadded * pmeGpu->common->pme_order;
    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();
    if (c_usePadding)
    {
        // TODO: could clear only the padding and not the whole thing, but this is a test-exclusive code anyway
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_gridlineIndices, 0,
                               pmeGpu->nAtomsAlloc * DIM, pmeGpu->archSpecific->pmeStream);
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_dtheta, 0,
                               pmeGpu->nAtomsAlloc * pmeGpu->common->pme_order * DIM, pmeGpu->archSpecific->pmeStream);
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_theta, 0,
                               pmeGpu->nAtomsAlloc * pmeGpu->common->pme_order * DIM, pmeGpu->archSpecific->pmeStream);
    }
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_dtheta, pmeGpu->staging.h_dtheta,
                       0, splinesCount,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_theta, pmeGpu->staging.h_theta,
                       0, splinesCount,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_gridlineIndices, pmeGpu->staging.h_gridlineIndices,
                       0, kernelParamsPtr->atoms.nAtoms * DIM,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
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

    pmeGpu->maxGridWidthX = pmeGpu->deviceInfo->prop.maxGridSize[0];

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

int getSplineParamFullIndex(int order, int splineIndex, int dimIndex, int warpIndex, int atomWarpIndex)
{
    if (order != 4)
    {
        throw order;
    }
    constexpr int fixedOrder = 4;
    GMX_UNUSED_VALUE(fixedOrder);
    const int     indexBase  = getSplineParamIndexBase<fixedOrder>(warpIndex, atomWarpIndex);
    return getSplineParamIndex<fixedOrder>(indexBase, dimIndex, splineIndex);
}
