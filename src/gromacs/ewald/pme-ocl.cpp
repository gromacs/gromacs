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

//#include "gromacs/gpu_utils/cudautils.cuh"
//#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/gpu_utils/ocl_compiler.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme-types-ocl.h"
#include "pme-3dfft-ocl.h"
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
    gpuStreamSynchronize(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_alloc_energy_virial(const PmeGpu *pmeGpu)
{
    const size_t energyAndVirialSize = c_virialAndEnergyCount * sizeof(float);
    allocateDeviceBuffer(&pmeGpu->kernelParams->constants.d_virialAndEnergy, c_virialAndEnergyCount, pmeGpu->archSpecific->persistent->context);
    ocl_pmalloc((void **)&pmeGpu->staging.h_virialAndEnergy, energyAndVirialSize);
}

void pme_gpu_free_energy_virial(PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->constants.d_virialAndEnergy);
    ocl_pfree(pmeGpu->staging.h_virialAndEnergy);
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
                           &pmeGpu->archSpecific->splineValuesSize, &pmeGpu->archSpecific->splineValuesSizeAlloc, pmeGpu->archSpecific->persistent->context);
    if (shouldRealloc)
    {
        /* Reallocate the host buffer */
        ocl_pfree(pmeGpu->staging.h_splineModuli);
        ocl_pmalloc((void **)&pmeGpu->staging.h_splineModuli, newSplineValuesSize * sizeof(float));
    }
    for (int i = 0; i < DIM; i++)
    {
        memcpy(pmeGpu->staging.h_splineModuli + splineValuesOffset[i], pmeGpu->common->bsp_mod[i].data(), pmeGpu->common->bsp_mod[i].size() * sizeof(float));
    }
    /* TODO: pin original buffer instead! */
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli, pmeGpu->staging.h_splineModuli, 0,
                newSplineValuesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_free_bspline_values(const PmeGpu *pmeGpu)
{
    ocl_pfree(pmeGpu->staging.h_splineModuli);
    freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli);
}

void pme_gpu_realloc_forces(PmeGpu *pmeGpu)
{
    const size_t newForcesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces, newForcesSize,
                           &pmeGpu->archSpecific->forcesSize, &pmeGpu->archSpecific->forcesSizeAlloc, pmeGpu->archSpecific->persistent->context);
    pmeGpu->staging.h_forces.reserve(pmeGpu->nAtomsAlloc);
    pmeGpu->staging.h_forces.resize(pmeGpu->kernelParams->atoms.nAtoms);
}

void pme_gpu_free_forces(const PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces);
}

template <typename Thing>
inline float *FIXMEcast(Thing *thing)
{
    return (float *)thing;
}

void pme_gpu_copy_input_forces(PmeGpu *pmeGpu)
{
    const size_t forcesSize = DIM * pmeGpu->kernelParams->atoms.nAtoms;
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces, FIXMEcast(pmeGpu->staging.h_forces.data()), 0,
                       forcesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_output_forces(PmeGpu *pmeGpu)
{
    const size_t forcesSize   = DIM * pmeGpu->kernelParams->atoms.nAtoms;
    GMX_ASSERT(forcesSize > 0, "Bad number of atoms in PME GPU");
    copyFromDeviceBuffer(FIXMEcast(pmeGpu->staging.h_forces.data()), &pmeGpu->kernelParams->atoms.d_forces, 0,
                         forcesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_realloc_coordinates(const PmeGpu *pmeGpu)
{
    const size_t newCoordinatesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newCoordinatesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates, newCoordinatesSize,
                           &pmeGpu->archSpecific->coordinatesSize, &pmeGpu->archSpecific->coordinatesSizeAlloc, pmeGpu->archSpecific->persistent->context);
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
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates, FIXMEcast(h_coordinates), 0,
                pmeGpu->kernelParams->atoms.nAtoms * DIM, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
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
                           &pmeGpu->archSpecific->coefficientsSize, &pmeGpu->archSpecific->coefficientsSizeAlloc, pmeGpu->archSpecific->persistent->context);
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients, FIXMEcast(h_coefficients), 0,
                pmeGpu->kernelParams->atoms.nAtoms, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
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
                           &currentSizeTemp, &currentSizeTempAlloc, pmeGpu->archSpecific->persistent->context);
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta, newSplineDataSize,
                           &pmeGpu->archSpecific->splineDataSize, &pmeGpu->archSpecific->splineDataSizeAlloc, pmeGpu->archSpecific->persistent->context);
    // the host side reallocation
    if (shouldRealloc)
    {
        ocl_pfree(pmeGpu->staging.h_theta);
        ocl_pmalloc((void **)&pmeGpu->staging.h_theta, newSplineDataSize * sizeof(float));
        ocl_pfree(pmeGpu->staging.h_dtheta);
        ocl_pmalloc((void **)&pmeGpu->staging.h_dtheta, newSplineDataSize * sizeof(float));
    }
}

void pme_gpu_free_spline_data(const PmeGpu *pmeGpu)
{
    /* Two arrays of the same size */
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_theta);
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta);
    ocl_pfree(pmeGpu->staging.h_theta);
    ocl_pfree(pmeGpu->staging.h_dtheta);
}

void pme_gpu_realloc_grid_indices(const PmeGpu *pmeGpu)
{
    const size_t newIndicesSize = DIM * pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newIndicesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices, newIndicesSize,
                           &pmeGpu->archSpecific->gridlineIndicesSize, &pmeGpu->archSpecific->gridlineIndicesSizeAlloc, pmeGpu->archSpecific->persistent->context);
    ocl_pfree(pmeGpu->staging.h_gridlineIndices);
    ocl_pmalloc((void **)&pmeGpu->staging.h_gridlineIndices, newIndicesSize * sizeof(int));
}

void pme_gpu_free_grid_indices(const PmeGpu *pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices);
    ocl_pfree(pmeGpu->staging.h_gridlineIndices);
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
                               &pmeGpu->archSpecific->complexGridSize, &pmeGpu->archSpecific->complexGridSizeAlloc, pmeGpu->archSpecific->persistent->context);
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid, newRealGridSize,
                               &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc, pmeGpu->archSpecific->persistent->context);
    }
    else
    {
        /* A single buffer so that any grid will fit */
        const int newGridsSize = std::max(newRealGridSize, newComplexGridSize);
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid, newGridsSize,
                               &pmeGpu->archSpecific->realGridSize, &pmeGpu->archSpecific->realGridSizeAlloc, pmeGpu->archSpecific->persistent->context);
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

    initParamLookupTable(&kernelParamsPtr->grid.d_fractShiftsTable,
                         kernelParamsPtr->fractShiftsTableTexture,
                         pmeGpu->common->fsh.data(),
                         newFractShiftsSize,
                         pmeGpu->deviceInfo,
                         pmeGpu->archSpecific->persistent->context,
			 pmeGpu->archSpecific->pmeStream);

    initParamLookupTable(&kernelParamsPtr->grid.d_gridlineIndicesTable,
                         kernelParamsPtr->gridlineIndicesTableTexture,
                         pmeGpu->common->nn.data(),
                         newFractShiftsSize,
                         pmeGpu->deviceInfo,
                         pmeGpu->archSpecific->persistent->context,
			 pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_free_fract_shifts(const PmeGpu *pmeGpu)
{
    auto *kernelParamsPtr = pmeGpu->kernelParams.get();
    destroyParamLookupTable(&kernelParamsPtr->grid.d_fractShiftsTable,
                            kernelParamsPtr->fractShiftsTableTexture,
                            pmeGpu->deviceInfo);
    destroyParamLookupTable(&kernelParamsPtr->grid.d_gridlineIndicesTable,
                            kernelParamsPtr->gridlineIndicesTableTexture,
                            pmeGpu->deviceInfo);
}

bool pme_gpu_stream_query(const PmeGpu *pmeGpu)
{
    return haveStreamTasksCompleted(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_input_gather_grid(const PmeGpu *pmeGpu, float *h_grid)
{
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid, h_grid, 0,
                       pmeGpu->archSpecific->realGridSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_output_spread_grid(const PmeGpu *pmeGpu, float *h_grid)
{
    copyFromDeviceBuffer(h_grid, &pmeGpu->kernelParams->grid.d_realGrid, 0, pmeGpu->archSpecific->realGridSize,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    //FIXME just pass it directly into API instead?
    pmeGpu->archSpecific->syncSpreadGridD2H.markSyncPoint(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_output_spread_atom_data(PmeGpu *pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesSize     = DIM * nAtomsPadded * pmeGpu->common->pme_order;
    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();

    copyFromDeviceBuffer(pmeGpu->staging.h_dtheta, &kernelParamsPtr->atoms.d_dtheta, 0,
                         splinesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_theta, &kernelParamsPtr->atoms.d_theta, 0,
                         splinesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_gridlineIndices, &kernelParamsPtr->atoms.d_gridlineIndices, 0,
                         kernelParamsPtr->atoms.nAtoms * DIM, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_input_gather_atom_data(const PmeGpu *pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesSize     = DIM * nAtomsPadded * pmeGpu->common->pme_order;
    auto        *kernelParamsPtr = pmeGpu->kernelParams.get();
    if (c_usePadding)
    {
        const size_t splineValuesPerAtom = pmeGpu->common->pme_order * DIM;
        // TODO: could clear only the padding and not the whole thing, but this is a test-exclusive code anyway
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_gridlineIndices, 0,
                                pmeGpu->nAtomsAlloc * DIM, pmeGpu->archSpecific->pmeStream);
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_dtheta, 0,
                                pmeGpu->nAtomsAlloc * splineValuesPerAtom, pmeGpu->archSpecific->pmeStream);
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_theta, 0,
                                pmeGpu->nAtomsAlloc * splineValuesPerAtom, pmeGpu->archSpecific->pmeStream);
    }
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_dtheta, pmeGpu->staging.h_dtheta, 0,
                       splinesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_theta, pmeGpu->staging.h_theta, 0,
                       splinesSize, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_gridlineIndices, pmeGpu->staging.h_gridlineIndices, 0,
                kernelParamsPtr->atoms.nAtoms * DIM, pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_sync_spread_grid(const PmeGpu *pmeGpu)
{
    pmeGpu->archSpecific->syncSpreadGridD2H.waitForSyncPoint();
}

void pme_gpu_init_internal(PmeGpu *pmeGpu, PmePersistentDataHandle persistent)
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
    pmeGpu->archSpecific->useTiming = (getenv("GMX_ENABLE_GPU_TIMING") != nullptr); //FIXME: DISABLE for OpenCL?


   pmeGpu->archSpecific->persistent = persistent;
    if (!pmeGpu->archSpecific->persistent)
    pmeGpu->archSpecific->persistent = std::make_shared<PmeGpuPersistentData>(pmeGpu);

    
    /* Creating a GPU stream */

    // TODO wrapper;
    // TODO priorities/out of order? CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
#if GMX_GPU == GMX_GPU_OPENCL
    cl_command_queue_properties queueProperties = pmeGpu->archSpecific->useTiming ? CL_QUEUE_PROFILING_ENABLE : 0;

    //FIXME copypaste
    cl_platform_id            platform_id;
    cl_device_id              device_id;
    cl_int                    clError;

    platform_id      = pmeGpu->deviceInfo->ocl_gpu_id.ocl_platform_id;
    device_id        = pmeGpu->deviceInfo->ocl_gpu_id.ocl_device_id;

    
    /* local/non-local GPU streams */
    pmeGpu->archSpecific->pmeStream = clCreateCommandQueue(pmeGpu->archSpecific->persistent->context,
                                                           device_id, queueProperties, &clError);
    GMX_RELEASE_ASSERT(clError == CL_SUCCESS, "Failed to create command queue");
#endif
#if GMX_GPU == GMX_GPU_CUDA
    cudaError_t stat;
    int         highest_priority, lowest_priority;
    stat = cudaDeviceGetStreamPriorityRange(&lowest_priority, &highest_priority);
    CU_RET_ERR(stat, "PME cudaDeviceGetStreamPriorityRange failed");
    stat = cudaStreamCreateWithPriority(&pmeGpu->archSpecific->pmeStream,
                                        cudaStreamDefault, //cudaStreamNonBlocking,
                                        highest_priority);
    CU_RET_ERR(stat, "cudaStreamCreateWithPriority on the PME stream failed");
#endif
}

void pme_gpu_destroy_specific(const PmeGpu *pmeGpu)
{
    //FIXME do we care abotu errors here at all?
    //FIXME retain all the stuff for the unit tests???

    /*
    */

    /* Free command queues */
    cl_int clError = clReleaseCommandQueue(pmeGpu->archSpecific->pmeStream);
    //throwUponFailure(clError);
    //GMX_RELEASE_ASSERT(clError == CL_SUCCESS, "PME stream destruction error");

    //FIXME
    /*
    free_gpu_device_runtime_data(nb->dev_rundata);
    sfree(nb->dev_rundata);
    */
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

//FIXME: these guys are here because they need to be instantiated for pme.cpp calls

void pme_gpu_destroy(PmeGpu *pmeGpu)
{
    /* Free lots of data */
    pme_gpu_free_energy_virial(pmeGpu);
    pme_gpu_free_bspline_values(pmeGpu);
    pme_gpu_free_forces(pmeGpu);
    pme_gpu_free_coordinates(pmeGpu);
    pme_gpu_free_coefficients(pmeGpu);
    pme_gpu_free_spline_data(pmeGpu);
    pme_gpu_free_grid_indices(pmeGpu);
    pme_gpu_free_fract_shifts(pmeGpu);
    pme_gpu_free_grids(pmeGpu);

    pme_gpu_destroy_3dfft(pmeGpu);

    /* Free the GPU-framework specific data last */
    pme_gpu_destroy_specific(pmeGpu);

    delete pmeGpu;
}

void pme_gpu_reinit(gmx_pme_t *pme, gmx_device_info_t *gpuInfo, PmePersistentDataHandle persistent)
{
    if (!pme_gpu_active(pme))
    {
        return;
    }

    if (!pme->gpu)
    {
        /* First-time initialization */
      pme_gpu_init(pme, gpuInfo, persistent);
    }
    else
    {
        /* After this call nothing in the GPU code should refer to the gmx_pme_t *pme itself - until the next pme_gpu_reinit */
        pme_gpu_copy_common_data_from(pme);
    }
    /* GPU FFT will only get used for a single rank.*/
    pme->gpu->settings.performGPUFFT   = (pme->gpu->common->runMode == PmeRunMode::GPU) && !pme_gpu_uses_dd(pme->gpu);
    pme->gpu->settings.performGPUSolve = (pme->gpu->common->runMode == PmeRunMode::GPU);

    /* Reinit active timers */
    pme_gpu_reinit_timings(pme->gpu);

    pme_gpu_reinit_grids(pme->gpu);
    pme_gpu_reinit_computation(pme->gpu);
    /* Clear the previous box - doesn't hurt, and forces the PME CPU recipbox
     * update for mixed mode on grid switch. TODO: use shared recipbox field.
     */
    std::memset(pme->gpu->common->previousBox, 0, sizeof(pme->gpu->common->previousBox));
}

void pme_gpu_reinit_atoms(PmeGpu *pmeGpu, const int nAtoms, const real *charges)
{
    auto      *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->atoms.nAtoms = nAtoms;
    const int  alignment = pme_gpu_get_atom_data_alignment(pmeGpu);
    pmeGpu->nAtomsPadded = ((nAtoms + alignment - 1) / alignment) * alignment;
    const int  nAtomsAlloc   = c_usePadding ? pmeGpu->nAtomsPadded : nAtoms;
    const bool haveToRealloc = (pmeGpu->nAtomsAlloc < nAtomsAlloc); /* This check might be redundant, but is logical */
    pmeGpu->nAtomsAlloc = nAtomsAlloc;

#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision supported");
    GMX_UNUSED_VALUE(charges);
#else
    pme_gpu_realloc_and_copy_input_coefficients(pmeGpu, reinterpret_cast<const float *>(charges));
    /* Could also be checked for haveToRealloc, but the copy always needs to be performed */
#endif

    if (haveToRealloc)
    {
        pme_gpu_realloc_coordinates(pmeGpu);
        pme_gpu_realloc_forces(pmeGpu);
        pme_gpu_realloc_spline_data(pmeGpu);
        pme_gpu_realloc_grid_indices(pmeGpu);
    }
}

void pme_gpu_get_real_grid_sizes(const PmeGpu *pmeGpu, gmx::IVec *gridSize, gmx::IVec *paddedGridSize)
{
    GMX_ASSERT(gridSize != nullptr, "");
    GMX_ASSERT(paddedGridSize != nullptr, "");
    GMX_ASSERT(pmeGpu != nullptr, "");
    auto *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    for (int i = 0; i < DIM; i++)
    {
        (*gridSize)[i]       = kernelParamsPtr->grid.realGridSize[i];
        (*paddedGridSize)[i] = kernelParamsPtr->grid.realGridSizePadded[i];
    }
}

void pme_gpu_transform_spline_atom_data(const PmeGpu *pmeGpu, const pme_atomcomm_t *atc,
                                        PmeSplineDataType type, int dimIndex, PmeLayoutTransform transform)
{
    // The GPU atom spline data is laid out in a different way currently than the CPU one.
    // This function converts the data from GPU to CPU layout (in the host memory).
    // It is only intended for testing purposes so far.
    // Ideally we should use similar layouts on CPU and GPU if we care about mixed modes and their performance
    // (e.g. spreading on GPU, gathering on CPU).
    GMX_RELEASE_ASSERT(atc->nthread == 1, "Only the serial PME data layout is supported");
    const uintmax_t threadIndex  = 0;
    const auto      atomCount    = pme_gpu_get_kernel_params_base_ptr(pmeGpu)->atoms.nAtoms;
    const auto      atomsPerWarp = pme_gpu_get_atoms_per_warp(pmeGpu);
    const auto      pmeOrder     = pmeGpu->common->pme_order;

    real           *cpuSplineBuffer;
    float          *h_splineBuffer;
    switch (type)
    {
        case PmeSplineDataType::Values:
            cpuSplineBuffer = atc->spline[threadIndex].theta[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_theta;
            break;

        case PmeSplineDataType::Derivatives:
            cpuSplineBuffer = atc->spline[threadIndex].dtheta[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_dtheta;
            break;

        default:
            GMX_THROW(gmx::InternalError("Unknown spline data type"));
    }

    for (auto atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        auto atomWarpIndex = atomIndex % atomsPerWarp;
        auto warpIndex     = atomIndex / atomsPerWarp;
       
        for (auto orderIndex = 0; orderIndex < pmeOrder; orderIndex++)
        {
            const auto gpuValueIndex = ((pmeOrder * warpIndex + orderIndex) * DIM + dimIndex) * atomsPerWarp + atomWarpIndex;
            const auto cpuValueIndex = atomIndex * pmeOrder + orderIndex;
            GMX_ASSERT(cpuValueIndex < atomCount * pmeOrder, "Atom spline data index out of bounds (while transforming GPU data layout for host)");
            switch (transform)
            {
                case PmeLayoutTransform::GpuToHost:
                    cpuSplineBuffer[cpuValueIndex] = h_splineBuffer[gpuValueIndex];
                    break;

                case PmeLayoutTransform::HostToGpu:
                    h_splineBuffer[gpuValueIndex] = cpuSplineBuffer[cpuValueIndex];
                    break;

                default:
                    GMX_THROW(gmx::InternalError("Unknown layout transform"));
            }
        }
    }
}

#include "gromacs/utility/arrayref.h"

gmx::ArrayRef<const gmx::IVec> pmeGpuGetGridlineIndices(const PmeGpu *pmeGpu)
{
    return arrayRefFromArray(reinterpret_cast<gmx::IVec *>(pmeGpu->staging.h_gridlineIndices), pmeGpu->kernelParams->atoms.nAtoms);
}

void pmeGpuSetGridlineIndices(PmeGpu *pmeGpu, gmx::ArrayRef<const gmx::IVec> gridlineIndices)
{
    memcpy(pmeGpu->staging.h_gridlineIndices, gridlineIndices.data(), gridlineIndices.size() * sizeof(gridlineIndices[0]));
}

PmePersistentDataHandle pmeGpuAcquirePersistentData(PmeGpu *pmeGpu)
{
    GMX_ASSERT(pmeGpu != nullptr, "PME GPU must have been initialized");
    return pmeGpu->archSpecific->persistent;
}
