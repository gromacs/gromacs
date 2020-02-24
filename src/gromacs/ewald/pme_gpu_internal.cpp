/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020, by the GROMACS development team, led by
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
 *
 * \brief This file contains internal function implementations
 * for performing the PME calculations on GPU.
 *
 * Note that this file is compiled as regular C++ source in OpenCL builds, but
 * it is treated as CUDA source in CUDA-enabled GPU builds.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "pme_gpu_internal.h"

#include "config.h"

#include <list>
#include <memory>
#include <string>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#if GMX_GPU == GMX_GPU_CUDA
#    include "gromacs/gpu_utils/pmalloc_cuda.h"

#    include "pme.cuh"
#elif GMX_GPU == GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gmxopencl.h"
#endif

#include "gromacs/ewald/pme.h"

#include "pme_gpu_3dfft.h"
#include "pme_gpu_constants.h"
#include "pme_gpu_program_impl.h"
#include "pme_gpu_timings.h"
#include "pme_gpu_types.h"
#include "pme_gpu_types_host.h"
#include "pme_gpu_types_host_impl.h"
#include "pme_gpu_utils.h"
#include "pme_grid.h"
#include "pme_internal.h"
#include "pme_solve.h"

/*! \brief
 * CUDA only
 * Atom limit above which it is advantageous to turn on the
 * recalcuating of the splines in the gather and using less threads per atom in the spline and spread
 */
constexpr int c_pmeGpuPerformanceAtomLimit = 23000;

/*! \internal \brief
 * Wrapper for getting a pointer to the plain C++ part of the GPU kernel parameters structure.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 * \returns The pointer to the kernel parameters.
 */
static PmeGpuKernelParamsBase* pme_gpu_get_kernel_params_base_ptr(const PmeGpu* pmeGpu)
{
    // reinterpret_cast is needed because the derived CUDA structure is not known in this file
    auto* kernelParamsPtr = reinterpret_cast<PmeGpuKernelParamsBase*>(pmeGpu->kernelParams.get());
    return kernelParamsPtr;
}

int pme_gpu_get_atom_data_alignment(const PmeGpu* /*unused*/)
{
    // TODO: this can be simplified, as c_pmeAtomDataAlignment is now constant
    if (c_usePadding)
    {
        return c_pmeAtomDataAlignment;
    }
    else
    {
        return 0;
    }
}

int pme_gpu_get_atoms_per_warp(const PmeGpu* pmeGpu)
{
    if (pmeGpu->settings.useOrderThreadsPerAtom)
    {
        return pmeGpu->programHandle_->impl_->warpSize / c_pmeSpreadGatherThreadsPerAtom4ThPerAtom;
    }
    else
    {
        return pmeGpu->programHandle_->impl_->warpSize / c_pmeSpreadGatherThreadsPerAtom;
    }
}

void pme_gpu_synchronize(const PmeGpu* pmeGpu)
{
    gpuStreamSynchronize(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_alloc_energy_virial(PmeGpu* pmeGpu)
{
    const size_t energyAndVirialSize = c_virialAndEnergyCount * sizeof(float);
    allocateDeviceBuffer(&pmeGpu->kernelParams->constants.d_virialAndEnergy, c_virialAndEnergyCount,
                         pmeGpu->archSpecific->context);
    pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_virialAndEnergy), energyAndVirialSize);
}

void pme_gpu_free_energy_virial(PmeGpu* pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->constants.d_virialAndEnergy);
    pfree(pmeGpu->staging.h_virialAndEnergy);
    pmeGpu->staging.h_virialAndEnergy = nullptr;
}

void pme_gpu_clear_energy_virial(const PmeGpu* pmeGpu)
{
    clearDeviceBufferAsync(&pmeGpu->kernelParams->constants.d_virialAndEnergy, 0,
                           c_virialAndEnergyCount, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_realloc_and_copy_bspline_values(PmeGpu* pmeGpu)
{
    const int splineValuesOffset[DIM] = { 0, pmeGpu->kernelParams->grid.realGridSize[XX],
                                          pmeGpu->kernelParams->grid.realGridSize[XX]
                                                  + pmeGpu->kernelParams->grid.realGridSize[YY] };
    memcpy(&pmeGpu->kernelParams->grid.splineValuesOffset, &splineValuesOffset, sizeof(splineValuesOffset));

    const int newSplineValuesSize = pmeGpu->kernelParams->grid.realGridSize[XX]
                                    + pmeGpu->kernelParams->grid.realGridSize[YY]
                                    + pmeGpu->kernelParams->grid.realGridSize[ZZ];
    const bool shouldRealloc = (newSplineValuesSize > pmeGpu->archSpecific->splineValuesSize);
    reallocateDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli, newSplineValuesSize,
                           &pmeGpu->archSpecific->splineValuesSize,
                           &pmeGpu->archSpecific->splineValuesSizeAlloc, pmeGpu->archSpecific->context);
    if (shouldRealloc)
    {
        /* Reallocate the host buffer */
        pfree(pmeGpu->staging.h_splineModuli);
        pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_splineModuli),
                newSplineValuesSize * sizeof(float));
    }
    for (int i = 0; i < DIM; i++)
    {
        memcpy(pmeGpu->staging.h_splineModuli + splineValuesOffset[i],
               pmeGpu->common->bsp_mod[i].data(), pmeGpu->common->bsp_mod[i].size() * sizeof(float));
    }
    /* TODO: pin original buffer instead! */
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli, pmeGpu->staging.h_splineModuli,
                       0, newSplineValuesSize, pmeGpu->archSpecific->pmeStream,
                       pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_free_bspline_values(const PmeGpu* pmeGpu)
{
    pfree(pmeGpu->staging.h_splineModuli);
    freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli);
}

void pme_gpu_realloc_forces(PmeGpu* pmeGpu)
{
    const size_t newForcesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces, newForcesSize,
                           &pmeGpu->archSpecific->forcesSize,
                           &pmeGpu->archSpecific->forcesSizeAlloc, pmeGpu->archSpecific->context);
    pmeGpu->staging.h_forces.reserveWithPadding(pmeGpu->nAtomsAlloc);
    pmeGpu->staging.h_forces.resizeWithPadding(pmeGpu->kernelParams->atoms.nAtoms);
}

void pme_gpu_free_forces(const PmeGpu* pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces);
}

void pme_gpu_copy_input_forces(PmeGpu* pmeGpu)
{
    GMX_ASSERT(pmeGpu->kernelParams->atoms.nAtoms > 0, "Bad number of atoms in PME GPU");
    float* h_forcesFloat = reinterpret_cast<float*>(pmeGpu->staging.h_forces.data());
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces, h_forcesFloat, 0,
                       DIM * pmeGpu->kernelParams->atoms.nAtoms, pmeGpu->archSpecific->pmeStream,
                       pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_output_forces(PmeGpu* pmeGpu)
{
    GMX_ASSERT(pmeGpu->kernelParams->atoms.nAtoms > 0, "Bad number of atoms in PME GPU");
    float* h_forcesFloat = reinterpret_cast<float*>(pmeGpu->staging.h_forces.data());
    copyFromDeviceBuffer(h_forcesFloat, &pmeGpu->kernelParams->atoms.d_forces, 0,
                         DIM * pmeGpu->kernelParams->atoms.nAtoms, pmeGpu->archSpecific->pmeStream,
                         pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_realloc_coordinates(const PmeGpu* pmeGpu)
{
    const size_t newCoordinatesSize = pmeGpu->nAtomsAlloc * DIM;
    GMX_ASSERT(newCoordinatesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates, newCoordinatesSize,
                           &pmeGpu->archSpecific->coordinatesSize,
                           &pmeGpu->archSpecific->coordinatesSizeAlloc, pmeGpu->archSpecific->context);
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

void pme_gpu_free_coordinates(const PmeGpu* pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coordinates);
}

void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu* pmeGpu, const float* h_coefficients)
{
    GMX_ASSERT(h_coefficients, "Bad host-side charge buffer in PME GPU");
    const size_t newCoefficientsSize = pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newCoefficientsSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients, newCoefficientsSize,
                           &pmeGpu->archSpecific->coefficientsSize,
                           &pmeGpu->archSpecific->coefficientsSizeAlloc, pmeGpu->archSpecific->context);
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients,
                       const_cast<float*>(h_coefficients), 0, pmeGpu->kernelParams->atoms.nAtoms,
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

void pme_gpu_free_coefficients(const PmeGpu* pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients);
}

void pme_gpu_realloc_spline_data(PmeGpu* pmeGpu)
{
    const int    order        = pmeGpu->common->pme_order;
    const int    alignment    = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const int    newSplineDataSize = DIM * order * nAtomsPadded;
    GMX_ASSERT(newSplineDataSize > 0, "Bad number of atoms in PME GPU");
    /* Two arrays of the same size */
    const bool shouldRealloc        = (newSplineDataSize > pmeGpu->archSpecific->splineDataSize);
    int        currentSizeTemp      = pmeGpu->archSpecific->splineDataSize;
    int        currentSizeTempAlloc = pmeGpu->archSpecific->splineDataSizeAlloc;
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_theta, newSplineDataSize,
                           &currentSizeTemp, &currentSizeTempAlloc, pmeGpu->archSpecific->context);
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta, newSplineDataSize,
                           &pmeGpu->archSpecific->splineDataSize,
                           &pmeGpu->archSpecific->splineDataSizeAlloc, pmeGpu->archSpecific->context);
    // the host side reallocation
    if (shouldRealloc)
    {
        pfree(pmeGpu->staging.h_theta);
        pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_theta), newSplineDataSize * sizeof(float));
        pfree(pmeGpu->staging.h_dtheta);
        pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_dtheta), newSplineDataSize * sizeof(float));
    }
}

void pme_gpu_free_spline_data(const PmeGpu* pmeGpu)
{
    /* Two arrays of the same size */
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_theta);
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta);
    pfree(pmeGpu->staging.h_theta);
    pfree(pmeGpu->staging.h_dtheta);
}

void pme_gpu_realloc_grid_indices(PmeGpu* pmeGpu)
{
    const size_t newIndicesSize = DIM * pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newIndicesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices, newIndicesSize,
                           &pmeGpu->archSpecific->gridlineIndicesSize,
                           &pmeGpu->archSpecific->gridlineIndicesSizeAlloc, pmeGpu->archSpecific->context);
    pfree(pmeGpu->staging.h_gridlineIndices);
    pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_gridlineIndices), newIndicesSize * sizeof(int));
}

void pme_gpu_free_grid_indices(const PmeGpu* pmeGpu)
{
    freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices);
    pfree(pmeGpu->staging.h_gridlineIndices);
}

void pme_gpu_realloc_grids(PmeGpu* pmeGpu)
{
    auto*     kernelParamsPtr = pmeGpu->kernelParams.get();
    const int newRealGridSize = kernelParamsPtr->grid.realGridSizePadded[XX]
                                * kernelParamsPtr->grid.realGridSizePadded[YY]
                                * kernelParamsPtr->grid.realGridSizePadded[ZZ];
    const int newComplexGridSize = kernelParamsPtr->grid.complexGridSizePadded[XX]
                                   * kernelParamsPtr->grid.complexGridSizePadded[YY]
                                   * kernelParamsPtr->grid.complexGridSizePadded[ZZ] * 2;
    // Multiplied by 2 because we count complex grid size for complex numbers, but all allocations/pointers are float
    if (pmeGpu->archSpecific->performOutOfPlaceFFT)
    {
        /* 2 separate grids */
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_fourierGrid, newComplexGridSize,
                               &pmeGpu->archSpecific->complexGridSize,
                               &pmeGpu->archSpecific->complexGridSizeAlloc, pmeGpu->archSpecific->context);
        reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid, newRealGridSize,
                               &pmeGpu->archSpecific->realGridSize,
                               &pmeGpu->archSpecific->realGridSizeAlloc, pmeGpu->archSpecific->context);
    }
    else
    {
        /* A single buffer so that any grid will fit */
        const int newGridsSize = std::max(newRealGridSize, newComplexGridSize);
        reallocateDeviceBuffer(
                &kernelParamsPtr->grid.d_realGrid, newGridsSize, &pmeGpu->archSpecific->realGridSize,
                &pmeGpu->archSpecific->realGridSizeAlloc, pmeGpu->archSpecific->context);
        kernelParamsPtr->grid.d_fourierGrid   = kernelParamsPtr->grid.d_realGrid;
        pmeGpu->archSpecific->complexGridSize = pmeGpu->archSpecific->realGridSize;
        // the size might get used later for copying the grid
    }
}

void pme_gpu_free_grids(const PmeGpu* pmeGpu)
{
    if (pmeGpu->archSpecific->performOutOfPlaceFFT)
    {
        freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_fourierGrid);
    }
    freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid);
}

void pme_gpu_clear_grids(const PmeGpu* pmeGpu)
{
    clearDeviceBufferAsync(&pmeGpu->kernelParams->grid.d_realGrid, 0,
                           pmeGpu->archSpecific->realGridSize, pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_realloc_and_copy_fract_shifts(PmeGpu* pmeGpu)
{
    pme_gpu_free_fract_shifts(pmeGpu);

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    const int nx                  = kernelParamsPtr->grid.realGridSize[XX];
    const int ny                  = kernelParamsPtr->grid.realGridSize[YY];
    const int nz                  = kernelParamsPtr->grid.realGridSize[ZZ];
    const int cellCount           = c_pmeNeighborUnitcellCount;
    const int gridDataOffset[DIM] = { 0, cellCount * nx, cellCount * (nx + ny) };

    memcpy(kernelParamsPtr->grid.tablesOffsets, &gridDataOffset, sizeof(gridDataOffset));

    const int newFractShiftsSize = cellCount * (nx + ny + nz);

#if GMX_GPU == GMX_GPU_CUDA
    initParamLookupTable(kernelParamsPtr->grid.d_fractShiftsTable, kernelParamsPtr->fractShiftsTableTexture,
                         pmeGpu->common->fsh.data(), newFractShiftsSize);

    initParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                         kernelParamsPtr->gridlineIndicesTableTexture, pmeGpu->common->nn.data(),
                         newFractShiftsSize);
#elif GMX_GPU == GMX_GPU_OPENCL
    // No dedicated texture routines....
    allocateDeviceBuffer(&kernelParamsPtr->grid.d_fractShiftsTable, newFractShiftsSize,
                         pmeGpu->archSpecific->context);
    allocateDeviceBuffer(&kernelParamsPtr->grid.d_gridlineIndicesTable, newFractShiftsSize,
                         pmeGpu->archSpecific->context);
    copyToDeviceBuffer(&kernelParamsPtr->grid.d_fractShiftsTable, pmeGpu->common->fsh.data(), 0,
                       newFractShiftsSize, pmeGpu->archSpecific->pmeStream,
                       GpuApiCallBehavior::Async, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->grid.d_gridlineIndicesTable, pmeGpu->common->nn.data(), 0,
                       newFractShiftsSize, pmeGpu->archSpecific->pmeStream,
                       GpuApiCallBehavior::Async, nullptr);
#endif
}

void pme_gpu_free_fract_shifts(const PmeGpu* pmeGpu)
{
    auto* kernelParamsPtr = pmeGpu->kernelParams.get();
#if GMX_GPU == GMX_GPU_CUDA
    destroyParamLookupTable(kernelParamsPtr->grid.d_fractShiftsTable,
                            kernelParamsPtr->fractShiftsTableTexture);
    destroyParamLookupTable(kernelParamsPtr->grid.d_gridlineIndicesTable,
                            kernelParamsPtr->gridlineIndicesTableTexture);
#elif GMX_GPU == GMX_GPU_OPENCL
    freeDeviceBuffer(&kernelParamsPtr->grid.d_fractShiftsTable);
    freeDeviceBuffer(&kernelParamsPtr->grid.d_gridlineIndicesTable);
#endif
}

bool pme_gpu_stream_query(const PmeGpu* pmeGpu)
{
    return haveStreamTasksCompleted(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_input_gather_grid(const PmeGpu* pmeGpu, float* h_grid)
{
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid, h_grid, 0, pmeGpu->archSpecific->realGridSize,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_output_spread_grid(const PmeGpu* pmeGpu, float* h_grid)
{
    copyFromDeviceBuffer(h_grid, &pmeGpu->kernelParams->grid.d_realGrid, 0,
                         pmeGpu->archSpecific->realGridSize, pmeGpu->archSpecific->pmeStream,
                         pmeGpu->settings.transferKind, nullptr);
    pmeGpu->archSpecific->syncSpreadGridD2H.markEvent(pmeGpu->archSpecific->pmeStream);
}

void pme_gpu_copy_output_spread_atom_data(const PmeGpu* pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesCount    = DIM * nAtomsPadded * pmeGpu->common->pme_order;
    auto*        kernelParamsPtr = pmeGpu->kernelParams.get();
    copyFromDeviceBuffer(pmeGpu->staging.h_dtheta, &kernelParamsPtr->atoms.d_dtheta, 0, splinesCount,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_theta, &kernelParamsPtr->atoms.d_theta, 0, splinesCount,
                         pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_gridlineIndices, &kernelParamsPtr->atoms.d_gridlineIndices,
                         0, kernelParamsPtr->atoms.nAtoms * DIM, pmeGpu->archSpecific->pmeStream,
                         pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_copy_input_gather_atom_data(const PmeGpu* pmeGpu)
{
    const int    alignment       = pme_gpu_get_atoms_per_warp(pmeGpu);
    const size_t nAtomsPadded    = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
    const size_t splinesCount    = DIM * nAtomsPadded * pmeGpu->common->pme_order;
    auto*        kernelParamsPtr = pmeGpu->kernelParams.get();
    if (c_usePadding)
    {
        // TODO: could clear only the padding and not the whole thing, but this is a test-exclusive code anyway
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_gridlineIndices, 0,
                               pmeGpu->nAtomsAlloc * DIM, pmeGpu->archSpecific->pmeStream);
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_dtheta, 0,
                               pmeGpu->nAtomsAlloc * pmeGpu->common->pme_order * DIM,
                               pmeGpu->archSpecific->pmeStream);
        clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_theta, 0,
                               pmeGpu->nAtomsAlloc * pmeGpu->common->pme_order * DIM,
                               pmeGpu->archSpecific->pmeStream);
    }
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_dtheta, pmeGpu->staging.h_dtheta, 0, splinesCount,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_theta, pmeGpu->staging.h_theta, 0, splinesCount,
                       pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_gridlineIndices, pmeGpu->staging.h_gridlineIndices,
                       0, kernelParamsPtr->atoms.nAtoms * DIM, pmeGpu->archSpecific->pmeStream,
                       pmeGpu->settings.transferKind, nullptr);
}

void pme_gpu_sync_spread_grid(const PmeGpu* pmeGpu)
{
    pmeGpu->archSpecific->syncSpreadGridD2H.waitForEvent();
}

void pme_gpu_init_internal(PmeGpu* pmeGpu)
{
#if GMX_GPU == GMX_GPU_CUDA
    // Prepare to use the device that this PME task was assigned earlier.
    // Other entities, such as CUDA timing events, are known to implicitly use the device context.
    CU_RET_ERR(cudaSetDevice(pmeGpu->deviceInfo->id), "Switching to PME CUDA device");
#endif

    /* Allocate the target-specific structures */
    pmeGpu->archSpecific.reset(new PmeGpuSpecific());
    pmeGpu->kernelParams.reset(new PmeGpuKernelParams());

    pmeGpu->archSpecific->performOutOfPlaceFFT = true;
    /* This should give better performance, according to the cuFFT documentation.
     * The performance seems to be the same though.
     * TODO: PME could also try to pick up nice grid sizes (with factors of 2, 3, 5, 7).
     */

    // TODO: this is just a convenient reuse because programHandle_ currently is in charge of creating context
    pmeGpu->archSpecific->context = pmeGpu->programHandle_->impl_->context;

    // timing enabling - TODO put this in gpu_utils (even though generally this is just option handling?) and reuse in NB
    if (GMX_GPU == GMX_GPU_CUDA)
    {
        /* WARNING: CUDA timings are incorrect with multiple streams.
         *          This is the main reason why they are disabled by default.
         */
        // TODO: Consider turning on by default when we can detect nr of streams.
        pmeGpu->archSpecific->useTiming = (getenv("GMX_ENABLE_GPU_TIMING") != nullptr);
    }
    else if (GMX_GPU == GMX_GPU_OPENCL)
    {
        pmeGpu->archSpecific->useTiming = (getenv("GMX_DISABLE_GPU_TIMING") == nullptr);
    }

#if GMX_GPU == GMX_GPU_CUDA
    pmeGpu->maxGridWidthX = pmeGpu->deviceInfo->prop.maxGridSize[0];
#elif GMX_GPU == GMX_GPU_OPENCL
    pmeGpu->maxGridWidthX = INT32_MAX / 2;
    // TODO: is there no really global work size limit in OpenCL?
#endif

    /* Creating a PME GPU stream:
     * - default high priority with CUDA
     * - no priorities implemented yet with OpenCL; see #2532
     */
#if GMX_GPU == GMX_GPU_CUDA
    cudaError_t stat;
    int         highest_priority, lowest_priority;
    stat = cudaDeviceGetStreamPriorityRange(&lowest_priority, &highest_priority);
    CU_RET_ERR(stat, "PME cudaDeviceGetStreamPriorityRange failed");
    stat = cudaStreamCreateWithPriority(&pmeGpu->archSpecific->pmeStream,
                                        cudaStreamDefault, // cudaStreamNonBlocking,
                                        highest_priority);
    CU_RET_ERR(stat, "cudaStreamCreateWithPriority on the PME stream failed");
#elif GMX_GPU == GMX_GPU_OPENCL
    cl_command_queue_properties queueProperties =
            pmeGpu->archSpecific->useTiming ? CL_QUEUE_PROFILING_ENABLE : 0;
    cl_device_id device_id = pmeGpu->deviceInfo->ocl_gpu_id.ocl_device_id;
    cl_int       clError;
    pmeGpu->archSpecific->pmeStream =
            clCreateCommandQueue(pmeGpu->archSpecific->context, device_id, queueProperties, &clError);
    if (clError != CL_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("Failed to create PME command queue"));
    }
#endif
}

void pme_gpu_destroy_specific(const PmeGpu* pmeGpu)
{
#if GMX_GPU == GMX_GPU_CUDA
    /* Destroy the CUDA stream */
    cudaError_t stat = cudaStreamDestroy(pmeGpu->archSpecific->pmeStream);
    CU_RET_ERR(stat, "PME cudaStreamDestroy error");
#elif GMX_GPU == GMX_GPU_OPENCL
    cl_int clError = clReleaseCommandQueue(pmeGpu->archSpecific->pmeStream);
    if (clError != CL_SUCCESS)
    {
        gmx_warning("Failed to destroy PME command queue");
    }
#endif
}

void pme_gpu_reinit_3dfft(const PmeGpu* pmeGpu)
{
    if (pme_gpu_performs_FFT(pmeGpu))
    {
        pmeGpu->archSpecific->fftSetup.resize(0);
        for (int i = 0; i < pmeGpu->common->ngrids; i++)
        {
            pmeGpu->archSpecific->fftSetup.push_back(std::make_unique<GpuParallel3dFft>(pmeGpu));
        }
    }
}

void pme_gpu_destroy_3dfft(const PmeGpu* pmeGpu)
{
    pmeGpu->archSpecific->fftSetup.resize(0);
}

int getSplineParamFullIndex(int order, int splineIndex, int dimIndex, int atomIndex, int atomsPerWarp)
{
    if (order != c_pmeGpuOrder)
    {
        throw order;
    }
    constexpr int fixedOrder = c_pmeGpuOrder;
    GMX_UNUSED_VALUE(fixedOrder);

    const int atomWarpIndex = atomIndex % atomsPerWarp;
    const int warpIndex     = atomIndex / atomsPerWarp;
    int       indexBase, result;
    switch (atomsPerWarp)
    {
        case 1:
            indexBase = getSplineParamIndexBase<fixedOrder, 1>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 1>(indexBase, dimIndex, splineIndex);
            break;

        case 2:
            indexBase = getSplineParamIndexBase<fixedOrder, 2>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 2>(indexBase, dimIndex, splineIndex);
            break;

        case 4:
            indexBase = getSplineParamIndexBase<fixedOrder, 4>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 4>(indexBase, dimIndex, splineIndex);
            break;

        case 8:
            indexBase = getSplineParamIndexBase<fixedOrder, 8>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 8>(indexBase, dimIndex, splineIndex);
            break;

        default:
            GMX_THROW(gmx::NotImplementedError(
                    gmx::formatString("Test function call not unrolled for atomsPerWarp = %d in "
                                      "getSplineParamFullIndex",
                                      atomsPerWarp)));
    }
    return result;
}

void pme_gpu_getEnergyAndVirial(const gmx_pme_t& pme, PmeOutput* output)
{
    const PmeGpu* pmeGpu = pme.gpu;
    for (int j = 0; j < c_virialAndEnergyCount; j++)
    {
        GMX_ASSERT(std::isfinite(pmeGpu->staging.h_virialAndEnergy[j]),
                   "PME GPU produces incorrect energy/virial.");
    }

    unsigned int j                 = 0;
    output->coulombVirial_[XX][XX] = 0.25F * pmeGpu->staging.h_virialAndEnergy[j++];
    output->coulombVirial_[YY][YY] = 0.25F * pmeGpu->staging.h_virialAndEnergy[j++];
    output->coulombVirial_[ZZ][ZZ] = 0.25F * pmeGpu->staging.h_virialAndEnergy[j++];
    output->coulombVirial_[XX][YY] = output->coulombVirial_[YY][XX] =
            0.25F * pmeGpu->staging.h_virialAndEnergy[j++];
    output->coulombVirial_[XX][ZZ] = output->coulombVirial_[ZZ][XX] =
            0.25F * pmeGpu->staging.h_virialAndEnergy[j++];
    output->coulombVirial_[YY][ZZ] = output->coulombVirial_[ZZ][YY] =
            0.25F * pmeGpu->staging.h_virialAndEnergy[j++];
    output->coulombEnergy_ = 0.5F * pmeGpu->staging.h_virialAndEnergy[j++];
}

/*! \brief Sets the force-related members in \p output
 *
 * \param[in]   pmeGpu      PME GPU data structure
 * \param[out]  output      Pointer to PME output data structure
 */
static void pme_gpu_getForceOutput(PmeGpu* pmeGpu, PmeOutput* output)
{
    output->haveForceOutput_ = !pmeGpu->settings.useGpuForceReduction;
    if (output->haveForceOutput_)
    {
        output->forces_ = pmeGpu->staging.h_forces;
    }
}

PmeOutput pme_gpu_getOutput(const gmx_pme_t& pme, const int flags)
{
    PmeGpu*    pmeGpu                      = pme.gpu;
    const bool haveComputedEnergyAndVirial = (flags & GMX_PME_CALC_ENER_VIR) != 0;

    PmeOutput output;

    pme_gpu_getForceOutput(pmeGpu, &output);

    // The caller knows from the flags that the energy and the virial are not usable
    // on the else branch
    if (haveComputedEnergyAndVirial)
    {
        if (pme_gpu_performs_solve(pmeGpu))
        {
            pme_gpu_getEnergyAndVirial(pme, &output);
        }
        else
        {
            get_pme_ener_vir_q(pme.solve_work, pme.nthread, &output);
        }
    }
    return output;
}

void pme_gpu_update_input_box(PmeGpu gmx_unused* pmeGpu, const matrix gmx_unused box)
{
#if GMX_DOUBLE
    GMX_THROW(gmx::NotImplementedError("PME is implemented for single-precision only on GPU"));
#else
    matrix scaledBox;
    pmeGpu->common->boxScaler->scaleBox(box, scaledBox);
    auto* kernelParamsPtr              = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->current.boxVolume = scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ];
    GMX_ASSERT(kernelParamsPtr->current.boxVolume != 0.0F, "Zero volume of the unit cell");
    matrix recipBox;
    gmx::invertBoxMatrix(scaledBox, recipBox);

    /* The GPU recipBox is transposed as compared to the CPU recipBox.
     * Spread uses matrix columns (while solve and gather use rows).
     * There is no particular reason for this; it might be further rethought/optimized for better access patterns.
     */
    const real newRecipBox[DIM][DIM] = { { recipBox[XX][XX], recipBox[YY][XX], recipBox[ZZ][XX] },
                                         { 0.0, recipBox[YY][YY], recipBox[ZZ][YY] },
                                         { 0.0, 0.0, recipBox[ZZ][ZZ] } };
    memcpy(kernelParamsPtr->current.recipBox, newRecipBox, sizeof(matrix));
#endif
}

/*! \brief \libinternal
 * (Re-)initializes all the PME GPU data related to the grid size and cut-off.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
static void pme_gpu_reinit_grids(PmeGpu* pmeGpu)
{
    auto* kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->grid.ewaldFactor =
            (M_PI * M_PI) / (pmeGpu->common->ewaldcoeff_q * pmeGpu->common->ewaldcoeff_q);

    /* The grid size variants */
    for (int i = 0; i < DIM; i++)
    {
        kernelParamsPtr->grid.realGridSize[i] = pmeGpu->common->nk[i];
        kernelParamsPtr->grid.realGridSizeFP[i] =
                static_cast<float>(kernelParamsPtr->grid.realGridSize[i]);
        kernelParamsPtr->grid.realGridSizePadded[i] = kernelParamsPtr->grid.realGridSize[i];

        // The complex grid currently uses no padding;
        // if it starts to do so, then another test should be added for that
        kernelParamsPtr->grid.complexGridSize[i]       = kernelParamsPtr->grid.realGridSize[i];
        kernelParamsPtr->grid.complexGridSizePadded[i] = kernelParamsPtr->grid.realGridSize[i];
    }
    /* FFT: n real elements correspond to (n / 2 + 1) complex elements in minor dimension */
    if (!pme_gpu_performs_FFT(pmeGpu))
    {
        // This allows for GPU spreading grid and CPU fftgrid to have the same layout, so that we can copy the data directly
        kernelParamsPtr->grid.realGridSizePadded[ZZ] =
                (kernelParamsPtr->grid.realGridSize[ZZ] / 2 + 1) * 2;
    }

    /* GPU FFT: n real elements correspond to (n / 2 + 1) complex elements in minor dimension */
    kernelParamsPtr->grid.complexGridSize[ZZ] /= 2;
    kernelParamsPtr->grid.complexGridSize[ZZ]++;
    kernelParamsPtr->grid.complexGridSizePadded[ZZ] = kernelParamsPtr->grid.complexGridSize[ZZ];

    pme_gpu_realloc_and_copy_fract_shifts(pmeGpu);
    pme_gpu_realloc_and_copy_bspline_values(pmeGpu);
    pme_gpu_realloc_grids(pmeGpu);
    pme_gpu_reinit_3dfft(pmeGpu);
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
static void pme_gpu_copy_common_data_from(const gmx_pme_t* pme)
{
    /* TODO: Consider refactoring the CPU PME code to use the same structure,
     * so that this function becomes 2 lines */
    PmeGpu* pmeGpu               = pme->gpu;
    pmeGpu->common->ngrids       = pme->ngrids;
    pmeGpu->common->epsilon_r    = pme->epsilon_r;
    pmeGpu->common->ewaldcoeff_q = pme->ewaldcoeff_q;
    pmeGpu->common->nk[XX]       = pme->nkx;
    pmeGpu->common->nk[YY]       = pme->nky;
    pmeGpu->common->nk[ZZ]       = pme->nkz;
    pmeGpu->common->pme_order    = pme->pme_order;
    if (pmeGpu->common->pme_order != c_pmeGpuOrder)
    {
        GMX_THROW(gmx::NotImplementedError("pme_order != 4 is not implemented!"));
    }
    for (int i = 0; i < DIM; i++)
    {
        pmeGpu->common->bsp_mod[i].assign(pme->bsp_mod[i], pme->bsp_mod[i] + pmeGpu->common->nk[i]);
    }
    const int cellCount = c_pmeNeighborUnitcellCount;
    pmeGpu->common->fsh.resize(0);
    pmeGpu->common->fsh.insert(pmeGpu->common->fsh.end(), pme->fshx, pme->fshx + cellCount * pme->nkx);
    pmeGpu->common->fsh.insert(pmeGpu->common->fsh.end(), pme->fshy, pme->fshy + cellCount * pme->nky);
    pmeGpu->common->fsh.insert(pmeGpu->common->fsh.end(), pme->fshz, pme->fshz + cellCount * pme->nkz);
    pmeGpu->common->nn.resize(0);
    pmeGpu->common->nn.insert(pmeGpu->common->nn.end(), pme->nnx, pme->nnx + cellCount * pme->nkx);
    pmeGpu->common->nn.insert(pmeGpu->common->nn.end(), pme->nny, pme->nny + cellCount * pme->nky);
    pmeGpu->common->nn.insert(pmeGpu->common->nn.end(), pme->nnz, pme->nnz + cellCount * pme->nkz);
    pmeGpu->common->runMode       = pme->runMode;
    pmeGpu->common->isRankPmeOnly = !pme->bPPnode;
    pmeGpu->common->boxScaler     = pme->boxScaler;
}

/*! \libinternal \brief
 * uses heuristics to select the best performing PME gather and scatter kernels
 *
 * \param[in,out] pmeGpu         The PME GPU structure.
 */
static void pme_gpu_select_best_performing_pme_spreadgather_kernels(PmeGpu* pmeGpu)
{
    if (pmeGpu->kernelParams->atoms.nAtoms > c_pmeGpuPerformanceAtomLimit && (GMX_GPU == GMX_GPU_CUDA))
    {
        pmeGpu->settings.useOrderThreadsPerAtom = true;
        pmeGpu->settings.recalculateSplines     = true;
    }
    else
    {
        pmeGpu->settings.useOrderThreadsPerAtom = false;
        pmeGpu->settings.recalculateSplines     = false;
    }
}


/*! \libinternal \brief
 * Initializes the PME GPU data at the beginning of the run.
 * TODO: this should become PmeGpu::PmeGpu()
 *
 * \param[in,out] pme            The PME structure.
 * \param[in,out] gpuInfo        The GPU information structure.
 * \param[in]     pmeGpuProgram  The handle to the program/kernel data created outside (e.g. in unit tests/runner)
 */
static void pme_gpu_init(gmx_pme_t* pme, const gmx_device_info_t* gpuInfo, PmeGpuProgramHandle pmeGpuProgram)
{
    pme->gpu       = new PmeGpu();
    PmeGpu* pmeGpu = pme->gpu;
    changePinningPolicy(&pmeGpu->staging.h_forces, pme_get_pinning_policy());
    pmeGpu->common = std::make_shared<PmeShared>();

    /* These settings are set here for the whole run; dynamic ones are set in pme_gpu_reinit() */
    /* A convenience variable. */
    pmeGpu->settings.useDecomposition = (pme->nnodes == 1);
    /* TODO: CPU gather with GPU spread is broken due to different theta/dtheta layout. */
    pmeGpu->settings.performGPUGather = true;
    // By default GPU-side reduction is off (explicitly set here for tests, otherwise reset per-step)
    pmeGpu->settings.useGpuForceReduction = false;

    pme_gpu_set_testing(pmeGpu, false);

    pmeGpu->deviceInfo = gpuInfo;
    GMX_ASSERT(pmeGpuProgram != nullptr, "GPU kernels must be already compiled");
    pmeGpu->programHandle_ = pmeGpuProgram;

    pmeGpu->initializedClfftLibrary_ = std::make_unique<gmx::ClfftInitializer>();

    pme_gpu_init_internal(pmeGpu);
    pme_gpu_alloc_energy_virial(pmeGpu);

    pme_gpu_copy_common_data_from(pme);

    GMX_ASSERT(pmeGpu->common->epsilon_r != 0.0F, "PME GPU: bad electrostatic coefficient");

    auto* kernelParamsPtr               = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->constants.elFactor = ONE_4PI_EPS0 / pmeGpu->common->epsilon_r;
}

void pme_gpu_transform_spline_atom_data(const PmeGpu*      pmeGpu,
                                        const PmeAtomComm* atc,
                                        PmeSplineDataType  type,
                                        int                dimIndex,
                                        PmeLayoutTransform transform)
{
    // The GPU atom spline data is laid out in a different way currently than the CPU one.
    // This function converts the data from GPU to CPU layout (in the host memory).
    // It is only intended for testing purposes so far.
    // Ideally we should use similar layouts on CPU and GPU if we care about mixed modes and their
    // performance (e.g. spreading on GPU, gathering on CPU).
    GMX_RELEASE_ASSERT(atc->nthread == 1, "Only the serial PME data layout is supported");
    const uintmax_t threadIndex  = 0;
    const auto      atomCount    = pme_gpu_get_kernel_params_base_ptr(pmeGpu)->atoms.nAtoms;
    const auto      atomsPerWarp = pme_gpu_get_atoms_per_warp(pmeGpu);
    const auto      pmeOrder     = pmeGpu->common->pme_order;
    GMX_ASSERT(pmeOrder == c_pmeGpuOrder, "Only PME order 4 is implemented");

    real*  cpuSplineBuffer;
    float* h_splineBuffer;
    switch (type)
    {
        case PmeSplineDataType::Values:
            cpuSplineBuffer = atc->spline[threadIndex].theta.coefficients[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_theta;
            break;

        case PmeSplineDataType::Derivatives:
            cpuSplineBuffer = atc->spline[threadIndex].dtheta.coefficients[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_dtheta;
            break;

        default: GMX_THROW(gmx::InternalError("Unknown spline data type"));
    }

    for (auto atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        for (auto orderIndex = 0; orderIndex < pmeOrder; orderIndex++)
        {
            const auto gpuValueIndex =
                    getSplineParamFullIndex(pmeOrder, orderIndex, dimIndex, atomIndex, atomsPerWarp);
            const auto cpuValueIndex = atomIndex * pmeOrder + orderIndex;
            GMX_ASSERT(cpuValueIndex < atomCount * pmeOrder,
                       "Atom spline data index out of bounds (while transforming GPU data layout "
                       "for host)");
            switch (transform)
            {
                case PmeLayoutTransform::GpuToHost:
                    cpuSplineBuffer[cpuValueIndex] = h_splineBuffer[gpuValueIndex];
                    break;

                case PmeLayoutTransform::HostToGpu:
                    h_splineBuffer[gpuValueIndex] = cpuSplineBuffer[cpuValueIndex];
                    break;

                default: GMX_THROW(gmx::InternalError("Unknown layout transform"));
            }
        }
    }
}

void pme_gpu_get_real_grid_sizes(const PmeGpu* pmeGpu, gmx::IVec* gridSize, gmx::IVec* paddedGridSize)
{
    GMX_ASSERT(gridSize != nullptr, "");
    GMX_ASSERT(paddedGridSize != nullptr, "");
    GMX_ASSERT(pmeGpu != nullptr, "");
    auto* kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    for (int i = 0; i < DIM; i++)
    {
        (*gridSize)[i]       = kernelParamsPtr->grid.realGridSize[i];
        (*paddedGridSize)[i] = kernelParamsPtr->grid.realGridSizePadded[i];
    }
}

void pme_gpu_reinit(gmx_pme_t* pme, const gmx_device_info_t* gpuInfo, PmeGpuProgramHandle pmeGpuProgram)
{
    if (!pme_gpu_active(pme))
    {
        return;
    }

    if (!pme->gpu)
    {
        /* First-time initialization */
        pme_gpu_init(pme, gpuInfo, pmeGpuProgram);
    }
    else
    {
        /* After this call nothing in the GPU code should refer to the gmx_pme_t *pme itself - until the next pme_gpu_reinit */
        pme_gpu_copy_common_data_from(pme);
    }
    /* GPU FFT will only get used for a single rank.*/
    pme->gpu->settings.performGPUFFT =
            (pme->gpu->common->runMode == PmeRunMode::GPU) && !pme_gpu_uses_dd(pme->gpu);
    pme->gpu->settings.performGPUSolve = (pme->gpu->common->runMode == PmeRunMode::GPU);

    /* Reinit active timers */
    pme_gpu_reinit_timings(pme->gpu);

    pme_gpu_reinit_grids(pme->gpu);
    // Note: if timing the reinit launch overhead becomes more relevant
    // (e.g. with regulat PP-PME re-balancing), we should pass wcycle here.
    pme_gpu_reinit_computation(pme, nullptr);
    /* Clear the previous box - doesn't hurt, and forces the PME CPU recipbox
     * update for mixed mode on grid switch. TODO: use shared recipbox field.
     */
    std::memset(pme->gpu->common->previousBox, 0, sizeof(pme->gpu->common->previousBox));
}

void pme_gpu_destroy(PmeGpu* pmeGpu)
{
    /* Free lots of data */
    pme_gpu_free_energy_virial(pmeGpu);
    pme_gpu_free_bspline_values(pmeGpu);
    pme_gpu_free_forces(pmeGpu);
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

void pme_gpu_reinit_atoms(PmeGpu* pmeGpu, const int nAtoms, const real* charges)
{
    auto* kernelParamsPtr         = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->atoms.nAtoms = nAtoms;
    const int alignment           = pme_gpu_get_atom_data_alignment(pmeGpu);
    pmeGpu->nAtomsPadded          = ((nAtoms + alignment - 1) / alignment) * alignment;
    const int  nAtomsAlloc        = c_usePadding ? pmeGpu->nAtomsPadded : nAtoms;
    const bool haveToRealloc =
            (pmeGpu->nAtomsAlloc < nAtomsAlloc); /* This check might be redundant, but is logical */
    pmeGpu->nAtomsAlloc = nAtomsAlloc;

#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision supported");
    GMX_UNUSED_VALUE(charges);
#else
    pme_gpu_realloc_and_copy_input_coefficients(pmeGpu, reinterpret_cast<const float*>(charges));
    /* Could also be checked for haveToRealloc, but the copy always needs to be performed */
#endif

    if (haveToRealloc)
    {
        pme_gpu_realloc_forces(pmeGpu);
        pme_gpu_realloc_spline_data(pmeGpu);
        pme_gpu_realloc_grid_indices(pmeGpu);
    }
    pme_gpu_select_best_performing_pme_spreadgather_kernels(pmeGpu);
}

void pme_gpu_3dfft(const PmeGpu* pmeGpu, gmx_fft_direction dir, int grid_index)
{
    int timerId = (dir == GMX_FFT_REAL_TO_COMPLEX) ? gtPME_FFT_R2C : gtPME_FFT_C2R;

    pme_gpu_start_timing(pmeGpu, timerId);
    pmeGpu->archSpecific->fftSetup[grid_index]->perform3dFft(
            dir, pme_gpu_fetch_timing_event(pmeGpu, timerId));
    pme_gpu_stop_timing(pmeGpu, timerId);
}

/*! \brief
 * Given possibly large \p blockCount, returns a compact 1D or 2D grid for kernel scheduling,
 * to minimize number of unused blocks.
 */
std::pair<int, int> inline pmeGpuCreateGrid(const PmeGpu* pmeGpu, int blockCount)
{
    // How many maximum widths in X do we need (hopefully just one)
    const int minRowCount = (blockCount + pmeGpu->maxGridWidthX - 1) / pmeGpu->maxGridWidthX;
    // Trying to make things even
    const int colCount = (blockCount + minRowCount - 1) / minRowCount;
    GMX_ASSERT((colCount * minRowCount - blockCount) >= 0, "pmeGpuCreateGrid: totally wrong");
    GMX_ASSERT((colCount * minRowCount - blockCount) < minRowCount,
               "pmeGpuCreateGrid: excessive blocks");
    return std::pair<int, int>(colCount, minRowCount);
}

/*! \brief
 * Returns a pointer to appropriate spline and spread kernel based on the input bool values
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  useOrderThreadsPerAtom   bool controlling if we should use order or order*order threads per atom
 * \param[in]  writeSplinesToGlobal     bool controlling if we should write spline data to global memory
 *
 * \return Pointer to CUDA kernel
 */
static auto selectSplineAndSpreadKernelPtr(const PmeGpu* pmeGpu, bool useOrderThreadsPerAtom, bool writeSplinesToGlobal)
{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (writeSplinesToGlobal)
    {
        if (useOrderThreadsPerAtom)
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelWriteSplinesThPerAtom4;
        }
        else
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelWriteSplines;
        }
    }
    else
    {
        if (useOrderThreadsPerAtom)
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelThPerAtom4;
        }
        else
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernel;
        }
    }

    return kernelPtr;
}

/*! \brief
 * Returns a pointer to appropriate spline kernel based on the input bool values
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  useOrderThreadsPerAtom   bool controlling if we should use order or order*order threads per atom
 * \param[in]  writeSplinesToGlobal     bool controlling if we should write spline data to global memory
 *
 * \return Pointer to CUDA kernel
 */
static auto selectSplineKernelPtr(const PmeGpu* pmeGpu, bool useOrderThreadsPerAtom, bool gmx_unused writeSplinesToGlobal)
{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    GMX_ASSERT(
            writeSplinesToGlobal,
            "Spline data should always be written to global memory when just calculating splines");

    if (useOrderThreadsPerAtom)
    {
        kernelPtr = pmeGpu->programHandle_->impl_->splineKernelThPerAtom4;
    }
    else
    {
        kernelPtr = pmeGpu->programHandle_->impl_->splineKernel;
    }
    return kernelPtr;
}

/*! \brief
 * Returns a pointer to appropriate spread kernel based on the input bool values
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  useOrderThreadsPerAtom   bool controlling if we should use order or order*order threads per atom
 * \param[in]  writeSplinesToGlobal     bool controlling if we should write spline data to global memory
 *
 * \return Pointer to CUDA kernel
 */
static auto selectSpreadKernelPtr(const PmeGpu* pmeGpu, bool useOrderThreadsPerAtom, bool writeSplinesToGlobal)
{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (writeSplinesToGlobal)
    {
        if (useOrderThreadsPerAtom)
        {
            kernelPtr = pmeGpu->programHandle_->impl_->spreadKernelThPerAtom4;
        }
        else
        {
            kernelPtr = pmeGpu->programHandle_->impl_->spreadKernel;
        }
    }
    else
    {
        /* if we are not saving the spline data we need to recalculate it
           using the spline and spread Kernel */
        if (useOrderThreadsPerAtom)
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelThPerAtom4;
        }
        else
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernel;
        }
    }
    return kernelPtr;
}

void pme_gpu_spread(const PmeGpu*         pmeGpu,
                    GpuEventSynchronizer* xReadyOnDevice,
                    int gmx_unused gridIndex,
                    real*          h_grid,
                    bool           computeSplines,
                    bool           spreadCharges)
{
    GMX_ASSERT(computeSplines || spreadCharges,
               "PME spline/spread kernel has invalid input (nothing to do)");
    const auto* kernelParamsPtr = pmeGpu->kernelParams.get();
    GMX_ASSERT(kernelParamsPtr->atoms.nAtoms > 0, "No atom data in PME GPU spread");

    const size_t blockSize = pmeGpu->programHandle_->impl_->spreadWorkGroupSize;

    const int order = pmeGpu->common->pme_order;
    GMX_ASSERT(order == c_pmeGpuOrder, "Only PME order 4 is implemented");
    const bool writeGlobal            = pmeGpu->settings.copyAllOutputs;
    const bool useOrderThreadsPerAtom = pmeGpu->settings.useOrderThreadsPerAtom;
    const bool recalculateSplines     = pmeGpu->settings.recalculateSplines;
#if GMX_GPU == GMX_GPU_OPENCL
    GMX_ASSERT(!useOrderThreadsPerAtom, "Only 16 threads per atom supported in OpenCL");
    GMX_ASSERT(!recalculateSplines, "Recalculating splines not supported in OpenCL");
#endif
    const int atomsPerBlock = useOrderThreadsPerAtom ? blockSize / c_pmeSpreadGatherThreadsPerAtom4ThPerAtom
                                                     : blockSize / c_pmeSpreadGatherThreadsPerAtom;

    // TODO: pick smaller block size in runtime if needed
    // (e.g. on 660 Ti where 50% occupancy is ~25% faster than 100% occupancy with RNAse (~17.8k atoms))
    // If doing so, change atomsPerBlock in the kernels as well.
    // TODO: test varying block sizes on modern arch-s as well
    // TODO: also consider using cudaFuncSetCacheConfig() for preferring shared memory on older architectures
    //(for spline data mostly)
    GMX_ASSERT(!c_usePadding || !(c_pmeAtomDataAlignment % atomsPerBlock),
               "inconsistent atom data padding vs. spreading block size");

    // Ensure that coordinates are ready on the device before launching spread;
    // only needed with CUDA on PP+PME ranks, not on separate PME ranks, in unit tests
    // nor in OpenCL as these cases use a single stream (hence xReadyOnDevice == nullptr).
    GMX_ASSERT(xReadyOnDevice != nullptr || (GMX_GPU != GMX_GPU_CUDA)
                       || pmeGpu->common->isRankPmeOnly || pme_gpu_is_testing(pmeGpu),
               "Need a valid coordinate synchronizer on PP+PME ranks with CUDA.");
    if (xReadyOnDevice)
    {
        xReadyOnDevice->enqueueWaitEvent(pmeGpu->archSpecific->pmeStream);
    }

    const int blockCount = pmeGpu->nAtomsPadded / atomsPerBlock;
    auto      dimGrid    = pmeGpuCreateGrid(pmeGpu, blockCount);

    KernelLaunchConfig config;
    config.blockSize[0] = order;
    config.blockSize[1] = useOrderThreadsPerAtom ? 1 : order;
    config.blockSize[2] = atomsPerBlock;
    config.gridSize[0]  = dimGrid.first;
    config.gridSize[1]  = dimGrid.second;
    config.stream       = pmeGpu->archSpecific->pmeStream;

    int                                timingId;
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (computeSplines)
    {
        if (spreadCharges)
        {
            timingId  = gtPME_SPLINEANDSPREAD;
            kernelPtr = selectSplineAndSpreadKernelPtr(pmeGpu, useOrderThreadsPerAtom,
                                                       writeGlobal || (!recalculateSplines));
        }
        else
        {
            timingId  = gtPME_SPLINE;
            kernelPtr = selectSplineKernelPtr(pmeGpu, useOrderThreadsPerAtom,
                                              writeGlobal || (!recalculateSplines));
        }
    }
    else
    {
        timingId  = gtPME_SPREAD;
        kernelPtr = selectSpreadKernelPtr(pmeGpu, useOrderThreadsPerAtom,
                                          writeGlobal || (!recalculateSplines));
    }


    pme_gpu_start_timing(pmeGpu, timingId);
    auto* timingEvent = pme_gpu_fetch_timing_event(pmeGpu, timingId);
#if c_canEmbedBuffers
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
    const auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, config, kernelParamsPtr, &kernelParamsPtr->atoms.d_theta,
            &kernelParamsPtr->atoms.d_dtheta, &kernelParamsPtr->atoms.d_gridlineIndices,
            &kernelParamsPtr->grid.d_realGrid, &kernelParamsPtr->grid.d_fractShiftsTable,
            &kernelParamsPtr->grid.d_gridlineIndicesTable, &kernelParamsPtr->atoms.d_coefficients,
            &kernelParamsPtr->atoms.d_coordinates);
#endif

    launchGpuKernel(kernelPtr, config, timingEvent, "PME spline/spread", kernelArgs);
    pme_gpu_stop_timing(pmeGpu, timingId);

    const bool copyBackGrid =
            spreadCharges && (pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu));
    if (copyBackGrid)
    {
        pme_gpu_copy_output_spread_grid(pmeGpu, h_grid);
    }
    const bool copyBackAtomData =
            computeSplines && (pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_gather(pmeGpu));
    if (copyBackAtomData)
    {
        pme_gpu_copy_output_spread_atom_data(pmeGpu);
    }
}

void pme_gpu_solve(const PmeGpu* pmeGpu, t_complex* h_grid, GridOrdering gridOrdering, bool computeEnergyAndVirial)
{
    const bool copyInputAndOutputGrid = pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu);

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    float* h_gridFloat = reinterpret_cast<float*>(h_grid);
    if (copyInputAndOutputGrid)
    {
        copyToDeviceBuffer(&kernelParamsPtr->grid.d_fourierGrid, h_gridFloat, 0,
                           pmeGpu->archSpecific->complexGridSize, pmeGpu->archSpecific->pmeStream,
                           pmeGpu->settings.transferKind, nullptr);
    }

    int majorDim = -1, middleDim = -1, minorDim = -1;
    switch (gridOrdering)
    {
        case GridOrdering::YZX:
            majorDim  = YY;
            middleDim = ZZ;
            minorDim  = XX;
            break;

        case GridOrdering::XYZ:
            majorDim  = XX;
            middleDim = YY;
            minorDim  = ZZ;
            break;

        default: GMX_ASSERT(false, "Implement grid ordering here and below for the kernel launch");
    }

    const int maxBlockSize = pmeGpu->programHandle_->impl_->solveMaxWorkGroupSize;

    const int gridLineSize      = pmeGpu->kernelParams->grid.complexGridSize[minorDim];
    const int gridLinesPerBlock = std::max(maxBlockSize / gridLineSize, 1);
    const int blocksPerGridLine = (gridLineSize + maxBlockSize - 1) / maxBlockSize;
    int       cellsPerBlock;
    if (blocksPerGridLine == 1)
    {
        cellsPerBlock = gridLineSize * gridLinesPerBlock;
    }
    else
    {
        cellsPerBlock = (gridLineSize + blocksPerGridLine - 1) / blocksPerGridLine;
    }
    const int warpSize  = pmeGpu->programHandle_->impl_->warpSize;
    const int blockSize = (cellsPerBlock + warpSize - 1) / warpSize * warpSize;

    static_assert(GMX_GPU != GMX_GPU_CUDA || c_solveMaxWarpsPerBlock / 2 >= 4,
                  "The CUDA solve energy kernels needs at least 4 warps. "
                  "Here we launch at least half of the max warps.");

    KernelLaunchConfig config;
    config.blockSize[0] = blockSize;
    config.gridSize[0]  = blocksPerGridLine;
    // rounding up to full warps so that shuffle operations produce defined results
    config.gridSize[1] = (pmeGpu->kernelParams->grid.complexGridSize[middleDim] + gridLinesPerBlock - 1)
                         / gridLinesPerBlock;
    config.gridSize[2] = pmeGpu->kernelParams->grid.complexGridSize[majorDim];
    config.stream      = pmeGpu->archSpecific->pmeStream;

    int                                timingId  = gtPME_SOLVE;
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (gridOrdering == GridOrdering::YZX)
    {
        kernelPtr = computeEnergyAndVirial ? pmeGpu->programHandle_->impl_->solveYZXEnergyKernel
                                           : pmeGpu->programHandle_->impl_->solveYZXKernel;
    }
    else if (gridOrdering == GridOrdering::XYZ)
    {
        kernelPtr = computeEnergyAndVirial ? pmeGpu->programHandle_->impl_->solveXYZEnergyKernel
                                           : pmeGpu->programHandle_->impl_->solveXYZKernel;
    }

    pme_gpu_start_timing(pmeGpu, timingId);
    auto* timingEvent = pme_gpu_fetch_timing_event(pmeGpu, timingId);
#if c_canEmbedBuffers
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
    const auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, config, kernelParamsPtr, &kernelParamsPtr->grid.d_splineModuli,
            &kernelParamsPtr->constants.d_virialAndEnergy, &kernelParamsPtr->grid.d_fourierGrid);
#endif
    launchGpuKernel(kernelPtr, config, timingEvent, "PME solve", kernelArgs);
    pme_gpu_stop_timing(pmeGpu, timingId);

    if (computeEnergyAndVirial)
    {
        copyFromDeviceBuffer(pmeGpu->staging.h_virialAndEnergy,
                             &kernelParamsPtr->constants.d_virialAndEnergy, 0, c_virialAndEnergyCount,
                             pmeGpu->archSpecific->pmeStream, pmeGpu->settings.transferKind, nullptr);
    }

    if (copyInputAndOutputGrid)
    {
        copyFromDeviceBuffer(h_gridFloat, &kernelParamsPtr->grid.d_fourierGrid, 0,
                             pmeGpu->archSpecific->complexGridSize, pmeGpu->archSpecific->pmeStream,
                             pmeGpu->settings.transferKind, nullptr);
    }
}

/*! \brief
 * Returns a pointer to appropriate gather kernel based on the inputvalues
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  useOrderThreadsPerAtom   bool controlling if we should use order or order*order threads per atom
 * \param[in]  readSplinesFromGlobal    bool controlling if we should write spline data to global memory
 * \param[in]  forceTreatment           Controls if the forces from the gather should increment or replace the input forces.
 *
 * \return Pointer to CUDA kernel
 */
inline auto selectGatherKernelPtr(const PmeGpu*          pmeGpu,
                                  bool                   useOrderThreadsPerAtom,
                                  bool                   readSplinesFromGlobal,
                                  PmeForceOutputHandling forceTreatment)

{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;

    if (readSplinesFromGlobal)
    {
        if (useOrderThreadsPerAtom)
        {
            kernelPtr = (forceTreatment == PmeForceOutputHandling::Set)
                                ? pmeGpu->programHandle_->impl_->gatherKernelReadSplinesThPerAtom4
                                : pmeGpu->programHandle_->impl_->gatherReduceWithInputKernelReadSplinesThPerAtom4;
        }
        else
        {
            kernelPtr = (forceTreatment == PmeForceOutputHandling::Set)
                                ? pmeGpu->programHandle_->impl_->gatherKernelReadSplines
                                : pmeGpu->programHandle_->impl_->gatherReduceWithInputKernelReadSplines;
        }
    }
    else
    {
        if (useOrderThreadsPerAtom)
        {
            kernelPtr = (forceTreatment == PmeForceOutputHandling::Set)
                                ? pmeGpu->programHandle_->impl_->gatherKernelThPerAtom4
                                : pmeGpu->programHandle_->impl_->gatherReduceWithInputKernelThPerAtom4;
        }
        else
        {
            kernelPtr = (forceTreatment == PmeForceOutputHandling::Set)
                                ? pmeGpu->programHandle_->impl_->gatherKernel
                                : pmeGpu->programHandle_->impl_->gatherReduceWithInputKernel;
        }
    }
    return kernelPtr;
}


void pme_gpu_gather(PmeGpu* pmeGpu, PmeForceOutputHandling forceTreatment, const float* h_grid)
{
    /* Copying the input CPU forces for reduction */
    if (forceTreatment != PmeForceOutputHandling::Set)
    {
        pme_gpu_copy_input_forces(pmeGpu);
    }

    if (!pme_gpu_performs_FFT(pmeGpu) || pme_gpu_is_testing(pmeGpu))
    {
        pme_gpu_copy_input_gather_grid(pmeGpu, const_cast<float*>(h_grid));
    }

    if (pme_gpu_is_testing(pmeGpu))
    {
        pme_gpu_copy_input_gather_atom_data(pmeGpu);
    }

    /* Set if we have unit tests */
    const bool   readGlobal             = pmeGpu->settings.copyAllOutputs;
    const size_t blockSize              = pmeGpu->programHandle_->impl_->gatherWorkGroupSize;
    const bool   useOrderThreadsPerAtom = pmeGpu->settings.useOrderThreadsPerAtom;
    const bool   recalculateSplines     = pmeGpu->settings.recalculateSplines;
#if GMX_GPU == GMX_GPU_OPENCL
    GMX_ASSERT(!useOrderThreadsPerAtom, "Only 16 threads per atom supported in OpenCL");
    GMX_ASSERT(!recalculateSplines, "Recalculating splines not supported in OpenCL");
#endif
    const int atomsPerBlock = useOrderThreadsPerAtom ? blockSize / c_pmeSpreadGatherThreadsPerAtom4ThPerAtom
                                                     : blockSize / c_pmeSpreadGatherThreadsPerAtom;

    GMX_ASSERT(!c_usePadding || !(c_pmeAtomDataAlignment % atomsPerBlock),
               "inconsistent atom data padding vs. gathering block size");

    const int blockCount = pmeGpu->nAtomsPadded / atomsPerBlock;
    auto      dimGrid    = pmeGpuCreateGrid(pmeGpu, blockCount);

    const int order = pmeGpu->common->pme_order;
    GMX_ASSERT(order == c_pmeGpuOrder, "Only PME order 4 is implemented");

    KernelLaunchConfig config;
    config.blockSize[0] = order;
    config.blockSize[1] = useOrderThreadsPerAtom ? 1 : order;
    config.blockSize[2] = atomsPerBlock;
    config.gridSize[0]  = dimGrid.first;
    config.gridSize[1]  = dimGrid.second;
    config.stream       = pmeGpu->archSpecific->pmeStream;

    // TODO test different cache configs

    int                                timingId  = gtPME_GATHER;
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = selectGatherKernelPtr(
            pmeGpu, useOrderThreadsPerAtom, readGlobal || (!recalculateSplines), forceTreatment);
    // TODO design kernel selection getters and make PmeGpu a friend of PmeGpuProgramImpl

    pme_gpu_start_timing(pmeGpu, timingId);
    auto*       timingEvent     = pme_gpu_fetch_timing_event(pmeGpu, timingId);
    const auto* kernelParamsPtr = pmeGpu->kernelParams.get();
#if c_canEmbedBuffers
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
    const auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, config, kernelParamsPtr, &kernelParamsPtr->atoms.d_coefficients,
            &kernelParamsPtr->grid.d_realGrid, &kernelParamsPtr->atoms.d_theta,
            &kernelParamsPtr->atoms.d_dtheta, &kernelParamsPtr->atoms.d_gridlineIndices,
            &kernelParamsPtr->atoms.d_forces);
#endif
    launchGpuKernel(kernelPtr, config, timingEvent, "PME gather", kernelArgs);
    pme_gpu_stop_timing(pmeGpu, timingId);

    if (pmeGpu->settings.useGpuForceReduction)
    {
        pmeGpu->archSpecific->pmeForcesReady.markEvent(pmeGpu->archSpecific->pmeStream);
    }
    else
    {
        pme_gpu_copy_output_forces(pmeGpu);
    }
}

DeviceBuffer<float> pme_gpu_get_kernelparam_coordinates(const PmeGpu* pmeGpu)
{
    GMX_ASSERT(pmeGpu && pmeGpu->kernelParams,
               "PME GPU device buffer was requested in non-GPU build or before the GPU PME was "
               "initialized.");

    return pmeGpu->kernelParams->atoms.d_coordinates;
}

void* pme_gpu_get_kernelparam_forces(const PmeGpu* pmeGpu)
{
    if (pmeGpu && pmeGpu->kernelParams)
    {
        return pmeGpu->kernelParams->atoms.d_forces;
    }
    else
    {
        return nullptr;
    }
}

/*! \brief Check the validity of the device buffer.
 *
 * Checks if the buffer is not nullptr and, when possible, if it is big enough.
 *
 * \todo Split and move this function to gpu_utils.
 *
 * \param[in] buffer        Device buffer to be checked.
 * \param[in] requiredSize  Number of elements that the buffer will have to accommodate.
 *
 * \returns If the device buffer can be set.
 */
template<typename T>
static bool checkDeviceBuffer(gmx_unused DeviceBuffer<T> buffer, gmx_unused int requiredSize)
{
#if GMX_GPU == GMX_GPU_CUDA
    GMX_ASSERT(buffer != nullptr, "The device pointer is nullptr");
    return buffer != nullptr;
#elif GMX_GPU == GMX_GPU_OPENCL
    size_t size;
    int    retval = clGetMemObjectInfo(buffer, CL_MEM_SIZE, sizeof(size), &size, nullptr);
    GMX_ASSERT(retval == CL_SUCCESS,
               gmx::formatString("clGetMemObjectInfo failed with error code #%d", retval).c_str());
    GMX_ASSERT(static_cast<int>(size) >= requiredSize,
               "Number of atoms in device buffer is smaller then required size.");
    return retval == CL_SUCCESS && static_cast<int>(size) >= requiredSize;
#elif GMX_GPU == GMX_GPU_NONE
    GMX_ASSERT(false, "Setter for device-side coordinates was called in non-GPU build.");
    return false;
#endif
}

void pme_gpu_set_kernelparam_coordinates(const PmeGpu* pmeGpu, DeviceBuffer<float> d_x)
{
    GMX_ASSERT(pmeGpu && pmeGpu->kernelParams,
               "PME GPU device buffer can not be set in non-GPU builds or before the GPU PME was "
               "initialized.");

    GMX_ASSERT(checkDeviceBuffer(d_x, pmeGpu->kernelParams->atoms.nAtoms),
               "The device-side buffer can not be set.");

    pmeGpu->kernelParams->atoms.d_coordinates = d_x;
}

void* pme_gpu_get_stream(const PmeGpu* pmeGpu)
{
    if (pmeGpu)
    {
        return static_cast<void*>(&pmeGpu->archSpecific->pmeStream);
    }
    else
    {
        return nullptr;
    }
}

void* pme_gpu_get_context(const PmeGpu* pmeGpu)
{
    if (pmeGpu)
    {
        return static_cast<void*>(&pmeGpu->archSpecific->context);
    }
    else
    {
        return nullptr;
    }
}

GpuEventSynchronizer* pme_gpu_get_forces_ready_synchronizer(const PmeGpu* pmeGpu)
{
    if (pmeGpu && pmeGpu->kernelParams)
    {
        return &pmeGpu->archSpecific->pmeForcesReady;
    }
    else
    {
        return nullptr;
    }
}
