/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
#include "gromacs/fft/gpu_3dfft.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/pmalloc.h"
#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/syclutils.h"
#endif
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_coordinate_receiver_gpu.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#if GMX_GPU_CUDA
#    include "pme.cuh"
#endif

#include "pme_gpu_calculate_splines.h"
#include "pme_gpu_constants.h"
#include "pme_gpu_grid.h"
#include "pme_gpu_program_impl.h"
#include "pme_gpu_timings.h"
#include "pme_gpu_types.h"
#include "pme_gpu_types_host.h"
#include "pme_gpu_types_host_impl.h"
#include "pme_grid.h"
#include "pme_internal.h"
#include "pme_solve.h"

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

/*! \brief
 * Atom data block size (in terms of number of atoms).
 * This is the least common multiple of number of atoms processed by
 * a single block/workgroup of the spread and gather kernels.
 * The GPU atom data buffers must be padded, which means that
 * the numbers of atoms used for determining the size of the memory
 * allocation must be divisible by this.
 */
#if !GMX_GPU_SYCL
constexpr int c_pmeAtomDataBlockSize = 64;
#else
// Use more padding to support 64-wide warps and ThreadsPerAtom::Order
constexpr int c_pmeAtomDataBlockSize = 128;
#endif

int pme_gpu_get_atom_data_block_size()
{
    return c_pmeAtomDataBlockSize;
}

int pme_gpu_get_atoms_per_warp(const PmeGpu* pmeGpu)
{
    const int order = pmeGpu->common->pme_order;
    const int threadsPerAtom =
            (pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::Order ? order : order * order);
    return pmeGpu->programHandle_->warpSize() / threadsPerAtom;
}

void pme_gpu_synchronize(const PmeGpu* pmeGpu)
{
    pmeGpu->archSpecific->pmeStream_.synchronize();
}

void pme_gpu_alloc_energy_virial(PmeGpu* pmeGpu)
{
    const size_t energyAndVirialSize = c_virialAndEnergyCount * sizeof(float);

    GMX_ASSERT(
            pmeGpu->common->ngrids == 1 || pmeGpu->common->ngrids == 2,
            "Only one (normal Coulomb PME) or two (FEP coulomb PME) PME grids can be used on GPU");

    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        allocateDeviceBuffer(&pmeGpu->kernelParams->constants.d_virialAndEnergy[gridIndex],
                             c_virialAndEnergyCount,
                             pmeGpu->archSpecific->deviceContext_);
        pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_virialAndEnergy[gridIndex]), energyAndVirialSize);
    }
}

void pme_gpu_free_energy_virial(PmeGpu* pmeGpu)
{
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        freeDeviceBuffer(&pmeGpu->kernelParams->constants.d_virialAndEnergy[gridIndex]);
        pfree(pmeGpu->staging.h_virialAndEnergy[gridIndex]);
        pmeGpu->staging.h_virialAndEnergy[gridIndex] = nullptr;
    }
}

void pme_gpu_clear_energy_virial(const PmeGpu* pmeGpu, const bool gpuGraphWithSeparatePmeRank)
{
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        clearDeviceBufferAsync(&pmeGpu->kernelParams->constants.d_virialAndEnergy[gridIndex],
                               0,
                               c_virialAndEnergyCount,
                               pmeGpu->archSpecific->pmeStream_);
    }
    if (pmeGpu->settings.useGpuForceReduction && gpuGraphWithSeparatePmeRank)
    {
        // Mark forces ready event after this clearing, otherwise CUDA graph capture fails due to unjoined work
        pmeGpu->archSpecific->pmeForcesReady.markEvent(pmeGpu->archSpecific->pmeStream_);
    }
}

void pme_gpu_realloc_and_copy_bspline_values(PmeGpu* pmeGpu, const int gridIndex)
{
    GMX_ASSERT(
            pmeGpu->common->ngrids == 1 || pmeGpu->common->ngrids == 2,
            "Only one (normal Coulomb PME) or two (FEP coulomb PME) PME grids can be used on GPU");
    GMX_ASSERT(gridIndex < pmeGpu->common->ngrids,
               "Invalid combination of gridIndex and number of grids");

    const int splineValuesOffset[DIM] = { 0,
                                          pmeGpu->kernelParams->grid.realGridSize[XX],
                                          pmeGpu->kernelParams->grid.realGridSize[XX]
                                                  + pmeGpu->kernelParams->grid.realGridSize[YY] };
    memcpy(&pmeGpu->kernelParams->grid.splineValuesOffset, &splineValuesOffset, sizeof(splineValuesOffset));

    const int newSplineValuesSize = pmeGpu->kernelParams->grid.realGridSize[XX]
                                    + pmeGpu->kernelParams->grid.realGridSize[YY]
                                    + pmeGpu->kernelParams->grid.realGridSize[ZZ];
    const bool shouldRealloc = (newSplineValuesSize > pmeGpu->archSpecific->splineValuesSize[gridIndex]);
    reallocateDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli[gridIndex],
                           newSplineValuesSize,
                           &pmeGpu->archSpecific->splineValuesSize[gridIndex],
                           &pmeGpu->archSpecific->splineValuesCapacity[gridIndex],
                           pmeGpu->archSpecific->deviceContext_);
    if (shouldRealloc)
    {
        /* Reallocate the host buffer */
        pfree(pmeGpu->staging.h_splineModuli[gridIndex]);
        pmalloc(reinterpret_cast<void**>(&pmeGpu->staging.h_splineModuli[gridIndex]),
                newSplineValuesSize * sizeof(float));
    }
    for (int i = 0; i < DIM; i++)
    {
        memcpy(pmeGpu->staging.h_splineModuli[gridIndex] + splineValuesOffset[i],
               pmeGpu->common->bsp_mod[i].data(),
               pmeGpu->common->bsp_mod[i].size() * sizeof(float));
    }
    /* TODO: pin original buffer instead! */
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli[gridIndex],
                       pmeGpu->staging.h_splineModuli[gridIndex],
                       0,
                       newSplineValuesSize,
                       pmeGpu->archSpecific->pmeStream_,
                       pmeGpu->settings.transferKind,
                       nullptr);
}

void pme_gpu_free_bspline_values(const PmeGpu* pmeGpu)
{
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        pfree(pmeGpu->staging.h_splineModuli[gridIndex]);
        freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_splineModuli[gridIndex]);
    }
}

void pme_gpu_realloc_forces(PmeGpu* pmeGpu)
{
    const size_t newForcesSize =
            pmeGpu->useNvshmem ? pmeGpu->nvshmemParams->nAtomsAlloc_symmetric : pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newForcesSize > 0, "Bad number of atoms in PME GPU");
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces,
                           newForcesSize,
                           &pmeGpu->archSpecific->forcesSize,
                           &pmeGpu->archSpecific->forcesSizeAlloc,
                           pmeGpu->archSpecific->deviceContext_,
                           pmeGpu->useNvshmem);

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
    copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_forces,
                       pmeGpu->staging.h_forces.data(),
                       0,
                       pmeGpu->kernelParams->atoms.nAtoms,
                       pmeGpu->archSpecific->pmeStream_,
                       pmeGpu->settings.transferKind,
                       nullptr);
}

void pme_gpu_copy_output_forces(PmeGpu* pmeGpu)
{
    GMX_ASSERT(pmeGpu->kernelParams->atoms.nAtoms > 0, "Bad number of atoms in PME GPU");
    copyFromDeviceBuffer(pmeGpu->staging.h_forces.data(),
                         &pmeGpu->kernelParams->atoms.d_forces,
                         0,
                         pmeGpu->kernelParams->atoms.nAtoms,
                         pmeGpu->archSpecific->pmeStream_,
                         pmeGpu->settings.transferKind,
                         nullptr);
}

void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu* pmeGpu,
                                                 const float*  h_coefficients,
                                                 const int     gridIndex)
{
    const size_t newCoefficientsSize = pmeGpu->nAtomsAlloc;
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients[gridIndex],
                           newCoefficientsSize,
                           &pmeGpu->archSpecific->coefficientsSize[gridIndex],
                           &pmeGpu->archSpecific->coefficientsCapacity[gridIndex],
                           pmeGpu->archSpecific->deviceContext_);
    if (newCoefficientsSize != 0)
    {
        GMX_ASSERT(h_coefficients, "Bad host-side charge buffer in PME GPU");
        GMX_ASSERT(newCoefficientsSize > 0, "Bad number of atoms in PME GPU");
        copyToDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients[gridIndex],
                           h_coefficients,
                           0,
                           pmeGpu->kernelParams->atoms.nAtoms,
                           pmeGpu->archSpecific->pmeStream_,
                           pmeGpu->settings.transferKind,
                           nullptr);

        const size_t paddingIndex = pmeGpu->kernelParams->atoms.nAtoms;
        const size_t paddingCount = pmeGpu->nAtomsAlloc - paddingIndex;
        if (paddingCount > 0)
        {
            clearDeviceBufferAsync(&pmeGpu->kernelParams->atoms.d_coefficients[gridIndex],
                                   paddingIndex,
                                   paddingCount,
                                   pmeGpu->archSpecific->pmeStream_);
        }
    }
}

void pme_gpu_free_coefficients(const PmeGpu* pmeGpu)
{
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        freeDeviceBuffer(&pmeGpu->kernelParams->atoms.d_coefficients[gridIndex]);
    }
}

void pme_gpu_realloc_spline_data(PmeGpu* pmeGpu)
{
    const int order             = pmeGpu->common->pme_order;
    const int newSplineDataSize = DIM * order * pmeGpu->nAtomsAlloc;
    GMX_ASSERT(newSplineDataSize > 0, "Bad number of atoms in PME GPU");
    /* Two arrays of the same size */
    const bool shouldRealloc        = (newSplineDataSize > pmeGpu->archSpecific->splineDataSize);
    int        currentSizeTemp      = pmeGpu->archSpecific->splineDataSize;
    int        currentSizeTempAlloc = pmeGpu->archSpecific->splineDataSizeAlloc;
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_theta,
                           newSplineDataSize,
                           &currentSizeTemp,
                           &currentSizeTempAlloc,
                           pmeGpu->archSpecific->deviceContext_);
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_dtheta,
                           newSplineDataSize,
                           &pmeGpu->archSpecific->splineDataSize,
                           &pmeGpu->archSpecific->splineDataSizeAlloc,
                           pmeGpu->archSpecific->deviceContext_);

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
    reallocateDeviceBuffer(&pmeGpu->kernelParams->atoms.d_gridlineIndices,
                           newIndicesSize,
                           &pmeGpu->archSpecific->gridlineIndicesSize,
                           &pmeGpu->archSpecific->gridlineIndicesSizeAlloc,
                           pmeGpu->archSpecific->deviceContext_);
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
    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    const int newRealGridSize = kernelParamsPtr->grid.realGridSizePadded[XX]
                                * kernelParamsPtr->grid.realGridSizePadded[YY]
                                * kernelParamsPtr->grid.realGridSizePadded[ZZ];
    const int newComplexGridSize = kernelParamsPtr->grid.complexGridSizePadded[XX]
                                   * kernelParamsPtr->grid.complexGridSizePadded[YY]
                                   * kernelParamsPtr->grid.complexGridSizePadded[ZZ] * 2;
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        // Multiplied by 2 because we count complex grid size for complex numbers, but all allocations/pointers are float
        if (pmeGpu->archSpecific->performOutOfPlaceFFT)
        {
            /* 2 separate grids */
            reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid[gridIndex],
                                   newRealGridSize,
                                   &pmeGpu->archSpecific->realGridSize[gridIndex],
                                   &pmeGpu->archSpecific->realGridCapacity[gridIndex],
                                   pmeGpu->archSpecific->deviceContext_);

            pmeGpu->archSpecific->complexGridSize[gridIndex] = newComplexGridSize;
        }
        else
        {
            /* A single buffer so that any grid will fit */
            const int newGridsSize = std::max(newRealGridSize, newComplexGridSize);
            reallocateDeviceBuffer(&kernelParamsPtr->grid.d_realGrid[gridIndex],
                                   newGridsSize,
                                   &pmeGpu->archSpecific->realGridSize[gridIndex],
                                   &pmeGpu->archSpecific->realGridCapacity[gridIndex],
                                   pmeGpu->archSpecific->deviceContext_);
            pmeGpu->archSpecific->complexGridSize[gridIndex] =
                    pmeGpu->archSpecific->realGridSize[gridIndex];
            // the size might get used later for copying the grid
        }
    }
}

void pme_gpu_free_grids(const PmeGpu* pmeGpu)
{
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        freeDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex]);
    }
}

void pme_gpu_reinit_haloexchange(PmeGpu* pmeGpu)
{
    // allocate overlap buffers needed for PME grid halo exchanges
    int overlapX     = 0;
    int overlapY     = 0;
    int overlapDown  = 0;
    int overlapRight = 0;

    // calculate neighboring rank ids in X-dimension
    int rankX = pmeGpu->common->nodeidX;
    int sizeX = pmeGpu->common->nnodesX;
    int down  = (rankX + 1) % sizeX;
    int up    = (rankX + sizeX - 1) % sizeX;

    // calculate neighboring rank ids in Y-dimension
    int rankY = pmeGpu->common->nodeidY;
    int sizeY = pmeGpu->common->nnodesY;
    int right = (rankY + 1) % sizeY;
    int left  = (rankY + sizeY - 1) % sizeY;

    // local grid size in X and Y dimension
    int myGridX = pmeGpu->common->s2g0X[rankX + 1] - pmeGpu->common->s2g0X[rankX];
    int myGridY = pmeGpu->common->s2g0Y[rankY + 1] - pmeGpu->common->s2g0Y[rankY];

    // populate haloExchange structure
    pmeGpu->haloExchange->gridSizeX = myGridX;
    pmeGpu->haloExchange->gridSizeY = myGridY;

    // halo sizes for all neighbours.
    // This can be different from pmeGpu->common->gridHalo
    // when there are only 2 ranks along a dimension
    pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Up]     = 0;
    pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Down]   = 0;
    pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Center] = 0;

    pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Left]   = 0;
    pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Right]  = 0;
    pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Center] = 0;

    // populate rank ids of neighbours along X- and Y-dimension
    pmeGpu->haloExchange->ranksX[gmx::DirectionX::Up]     = up;
    pmeGpu->haloExchange->ranksX[gmx::DirectionX::Down]   = down;
    pmeGpu->haloExchange->ranksX[gmx::DirectionX::Center] = rankX;

    pmeGpu->haloExchange->ranksY[gmx::DirectionY::Left]   = left;
    pmeGpu->haloExchange->ranksY[gmx::DirectionY::Right]  = right;
    pmeGpu->haloExchange->ranksY[gmx::DirectionY::Center] = rankY;

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    if (pmeGpu->common->nnodesX > 1)
    {
        int downGrid = pmeGpu->common->s2g0X[down + 1] - pmeGpu->common->s2g0X[down];
        int upGrid   = pmeGpu->common->s2g0X[up + 1] - pmeGpu->common->s2g0X[up];

        overlapX = pmeGpu->common->gridHalo;

        // current implementation transfers from/to only immediate neighbours
        GMX_RELEASE_ASSERT(overlapX <= myGridX && overlapX <= downGrid && overlapX <= upGrid,
                           "Exchange supported only with immediate neighbor");

        overlapDown   = overlapX;
        int overlapUp = overlapX;

        // if only 2 PME ranks in X-domain and overlap width more than slab width
        // just transfer all grid points from neighbor
        if (down == up && 2 * overlapX >= downGrid)
        {
            overlapX    = myGridX;
            overlapDown = downGrid;
            overlapUp   = 0;
        }

        pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Up]     = overlapUp;
        pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Down]   = overlapDown;
        pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Center] = overlapX;

        int overlapSize = std::max(overlapX, overlapDown) * kernelParamsPtr->grid.realGridSizePadded[YY]
                          * kernelParamsPtr->grid.realGridSizePadded[ZZ];

        for (std::size_t i = 0; i < std::size_t(gmx::DirectionX::Center); ++i)
        {
            reallocateDeviceBuffer(&pmeGpu->haloExchange->d_sendGrids[i][gmx::DirectionY::Center],
                                   overlapSize,
                                   &pmeGpu->haloExchange->overlapSendSize[i][gmx::DirectionY::Center],
                                   &pmeGpu->haloExchange->overlapSendCapacity[i][gmx::DirectionY::Center],
                                   pmeGpu->archSpecific->deviceContext_);

            reallocateDeviceBuffer(&pmeGpu->haloExchange->d_recvGrids[i][gmx::DirectionY::Center],
                                   overlapSize,
                                   &pmeGpu->haloExchange->overlapRecvSize[i][gmx::DirectionY::Center],
                                   &pmeGpu->haloExchange->overlapRecvCapacity[i][gmx::DirectionY::Center],
                                   pmeGpu->archSpecific->deviceContext_);
        }
    }

    if (pmeGpu->common->nnodesY > 1)
    {
        int rightGrid = pmeGpu->common->s2g0Y[right + 1] - pmeGpu->common->s2g0Y[right];
        int leftGrid  = pmeGpu->common->s2g0Y[left + 1] - pmeGpu->common->s2g0Y[left];

        overlapY = pmeGpu->common->gridHalo;

        // current implementation transfers from/to only immediate neighbours
        GMX_RELEASE_ASSERT(overlapY <= myGridY && overlapY <= rightGrid && overlapY <= leftGrid,
                           "Exchange supported only with immediate neighbor");

        overlapRight    = overlapY;
        int overlapLeft = overlapY;

        // if only 2 PME ranks in Y-domain and overlap width more than slab width
        // just transfer all grid points from neighbor
        if (right == left && 2 * overlapY >= rightGrid)
        {
            overlapY     = myGridY;
            overlapRight = rightGrid;
            overlapLeft  = 0;
        }

        pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Left]   = overlapLeft;
        pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Right]  = overlapRight;
        pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Center] = overlapY;

        int overlapSize = std::max(overlapY, overlapRight) * kernelParamsPtr->grid.realGridSizePadded[XX]
                          * kernelParamsPtr->grid.realGridSizePadded[ZZ];

        for (std::size_t i = 0; i < std::size_t(gmx::DirectionY::Center); ++i)
        {
            reallocateDeviceBuffer(&pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][i],
                                   overlapSize,
                                   &pmeGpu->haloExchange->overlapSendSize[gmx::DirectionX::Center][i],
                                   &pmeGpu->haloExchange->overlapSendCapacity[gmx::DirectionX::Center][i],
                                   pmeGpu->archSpecific->deviceContext_);

            reallocateDeviceBuffer(&pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][i],
                                   overlapSize,
                                   &pmeGpu->haloExchange->overlapRecvSize[gmx::DirectionX::Center][i],
                                   &pmeGpu->haloExchange->overlapRecvCapacity[gmx::DirectionX::Center][i],
                                   pmeGpu->archSpecific->deviceContext_);
        }
    }

    if (pmeGpu->common->nnodesX > 1 && pmeGpu->common->nnodesY > 1)
    {
        int overlapSize = std::max(overlapX, overlapDown) * std::max(overlapY, overlapRight)
                          * kernelParamsPtr->grid.realGridSizePadded[ZZ];

        for (std::size_t i = 0; i < std::size_t(gmx::DirectionX::Center); ++i)
        {
            for (std::size_t j = 0; j < std::size_t(gmx::DirectionY::Center); ++j)
            {
                reallocateDeviceBuffer(&pmeGpu->haloExchange->d_sendGrids[i][j],
                                       overlapSize,
                                       &pmeGpu->haloExchange->overlapSendSize[i][j],
                                       &pmeGpu->haloExchange->overlapSendCapacity[i][j],
                                       pmeGpu->archSpecific->deviceContext_);

                reallocateDeviceBuffer(&pmeGpu->haloExchange->d_recvGrids[i][j],
                                       overlapSize,
                                       &pmeGpu->haloExchange->overlapRecvSize[i][j],
                                       &pmeGpu->haloExchange->overlapRecvCapacity[i][j],
                                       pmeGpu->archSpecific->deviceContext_);
            }
        }
    }
}

void pme_gpu_free_haloexchange(const PmeGpu* pmeGpu)
{
    for (std::size_t i = 0; i < std::size_t(gmx::DirectionX::Count); ++i)
    {
        for (std::size_t j = 0; j < std::size_t(gmx::DirectionY::Count); ++j)
        {
            freeDeviceBuffer(&pmeGpu->haloExchange->d_sendGrids[i][j]);
            freeDeviceBuffer(&pmeGpu->haloExchange->d_recvGrids[i][j]);
        }
    }
}

void pme_gpu_clear_grids(const PmeGpu* pmeGpu)
{
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        clearDeviceBufferAsync(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                               0,
                               pmeGpu->archSpecific->realGridSize[gridIndex],
                               pmeGpu->archSpecific->pmeStream_);
    }
}

void pme_gpu_realloc_and_copy_fract_shifts(PmeGpu* pmeGpu)
{
    pme_gpu_free_fract_shifts(pmeGpu);

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    const int nx                  = pmeGpu->common->nk[XX];
    const int ny                  = pmeGpu->common->nk[YY];
    const int nz                  = pmeGpu->common->nk[ZZ];
    const int cellCount           = c_pmeNeighborUnitcellCount;
    const int gridDataOffset[DIM] = { 0, cellCount * nx, cellCount * (nx + ny) };

    memcpy(kernelParamsPtr->grid.tablesOffsets, &gridDataOffset, sizeof(gridDataOffset));

    const int newFractShiftsSize = cellCount * (nx + ny + nz);

    initParamLookupTable(&kernelParamsPtr->grid.d_fractShiftsTable,
                         &kernelParamsPtr->fractShiftsTableTexture,
                         pmeGpu->common->fsh.data(),
                         newFractShiftsSize,
                         pmeGpu->archSpecific->deviceContext_);

    initParamLookupTable(&kernelParamsPtr->grid.d_gridlineIndicesTable,
                         &kernelParamsPtr->gridlineIndicesTableTexture,
                         pmeGpu->common->nn.data(),
                         newFractShiftsSize,
                         pmeGpu->archSpecific->deviceContext_);
}

void pme_gpu_free_fract_shifts(const PmeGpu* pmeGpu)
{
    auto* kernelParamsPtr = pmeGpu->kernelParams.get();
#if GMX_GPU_CUDA
    destroyParamLookupTable(&kernelParamsPtr->grid.d_fractShiftsTable,
                            &kernelParamsPtr->fractShiftsTableTexture);
    destroyParamLookupTable(&kernelParamsPtr->grid.d_gridlineIndicesTable,
                            &kernelParamsPtr->gridlineIndicesTableTexture);
#elif GMX_GPU_OPENCL || GMX_GPU_SYCL
    freeDeviceBuffer(&kernelParamsPtr->grid.d_fractShiftsTable);
    freeDeviceBuffer(&kernelParamsPtr->grid.d_gridlineIndicesTable);
#endif
}

bool pme_gpu_stream_query(const PmeGpu* pmeGpu)
{
    return haveStreamTasksCompleted(pmeGpu->archSpecific->pmeStream_);
}

void pme_gpu_copy_input_gather_grid(const PmeGpu* pmeGpu, const float* h_grid, const int gridIndex)
{
    copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                       h_grid,
                       0,
                       pmeGpu->archSpecific->realGridSize[gridIndex],
                       pmeGpu->archSpecific->pmeStream_,
                       pmeGpu->settings.transferKind,
                       nullptr);
}

void pme_gpu_copy_output_spread_grid(const PmeGpu* pmeGpu, float* h_grid, const int gridIndex)
{
    copyFromDeviceBuffer(h_grid,
                         &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                         0,
                         pmeGpu->archSpecific->realGridSize[gridIndex],
                         pmeGpu->archSpecific->pmeStream_,
                         pmeGpu->settings.transferKind,
                         nullptr);
    pmeGpu->archSpecific->syncSpreadGridD2H.markEvent(pmeGpu->archSpecific->pmeStream_);
}

void pme_gpu_copy_output_spread_atom_data(const PmeGpu* pmeGpu)
{
    auto*     kernelParamsPtr = pmeGpu->kernelParams.get();
    const int nAtoms          = pmeGpu->kernelParams->atoms.nAtoms;

    copyFromDeviceBuffer(pmeGpu->staging.h_dtheta,
                         &kernelParamsPtr->atoms.d_dtheta,
                         0,
                         pmeGpu->archSpecific->splineCountActive,
                         pmeGpu->archSpecific->pmeStream_,
                         pmeGpu->settings.transferKind,
                         nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_theta,
                         &kernelParamsPtr->atoms.d_theta,
                         0,
                         pmeGpu->archSpecific->splineCountActive,
                         pmeGpu->archSpecific->pmeStream_,
                         pmeGpu->settings.transferKind,
                         nullptr);
    copyFromDeviceBuffer(pmeGpu->staging.h_gridlineIndices,
                         &kernelParamsPtr->atoms.d_gridlineIndices,
                         0,
                         nAtoms * DIM,
                         pmeGpu->archSpecific->pmeStream_,
                         pmeGpu->settings.transferKind,
                         nullptr);
}

void pme_gpu_copy_input_gather_atom_data(const PmeGpu* pmeGpu)
{
    const size_t splineDataSize  = pmeGpu->archSpecific->splineDataSize;
    auto*        kernelParamsPtr = pmeGpu->kernelParams.get();

    // TODO: could clear only the padding and not the whole thing, but this is a test-exclusive code anyway
    clearDeviceBufferAsync(&kernelParamsPtr->atoms.d_gridlineIndices,
                           0,
                           pmeGpu->nAtomsAlloc * DIM,
                           pmeGpu->archSpecific->pmeStream_);
    clearDeviceBufferAsync(
            &kernelParamsPtr->atoms.d_dtheta, 0, splineDataSize, pmeGpu->archSpecific->pmeStream_);
    clearDeviceBufferAsync(
            &kernelParamsPtr->atoms.d_theta, 0, splineDataSize, pmeGpu->archSpecific->pmeStream_);

    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_dtheta,
                       pmeGpu->staging.h_dtheta,
                       0,
                       splineDataSize,
                       pmeGpu->archSpecific->pmeStream_,
                       pmeGpu->settings.transferKind,
                       nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_theta,
                       pmeGpu->staging.h_theta,
                       0,
                       splineDataSize,
                       pmeGpu->archSpecific->pmeStream_,
                       pmeGpu->settings.transferKind,
                       nullptr);
    copyToDeviceBuffer(&kernelParamsPtr->atoms.d_gridlineIndices,
                       pmeGpu->staging.h_gridlineIndices,
                       0,
                       kernelParamsPtr->atoms.nAtoms * DIM,
                       pmeGpu->archSpecific->pmeStream_,
                       pmeGpu->settings.transferKind,
                       nullptr);
}

void pme_gpu_sync_spread_grid(const PmeGpu* pmeGpu)
{
    pmeGpu->archSpecific->syncSpreadGridD2H.waitForEvent();
}

/*! \brief Internal GPU initialization for PME.
 *
 * \param[in]  pmeGpu         GPU PME data.
 * \param[in]  deviceContext  GPU context.
 * \param[in]  deviceStream   GPU stream.
 */
static void pme_gpu_init_internal(PmeGpu* pmeGpu, const DeviceContext& deviceContext, const DeviceStream& deviceStream)
{
    /* Allocate the target-specific structures */
    pmeGpu->archSpecific.reset(new PmeGpuSpecific(deviceContext, deviceStream));
    pmeGpu->kernelParams.reset(new PmeGpuKernelParams());

    // Use in-place FFT with cuFFTMp or BBFFT.
    pmeGpu->archSpecific->performOutOfPlaceFFT =
            !((pmeGpu->settings.useDecomposition && GMX_USE_cuFFTMp) || GMX_GPU_FFT_BBFFT);

    /* This should give better performance, according to the cuFFT documentation.
     * The performance seems to be the same though.
     * TODO: PME could also try to pick up nice grid sizes (with factors of 2, 3, 5, 7).
     */

#if GMX_GPU_CUDA || GMX_GPU_SYCL
    pmeGpu->kernelParams->usePipeline       = char(false);
    pmeGpu->kernelParams->pipelineAtomStart = 0;
    pmeGpu->kernelParams->pipelineAtomEnd   = 0;
#endif
#if GMX_GPU_CUDA
    pmeGpu->maxGridWidthX = deviceContext.deviceInfo().prop.maxGridSize[0];
#else
    // Use this path for any non-CUDA GPU acceleration
    // TODO: is there no really global work size limit in OpenCL?
    pmeGpu->maxGridWidthX = INT32_MAX / 2;
#endif

    if (pmeGpu->settings.useDecomposition)
    {
        pmeGpu->haloExchange = std::make_unique<PmeGpuHaloExchange>();
        for (std::size_t i = 0; i < std::size_t(gmx::DirectionX::Count); ++i)
        {
            for (std::size_t j = 0; j < std::size_t(gmx::DirectionY::Count); ++j)
            {
                // initialize capacity so that realloc can be used
                pmeGpu->haloExchange->d_sendGrids[i][j] = nullptr;
                pmeGpu->haloExchange->d_recvGrids[i][j] = nullptr;

                pmeGpu->haloExchange->overlapSendSize[i][j] = 0;
                pmeGpu->haloExchange->overlapRecvSize[i][j] = 0;

                pmeGpu->haloExchange->overlapSendCapacity[i][j] = 0;
                pmeGpu->haloExchange->overlapRecvCapacity[i][j] = 0;
            }
        }
    }
}

/*! \internal \brief
 * Wrapper for getting FFT backend depending on Gromacs GPU backend compilation flags.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 * \returns FftBackend enum value.
 */
static gmx::FftBackend getFftBackend(const PmeGpu* pmeGpu)
{
    if (GMX_GPU_CUDA)
    {
        if (!pmeGpu->settings.useDecomposition)
        {
            return gmx::FftBackend::Cufft;
        }
        else
        {
            if (GMX_USE_cuFFTMp)
            {
                return gmx::FftBackend::CuFFTMp;
            }
            else if (GMX_USE_Heffte)
            {
                return gmx::FftBackend::HeFFTe_CUDA;
            }
            else
            {
                GMX_RELEASE_ASSERT(
                        false,
                        "Gromacs must be built with cuFFTMp or Heffte to enable GPU-based "
                        "PME decomposition on CUDA-compatible GPUs");
                return gmx::FftBackend::Count;
            }
        }
    }
    else if (GMX_GPU_OPENCL)
    {
        if (GMX_GPU_FFT_VKFFT)
        {
            return gmx::FftBackend::OclVkfft;
        }
        else
        {
            GMX_RELEASE_ASSERT(GMX_GPU_FFT_CLFFT, "Only clFFT and VkFFT are supported with OpenCL");
            return gmx::FftBackend::Ocl;
        }
    }
    else if (GMX_GPU_SYCL)
    {
        if (GMX_GPU_FFT_MKL)
        {
            if (!pmeGpu->settings.useDecomposition)
            {
                return gmx::FftBackend::SyclMkl;
            }
            else if (GMX_USE_Heffte)
            {
                return gmx::FftBackend::HeFFTe_Sycl_OneMkl;
            }
            else
            {
                GMX_THROW(gmx::NotImplementedError(
                        "GROMACS must be built with HeFFTe to enable fully GPU-offloaded "
                        "PME decomposition on oneAPI-compatible GPUs"));
            }
        }
        else if (GMX_GPU_FFT_BBFFT)
        {
            return gmx::FftBackend::SyclBbfft;
        }
        else if (GMX_GPU_FFT_ROCFFT)
        {
            if (!pmeGpu->settings.useDecomposition)
            {
                return gmx::FftBackend::SyclRocfft;
            }
            else if (GMX_USE_Heffte)
            {
                return gmx::FftBackend::HeFFTe_Sycl_Rocfft;
            }
            else
            {
                GMX_THROW(gmx::NotImplementedError(
                        "GROMACS must be built with HeFFTe to enable fully GPU-offloaded "
                        "PME decomposition on ROCm-compatible GPUs"));
            }
        }
        else if (GMX_GPU_FFT_CUFFT)
        {
            if (GMX_USE_Heffte)
            {
                if (pmeGpu->settings.useDecomposition)
                {
                    return gmx::FftBackend::HeFFTe_Sycl_cuFFT;
                }
                else
                {
                    GMX_THROW(gmx::NotImplementedError(
                            "GROMACS can only do multi-GPU FFT in SYCL+cuFFT+HeFFTe build"));
                }
            }
            else
            {
                GMX_THROW(gmx::NotImplementedError("GROMACS does not support cuFFT in SYCL build"));
            }
        }
        else if (GMX_GPU_FFT_VKFFT)
        {
            return gmx::FftBackend::SyclVkfft;
        }
        else
        {
            return gmx::FftBackend::Sycl;
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Unknown GPU backend");
        return gmx::FftBackend::Count;
    }
}

void pme_gpu_reinit_3dfft(const PmeGpu* pmeGpu)
{
    if (pme_gpu_settings(pmeGpu).performGPUFFT)
    {
        pmeGpu->archSpecific->fftSetup.resize(0);
        const bool       performOutOfPlaceFFT = pmeGpu->archSpecific->performOutOfPlaceFFT;
        MPI_Comm         comm                 = pmeGpu->common->mpiComm;
        std::vector<int> gridSizesInXForEachRank(pmeGpu->common->nnodesX);
        std::vector<int> gridSizesInYForEachRank(pmeGpu->common->nnodesY);

        for (int i = 0; i < pmeGpu->common->nnodesX; ++i)
        {
            gridSizesInXForEachRank[i] = pmeGpu->common->s2g0X[i + 1] - pmeGpu->common->s2g0X[i];
        }
        for (int i = 0; i < pmeGpu->common->nnodesY; ++i)
        {
            gridSizesInYForEachRank[i] = pmeGpu->common->s2g0Y[i + 1] - pmeGpu->common->s2g0Y[i];
        }

        const bool allocateRealGrid = pmeGpu->settings.useDecomposition;

        const gmx::FftBackend backend = getFftBackend(pmeGpu);

        PmeGpuGridParams& grid = pme_gpu_get_kernel_params_base_ptr(pmeGpu)->grid;
        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            const bool useDecomposition = pmeGpu->settings.useDecomposition;

            // grid needs to be alloacted only with decomposition
            GMX_RELEASE_ASSERT(
                    allocateRealGrid == useDecomposition,
                    "Separate FFT real grid needs to be allocated only with decomposition");

            if (!useDecomposition)
            {
                memcpy(pmeGpu->archSpecific->localRealGridSize, grid.realGridSize, DIM * sizeof(int));
                memcpy(pmeGpu->archSpecific->localRealGridSizePadded,
                       grid.realGridSizePadded,
                       DIM * sizeof(int));
                memcpy(grid.localComplexGridSizePadded, grid.complexGridSizePadded, DIM * sizeof(int));

                // PME grid is same as FFT real grid in case of no decomposition
                pmeGpu->archSpecific->d_fftRealGrid[gridIndex] = grid.d_realGrid[gridIndex];
            }

            pmeGpu->archSpecific->fftSetup.push_back(std::make_unique<gmx::Gpu3dFft>(
                    backend,
                    allocateRealGrid,
                    comm,
                    gridSizesInXForEachRank,
                    gridSizesInYForEachRank,
                    grid.realGridSize[ZZ],
                    performOutOfPlaceFFT,
                    pmeGpu->archSpecific->deviceContext_,
                    pmeGpu->archSpecific->pmeStream_,
                    pmeGpu->archSpecific->localRealGridSize,
                    pmeGpu->archSpecific->localRealGridSizePadded,
                    grid.localComplexGridSizePadded,
                    &(pmeGpu->archSpecific->d_fftRealGrid[gridIndex]),
                    &(grid.d_fftComplexGrid[gridIndex])));

            // no difference in padded and unpadded size
            memcpy(grid.localComplexGridSize, grid.localComplexGridSizePadded, DIM * sizeof(int));
        }
    }
    else
    {
        // Initialize fft complex grid and size.
        // These values needs to be initialized for unit tests which run pme_gpu_solve even in mixed
        // mode. In real world cases, pme_gpu_solve is never called in mixed mode.
        PmeGpuGridParams& grid = pme_gpu_get_kernel_params_base_ptr(pmeGpu)->grid;
        memcpy(grid.localComplexGridSizePadded, grid.complexGridSizePadded, DIM * sizeof(int));
        memcpy(grid.localComplexGridSize, grid.complexGridSize, DIM * sizeof(int));

        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            grid.d_fftComplexGrid[gridIndex] = grid.d_realGrid[gridIndex];
        }
    }
}

void pme_gpu_destroy_3dfft(const PmeGpu* pmeGpu)
{
    pmeGpu->archSpecific->fftSetup.resize(0);
}

void pme_gpu_getEnergyAndVirial(const gmx_pme_t& pme, const float lambda, PmeOutput* output)
{
    const PmeGpu* pmeGpu = pme.gpu;

    GMX_ASSERT(lambda == 1.0 || pmeGpu->common->ngrids == 2,
               "Invalid combination of lambda and number of grids");

    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        for (int j = 0; j < c_virialAndEnergyCount; j++)
        {
            GMX_ASSERT(std::isfinite(pmeGpu->staging.h_virialAndEnergy[gridIndex][j]),
                       "PME GPU produces incorrect energy/virial.");
        }
    }
    for (int dim1 = 0; dim1 < DIM; dim1++)
    {
        for (int dim2 = 0; dim2 < DIM; dim2++)
        {
            output->coulombVirial_[dim1][dim2] = 0;
        }
    }
    output->coulombEnergy_ = 0;
    float scale            = 1.0;
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        if (pmeGpu->common->ngrids == 2)
        {
            scale = gridIndex == 0 ? (1.0 - lambda) : lambda;
        }
        output->coulombVirial_[XX][XX] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][0];
        output->coulombVirial_[YY][YY] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][1];
        output->coulombVirial_[ZZ][ZZ] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][2];
        output->coulombVirial_[XX][YY] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][3];
        output->coulombVirial_[YY][XX] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][3];
        output->coulombVirial_[XX][ZZ] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][4];
        output->coulombVirial_[ZZ][XX] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][4];
        output->coulombVirial_[YY][ZZ] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][5];
        output->coulombVirial_[ZZ][YY] +=
                scale * 0.25F * pmeGpu->staging.h_virialAndEnergy[gridIndex][5];
        output->coulombEnergy_ += scale * 0.5F * pmeGpu->staging.h_virialAndEnergy[gridIndex][6];
    }
    if (pmeGpu->common->ngrids > 1)
    {
        output->coulombDvdl_ = 0.5F
                               * (pmeGpu->staging.h_virialAndEnergy[FEP_STATE_B][6]
                                  - pmeGpu->staging.h_virialAndEnergy[FEP_STATE_A][6]);
    }
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

PmeOutput pme_gpu_getOutput(const gmx_pme_t& pme, const bool computeEnergyAndVirial, const real lambdaQ)
{
    PmeGpu* pmeGpu = pme.gpu;

    PmeOutput output;

    pme_gpu_getForceOutput(pmeGpu, &output);

    if (computeEnergyAndVirial)
    {
        if (pme_gpu_settings(pmeGpu).performGPUSolve)
        {
            pme_gpu_getEnergyAndVirial(pme, lambdaQ, &output);
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

    GMX_ASSERT(
            pmeGpu->common->ngrids == 1 || pmeGpu->common->ngrids == 2,
            "Only one (normal Coulomb PME) or two (FEP coulomb PME) PME grids can be used on GPU");

    kernelParamsPtr->grid.ewaldFactor =
            (M_PI * M_PI) / (pmeGpu->common->ewaldcoeff_q * pmeGpu->common->ewaldcoeff_q);
    /* The grid size variants */
    for (int i = 0; i < DIM; i++)
    {
        // allocate only decomposed grid size on device
        kernelParamsPtr->grid.realGridSize[i] = pmeGpu->common->nk[i];

        // realGridSizeFP is expected to be full-grid size in kernels
        kernelParamsPtr->grid.realGridSizeFP[i]     = static_cast<float>(pmeGpu->common->nk[i]);
        kernelParamsPtr->grid.realGridSizePadded[i] = kernelParamsPtr->grid.realGridSize[i];

        // The complex grid currently uses no padding;
        // if it starts to do so, then another test should be added for that
        kernelParamsPtr->grid.complexGridSize[i]       = kernelParamsPtr->grid.realGridSize[i];
        kernelParamsPtr->grid.complexGridSizePadded[i] = kernelParamsPtr->grid.realGridSize[i];
    }
    /* FFT: n real elements correspond to (n / 2 + 1) complex elements in minor dimension */
    // NOLINTNEXTLINE(misc-redundant-expression)
    const bool needWorkaroundForOneMkl = (GMX_GPU_SYCL != 0) && (GMX_SYCL_DPCPP != 0)
                                         && pme_gpu_settings(pmeGpu).performGPUFFT; // Issue #4219.
    if (!pme_gpu_settings(pmeGpu).performGPUFFT || needWorkaroundForOneMkl)
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
    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        pme_gpu_realloc_and_copy_bspline_values(pmeGpu, gridIndex);
    }

    pme_gpu_realloc_grids(pmeGpu);
    if (pmeGpu->settings.useDecomposition)
    {
        pme_gpu_reinit_haloexchange(pmeGpu);
    }
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
    pmeGpu->common->ngrids       = pme->bFEP_q ? 2 : 1;
    pmeGpu->common->epsilon_r    = pme->epsilon_r;
    pmeGpu->common->ewaldcoeff_q = pme->ewaldcoeff_q;
    pmeGpu->common->nk[XX]       = pme->nkx;
    pmeGpu->common->nk[YY]       = pme->nky;
    pmeGpu->common->nk[ZZ]       = pme->nkz;
    pmeGpu->common->pme_order    = pme->pme_order;
    // overlap region is needed only in case of decomposition
    pmeGpu->common->pmegridNk[XX] = pme->nnodes_major > 1 ? pme->pmegrid_nx : pme->nkx;
    pmeGpu->common->pmegridNk[YY] = pme->nnodes_minor > 1 ? pme->pmegrid_ny : pme->nky;
    pmeGpu->common->pmegridNk[ZZ] = pme->nkz;
    pmeGpu->common->ndecompdim    = pme->ndecompdim;
    pmeGpu->common->nodeid        = pme->nodeid;
    pmeGpu->common->nodeidX       = pme->nodeid_major;
    pmeGpu->common->nodeidY       = pme->nodeid_minor;
    pmeGpu->common->nnodes        = pme->nnodes;
    pmeGpu->common->nnodesX       = pme->nnodes_major;
    pmeGpu->common->nnodesY       = pme->nnodes_minor;
    pmeGpu->common->s2g0X         = pme->overlap[0].s2g0;
    pmeGpu->common->s2g1X         = pme->overlap[0].s2g1;
    pmeGpu->common->s2g0Y         = pme->overlap[1].s2g0;
    pmeGpu->common->s2g1Y         = pme->overlap[1].s2g1;
    pmeGpu->common->gridHalo      = pme->pmeGpuGridHalo;

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
    pmeGpu->common->boxScaler     = pme->boxScaler.get();
    pmeGpu->common->mpiCommX      = pme->mpi_comm_d[0];
    pmeGpu->common->mpiCommY      = pme->mpi_comm_d[1];
    pmeGpu->common->mpiComm       = pme->mpi_comm;
}

/*! \libinternal \brief
 * uses heuristics to select the best performing PME gather and scatter kernels
 *
 * \param[in,out] pmeGpu         The PME GPU structure.
 */
static void pme_gpu_select_best_performing_pme_spreadgather_kernels(PmeGpu* pmeGpu)
{
    if (((GMX_GPU_CUDA != 0) || (GMX_GPU_SYCL != 0))
        && pmeGpu->kernelParams->atoms.nAtoms > pmeGpu->minParticleCountToRecalculateSplines)
    {
        pmeGpu->settings.threadsPerAtom     = ThreadsPerAtom::Order;
        pmeGpu->settings.recalculateSplines = true;
    }
    else
    {
        pmeGpu->settings.threadsPerAtom     = ThreadsPerAtom::OrderSquared;
        pmeGpu->settings.recalculateSplines = false;
    }
}


/*! \libinternal \brief
 * Initializes the PME GPU data at the beginning of the run.
 * TODO: this should become PmeGpu::PmeGpu()
 *
 * \param[in,out] pme            The PME structure.
 * \param[in]     deviceContext  The GPU context.
 * \param[in]     deviceStream   The GPU stream.
 * \param[in,out] pmeGpuProgram  The handle to the program/kernel data created outside (e.g. in unit tests/runner)
 */
static void pme_gpu_init(gmx_pme_t*           pme,
                         const DeviceContext& deviceContext,
                         const DeviceStream&  deviceStream,
                         const PmeGpuProgram* pmeGpuProgram)
{
    pme->gpu       = new PmeGpu();
    PmeGpu* pmeGpu = pme->gpu;
    changePinningPolicy(&pmeGpu->staging.h_forces, pme_get_pinning_policy());
    pmeGpu->common = std::make_shared<PmeShared>();

    /* These settings are set here for the whole run; dynamic ones are set in pme_gpu_reinit() */
    /* A convenience variable. */
    pmeGpu->settings.useDecomposition = (pme->nnodes != 1);
    /* TODO: CPU gather with GPU spread is broken due to different theta/dtheta layout. */
    pmeGpu->settings.performGPUGather = true;
    // By default GPU-side reduction is off (explicitly set here for tests, otherwise reset per-step)
    pmeGpu->settings.useGpuForceReduction = false;

    pme_gpu_set_testing(pmeGpu, false);

    GMX_ASSERT(pmeGpuProgram != nullptr, "GPU kernels must be already compiled");
    pmeGpu->programHandle_ = pmeGpuProgram;

    pmeGpu->initializedClfftLibrary_ = std::make_unique<gmx::ClfftInitializer>();

    pme_gpu_init_internal(pmeGpu, deviceContext, deviceStream);

    pme_gpu_copy_common_data_from(pme);
    pme_gpu_alloc_energy_virial(pmeGpu);

    GMX_ASSERT(pmeGpu->common->epsilon_r != 0.0F, "PME GPU: bad electrostatic coefficient");

    auto* kernelParamsPtr               = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->constants.elFactor = gmx::c_one4PiEps0 / pmeGpu->common->epsilon_r;
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

void pme_gpu_reinit(gmx_pme_t*           pme,
                    const DeviceContext* deviceContext,
                    const DeviceStream*  deviceStream,
                    const PmeGpuProgram* pmeGpuProgram,
                    const bool           useMdGpuGraph)
{
    GMX_ASSERT(pme != nullptr, "Need valid PME object");

    if (!pme->gpu)
    {
        GMX_RELEASE_ASSERT(deviceContext != nullptr,
                           "Device context can not be nullptr when setting up PME on GPU.");
        GMX_RELEASE_ASSERT(deviceStream != nullptr,
                           "Device stream can not be nullptr when setting up PME on GPU.");
        /* First-time initialization */
        pme_gpu_init(pme, *deviceContext, *deviceStream, pmeGpuProgram);
    }
    else
    {
        /* After this call nothing in the GPU code should refer to the gmx_pme_t *pme itself - until the next pme_gpu_reinit */
        pme_gpu_copy_common_data_from(pme);
    }

    pme->gpu->settings.performGPUFFT   = (pme->gpu->common->runMode == PmeRunMode::GPU);
    pme->gpu->settings.performGPUSolve = (pme->gpu->common->runMode == PmeRunMode::GPU);

    /* Reinit active timers */
    pme_gpu_reinit_timings(pme->gpu);

    pme_gpu_reinit_grids(pme->gpu);
    // Note: if timing the reinit launch overhead becomes more relevant
    // (e.g. with regulat PP-PME re-balancing), we should pass wcycle here.
    pme_gpu_reinit_computation(pme, useMdGpuGraph, nullptr);
    /* Clear the previous box - doesn't hurt, and forces the PME CPU recipbox
     * update for mixed mode on grid switch. TODO: use shared recipbox field.
     */
    std::memset(pme->gpu->common->previousBox, 0, sizeof(pme->gpu->common->previousBox));
}

void pme_gpu_destroy(PmeGpu* pmeGpu)
{
    // Wait for all the tasks to complete before freeing the memory. See #4519.
    pmeGpu->archSpecific->pmeStream_.synchronize();

    /* Free lots of data */
    pme_gpu_free_energy_virial(pmeGpu);
    pme_gpu_free_bspline_values(pmeGpu);
    pme_gpu_free_forces(pmeGpu);
    pme_gpu_free_coefficients(pmeGpu);
    pme_gpu_free_spline_data(pmeGpu);
    pme_gpu_free_grid_indices(pmeGpu);
    pme_gpu_free_fract_shifts(pmeGpu);
    pme_gpu_free_grids(pmeGpu);
    if (pmeGpu->settings.useDecomposition)
    {
        pme_gpu_free_haloexchange(pmeGpu);
    }
    pme_gpu_destroy_3dfft(pmeGpu);

    delete pmeGpu;
}

void pme_gpu_reinit_atoms(PmeGpu* pmeGpu, const int nAtoms, const real* chargesA, const real* chargesB)
{
    auto* kernelParamsPtr         = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->atoms.nAtoms = nAtoms;
    const int  blockSize          = pme_gpu_get_atom_data_block_size();
    const int  nAtomsNewPadded    = ((nAtoms + blockSize - 1) / blockSize) * blockSize;
    const bool haveToRealloc      = (pmeGpu->nAtomsAlloc < nAtomsNewPadded);
    pmeGpu->nAtomsAlloc           = nAtomsNewPadded;

    const auto atomsPerWarp                 = pme_gpu_get_atoms_per_warp(pmeGpu);
    const int  nWarps                       = ((nAtoms + atomsPerWarp - 1) / atomsPerWarp);
    pmeGpu->archSpecific->splineCountActive = DIM * nWarps * atomsPerWarp * pmeGpu->common->pme_order;

    if (pmeGpu->useNvshmem)
    {
        // find the max nAtomsAlloc among all the ranks for symmetric forces buffer allocation.
#if GMX_MPI
        MPI_Allreduce(
                &pmeGpu->nAtomsAlloc, &pmeGpu->nvshmemParams->nAtomsAlloc_symmetric, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

        int myRank     = -1;
        int numPpRanks = pmeGpu->nvshmemParams->ppRanksRef.size();
#if GMX_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Bcast(&numPpRanks, 1, MPI_INT, myRank, MPI_COMM_WORLD);
#endif
        // symmetric buffer allocation used for synchronization purpose
        // 1 to be used to signal PME to PP rank of put, and
        // numPpRanks is intended to be used for each PP rank buffer consumption completion
        // signal to PME to allow to produce it again.
        reallocateDeviceBuffer(&pmeGpu->kernelParams->forcesReadyNvshmemFlags,
                               1 + numPpRanks,
                               &pmeGpu->nvshmemParams->forcesReadyNvshmemFlagsSize,
                               &pmeGpu->nvshmemParams->forcesReadyNvshmemFlagsSizeAlloc,
                               pmeGpu->archSpecific->deviceContext_,
                               true);

        pme_gpu_realloc_forces(pmeGpu);

        /* allocate the ppRanksFInfo host buffer */
        kernelParamsPtr->ppRanksInfoSize = pmeGpu->nvshmemParams->ppRanksRef.size();

        if (pmeGpu->nvshmemParams->ppRanksFInfo == nullptr)
        {
            pmalloc(reinterpret_cast<void**>(&pmeGpu->nvshmemParams->ppRanksFInfo),
                    kernelParamsPtr->ppRanksInfoSize * sizeof(PpRanksSendFInfo));
        }

        // prepare the ppRanksFInfo struct for sending it to gpu.
        int receiverIndex = 0;
        for (const auto& receiver : pmeGpu->nvshmemParams->ppRanksRef)
        {
            auto& ppRankFInfo_prev = pmeGpu->nvshmemParams->ppRanksFInfo[receiverIndex - 1];
            int   startIndex =
                    receiverIndex ? ppRankFInfo_prev.startAtomOffset + ppRankFInfo_prev.numAtoms : 0;

            pmeGpu->nvshmemParams->ppRanksFInfo[receiverIndex] = { receiver.rankId, receiver.numAtoms, startIndex };
            receiverIndex++;
        }

        reallocateDeviceBuffer(&kernelParamsPtr->ppRanksInfo,
                               kernelParamsPtr->ppRanksInfoSize,
                               &pmeGpu->nvshmemParams->ppRanksFInfoSize,
                               &pmeGpu->nvshmemParams->ppRanksFInfoSizeAlloc,
                               pmeGpu->archSpecific->deviceContext_);
        copyToDeviceBuffer(&kernelParamsPtr->ppRanksInfo,
                           pmeGpu->nvshmemParams->ppRanksFInfo,
                           0,
                           kernelParamsPtr->ppRanksInfoSize,
                           pmeGpu->archSpecific->pmeStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);
        reallocateDeviceBuffer(&kernelParamsPtr->lastProcessedBlockPerPpRank,
                               kernelParamsPtr->ppRanksInfoSize,
                               &pmeGpu->nvshmemParams->lastProcessedBlockPerPpRankSize,
                               &pmeGpu->nvshmemParams->lastProcessedBlockPerPpRankSizeAlloc,
                               pmeGpu->archSpecific->deviceContext_);
        clearDeviceBufferAsync(&kernelParamsPtr->lastProcessedBlockPerPpRank,
                               0,
                               kernelParamsPtr->ppRanksInfoSize,
                               pmeGpu->archSpecific->pmeStream_);
    }
    else if (haveToRealloc)
    {
        pme_gpu_realloc_forces(pmeGpu);
    }

#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision supported");
    GMX_UNUSED_VALUE(charges);
#else
    int gridIndex = 0;
    /* Could also be checked for haveToRealloc, but the copy always needs to be performed */
    pme_gpu_realloc_and_copy_input_coefficients(pmeGpu, reinterpret_cast<const float*>(chargesA), gridIndex);
    gridIndex++;
    if (chargesB != nullptr)
    {
        pme_gpu_realloc_and_copy_input_coefficients(
                pmeGpu, reinterpret_cast<const float*>(chargesB), gridIndex);
    }
    else
    {
        /* Fill the second set of coefficients with chargesA as well to be able to avoid
         * conditionals in the GPU kernels */
        /* FIXME: This should be avoided by making a separate templated version of the
         * relevant kernel(s) (probably only pme_gather_kernel). That would require a
         * reduction of the current number of templated parameters of that kernel. */
        if (pmeGpu->common->ngrids > 1)
        {
            pme_gpu_realloc_and_copy_input_coefficients(
                    pmeGpu, reinterpret_cast<const float*>(chargesA), gridIndex);
        }
    }
#endif

    if (haveToRealloc)
    {
        pme_gpu_realloc_spline_data(pmeGpu);
        pme_gpu_realloc_grid_indices(pmeGpu);
    }
    else if (pmeGpu->staging.h_forces.size() != pmeGpu->kernelParams->atoms.nAtoms)
    {
        // re-alloc not needed but resizing is needed if nAtoms changed
        pmeGpu->staging.h_forces.resizeWithPadding(pmeGpu->kernelParams->atoms.nAtoms);
    }
    pme_gpu_select_best_performing_pme_spreadgather_kernels(pmeGpu);
}

/*! \internal \brief
 * Returns raw timing event from the corresponding GpuRegionTimer (if timings are enabled).
 * In CUDA result can be nullptr stub, per GpuRegionTimer implementation.
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \param[in] pmeStageId     The PME GPU stage gtPME_ index from the enum in src/gromacs/timing/gpu_timing.h
 */
static CommandEvent* pme_gpu_fetch_timing_event(const PmeGpu* pmeGpu, PmeStage pmeStageId)
{
    CommandEvent* timingEvent = nullptr;
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        GMX_ASSERT(pmeStageId < PmeStage::Count, "Wrong PME GPU timing event index");
        timingEvent = pmeGpu->archSpecific->timingEvents[pmeStageId].fetchNextEvent();
    }
    return timingEvent;
}

void pme_gpu_3dfft(const PmeGpu* pmeGpu, gmx_fft_direction dir, const int grid_index)
{
    PmeStage timerId = (dir == GMX_FFT_REAL_TO_COMPLEX) ? PmeStage::FftTransformR2C
                                                        : PmeStage::FftTransformC2R;

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
 * \param[in]  threadsPerAtom           Controls whether we should use order or order*order threads per atom
 * \param[in]  writeSplinesToGlobal     bool controlling if we should write spline data to global memory
 * \param[in]  numGrids                 Number of grids to use. numGrids == 2 if Coulomb is perturbed.
 *
 * \return Pointer to CUDA kernel
 */
static auto selectSplineAndSpreadKernelPtr(const PmeGpu*  pmeGpu,
                                           ThreadsPerAtom threadsPerAtom,
                                           bool           writeSplinesToGlobal,
                                           const int      numGrids)
{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (writeSplinesToGlobal)
    {
        if (threadsPerAtom == ThreadsPerAtom::Order)
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelWriteSplinesThPerAtom4Dual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelWriteSplinesThPerAtom4Single;
            }
        }
        else
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelWriteSplinesDual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelWriteSplinesSingle;
            }
        }
    }
    else
    {
        if (threadsPerAtom == ThreadsPerAtom::Order)
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelThPerAtom4Dual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelThPerAtom4Single;
            }
        }
        else
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelDual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelSingle;
            }
        }
    }

    return kernelPtr;
}

/*! \brief
 * Returns a pointer to appropriate spline kernel based on the input bool values
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  threadsPerAtom           Controls whether we should use order or order*order threads per atom
 * \param[in]  writeSplinesToGlobal     bool controlling if we should write spline data to global memory
 * \param[in]  numGrids                 Number of grids to use. numGrids == 2 if Coulomb is perturbed.
 *
 * \return Pointer to CUDA kernel
 */
static auto selectSplineKernelPtr(const PmeGpu*   pmeGpu,
                                  ThreadsPerAtom  threadsPerAtom,
                                  bool gmx_unused writeSplinesToGlobal,
                                  const int       numGrids)
{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    GMX_ASSERT(
            writeSplinesToGlobal,
            "Spline data should always be written to global memory when just calculating splines");

    if (threadsPerAtom == ThreadsPerAtom::Order)
    {
        if (numGrids == 2)
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineKernelThPerAtom4Dual;
        }
        else
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineKernelThPerAtom4Single;
        }
    }
    else
    {
        if (numGrids == 2)
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineKernelDual;
        }
        else
        {
            kernelPtr = pmeGpu->programHandle_->impl_->splineKernelSingle;
        }
    }
    return kernelPtr;
}

/*! \brief
 * Returns a pointer to appropriate spread kernel based on the input bool values
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  threadsPerAtom           Controls whether we should use order or order*order threads per atom
 * \param[in]  writeSplinesToGlobal     bool controlling if we should write spline data to global memory
 * \param[in]  numGrids                 Number of grids to use. numGrids == 2 if Coulomb is perturbed.
 *
 * \return Pointer to CUDA kernel
 */
static auto selectSpreadKernelPtr(const PmeGpu*  pmeGpu,
                                  ThreadsPerAtom threadsPerAtom,
                                  bool           writeSplinesToGlobal,
                                  const int      numGrids)
{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (writeSplinesToGlobal)
    {
        if (threadsPerAtom == ThreadsPerAtom::Order)
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->spreadKernelThPerAtom4Dual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->spreadKernelThPerAtom4Single;
            }
        }
        else
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->spreadKernelDual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->spreadKernelSingle;
            }
        }
    }
    else
    {
        /* if we are not saving the spline data we need to recalculate it
           using the spline and spread Kernel */
        if (threadsPerAtom == ThreadsPerAtom::Order)
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelThPerAtom4Dual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelThPerAtom4Single;
            }
        }
        else
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelDual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->splineAndSpreadKernelSingle;
            }
        }
    }
    return kernelPtr;
}

/*! \brief
 * Manages synchronization with remote GPU's PP coordinate sender, for a stage of the communication operation.
 * For thread-MPI, an event associated with a stage of the operation is enqueued to the GPU stream that will be used by the consumer.
 * For lib-MPI, the CPU task executing this function will wait for a stage to be completed.
 * In each case, the rank of the sender associated with the corresponding stage is returned.
 *
 * \param[in]  pmeGpu                    The PME GPU structure.
 * \param[in]  pmeCoordinateReceiverGpu  The PME coordinate reciever GPU object
 * \param[in]  usePipeline               Whether pipelining is in use for PME-PP communication
 * \param[in]  pipelineStage             Stage of the PME-PP communication pipeline
 *
 * \return Rank of remote sender associated with this stage
 */
static int manageSyncWithPpCoordinateSenderGpu(const PmeGpu*                  pmeGpu,
                                               gmx::PmeCoordinateReceiverGpu* pmeCoordinateReceiverGpu,
                                               bool                           usePipeline   = false,
                                               int                            pipelineStage = 0)
{
    int senderRank;
    if (GMX_THREAD_MPI)
    {
        GpuEventSynchronizer* event;
        std::tie(senderRank, event) =
                pmeCoordinateReceiverGpu->receivePpCoordinateSendEvent(pipelineStage);
        if (usePipeline)
        {
            event->enqueueWaitEvent(*(pmeCoordinateReceiverGpu->ppCommStream(senderRank)));
        }
        else
        {
            event->enqueueWaitEvent(pmeGpu->archSpecific->pmeStream_);
        }
    }
    else
    {
        senderRank = pmeCoordinateReceiverGpu->waitForCoordinatesFromAnyPpRank();
    }
    return senderRank;
}

void pme_gpu_spread(const PmeGpu*                  pmeGpu,
                    GpuEventSynchronizer*          xReadyOnDevice,
                    real**                         h_grids,
                    gmx_parallel_3dfft_t*          fftSetup,
                    bool                           computeSplines,
                    bool                           spreadCharges,
                    const real                     lambda,
                    const bool                     useGpuDirectComm,
                    gmx::PmeCoordinateReceiverGpu* pmeCoordinateReceiverGpu,
                    const bool                     useMdGpuGraph,
                    gmx_wallcycle*                 wcycle)
{
    GMX_ASSERT(
            pmeGpu->common->ngrids == 1 || pmeGpu->common->ngrids == 2,
            "Only one (normal Coulomb PME) or two (FEP coulomb PME) PME grids can be used on GPU");

    GMX_ASSERT(computeSplines || spreadCharges,
               "PME spline/spread kernel has invalid input (nothing to do)");
    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    const size_t blockSize = pmeGpu->programHandle_->impl_->spreadWorkGroupSize;

    const int order = pmeGpu->common->pme_order;
    GMX_ASSERT(order == c_pmeGpuOrder, "Only PME order 4 is implemented");
    const bool writeGlobal = pmeGpu->settings.copyAllOutputs;
    const int  threadsPerAtom =
            (pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::Order ? order : order * order);
    const bool recalculateSplines = pmeGpu->settings.recalculateSplines;

    GMX_ASSERT(!GMX_GPU_OPENCL || pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::OrderSquared,
               "Only 16 threads per atom supported in OpenCL");
    GMX_ASSERT(!GMX_GPU_OPENCL || !recalculateSplines,
               "Recalculating splines not supported in OpenCL");

    const int atomsPerBlock = blockSize / threadsPerAtom;

    // TODO: pick smaller block size in runtime if needed
    // (e.g. on 660 Ti where 50% occupancy is ~25% faster than 100% occupancy with RNAse (~17.8k atoms))
    // If doing so, change atomsPerBlock in the kernels as well.
    // TODO: test varying block sizes on modern arch-s as well
    // TODO: also consider using cudaFuncSetCacheConfig() for preferring shared memory on older architectures
    //(for spline data mostly)
    GMX_ASSERT(!(c_pmeAtomDataBlockSize % atomsPerBlock),
               "inconsistent atom data padding vs. spreading block size");

    // Ensure that coordinates are ready on the device before launching spread;
    // only needed on PP+PME ranks, not on separate PME ranks, in unit tests
    // as these cases use a single stream (hence xReadyOnDevice == nullptr).
    GMX_ASSERT(xReadyOnDevice != nullptr || pmeGpu->common->isRankPmeOnly
                       || pme_gpu_settings(pmeGpu).copyAllOutputs,
               "Need a valid coordinate synchronizer on PP+PME ranks with CUDA.");

    if (xReadyOnDevice)
    {
        xReadyOnDevice->enqueueWaitEvent(pmeGpu->archSpecific->pmeStream_);
    }

    // launch spread only if nAtoms > 0
    if (kernelParamsPtr->atoms.nAtoms > 0)
    {
        const int blockCount = pmeGpu->nAtomsAlloc / atomsPerBlock;
        auto      dimGrid    = pmeGpuCreateGrid(pmeGpu, blockCount);

        if (pmeGpu->common->ngrids == 1)
        {
            kernelParamsPtr->current.scale = 1.0;
        }
        else
        {
            kernelParamsPtr->current.scale = 1.0 - lambda;
        }

        KernelLaunchConfig config;
        config.blockSize[0] = order;
        config.blockSize[1] = (pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::Order ? 1 : order);
        config.blockSize[2] = atomsPerBlock;
        config.gridSize[0]  = dimGrid.first;
        config.gridSize[1]  = dimGrid.second;

        PmeStage                           timingId;
        PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
        const bool writeGlobalOrSaveSplines          = writeGlobal || (!recalculateSplines);
        if (computeSplines)
        {
            if (spreadCharges)
            {
                timingId  = PmeStage::SplineAndSpread;
                kernelPtr = selectSplineAndSpreadKernelPtr(pmeGpu,
                                                           pmeGpu->settings.threadsPerAtom,
                                                           writeGlobalOrSaveSplines,
                                                           pmeGpu->common->ngrids);
            }
            else
            {
                timingId  = PmeStage::Spline;
                kernelPtr = selectSplineKernelPtr(pmeGpu,
                                                  pmeGpu->settings.threadsPerAtom,
                                                  writeGlobalOrSaveSplines,
                                                  pmeGpu->common->ngrids);
            }
        }
        else
        {
            timingId  = PmeStage::Spread;
            kernelPtr = selectSpreadKernelPtr(pmeGpu,
                                              pmeGpu->settings.threadsPerAtom,
                                              writeGlobalOrSaveSplines,
                                              pmeGpu->common->ngrids);
        }


        pme_gpu_start_timing(pmeGpu, timingId);
        auto* timingEvent = pme_gpu_fetch_timing_event(pmeGpu, timingId);

        kernelParamsPtr->usePipeline = char(computeSplines && spreadCharges && useGpuDirectComm
                                            && (pmeCoordinateReceiverGpu->ppCommNumSenderRanks() > 1)
                                            && !writeGlobalOrSaveSplines);
        if (kernelParamsPtr->usePipeline != 0)
        {
            const int numStagesInPipeline = pmeCoordinateReceiverGpu->ppCommNumSenderRanks();

            GpuEventSynchronizer* gridsReadyForSpread = &pmeGpu->archSpecific->pmeGridsReadyForSpread;
            // Sync on grid zeroing is required except when GPU graphs are in use,
            // In which case the sync is already present through the zeroing being
            // explicitly included in the graph
            if (!useMdGpuGraph)
            {
                gridsReadyForSpread->markEvent(pmeGpu->archSpecific->pmeStream_);
                gridsReadyForSpread->setConsumptionLimits(numStagesInPipeline, numStagesInPipeline);
            }

            for (int i = 0; i < numStagesInPipeline; i++)
            {
                wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmePPRecvX);
                int senderRank = manageSyncWithPpCoordinateSenderGpu(
                        pmeGpu, pmeCoordinateReceiverGpu, kernelParamsPtr->usePipeline != 0, i);
                wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmePPRecvX);

                wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

                DeviceStream* launchStream = pmeCoordinateReceiverGpu->ppCommStream(senderRank);
                if (!useMdGpuGraph)
                {
                    gridsReadyForSpread->enqueueWaitEvent(*launchStream);
                }
                // set kernel configuration options specific to this stage of the pipeline
                std::tie(kernelParamsPtr->pipelineAtomStart, kernelParamsPtr->pipelineAtomEnd) =
                        pmeCoordinateReceiverGpu->ppCommAtomRange(senderRank);
                const int blockCount = static_cast<int>(std::ceil(
                        static_cast<float>(kernelParamsPtr->pipelineAtomEnd - kernelParamsPtr->pipelineAtomStart)
                        / atomsPerBlock));
                auto      dimGrid    = pmeGpuCreateGrid(pmeGpu, blockCount);
                config.gridSize[0]   = dimGrid.first;
                config.gridSize[1]   = dimGrid.second;


#if c_canEmbedBuffers
                const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
                const auto kernelArgs =
                        prepareGpuKernelArguments(kernelPtr,
                                                  config,
                                                  kernelParamsPtr,
                                                  &kernelParamsPtr->atoms.d_theta,
                                                  &kernelParamsPtr->atoms.d_dtheta,
                                                  &kernelParamsPtr->atoms.d_gridlineIndices,
                                                  &kernelParamsPtr->grid.d_realGrid[FEP_STATE_A],
                                                  &kernelParamsPtr->grid.d_realGrid[FEP_STATE_B],
                                                  &kernelParamsPtr->grid.d_fractShiftsTable,
                                                  &kernelParamsPtr->grid.d_gridlineIndicesTable,
                                                  &kernelParamsPtr->atoms.d_coefficients[FEP_STATE_A],
                                                  &kernelParamsPtr->atoms.d_coefficients[FEP_STATE_B],
                                                  &kernelParamsPtr->atoms.d_coordinates);
#endif

                launchGpuKernel(kernelPtr, config, *launchStream, timingEvent, "PME spline/spread", kernelArgs);
                wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
            }
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);
            // Set dependencies for PME stream on all pipeline streams
            for (int i = 0; i < pmeCoordinateReceiverGpu->ppCommNumSenderRanks(); i++)
            {
                pmeCoordinateReceiverGpu->insertAsDependencyIntoStream(i, pmeGpu->archSpecific->pmeStream_);
            }
            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }
        else // pipelining is not in use
        {
            if (useGpuDirectComm) // Sync all PME-PP communications to PME stream
            {
                wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmePPRecvX);
                for (int i = 0; i < pmeCoordinateReceiverGpu->ppCommNumSenderRanks(); i++)
                {
                    manageSyncWithPpCoordinateSenderGpu(pmeGpu, pmeCoordinateReceiverGpu);
                }
                wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmePPRecvX);
            }

            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

#if c_canEmbedBuffers
            const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
            const auto kernelArgs =
                    prepareGpuKernelArguments(kernelPtr,
                                              config,
                                              kernelParamsPtr,
                                              &kernelParamsPtr->atoms.d_theta,
                                              &kernelParamsPtr->atoms.d_dtheta,
                                              &kernelParamsPtr->atoms.d_gridlineIndices,
                                              &kernelParamsPtr->grid.d_realGrid[FEP_STATE_A],
                                              &kernelParamsPtr->grid.d_realGrid[FEP_STATE_B],
                                              &kernelParamsPtr->grid.d_fractShiftsTable,
                                              &kernelParamsPtr->grid.d_gridlineIndicesTable,
                                              &kernelParamsPtr->atoms.d_coefficients[FEP_STATE_A],
                                              &kernelParamsPtr->atoms.d_coefficients[FEP_STATE_B],
                                              &kernelParamsPtr->atoms.d_coordinates);
#endif

            launchGpuKernel(kernelPtr,
                            config,
                            pmeGpu->archSpecific->pmeStream_,
                            timingEvent,
                            "PME spline/spread",
                            kernelArgs);

            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }

        pme_gpu_stop_timing(pmeGpu, timingId);
    }

    const auto& settings = pmeGpu->settings;

    // halo exchange
    if (settings.useDecomposition)
    {
        pmeGpuGridHaloExchange(pmeGpu, wcycle);
    }

    wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

    // full PME GPU decomposition
    const bool convertPmeToFftGridOnGpu = settings.performGPUFFT && settings.useDecomposition;
    if (convertPmeToFftGridOnGpu)
    {
        // non-contiguous data - need to run kernel
        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            convertPmeGridToFftGrid<true>(pmeGpu, &pmeGpu->archSpecific->d_fftRealGrid[gridIndex], gridIndex);
        }
    }

    const bool copyBackGrid = spreadCharges && (!settings.performGPUFFT || settings.copyAllOutputs);
    if (copyBackGrid)
    {
        if (settings.useDecomposition) // mixed mode - multiple PME rank
        {
            // non-contiguous data - need to run kernel
            for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
            {
                float* h_grid = h_grids[gridIndex];

                convertPmeGridToFftGrid<true>(pmeGpu, h_grid, fftSetup, gridIndex);
            }
        }
        else // mixed mode - single PME rank
        {
            for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
            {
                float* h_grid = h_grids[gridIndex];
                pme_gpu_copy_output_spread_grid(pmeGpu, h_grid, gridIndex);
            }
        }
    }
    const bool copyBackAtomData =
            computeSplines && (!settings.performGPUGather || settings.copyAllOutputs);
    if (copyBackAtomData)
    {
        pme_gpu_copy_output_spread_atom_data(pmeGpu);
    }

    wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
}

void pme_gpu_solve(const PmeGpu* pmeGpu,
                   const int     gridIndex,
                   t_complex*    h_grid,
                   GridOrdering  gridOrdering,
                   bool          computeEnergyAndVirial)
{
    GMX_ASSERT(
            pmeGpu->common->ngrids == 1 || pmeGpu->common->ngrids == 2,
            "Only one (normal Coulomb PME) or two (FEP coulomb PME) PME grids can be used on GPU");
    GMX_ASSERT(gridIndex < pmeGpu->common->ngrids,
               "Invalid combination of gridIndex and number of grids");

    const auto& settings               = pmeGpu->settings;
    const bool  copyInputAndOutputGrid = !settings.performGPUFFT || settings.copyAllOutputs;

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();

    float*    h_gridFloat              = reinterpret_cast<float*>(h_grid);
    const int localComplexGridElements = kernelParamsPtr->grid.localComplexGridSizePadded[XX]
                                         * kernelParamsPtr->grid.localComplexGridSizePadded[YY]
                                         * kernelParamsPtr->grid.localComplexGridSizePadded[ZZ] * 2;
    if (copyInputAndOutputGrid)
    {
        copyToDeviceBuffer(&kernelParamsPtr->grid.d_fftComplexGrid[gridIndex],
                           h_gridFloat,
                           0,
                           localComplexGridElements,
                           pmeGpu->archSpecific->pmeStream_,
                           pmeGpu->settings.transferKind,
                           nullptr);
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

    const int gridLineSize      = pmeGpu->kernelParams->grid.localComplexGridSize[minorDim];
    const int gridLinesPerBlock = std::max(maxBlockSize / gridLineSize, 1);
    const int blocksPerGridLine = (gridLineSize + maxBlockSize - 1) / maxBlockSize;
    int       cellsPerBlock;

    if (pmeGpu->common->nnodesY > 1
        && pmeGpu->kernelParams->grid.complexGridSize[ZZ] >= pmeGpu->common->nnodes)
    {
        pmeGpu->kernelParams->grid.kOffsets[XX] = 0;
        pmeGpu->kernelParams->grid.kOffsets[YY] = 0;
        pmeGpu->kernelParams->grid.kOffsets[ZZ] = pmeGpu->common->nodeid
                                                  * pmeGpu->kernelParams->grid.complexGridSize[ZZ]
                                                  / pmeGpu->common->nnodes;
    }
    else
    {
        pmeGpu->kernelParams->grid.kOffsets[XX] = 0;
        pmeGpu->kernelParams->grid.kOffsets[YY] = pmeGpu->common->nodeidX
                                                  * pmeGpu->kernelParams->grid.complexGridSize[YY]
                                                  / pmeGpu->common->nnodesX;
        pmeGpu->kernelParams->grid.kOffsets[ZZ] = pmeGpu->common->nodeidY
                                                  * pmeGpu->kernelParams->grid.complexGridSize[ZZ]
                                                  / pmeGpu->common->nnodesY;
    }

    if (blocksPerGridLine == 1)
    {
        cellsPerBlock = gridLineSize * gridLinesPerBlock;
    }
    else
    {
        cellsPerBlock = (gridLineSize + blocksPerGridLine - 1) / blocksPerGridLine;
    }
    const int warpSize  = pmeGpu->programHandle_->warpSize();
    const int blockSize = (cellsPerBlock + warpSize - 1) / warpSize * warpSize;

    static_assert(!GMX_GPU_CUDA || c_solveMaxWarpsPerBlock / 2 >= 4,
                  "The CUDA solve energy kernels needs at least 4 warps. "
                  "Here we launch at least half of the max warps.");

    KernelLaunchConfig config;
    config.blockSize[0] = blockSize;
    config.gridSize[0]  = blocksPerGridLine;
    // rounding up to full warps so that shuffle operations produce defined results
    config.gridSize[1] =
            (pmeGpu->kernelParams->grid.localComplexGridSize[middleDim] + gridLinesPerBlock - 1)
            / gridLinesPerBlock;

    config.gridSize[2] = pmeGpu->kernelParams->grid.localComplexGridSize[majorDim];

    PmeStage                           timingId  = PmeStage::Solve;
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;
    if (gridOrdering == GridOrdering::YZX)
    {
        if (gridIndex == 0)
        {
            kernelPtr = computeEnergyAndVirial ? pmeGpu->programHandle_->impl_->solveYZXEnergyKernelA
                                               : pmeGpu->programHandle_->impl_->solveYZXKernelA;
        }
        else
        {
            kernelPtr = computeEnergyAndVirial ? pmeGpu->programHandle_->impl_->solveYZXEnergyKernelB
                                               : pmeGpu->programHandle_->impl_->solveYZXKernelB;
        }
    }
    else if (gridOrdering == GridOrdering::XYZ)
    {
        if (gridIndex == 0)
        {
            kernelPtr = computeEnergyAndVirial ? pmeGpu->programHandle_->impl_->solveXYZEnergyKernelA
                                               : pmeGpu->programHandle_->impl_->solveXYZKernelA;
        }
        else
        {
            kernelPtr = computeEnergyAndVirial ? pmeGpu->programHandle_->impl_->solveXYZEnergyKernelB
                                               : pmeGpu->programHandle_->impl_->solveXYZKernelB;
        }
    }

    pme_gpu_start_timing(pmeGpu, timingId);
    auto* timingEvent = pme_gpu_fetch_timing_event(pmeGpu, timingId);
#if c_canEmbedBuffers
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
    const auto kernelArgs =
            prepareGpuKernelArguments(kernelPtr,
                                      config,
                                      kernelParamsPtr,
                                      &kernelParamsPtr->grid.d_splineModuli[gridIndex],
                                      &kernelParamsPtr->constants.d_virialAndEnergy[gridIndex],
                                      &kernelParamsPtr->grid.d_fftComplexGrid[gridIndex]);
#endif
    launchGpuKernel(kernelPtr, config, pmeGpu->archSpecific->pmeStream_, timingEvent, "PME solve", kernelArgs);
    pme_gpu_stop_timing(pmeGpu, timingId);

    if (computeEnergyAndVirial)
    {
        copyFromDeviceBuffer(pmeGpu->staging.h_virialAndEnergy[gridIndex],
                             &kernelParamsPtr->constants.d_virialAndEnergy[gridIndex],
                             0,
                             c_virialAndEnergyCount,
                             pmeGpu->archSpecific->pmeStream_,
                             pmeGpu->settings.transferKind,
                             nullptr);
    }

    if (copyInputAndOutputGrid)
    {
        copyFromDeviceBuffer(h_gridFloat,
                             &kernelParamsPtr->grid.d_fftComplexGrid[gridIndex],
                             0,
                             localComplexGridElements,
                             pmeGpu->archSpecific->pmeStream_,
                             pmeGpu->settings.transferKind,
                             nullptr);
    }
}

/*! \brief
 * Returns a pointer to appropriate gather kernel based on the inputvalues
 *
 * \param[in]  pmeGpu                   The PME GPU structure.
 * \param[in]  threadsPerAtom           Controls whether we should use order or order*order threads per atom
 * \param[in]  readSplinesFromGlobal    bool controlling if we should write spline data to global memory
 * \param[in]  numGrids                 Number of grids to use. numGrids == 2 if Coulomb is perturbed.
 *
 * \return Pointer to CUDA kernel
 */
inline auto selectGatherKernelPtr(const PmeGpu*  pmeGpu,
                                  ThreadsPerAtom threadsPerAtom,
                                  bool           readSplinesFromGlobal,
                                  const int      numGrids)

{
    PmeGpuProgramImpl::PmeKernelHandle kernelPtr = nullptr;

    if (readSplinesFromGlobal)
    {
        if (threadsPerAtom == ThreadsPerAtom::Order)
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelReadSplinesThPerAtom4Dual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelReadSplinesThPerAtom4Single;
            }
        }
        else
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelReadSplinesDual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelReadSplinesSingle;
            }
        }
    }
    else
    {
        if (threadsPerAtom == ThreadsPerAtom::Order)
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelThPerAtom4Dual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelThPerAtom4Single;
            }
        }
        else
        {
            if (numGrids == 2)
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelDual;
            }
            else
            {
                kernelPtr = pmeGpu->programHandle_->impl_->gatherKernelSingle;
            }
        }
    }
    return kernelPtr;
}

void pme_gpu_gather(PmeGpu*               pmeGpu,
                    real**                h_grids,
                    gmx_parallel_3dfft_t* fftSetup,
                    const float           lambda,
                    gmx_wallcycle*        wcycle,
                    bool                  computeVirial)
{
    GMX_ASSERT(
            pmeGpu->common->ngrids == 1 || pmeGpu->common->ngrids == 2,
            "Only one (normal Coulomb PME) or two (FEP coulomb PME) PME grids can be used on GPU");

    const auto& settings        = pmeGpu->settings;
    auto*       kernelParamsPtr = pmeGpu->kernelParams.get();

    wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

    // full PME GPU decomposition
    const bool convertFftToPmeGridOnGpu = settings.performGPUFFT && settings.useDecomposition;
    if (convertFftToPmeGridOnGpu)
    {
        // non-contiguous data - need to run kernel
        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            convertPmeGridToFftGrid<false>(
                    pmeGpu, &pmeGpu->archSpecific->d_fftRealGrid[gridIndex], gridIndex);
        }
    }

    const bool copyBackGrid = (!settings.performGPUFFT || settings.copyAllOutputs);

    if (copyBackGrid)
    {
        if (settings.useDecomposition) // mixed mode - multple PME rank
        {
            // non-contiguous data - need to run kernel
            for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
            {
                float* h_grid = h_grids[gridIndex];
                convertPmeGridToFftGrid<false>(pmeGpu, h_grid, fftSetup, gridIndex);
            }
        }
        else // mixed mode - single PME rank
        {
            for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
            {
                const float* h_grid = h_grids[gridIndex];
                pme_gpu_copy_input_gather_grid(pmeGpu, h_grid, gridIndex);
            }
        }
    }

    if (settings.copyAllOutputs)
    {
        pme_gpu_copy_input_gather_atom_data(pmeGpu);
    }

    wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);

    // reverse halo exchange
    if (settings.useDecomposition)
    {
        pmeGpuGridHaloExchangeReverse(pmeGpu, wcycle);
    }

    wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

    /* Set if we have unit tests */
    const bool   readGlobal = pmeGpu->settings.copyAllOutputs;
    const size_t blockSize  = pmeGpu->programHandle_->impl_->gatherWorkGroupSize;
    const int    order      = pmeGpu->common->pme_order;
    GMX_ASSERT(order == c_pmeGpuOrder, "Only PME order 4 is implemented");
    const int threadsPerAtom =
            (pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::Order ? order : order * order);
    const bool recalculateSplines = pmeGpu->settings.recalculateSplines;

    GMX_ASSERT(!GMX_GPU_OPENCL || pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::OrderSquared,
               "Only 16 threads per atom supported in OpenCL");
    GMX_ASSERT(!GMX_GPU_OPENCL || !recalculateSplines,
               "Recalculating splines not supported in OpenCL");

    const int atomsPerBlock = blockSize / threadsPerAtom;

    GMX_ASSERT(!(c_pmeAtomDataBlockSize % atomsPerBlock),
               "inconsistent atom data padding vs. gathering block size");

    // launch gather only if nAtoms > 0
    if (kernelParamsPtr->atoms.nAtoms > 0)
    {
        const int blockCount = pmeGpu->nAtomsAlloc / atomsPerBlock;
        auto      dimGrid    = pmeGpuCreateGrid(pmeGpu, blockCount);

        KernelLaunchConfig config;
        config.blockSize[0] = order;
        config.blockSize[1] = (pmeGpu->settings.threadsPerAtom == ThreadsPerAtom::Order ? 1 : order);
        config.blockSize[2] = atomsPerBlock;
        config.gridSize[0]  = dimGrid.first;
        config.gridSize[1]  = dimGrid.second;

        // TODO test different cache configs

        PmeStage                           timingId = PmeStage::Gather;
        PmeGpuProgramImpl::PmeKernelHandle kernelPtr =
                selectGatherKernelPtr(pmeGpu,
                                      pmeGpu->settings.threadsPerAtom,
                                      readGlobal || (!recalculateSplines),
                                      pmeGpu->common->ngrids);
        // TODO design kernel selection getters and make PmeGpu a friend of PmeGpuProgramImpl

        pme_gpu_start_timing(pmeGpu, timingId);
        auto* timingEvent = pme_gpu_fetch_timing_event(pmeGpu, timingId);
        if (pmeGpu->common->ngrids == 1)
        {
            kernelParamsPtr->current.scale = 1.0;
        }
        else
        {
            kernelParamsPtr->current.scale = 1.0 - lambda;
        }

#if GMX_NVSHMEM
        kernelParamsPtr->isVirialStep = computeVirial;
        if (!computeVirial && pmeGpu->useNvshmem)
        {
            kernelParamsPtr->forcesReadyNvshmemFlagsCounter++;
        }
#else
        GMX_UNUSED_VALUE(computeVirial);
#endif

#if c_canEmbedBuffers
        const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, kernelParamsPtr);
#else
        const auto kernelArgs =
                prepareGpuKernelArguments(kernelPtr,
                                          config,
                                          kernelParamsPtr,
                                          &kernelParamsPtr->atoms.d_coefficients[FEP_STATE_A],
                                          &kernelParamsPtr->atoms.d_coefficients[FEP_STATE_B],
                                          &kernelParamsPtr->grid.d_realGrid[FEP_STATE_A],
                                          &kernelParamsPtr->grid.d_realGrid[FEP_STATE_B],
                                          &kernelParamsPtr->atoms.d_theta,
                                          &kernelParamsPtr->atoms.d_dtheta,
                                          &kernelParamsPtr->atoms.d_gridlineIndices,
                                          &kernelParamsPtr->atoms.d_forces);
#endif
        launchGpuKernel(
                kernelPtr, config, pmeGpu->archSpecific->pmeStream_, timingEvent, "PME gather", kernelArgs);
        if (!computeVirial && pmeGpu->useNvshmem)
        {
            KernelLaunchConfig config;
            config.blockSize[0] = 128;
            config.blockSize[1] = 1;
            config.blockSize[2] = 1;
            config.gridSize[0]  = 1;
            config.gridSize[1]  = 1;
            auto kernelPtr_     = pmeGpu->programHandle_->impl_->nvshmemSignalKern;

            const auto kernelArgs_ = prepareGpuKernelArguments(kernelPtr_, config, kernelParamsPtr);
            launchGpuKernel(
                    kernelPtr_, config, pmeGpu->archSpecific->pmeStream_, timingEvent, "PME gather", kernelArgs_);
        }
        pme_gpu_stop_timing(pmeGpu, timingId);
    }

    if (pmeGpu->settings.useGpuForceReduction)
    {
        pmeGpu->archSpecific->pmeForcesReady.markEvent(pmeGpu->archSpecific->pmeStream_);
    }
    else if (pmeGpu->kernelParams->atoms.nAtoms > 0)
    {
        pme_gpu_copy_output_forces(pmeGpu);
    }

    wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
}

DeviceBuffer<gmx::RVec> pme_gpu_get_kernelparam_forces(const PmeGpu* pmeGpu)
{
    if (pmeGpu && pmeGpu->kernelParams)
    {
        return pmeGpu->kernelParams->atoms.d_forces;
    }
    else
    {
        return DeviceBuffer<gmx::RVec>{};
    }
}

void pme_gpu_set_kernelparam_useNvshmem(const PmeGpu* pmeGpu, bool useNvshmem)
{
    GMX_ASSERT(pmeGpu && pmeGpu->kernelParams,
               "PME GPU NVSHMEM support can not be set in non-GPU builds or before the GPU PME was "
               "initialized.");

    if (pmeGpu && pmeGpu->kernelParams)
    {
        pmeGpu->kernelParams->useNvshmem = static_cast<int>(useNvshmem);
    }
}

void pme_gpu_set_kernelparam_coordinates(const PmeGpu* pmeGpu, DeviceBuffer<gmx::RVec> d_x)
{
    GMX_ASSERT(pmeGpu && pmeGpu->kernelParams,
               "PME GPU device buffer can not be set in non-GPU builds or before the GPU PME was "
               "initialized.");

    GMX_ASSERT(checkDeviceBuffer(d_x, pmeGpu->kernelParams->atoms.nAtoms),
               "The device-side buffer can not be set.");

    pmeGpu->kernelParams->atoms.d_coordinates = d_x;
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
