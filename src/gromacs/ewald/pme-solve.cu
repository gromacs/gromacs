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
 *  \brief Implements PME GPU Fourier grid solving in CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "pme.cuh"
#include "pme-timings.cuh"

/*! \brief
 * PME complex grid solver kernel function.
 *
 * \tparam[in] computeEnergyAndVirial   Tells if the reciprocal energy and virial should be computed.
 * \tparam[in] yzxGridOrdering          True if the solving should be performed on the complex grid with YZX dimension order;
 *
 *                                      False if the grid has XYZ dimension order.
 * \param[in]  kernelParams             Input PME CUDA data in constant memory.
 */
template<
    const bool computeEnergyAndVirial,
    const bool yzxGridOrdering
    >
//__launch_bounds__(PME_SOLVE_THREADS_PER_BLOCK, PME_MIN_BLOCKS_PER_MP) // TODO: figure out when and why this produces "invalid launch argument"
__global__ void pme_solve_kernel(const struct pme_gpu_cuda_kernel_params_t kernelParams)
{
    /* This kernel supports 2 different dimension orderings: YZX and XYZ */
    const int   majorDim  = yzxGridOrdering ? YY : XX;
    const int   middleDim = yzxGridOrdering ? ZZ : YY;
    const int   minorDim  = yzxGridOrdering ? XX : ZZ;

    /* Global memory pointers */
    const float * __restrict__ gm_splineValueMajor    = kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[majorDim];
    const float * __restrict__ gm_splineValueMiddle   = kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[middleDim];
    const float * __restrict__ gm_splineValueMinor    = kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[minorDim];
    float * __restrict__       gm_virialAndEnergy     = kernelParams.constants.d_virialAndEnergy;
    float2 * __restrict__      gm_grid                = (float2 *)kernelParams.grid.d_fourierGrid;

    const int                  localOffsetMinor = 0, localOffsetMajor = 0, localOffsetMiddle = 0; //unused
    const int                  localSizeMinor   = kernelParams.grid.complexGridSizePadded[minorDim];
    const int                  localSizeMiddle  = kernelParams.grid.complexGridSizePadded[middleDim];
    const int                  localCountMajor  = kernelParams.grid.complexGridSize[majorDim];
    const int                  localCountMiddle = kernelParams.grid.complexGridSize[middleDim];
    const int                  localCountMinor  = kernelParams.grid.complexGridSize[minorDim];

    // Each thread works on one cell of the Fourier space complex 3D grid (float2 * __restrict__ grid).
    // Each block handles PME_SOLVE_THREADS_PER_BLOCK cells -
    // depending on the grid contiguous dimension size,
    // that can range from a part of a single gridline to several complete gridlines.
    const int              threadLocalId = (threadIdx.y * blockDim.x) + threadIdx.x;

    const int              nMajor  = kernelParams.grid.complexGridSize[majorDim];
    const int              nMiddle = kernelParams.grid.complexGridSize[middleDim];
    const int              nMinor  = kernelParams.grid.complexGridSize[minorDim];

    int                    maxkMajor  = (nMajor + 1) / 2;                                      //X or Y
    int                    maxkMiddle = (nMiddle + 1) / 2;                                     //Y OR Z => only check for !YZX
    int                    maxkMinor  = (nMinor + 1) / 2;                                      //Z or X => only check for YZX

    float                  energy = 0.0f;
    float                  virxx  = 0.0f, virxy = 0.0f, virxz = 0.0f, viryy = 0.0f, viryz = 0.0f, virzz = 0.0f;

    const int              indexMinor  = blockIdx.x * blockDim.x + threadIdx.x;
    const int              indexMiddle = blockIdx.y * blockDim.y + threadIdx.y;
    const int              indexMajor  = blockIdx.z * blockDim.z + threadIdx.z;

    if ((indexMajor < localCountMajor) && (indexMiddle < localCountMiddle) && (indexMinor < localCountMinor))
    {
        /* The offset should be equal to the global thread index */
        const int   gridIndex     = (indexMajor * localSizeMiddle + indexMiddle) * localSizeMinor + indexMinor;
        float2     *globalGridPtr = gm_grid + gridIndex;
        //TODO reuse indexing function from tests

        const int   kMajor = indexMajor + localOffsetMajor;
        /* Checking either X in XYZ, or Y in YZX cases */
        const float mMajor = (kMajor < maxkMajor) ? kMajor : (kMajor - nMajor);

        const int   kMiddle = indexMiddle + localOffsetMiddle;
        float       mMiddle = kMiddle;
        /* Checking Y in XYZ case */
        if (!yzxGridOrdering)
        {
            mMiddle = (kMiddle < maxkMiddle) ? kMiddle : (kMiddle - nMiddle);
        }
        const int       kMinor        = localOffsetMinor + indexMinor;
        float           mMinor        = kMinor;
        /* We should skip the k-space point (0,0,0) */
        const bool      notZeroPoint  = (kMinor > 0) | (kMajor > 0) | (kMiddle > 0);

        /* Checking X in YZX case */
        if (yzxGridOrdering)
        {
            mMinor = (kMinor < maxkMinor) ? kMinor : (kMinor - nMinor);
        }

        float mX, mY, mZ;
        if (yzxGridOrdering)
        {
            mX = mMinor;
            mY = mMajor;
            mZ = mMiddle;
        }
        else
        {
            mX = mMajor;
            mY = mMiddle;
            mZ = mMinor;
        }

        /* 0.5 correction for corner points of a minor dimension */
        float corner_fac = 1.0f;
        if (yzxGridOrdering)
        {
            if (kMiddle == 0 || kMiddle == maxkMiddle)
            {
                corner_fac = 0.5f;
            }
        }
        else
        {
            if (kMinor == 0 || kMinor == maxkMinor)
            {
                corner_fac = 0.5f;
            }
        }

        if (notZeroPoint)
        {
            const float mhxk       = mX * kernelParams.step.recipBox[XX][XX];
            const float mhyk       = mX * kernelParams.step.recipBox[XX][YY] + mY * kernelParams.step.recipBox[YY][YY];
            const float mhzk       = mX * kernelParams.step.recipBox[XX][ZZ] + mY * kernelParams.step.recipBox[YY][ZZ] + mZ * kernelParams.step.recipBox[ZZ][ZZ];

            const float m2k        = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
            assert(m2k != 0.0f);
            float       denom = m2k * float(M_PI) * kernelParams.step.boxVolume * gm_splineValueMajor[kMajor] * gm_splineValueMiddle[kMiddle] * gm_splineValueMinor[kMinor];
            assert(!isnan(denom));
            assert(denom != 0.0f);
            const float   tmp1   = expf(-kernelParams.grid.ewaldFactor * m2k);
            const float   etermk = kernelParams.constants.elFactor * tmp1 / denom;

            float2        gridValue    = *globalGridPtr;
            const float2  oldGridValue = gridValue;
            gridValue.x   *= etermk;
            gridValue.y   *= etermk;
            *globalGridPtr = gridValue;

            if (computeEnergyAndVirial)
            {
                const float tmp1k = 2.0f * (gridValue.x * oldGridValue.x + gridValue.y * oldGridValue.y);

                float       vfactor = (kernelParams.grid.ewaldFactor  + 1.0f / m2k) * 2.0f;
                float       ets2    = corner_fac * tmp1k;
                energy = ets2;

                float ets2vf  = ets2 * vfactor;

                virxx   = ets2vf * mhxk * mhxk - ets2;
                virxy   = ets2vf * mhxk * mhyk;
                virxz   = ets2vf * mhxk * mhzk;
                viryy   = ets2vf * mhyk * mhyk - ets2;
                viryz   = ets2vf * mhyk * mhzk;
                virzz   = ets2vf * mhzk * mhzk - ets2;
            }
        }
    }

    if (computeEnergyAndVirial)
    {
        /* The energy and virial reduction */

        /* TODO: there could be a shuffle reduction variant here as well!*/
        const int        blockSize = PME_SOLVE_ENERVIR_THREADS_PER_BLOCK;
        __shared__ float sm_virialAndEnergy[c_virialAndEnergyCount * blockSize];
        // 3.5k smem per block - a serious limiter!

        /*  a 7-thread reduction in shared memory inspired by reduce_force_j_generic */
        if (threadLocalId < blockSize)
        {
            sm_virialAndEnergy[threadLocalId + 0 * blockSize] = virxx;
            sm_virialAndEnergy[threadLocalId + 1 * blockSize] = viryy;
            sm_virialAndEnergy[threadLocalId + 2 * blockSize] = virzz;
            sm_virialAndEnergy[threadLocalId + 3 * blockSize] = virxy;
            sm_virialAndEnergy[threadLocalId + 4 * blockSize] = virxz;
            sm_virialAndEnergy[threadLocalId + 5 * blockSize] = viryz;
            sm_virialAndEnergy[threadLocalId + 6 * blockSize] = energy;
        }
        __syncthreads();

        /* Reducing every component to fit into warp_size */
        for (int s = blockSize >> 1; s >= warp_size; s >>= 1)
        {
#pragma unroll
            for (int i = 0; i < c_virialAndEnergyCount; i++)
            {
                if (threadLocalId < s) // split per threads?
                {
                    sm_virialAndEnergy[i * blockSize + threadLocalId] += sm_virialAndEnergy[i * blockSize + threadLocalId + s];
                }
            }
            __syncthreads();
        }

        const int threadsPerComponent    = warp_size / c_virialAndEnergyCount; // this is also the stride, will be 32 / 7 = 4
        const int contributionsPerThread = warp_size / threadsPerComponent;    // will be 32 / 4 = 8
        if (threadLocalId < c_virialAndEnergyCount * threadsPerComponent)
        {
            const int componentIndex        = threadLocalId / threadsPerComponent;
            const int threadComponentOffset = threadLocalId - componentIndex * threadsPerComponent;

            float     sum = 0.0f;
#pragma unroll
            for (int j = 0; j < contributionsPerThread; j++)
            {
                sum += sm_virialAndEnergy[componentIndex * blockSize + j * threadsPerComponent + threadComponentOffset];
            }
            atomicAdd(gm_virialAndEnergy + componentIndex, sum);
        }

        /* A naive reduction for debugging purposes */
        /*
           if (threadLocalId < blockSize)
           {
            sm_virialAndEnergy[sizing * threadLocalId + 0] = virxx;
            sm_virialAndEnergy[sizing * threadLocalId + 1] = viryy;
            sm_virialAndEnergy[sizing * threadLocalId + 2] = virzz;
            sm_virialAndEnergy[sizing * threadLocalId + 3] = virxy;
            sm_virialAndEnergy[sizing * threadLocalId + 4] = virxz;
            sm_virialAndEnergy[sizing * threadLocalId + 5] = viryz;
            sm_virialAndEnergy[sizing * threadLocalId + 6] = energy;
           }
           __syncthreads();
           #pragma unroll
           for (unsigned int stride = 1; stride < blockSize; stride <<= 1)
           {
            if ((threadLocalId % (stride << 1) == 0))
            {
           #pragma unroll
                for (int i = 0; i < sizing; i++)
                    sm_virialAndEnergy[sizing * threadLocalId + i] += sm_virialAndEnergy[sizing * (threadLocalId + stride) + i];
            }
            __syncthreads();
           }
           if (threadLocalId < sizing)
           {
            atomicAdd(virialAndEnergyGlobal + threadLocalId, sm_virialAndEnergy[threadLocalId]);
           }
         */
    }
}

void pme_gpu_solve(const pme_gpu_t *pmeGpu, t_complex *h_grid,
                   bool computeEnergyAndVirial, bool yzxGridOrdering)
{
    /* do recip sum over local cells in grid */

    const bool   copyInputAndOutputGrid = pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu);

    cudaStream_t stream          = pmeGpu->archSpecific->pmeStream;
    const auto  *kernelParamsPtr = pmeGpu->kernelParams.get();

    float2      *d_grid = (float2 *)kernelParamsPtr->grid.d_fourierGrid;
    if (copyInputAndOutputGrid)
    {
        cu_copy_H2D_async(d_grid, h_grid, pmeGpu->archSpecific->complexGridSize * sizeof(float), stream);
    }

    const int   majorDim  = yzxGridOrdering ? YY : XX;
    const int   middleDim = yzxGridOrdering ? ZZ : YY;
    const int   minorDim  = yzxGridOrdering ? XX : ZZ;

    const int   maxBlockSize      = computeEnergyAndVirial ? PME_SOLVE_ENERVIR_THREADS_PER_BLOCK : PME_SOLVE_THREADS_PER_BLOCK;
    const int   gridLineSize      = pmeGpu->kernelParams->grid.complexGridSizePadded[minorDim];
    const int   gridLinesPerBlock = max(maxBlockSize / gridLineSize, 1);
    const int   blocksPerGridLine = (gridLineSize + maxBlockSize - 1) / maxBlockSize; // rounded up
    // Z-dimension is too small in CUDA limitations (64 on CC30?), so instead of major-middle-minor sizing we do minor-middle-major
    // TODO: maybe that is why YZX fails occasionally with launch bounds?
    dim3 threads((maxBlockSize + gridLinesPerBlock - 1) / gridLinesPerBlock, gridLinesPerBlock);
    const int blockSize = threads.x * threads.y * threads.z;
    GMX_ASSERT(blockSize >= maxBlockSize, "Wrong PME GPU solve launch parameters");
    dim3 blocks(blocksPerGridLine,
                (pmeGpu->kernelParams->grid.complexGridSize[middleDim] + gridLinesPerBlock - 1) / gridLinesPerBlock, // rounded up middle dimension block number
                pmeGpu->kernelParams->grid.complexGridSize[majorDim]);

    pme_gpu_start_timing(pmeGpu, gtPME_SOLVE);
    if (yzxGridOrdering)
    {
        if (computeEnergyAndVirial)
        {
            pme_solve_kernel<true, true> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
        else
        {
            pme_solve_kernel<false, true> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
    }
    else
    {
        if (computeEnergyAndVirial)
        {
            pme_solve_kernel<true, false> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
        else
        {
            pme_solve_kernel<false, false> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
    }
    CU_LAUNCH_ERR("pme_solve_kernel");
    pme_gpu_stop_timing(pmeGpu, gtPME_SOLVE);

    if (computeEnergyAndVirial)
    {
        cu_copy_D2H_async(pmeGpu->staging.h_virialAndEnergy, kernelParamsPtr->constants.d_virialAndEnergy,
                          c_virialAndEnergyCount * sizeof(float), stream);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncEnerVirD2H, stream);
        CU_RET_ERR(stat, "PME solve energy/virial event record failure");
    }

    if (copyInputAndOutputGrid)
    {
        cu_copy_D2H_async(h_grid, d_grid, pmeGpu->archSpecific->complexGridSize * sizeof(float), stream);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncSolveGridD2H, stream);
        CU_RET_ERR(stat, "PME solve grid sync event record failure");
    }
}
