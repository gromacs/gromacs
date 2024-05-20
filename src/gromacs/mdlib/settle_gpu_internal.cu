/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief CUDA-specific routines for the GPU implementation of SETTLE constraints algorithm.
 *
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settle_gpu_internal.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include <algorithm>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"

namespace gmx
{

//! Number of CUDA threads in a block
constexpr static int sc_threadsPerBlock = 256;

//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int sc_maxThreadsPerBlock = sc_threadsPerBlock;

/*! \brief SETTLE constraints kernel
 *
 * Each thread corresponds to a single constraints triangle (i.e. single water molecule).
 *
 * See original CPU version in settle.cpp
 *
 * \param [in]      numSettles       Number of constraints triangles (water molecules).
 * \param [in]      gm_settles       Indexes of three atoms in the constraints triangle. The field .x of int3
 *                                   data type corresponds to Oxygen, fields .y and .z are two hydrogen atoms.
 * \param [in]      pars             Parameters for the algorithm (i.e. masses, target distances, etc.).
 * \param [in]      gm_x             Coordinates of atoms before the timestep.
 * \param [in,out]  gm_x             Coordinates of atoms after the timestep (constrained coordinates will be
 *                                   saved here).
 * \param [in]      invdt            Reciprocal timestep.
 * \param [in]      gm_v             Velocities of the particles.
 * \param [in]      gm_virialScaled  Virial tensor.
 * \param [in]      pbcAiuc          Periodic boundary conditions data.
 */
template<bool updateVelocities, bool computeVirial>
__launch_bounds__(sc_maxThreadsPerBlock) __global__
        void settle_kernel(const int numSettles,
                           const WaterMolecule* __restrict__ gm_settles,
                           const SettleParameters pars,
                           const float3* __restrict__ gm_x,
                           float3* __restrict__ gm_xprime,
                           float invdt,
                           float3* __restrict__ gm_v,
                           float* __restrict__ gm_virialScaled,
                           const PbcAiuc pbcAiuc)
{
    /* ******************************************************************* */
    /*                                                                  ** */
    /*    Original code by Shuichi Miyamoto, last update Oct. 1, 1992   ** */
    /*                                                                  ** */
    /*    Algorithm changes by Berk Hess:                               ** */
    /*    2004-07-15 Convert COM to double precision to avoid drift     ** */
    /*    2006-10-16 Changed velocity update to use differences         ** */
    /*    2012-09-24 Use oxygen as reference instead of COM             ** */
    /*    2016-02    Complete rewrite of the code for SIMD              ** */
    /*    2020-06    Completely remove use of COM to minimize drift     ** */
    /*                                                                  ** */
    /*    Reference for the SETTLE algorithm                            ** */
    /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992).    ** */
    /*                                                                  ** */
    /* ******************************************************************* */

    constexpr float almost_zero = real(1e-12);

    extern __shared__ float sm_threadVirial[];

    int tid = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);

    if (tid < numSettles)
    {
        // These are the indexes of three atoms in a single 'water' molecule.
        // TODO Can be reduced to one integer if atoms are consecutive in memory.
        WaterMolecule indices = gm_settles[tid];

        float3 x_ow1 = gm_x[indices.ow1];
        float3 x_hw2 = gm_x[indices.hw2];
        float3 x_hw3 = gm_x[indices.hw3];

        float3 xprime_ow1 = gm_xprime[indices.ow1];
        float3 xprime_hw2 = gm_xprime[indices.hw2];
        float3 xprime_hw3 = gm_xprime[indices.hw3];

        float3 dist21 = pbcDxAiuc(pbcAiuc, x_hw2, x_ow1);
        float3 dist31 = pbcDxAiuc(pbcAiuc, x_hw3, x_ow1);
        float3 doh2   = pbcDxAiuc(pbcAiuc, xprime_hw2, xprime_ow1);

        float3 doh3 = pbcDxAiuc(pbcAiuc, xprime_hw3, xprime_ow1);

        float3 a1 = (-doh2 - doh3) * pars.wh;

        float3 b1 = doh2 + a1;

        float3 c1 = doh3 + a1;

        float xakszd = dist21.y * dist31.z - dist21.z * dist31.y;
        float yakszd = dist21.z * dist31.x - dist21.x * dist31.z;
        float zakszd = dist21.x * dist31.y - dist21.y * dist31.x;

        float xaksxd = a1.y * zakszd - a1.z * yakszd;
        float yaksxd = a1.z * xakszd - a1.x * zakszd;
        float zaksxd = a1.x * yakszd - a1.y * xakszd;

        float xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        float yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        float zaksyd = xakszd * yaksxd - yakszd * xaksxd;

        float axlng = rsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        float aylng = rsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        float azlng = rsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        // TODO {1,2,3} indexes should be swapped with {.x, .y, .z} components.
        //      This way, we will be able to use vector ops more.
        float3 trns1, trns2, trns3;

        trns1.x = xaksxd * axlng;
        trns2.x = yaksxd * axlng;
        trns3.x = zaksxd * axlng;

        trns1.y = xaksyd * aylng;
        trns2.y = yaksyd * aylng;
        trns3.y = zaksyd * aylng;

        trns1.z = xakszd * azlng;
        trns2.z = yakszd * azlng;
        trns3.z = zakszd * azlng;


        float2 b0d, c0d;

        b0d.x = trns1.x * dist21.x + trns2.x * dist21.y + trns3.x * dist21.z;
        b0d.y = trns1.y * dist21.x + trns2.y * dist21.y + trns3.y * dist21.z;

        c0d.x = trns1.x * dist31.x + trns2.x * dist31.y + trns3.x * dist31.z;
        c0d.y = trns1.y * dist31.x + trns2.y * dist31.y + trns3.y * dist31.z;

        float3 b1d, c1d;

        float a1d_z = trns1.z * a1.x + trns2.z * a1.y + trns3.z * a1.z;

        b1d.x = trns1.x * b1.x + trns2.x * b1.y + trns3.x * b1.z;
        b1d.y = trns1.y * b1.x + trns2.y * b1.y + trns3.y * b1.z;
        b1d.z = trns1.z * b1.x + trns2.z * b1.y + trns3.z * b1.z;

        c1d.x = trns1.x * c1.x + trns2.x * c1.y + trns3.x * c1.z;
        c1d.y = trns1.y * c1.x + trns2.y * c1.y + trns3.y * c1.z;
        c1d.z = trns1.z * c1.x + trns2.z * c1.y + trns3.z * c1.z;


        float sinphi = a1d_z * rsqrt(pars.ra * pars.ra);
        float tmp2   = 1.0F - sinphi * sinphi;

        if (almost_zero > tmp2)
        {
            tmp2 = almost_zero;
        }

        float tmp    = rsqrt(tmp2);
        float cosphi = tmp2 * tmp;
        float sinpsi = (b1d.z - c1d.z) * pars.irc2 * tmp;
        tmp2         = 1.0F - sinpsi * sinpsi;

        float cospsi = tmp2 * rsqrt(tmp2);

        float a2d_y = pars.ra * cosphi;
        float b2d_x = -pars.rc * cospsi;
        float t1    = -pars.rb * cosphi;
        float t2    = pars.rc * sinpsi * sinphi;
        float b2d_y = t1 - t2;
        float c2d_y = t1 + t2;

        /*     --- Step3  al,be,ga            --- */
        float alpha  = b2d_x * (b0d.x - c0d.x) + b0d.y * b2d_y + c0d.y * c2d_y;
        float beta   = b2d_x * (c0d.y - b0d.y) + b0d.x * b2d_y + c0d.x * c2d_y;
        float gamma  = b0d.x * b1d.y - b1d.x * b0d.y + c0d.x * c1d.y - c1d.x * c0d.y;
        float al2be2 = alpha * alpha + beta * beta;
        tmp2         = (al2be2 - gamma * gamma);
        float sinthe = (alpha * gamma - beta * tmp2 * rsqrt(tmp2)) * rsqrt(al2be2 * al2be2);

        /*  --- Step4  A3' --- */
        tmp2         = 1.0F - sinthe * sinthe;
        float costhe = tmp2 * rsqrt(tmp2);

        float3 a3d, b3d, c3d;

        a3d.x = -a2d_y * sinthe;
        a3d.y = a2d_y * costhe;
        a3d.z = a1d_z;
        b3d.x = b2d_x * costhe - b2d_y * sinthe;
        b3d.y = b2d_x * sinthe + b2d_y * costhe;
        b3d.z = b1d.z;
        c3d.x = -b2d_x * costhe - c2d_y * sinthe;
        c3d.y = -b2d_x * sinthe + c2d_y * costhe;
        c3d.z = c1d.z;

        /*    --- Step5  A3 --- */
        float3 a3, b3, c3;

        a3.x = trns1.x * a3d.x + trns1.y * a3d.y + trns1.z * a3d.z;
        a3.y = trns2.x * a3d.x + trns2.y * a3d.y + trns2.z * a3d.z;
        a3.z = trns3.x * a3d.x + trns3.y * a3d.y + trns3.z * a3d.z;

        b3.x = trns1.x * b3d.x + trns1.y * b3d.y + trns1.z * b3d.z;
        b3.y = trns2.x * b3d.x + trns2.y * b3d.y + trns2.z * b3d.z;
        b3.z = trns3.x * b3d.x + trns3.y * b3d.y + trns3.z * b3d.z;

        c3.x = trns1.x * c3d.x + trns1.y * c3d.y + trns1.z * c3d.z;
        c3.y = trns2.x * c3d.x + trns2.y * c3d.y + trns2.z * c3d.z;
        c3.z = trns3.x * c3d.x + trns3.y * c3d.y + trns3.z * c3d.z;


        /* Compute and store the corrected new coordinate */
        const float3 dxOw1 = a3 - a1;
        const float3 dxHw2 = b3 - b1;
        const float3 dxHw3 = c3 - c1;

        gm_xprime[indices.ow1] = xprime_ow1 + dxOw1;
        gm_xprime[indices.hw2] = xprime_hw2 + dxHw2;
        gm_xprime[indices.hw3] = xprime_hw3 + dxHw3;

        if (updateVelocities)
        {
            float3 v_ow1 = gm_v[indices.ow1];
            float3 v_hw2 = gm_v[indices.hw2];
            float3 v_hw3 = gm_v[indices.hw3];

            /* Add the position correction divided by dt to the velocity */
            v_ow1 = dxOw1 * invdt + v_ow1;
            v_hw2 = dxHw2 * invdt + v_hw2;
            v_hw3 = dxHw3 * invdt + v_hw3;

            gm_v[indices.ow1] = v_ow1;
            gm_v[indices.hw2] = v_hw2;
            gm_v[indices.hw3] = v_hw3;
        }

        if (computeVirial)
        {
            float3 mdb = pars.mH * dxHw2;
            float3 mdc = pars.mH * dxHw3;
            float3 mdo = pars.mO * dxOw1 + mdb + mdc;

            sm_threadVirial[0 * blockDim.x + threadIdx.x] =
                    -(x_ow1.x * mdo.x + dist21.x * mdb.x + dist31.x * mdc.x);
            sm_threadVirial[1 * blockDim.x + threadIdx.x] =
                    -(x_ow1.x * mdo.y + dist21.x * mdb.y + dist31.x * mdc.y);
            sm_threadVirial[2 * blockDim.x + threadIdx.x] =
                    -(x_ow1.x * mdo.z + dist21.x * mdb.z + dist31.x * mdc.z);
            sm_threadVirial[3 * blockDim.x + threadIdx.x] =
                    -(x_ow1.y * mdo.y + dist21.y * mdb.y + dist31.y * mdc.y);
            sm_threadVirial[4 * blockDim.x + threadIdx.x] =
                    -(x_ow1.y * mdo.z + dist21.y * mdb.z + dist31.y * mdc.z);
            sm_threadVirial[5 * blockDim.x + threadIdx.x] =
                    -(x_ow1.z * mdo.z + dist21.z * mdb.z + dist31.z * mdc.z);
        }
    }
    else
    {
        // Filling data for dummy threads with zeroes
        if (computeVirial)
        {
            for (int d = 0; d < 6; d++)
            {
                sm_threadVirial[d * blockDim.x + threadIdx.x] = 0.0F;
            }
        }
    }
    // Basic reduction for the values inside single thread block
    // TODO what follows should be separated out as a standard virial reduction subroutine
    if (computeVirial)
    {
        // This is to ensure that all threads saved the data before reduction starts
        __syncthreads();
        // This casts unsigned into signed integers to avoid clang warnings
        int tib       = static_cast<int>(threadIdx.x);
        int blockSize = static_cast<int>(blockDim.x);
        // Reduce up to one virial per thread block
        // All blocks are divided by half, the first half of threads sums
        // two virials. Then the first half is divided by two and the first half
        // of it sums two values... The procedure continues until only one thread left.
        // Only works if the threads per blocks is a power of two.
        for (int divideBy = 2; divideBy <= blockSize; divideBy *= 2)
        {
            int dividedAt = blockSize / divideBy;
            if (tib < dividedAt)
            {
                for (int d = 0; d < 6; d++)
                {
                    sm_threadVirial[d * blockSize + tib] +=
                            sm_threadVirial[d * blockSize + (tib + dividedAt)];
                }
            }
            if (dividedAt > warpSize / 2)
            {
                __syncthreads();
            }
            else
            {
                __syncwarp();
            }
        }
        // First 6 threads in the block add the 6 components of virial to the global memory address
        if (tib < 6)
        {
            atomicAdd(&(gm_virialScaled[tib]), sm_threadVirial[tib * blockSize]);
        }
    }
}

/*! \brief Select templated kernel.
 *
 * Returns pointer to a CUDA kernel based on provided booleans.
 *
 * \param[in] updateVelocities  If the velocities should be constrained.
 * \param[in] computeVirial     If virial should be updated.
 *
 * \return                      Pointer to CUDA kernel
 */
inline auto getSettleKernelPtr(const bool updateVelocities, const bool computeVirial)
{

    auto kernelPtr = settle_kernel<true, true>;
    if (updateVelocities && computeVirial)
    {
        kernelPtr = settle_kernel<true, true>;
    }
    else if (updateVelocities && !computeVirial)
    {
        kernelPtr = settle_kernel<true, false>;
    }
    else if (!updateVelocities && computeVirial)
    {
        kernelPtr = settle_kernel<false, true>;
    }
    else if (!updateVelocities && !computeVirial)
    {
        kernelPtr = settle_kernel<false, false>;
    }
    return kernelPtr;
}

void launchSettleGpuKernel(const int                          numSettles,
                           const DeviceBuffer<WaterMolecule>& d_atomIds,
                           const SettleParameters&            settleParameters,
                           const DeviceBuffer<Float3>&        d_x,
                           DeviceBuffer<Float3>               d_xp,
                           const bool                         updateVelocities,
                           DeviceBuffer<Float3>               d_v,
                           const real                         invdt,
                           const bool                         computeVirial,
                           DeviceBuffer<float>                virialScaled,
                           const PbcAiuc&                     pbcAiuc,
                           const DeviceStream&                deviceStream)
{
    static_assert(
            gmx::isPowerOfTwo(sc_threadsPerBlock),
            "Number of threads per block should be a power of two in order for reduction to work.");

    auto kernelPtr = getSettleKernelPtr(updateVelocities, computeVirial);

    KernelLaunchConfig config;
    config.blockSize[0] = sc_threadsPerBlock;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;
    config.gridSize[0]  = divideRoundUp(numSettles, sc_threadsPerBlock);
    config.gridSize[1]  = 1;
    config.gridSize[2]  = 1;

    // Shared memory is only used for virial reduction
    if (computeVirial)
    {
        config.sharedMemorySize = sc_threadsPerBlock * 6 * sizeof(float);
    }
    else
    {
        config.sharedMemorySize = 0;
    }

    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      config,
                                                      &numSettles,
                                                      &d_atomIds,
                                                      &settleParameters,
                                                      asFloat3Pointer(&d_x),
                                                      asFloat3Pointer(&d_xp),
                                                      &invdt,
                                                      asFloat3Pointer(&d_v),
                                                      &virialScaled,
                                                      &pbcAiuc);

    launchGpuKernel(kernelPtr,
                    config,
                    deviceStream,
                    nullptr,
                    "settle_kernel<updateVelocities, computeVirial>",
                    kernelArgs);
}

} // namespace gmx
