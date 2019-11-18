/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Implements SETTLE using CUDA
 *
 * This file contains implementation of SETTLE constraints algorithm
 * using CUDA, including class initialization, data-structures management
 * and GPU kernel.
 *
 * \note Management of CUDA stream and periodic boundary should be unified with LINCS
 *       and removed from here once constraints are fully integrated with update module.
 * \todo Reconsider naming to use "gpu" suffix instead of "cuda".
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settle_cuda.cuh"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"

namespace gmx
{

//! Number of CUDA threads in a block
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

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
 * \param [in]      pbcAiuc          Periodic boundary conditions data.
 * \param [in]      invdt            Reciprocal timestep.
 * \param [in]      gm_v             Velocities of the particles.
 * \param [in]      gm_virialScaled  Virial tensor.
 */
template<bool updateVelocities, bool computeVirial>
__launch_bounds__(c_maxThreadsPerBlock) __global__
        void settle_kernel(const int numSettles,
                           const int3* __restrict__ gm_settles,
                           const SettleParameters pars,
                           const float3* __restrict__ gm_x,
                           float3* __restrict__ gm_xprime,
                           const PbcAiuc pbcAiuc,
                           float         invdt,
                           float3* __restrict__ gm_v,
                           float* __restrict__ gm_virialScaled)
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
        int3 indices = gm_settles[tid];

        float3 x_ow1 = gm_x[indices.x];
        float3 x_hw2 = gm_x[indices.y];
        float3 x_hw3 = gm_x[indices.z];

        float3 xprime_ow1 = gm_xprime[indices.x];
        float3 xprime_hw2 = gm_xprime[indices.y];
        float3 xprime_hw3 = gm_xprime[indices.z];

        float3 dist21 = pbcDxAiuc(pbcAiuc, x_hw2, x_ow1);
        float3 dist31 = pbcDxAiuc(pbcAiuc, x_hw3, x_ow1);
        float3 doh2   = pbcDxAiuc(pbcAiuc, xprime_hw2, xprime_ow1);

        float3 sh_hw2 = xprime_hw2 - (xprime_ow1 + doh2);
        xprime_hw2    = xprime_hw2 - sh_hw2;

        float3 doh3 = pbcDxAiuc(pbcAiuc, xprime_hw3, xprime_ow1);

        float3 sh_hw3 = xprime_hw3 - (xprime_ow1 + doh3);
        xprime_hw3    = xprime_hw3 - sh_hw3;

        float3 a1  = (-doh2 - doh3) * pars.wh;
        float3 com = xprime_ow1 - a1;

        float3 b1 = xprime_hw2 - com;

        float3 c1 = xprime_hw3 - com;

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
        float tmp2   = 1.0f - sinphi * sinphi;

        if (almost_zero > tmp2)
        {
            tmp2 = almost_zero;
        }

        float tmp    = rsqrt(tmp2);
        float cosphi = tmp2 * tmp;
        float sinpsi = (b1d.z - c1d.z) * pars.irc2 * tmp;
        tmp2         = 1.0f - sinpsi * sinpsi;

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
        tmp2         = 1.0f - sinthe * sinthe;
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
        xprime_ow1 = com + a3;
        xprime_hw2 = com + b3 + sh_hw2;
        xprime_hw3 = com + c3 + sh_hw3;

        gm_xprime[indices.x] = xprime_ow1;
        gm_xprime[indices.y] = xprime_hw2;
        gm_xprime[indices.z] = xprime_hw3;


        if (updateVelocities || computeVirial)
        {

            float3 da = a3 - a1;
            float3 db = b3 - b1;
            float3 dc = c3 - c1;

            if (updateVelocities)
            {

                float3 v_ow1 = gm_v[indices.x];
                float3 v_hw2 = gm_v[indices.y];
                float3 v_hw3 = gm_v[indices.z];

                /* Add the position correction divided by dt to the velocity */
                v_ow1 = da * invdt + v_ow1;
                v_hw2 = db * invdt + v_hw2;
                v_hw3 = dc * invdt + v_hw3;

                gm_v[indices.x] = v_ow1;
                gm_v[indices.y] = v_hw2;
                gm_v[indices.z] = v_hw3;
            }

            if (computeVirial)
            {

                float3 mdb = pars.mH * db;
                float3 mdc = pars.mH * dc;
                float3 mdo = pars.mO * da + mdb + mdc;

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
    }
    else
    {
        // Filling data for dummy threads with zeroes
        if (computeVirial)
        {
            for (int d = 0; d < 6; d++)
            {
                sm_threadVirial[d * blockDim.x + threadIdx.x] = 0.0f;
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
        }
        // First 6 threads in the block add the 6 components of virial to the global memory address
        if (tib < 6)
        {
            atomicAdd(&(gm_virialScaled[tib]), sm_threadVirial[tib * blockSize]);
        }
    }

    return;
}

/*! \brief Select templated kernel.
 *
 * Returns pointer to a CUDA kernel based on provided booleans.
 *
 * \param[in] updateVelocities  If the velocities should be constrained.
 * \param[in] bCalcVir          If virial should be updated.
 *
 * \retrun                      Pointer to CUDA kernel
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

void SettleCuda::apply(const float3* d_x,
                       float3*       d_xp,
                       const bool    updateVelocities,
                       float3*       d_v,
                       const real    invdt,
                       const bool    computeVirial,
                       tensor        virialScaled)
{

    ensureNoPendingCudaError("In CUDA version SETTLE");

    // Early exit if no settles
    if (numSettles_ == 0)
    {
        return;
    }

    if (computeVirial)
    {
        // Fill with zeros so the values can be reduced to it
        // Only 6 values are needed because virial is symmetrical
        clearDeviceBufferAsync(&d_virialScaled_, 0, 6, commandStream_);
    }

    auto kernelPtr = getSettleKernelPtr(updateVelocities, computeVirial);

    KernelLaunchConfig config;
    config.blockSize[0] = c_threadsPerBlock;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;
    config.gridSize[0]  = (numSettles_ + c_threadsPerBlock - 1) / c_threadsPerBlock;
    config.gridSize[1]  = 1;
    config.gridSize[2]  = 1;
    // Shared memory is only used for virial reduction
    if (computeVirial)
    {
        config.sharedMemorySize = c_threadsPerBlock * 6 * sizeof(float);
    }
    else
    {
        config.sharedMemorySize = 0;
    }
    config.stream = commandStream_;

    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config, &numSettles_, &d_atomIds_,
                                                      &settleParameters_, &d_x, &d_xp, &pbcAiuc_,
                                                      &invdt, &d_v, &d_virialScaled_);

    launchGpuKernel(kernelPtr, config, nullptr, "settle_kernel<updateVelocities, computeVirial>", kernelArgs);

    if (computeVirial)
    {
        copyFromDeviceBuffer(h_virialScaled_.data(), &d_virialScaled_, 0, 6, commandStream_,
                             GpuApiCallBehavior::Sync, nullptr);

        // Mapping [XX, XY, XZ, YY, YZ, ZZ] internal format to a tensor object
        virialScaled[XX][XX] += h_virialScaled_[0];
        virialScaled[XX][YY] += h_virialScaled_[1];
        virialScaled[XX][ZZ] += h_virialScaled_[2];

        virialScaled[YY][XX] += h_virialScaled_[1];
        virialScaled[YY][YY] += h_virialScaled_[3];
        virialScaled[YY][ZZ] += h_virialScaled_[4];

        virialScaled[ZZ][XX] += h_virialScaled_[2];
        virialScaled[ZZ][YY] += h_virialScaled_[4];
        virialScaled[ZZ][ZZ] += h_virialScaled_[5];
    }

    return;
}

SettleCuda::SettleCuda(const gmx_mtop_t& mtop, CommandStream commandStream) :
    commandStream_(commandStream)
{
    static_assert(sizeof(real) == sizeof(float),
                  "Real numbers should be in single precision in GPU code.");
    static_assert(
            c_threadsPerBlock > 0 && ((c_threadsPerBlock & (c_threadsPerBlock - 1)) == 0),
            "Number of threads per block should be a power of two in order for reduction to work.");

    // This is to prevent the assertion failure for the systems without water
    int totalSettles = 0;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype.at(mt).ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        totalSettles += iatoms.size() / nral1;
    }
    if (totalSettles == 0)
    {
        return;
    }

    // TODO This should be lifted to a separate subroutine that gets the values of Oxygen and
    // Hydrogen masses, checks if they are consistent across the topology and if there is no more
    // than two values for each mass if the free energy perturbation is enabled. In later case,
    // masses may need to be updated on a regular basis (i.e. in set(...) method).
    // TODO Do the checks for FEP
    real mO = -1.0;
    real mH = -1.0;

    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype.at(mt).ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        for (unsigned i = 0; i < iatoms.size() / nral1; i++)
        {
            int3 settler;
            settler.x  = iatoms[i * nral1 + 1]; // Oxygen index
            settler.y  = iatoms[i * nral1 + 2]; // First hydrogen index
            settler.z  = iatoms[i * nral1 + 3]; // Second hydrogen index
            t_atom ow1 = mtop.moltype.at(mt).atoms.atom[settler.x];
            t_atom hw2 = mtop.moltype.at(mt).atoms.atom[settler.y];
            t_atom hw3 = mtop.moltype.at(mt).atoms.atom[settler.z];

            if (mO < 0)
            {
                mO = ow1.m;
            }
            if (mH < 0)
            {
                mH = hw2.m;
            }
            GMX_RELEASE_ASSERT(mO == ow1.m,
                               "Topology has different values for oxygen mass. Should be identical "
                               "in order to use SETTLE.");
            GMX_RELEASE_ASSERT(hw2.m == hw3.m && hw2.m == mH,
                               "Topology has different values for hydrogen mass. Should be "
                               "identical in order to use SETTLE.");
        }
    }

    GMX_RELEASE_ASSERT(mO > 0, "Could not find oxygen mass in the topology. Needed in SETTLE.");
    GMX_RELEASE_ASSERT(mH > 0, "Could not find hydrogen mass in the topology. Needed in SETTLE.");

    // TODO Very similar to SETTLE initialization on CPU. Should be lifted to a separate method
    // (one that gets dOH and dHH values and checks them for consistency)
    int settle_type = -1;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int       nral1           = 1 + NRAL(F_SETTLE);
        InteractionList interactionList = mtop.moltype.at(mt).ilist[F_SETTLE];
        for (int i = 0; i < interactionList.size(); i += nral1)
        {
            if (settle_type == -1)
            {
                settle_type = interactionList.iatoms[i];
            }
            else if (interactionList.iatoms[i] != settle_type)
            {
                gmx_fatal(FARGS,
                          "The [molecules] section of your topology specifies more than one block "
                          "of\n"
                          "a [moleculetype] with a [settles] block. Only one such is allowed.\n"
                          "If you are trying to partition your solvent into different *groups*\n"
                          "(e.g. for freezing, T-coupling, etc.), you are using the wrong "
                          "approach. Index\n"
                          "files specify groups. Otherwise, you may wish to change the least-used\n"
                          "block of molecules with SETTLE constraints into 3 normal constraints.");
            }
        }
    }

    GMX_RELEASE_ASSERT(settle_type >= 0, "settle_init called without settles");

    real dOH = mtop.ffparams.iparams[settle_type].settle.doh;
    real dHH = mtop.ffparams.iparams[settle_type].settle.dhh;

    initSettleParameters(&settleParameters_, mO, mH, dOH, dHH);

    allocateDeviceBuffer(&d_virialScaled_, 6, nullptr);
    h_virialScaled_.resize(6);
}

SettleCuda::~SettleCuda()
{
    // Early exit if there is no settles
    if (numSettles_ == 0)
    {
        return;
    }
    freeDeviceBuffer(&d_virialScaled_);
    if (numAtomIdsAlloc_ > 0)
    {
        freeDeviceBuffer(&d_atomIds_);
    }
}

void SettleCuda::set(const t_idef& idef, const t_mdatoms gmx_unused& md)
{
    const int nral1     = 1 + NRAL(F_SETTLE);
    t_ilist   il_settle = idef.il[F_SETTLE];
    t_iatom*  iatoms    = il_settle.iatoms;
    numSettles_         = il_settle.nr / nral1;

    reallocateDeviceBuffer(&d_atomIds_, numSettles_, &numAtomIds_, &numAtomIdsAlloc_, nullptr);
    h_atomIds_.resize(numSettles_);
    for (int i = 0; i < numSettles_; i++)
    {
        int3 settler;
        settler.x        = iatoms[i * nral1 + 1]; // Oxygen index
        settler.y        = iatoms[i * nral1 + 2]; // First hydrogen index
        settler.z        = iatoms[i * nral1 + 3]; // Second hydrogen index
        h_atomIds_.at(i) = settler;
    }
    copyToDeviceBuffer(&d_atomIds_, h_atomIds_.data(), 0, numSettles_, commandStream_,
                       GpuApiCallBehavior::Sync, nullptr);
}

void SettleCuda::setPbc(const t_pbc* pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
}

} // namespace gmx
