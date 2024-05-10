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
 * \brief Implements LINCS kernels using CUDA
 *
 * This file contains CUDA kernels of LINCS constraints algorithm.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "lincs_gpu_internal.h"

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/mdlib/lincs_gpu.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"

namespace gmx
{

//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Main kernel for LINCS constraints.
 *
 * See Hess et al., J. Comput. Chem. 18: 1463-1472 (1997) for the description of the algorithm.
 *
 * In GPU version, one thread is responsible for all computations for one constraint. The blocks are
 * filled in a way that no constraint is coupled to the constraint from the next block. This is achieved
 * by moving active threads to the next block, if the correspondent group of coupled constraints is to big
 * to fit the current thread block. This may leave some 'dummy' threads in the end of the thread block, i.e.
 * threads that are not required to do actual work. Since constraints from different blocks are not coupled,
 * there is no need to synchronize across the device. However, extensive communication in a thread block
 * are still needed.
 *
 * \todo Reduce synchronization overhead. Some ideas are:
 *        1. Consider going to warp-level synchronization for the coupled constraints.
 *        2. Move more data to local/shared memory and try to get rid of atomic operations (at least on
 *           the device level).
 *        3. Use analytical solution for matrix A inversion.
 *        4. Introduce mapping of thread id to both single constraint and single atom, thus designating
 *           Nth threads to deal with Nat <= Nth coupled atoms and Nc <= Nth coupled constraints.
 *       See Issue #2885 for details (https://gitlab.com/gromacs/gromacs/-/issues/2885)
 * \todo The use of __restrict__  for gm_xp and gm_v causes failure, probably because of the atomic
         operations. Investigate this issue further.
 *
 * \param[in,out] kernelParams  All parameters and pointers for the kernel condensed in single struct.
 * \param[in]     invdt         Inverse timestep (needed to update velocities).
 */
template<bool updateVelocities, bool computeVirial>
__launch_bounds__(c_maxThreadsPerBlock) __global__
        void lincs_kernel(LincsGpuKernelParameters kernelParams,
                          const float3* __restrict__ gm_x,
                          float3*     gm_xp,
                          float3*     gm_v,
                          const float invdt)
{
    const PbcAiuc pbcAiuc                                 = kernelParams.pbcAiuc;
    const int     numConstraintsThreads                   = kernelParams.numConstraintsThreads;
    const int     numIterations                           = kernelParams.numIterations;
    const int     expansionOrder                          = kernelParams.expansionOrder;
    const AtomPair* __restrict__ gm_constraints           = kernelParams.d_constraints;
    const float* __restrict__ gm_constraintsTargetLengths = kernelParams.d_constraintsTargetLengths;
    const int* __restrict__ gm_coupledConstraintsCounts   = kernelParams.d_coupledConstraintsCounts;
    const int* __restrict__ gm_coupledConstraintsIndices = kernelParams.d_coupledConstraintsIndices;
    const float* __restrict__ gm_massFactors             = kernelParams.d_massFactors;
    float* __restrict__ gm_matrixA                       = kernelParams.d_matrixA;
    const float* __restrict__ gm_inverseMasses           = kernelParams.d_inverseMasses;
    float* __restrict__ gm_virialScaled                  = kernelParams.d_virialScaled;

    const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

    // numConstraintsThreads should be a integer multiple of blockSize (numConstraintsThreads = numBlocks*blockSize).
    // This is to ensure proper synchronizations and reduction. All array are padded to the required size.
    assert(threadIndex < numConstraintsThreads);

    // Vectors connecting constrained atoms before algorithm was applied.
    // Needed to construct constrain matrix A
    extern __shared__ float3 sm_r[];

    AtomPair pair = gm_constraints[threadIndex];
    int      i    = pair.i;
    int      j    = pair.j;

    // Mass-scaled Lagrange multiplier
    float lagrangeScaled = 0.0F;

    float targetLength;
    float inverseMassi;
    float inverseMassj;
    float sqrtReducedMass;

    float3 xi;
    float3 xj;
    float3 rc;

    // i == -1 indicates dummy constraint at the end of the thread block.
    bool isDummyThread = (i == -1);

    // Everything computed for these dummies will be equal to zero
    if (isDummyThread)
    {
        targetLength    = 0.0F;
        inverseMassi    = 0.0F;
        inverseMassj    = 0.0F;
        sqrtReducedMass = 0.0F;

        xi = make_float3(0.0F, 0.0F, 0.0F);
        xj = make_float3(0.0F, 0.0F, 0.0F);
        rc = make_float3(0.0F, 0.0F, 0.0F);
    }
    else
    {
        // Collecting data
        targetLength    = gm_constraintsTargetLengths[threadIndex];
        inverseMassi    = gm_inverseMasses[i];
        inverseMassj    = gm_inverseMasses[j];
        sqrtReducedMass = rsqrt(inverseMassi + inverseMassj);

        xi = gm_x[i];
        xj = gm_x[j];

        float3 dx = pbcDxAiuc(pbcAiuc, xi, xj);

        float rlen = rsqrtf(dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
        rc         = rlen * dx;
    }

    sm_r[threadIdx.x] = rc;
    // Make sure that all r's are saved into shared memory
    // before they are accessed in the loop below
    __syncthreads();

    /*
     * Constructing LINCS matrix (A)
     */

    // Only non-zero values are saved (for coupled constraints)
    int coupledConstraintsCount = gm_coupledConstraintsCounts[threadIndex];
    for (int n = 0; n < coupledConstraintsCount; n++)
    {
        int index = n * numConstraintsThreads + threadIndex;
        int c1    = gm_coupledConstraintsIndices[index];

        float3 rc1        = sm_r[c1];
        gm_matrixA[index] = gm_massFactors[index] * (rc.x * rc1.x + rc.y * rc1.y + rc.z * rc1.z);
    }

    // Skipping in dummy threads
    if (!isDummyThread)
    {
        xi = gm_xp[i];
        xj = gm_xp[j];
    }

    float3 dx = pbcDxAiuc(pbcAiuc, xi, xj);

    float sol = sqrtReducedMass * ((rc.x * dx.x + rc.y * dx.y + rc.z * dx.z) - targetLength);

    /*
     *  Inverse matrix using a set of expansionOrder matrix multiplications
     */

    // Make sure that we don't overwrite the sm_r[..] array.
    __syncthreads();

    // This will use the same memory space as sm_r, which is no longer needed.
    extern __shared__ float sm_rhs[];
    // Save current right-hand-side vector in the shared memory
    sm_rhs[threadIdx.x] = sol;

    for (int rec = 0; rec < expansionOrder; rec++)
    {
        // Making sure that all sm_rhs are saved before they are accessed in a loop below
        __syncthreads();
        float mvb = 0.0F;

        for (int n = 0; n < coupledConstraintsCount; n++)
        {
            int index = n * numConstraintsThreads + threadIndex;
            int c1    = gm_coupledConstraintsIndices[index];
            // Convolute current right-hand-side with A
            // Different, non overlapping parts of sm_rhs[..] are read during odd and even iterations
            mvb = mvb + gm_matrixA[index] * sm_rhs[c1 + blockDim.x * (rec % 2)];
        }
        // 'Switch' rhs vectors, save current result
        // These values will be accessed in the loop above during the next iteration.
        sm_rhs[threadIdx.x + blockDim.x * ((rec + 1) % 2)] = mvb;
        sol                                                = sol + mvb;
    }

    // Current mass-scaled Lagrange multipliers
    lagrangeScaled = sqrtReducedMass * sol;

    // Save updated coordinates before correction for the rotational lengthening
    float3 tmp = rc * lagrangeScaled;

    // Writing for all but dummy constraints
    if (!isDummyThread)
    {
        atomicAdd(&gm_xp[i], -tmp * inverseMassi);
        atomicAdd(&gm_xp[j], tmp * inverseMassj);
    }

    /*
     *  Correction for centripetal effects
     */
    for (int iter = 0; iter < numIterations; iter++)
    {
        // Make sure that all xp's are saved: atomic operation calls before are
        // communicating current xp[..] values across thread block.
        __syncthreads();

        if (!isDummyThread)
        {
            xi = gm_xp[i];
            xj = gm_xp[j];
        }

        float3 dx = pbcDxAiuc(pbcAiuc, xi, xj);

        float len2  = targetLength * targetLength;
        float dlen2 = 2.0F * len2 - norm2(dx);

        // TODO A little bit more effective but slightly less readable version of the below would be:
        //      float proj = sqrtReducedMass*(targetLength - (dlen2 > 0.0f ? 1.0f : 0.0f)*dlen2*rsqrt(dlen2));
        float proj;
        if (dlen2 > 0.0F)
        {
            proj = sqrtReducedMass * (targetLength - dlen2 * rsqrt(dlen2));
        }
        else
        {
            proj = sqrtReducedMass * targetLength;
        }

        sm_rhs[threadIdx.x] = proj;
        float sol           = proj;

        /*
         * Same matrix inversion as above is used for updated data
         */
        for (int rec = 0; rec < expansionOrder; rec++)
        {
            // Make sure that all elements of rhs are saved into shared memory
            __syncthreads();
            float mvb = 0;

            for (int n = 0; n < coupledConstraintsCount; n++)
            {
                int index = n * numConstraintsThreads + threadIndex;
                int c1    = gm_coupledConstraintsIndices[index];

                mvb = mvb + gm_matrixA[index] * sm_rhs[c1 + blockDim.x * (rec % 2)];
            }
            sm_rhs[threadIdx.x + blockDim.x * ((rec + 1) % 2)] = mvb;
            sol                                                = sol + mvb;
        }

        // Add corrections to Lagrange multipliers
        float sqrtmu_sol = sqrtReducedMass * sol;
        lagrangeScaled += sqrtmu_sol;

        // Save updated coordinates for the next iteration
        // Dummy constraints are skipped
        if (!isDummyThread)
        {
            float3 tmp = rc * sqrtmu_sol;
            atomicAdd(&gm_xp[i], -tmp * inverseMassi);
            atomicAdd(&gm_xp[j], tmp * inverseMassj);
        }
    }

    // Updating particle velocities for all but dummy threads
    if (updateVelocities && !isDummyThread)
    {
        float3 tmp = rc * invdt * lagrangeScaled;
        atomicAdd(&gm_v[i], -tmp * inverseMassi);
        atomicAdd(&gm_v[j], tmp * inverseMassj);
    }


    if (computeVirial)
    {
        // Virial is computed from Lagrange multiplier (lagrangeScaled), target constrain length
        // (targetLength) and the normalized vector connecting constrained atoms before
        // the algorithm was applied (rc). The evaluation of virial in each thread is
        // followed by basic reduction for the values inside single thread block.
        // Then, the values are reduced across grid by atomicAdd(...).
        //
        // TODO Shuffle reduction.
        // TODO Should be unified and/or done once when virial is actually needed.
        // TODO Recursive version that removes atomicAdd(...)'s entirely is needed. Ideally,
        //      one that works for any datatype.

        // Save virial for each thread into the shared memory. Tensor is symmetrical, hence only
        // 6 values are saved. Dummy threads will have zeroes in their virial: targetLength,
        // lagrangeScaled and rc are all set to zero for them in the beginning of the kernel.
        // The sm_threadVirial[..] will overlap with the sm_r[..] and sm_rhs[..], but the latter
        // two are no longer in use, which we make sure by waiting for all threads in block.
        __syncthreads();
        extern __shared__ float sm_threadVirial[];
        float                   mult                  = targetLength * lagrangeScaled;
        sm_threadVirial[0 * blockDim.x + threadIdx.x] = mult * rc.x * rc.x;
        sm_threadVirial[1 * blockDim.x + threadIdx.x] = mult * rc.x * rc.y;
        sm_threadVirial[2 * blockDim.x + threadIdx.x] = mult * rc.x * rc.z;
        sm_threadVirial[3 * blockDim.x + threadIdx.x] = mult * rc.y * rc.y;
        sm_threadVirial[4 * blockDim.x + threadIdx.x] = mult * rc.y * rc.z;
        sm_threadVirial[5 * blockDim.x + threadIdx.x] = mult * rc.z * rc.z;

        __syncthreads();

        // Reduce up to one virial per thread block. All blocks are divided by half, the first
        // half of threads sums two virials. Then the first half is divided by two and the first
        // half of it sums two values. This procedure is repeated until only one thread is left.
        // Only works if the threads per blocks is a power of two (hence static_assert
        // in the beginning of the kernel).
        for (int divideBy = 2; divideBy <= static_cast<int>(blockDim.x); divideBy *= 2)
        {
            int dividedAt = blockDim.x / divideBy;
            if (static_cast<int>(threadIdx.x) < dividedAt)
            {
                for (int d = 0; d < 6; d++)
                {
                    sm_threadVirial[d * blockDim.x + threadIdx.x] +=
                            sm_threadVirial[d * blockDim.x + (threadIdx.x + dividedAt)];
                }
            }
            // Syncronize if not within one warp
            if (dividedAt > warpSize / 2)
            {
                __syncthreads();
            }
            else
            {
                __syncwarp();
            }
        }
        // First 6 threads in the block add the results of 6 tensor components to the global memory address.
        if (threadIdx.x < 6)
        {
            atomicAdd(&(gm_virialScaled[threadIdx.x]), sm_threadVirial[threadIdx.x * blockDim.x]);
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
inline auto getLincsKernelPtr(const bool updateVelocities, const bool computeVirial)
{

    auto kernelPtr = lincs_kernel<true, true>;
    if (updateVelocities && computeVirial)
    {
        kernelPtr = lincs_kernel<true, true>;
    }
    else if (updateVelocities && !computeVirial)
    {
        kernelPtr = lincs_kernel<true, false>;
    }
    else if (!updateVelocities && computeVirial)
    {
        kernelPtr = lincs_kernel<false, true>;
    }
    else if (!updateVelocities && !computeVirial)
    {
        kernelPtr = lincs_kernel<false, false>;
    }
    return kernelPtr;
}

void launchLincsGpuKernel(LincsGpuKernelParameters*   kernelParams,
                          const DeviceBuffer<Float3>& d_x,
                          DeviceBuffer<Float3>        d_xp,
                          const bool                  updateVelocities,
                          DeviceBuffer<Float3>        d_v,
                          const real                  invdt,
                          const bool                  computeVirial,
                          const DeviceStream&         deviceStream)
{

    auto kernelPtr = getLincsKernelPtr(updateVelocities, computeVirial);

    KernelLaunchConfig config;
    config.blockSize[0] = c_threadsPerBlock;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;
    config.gridSize[0]  = divideRoundUp(kernelParams->numConstraintsThreads, c_threadsPerBlock);
    config.gridSize[1]  = 1;
    config.gridSize[2]  = 1;

    // Shared memory is used to store:
    // -- Current coordinates (3 floats per thread)
    // -- Right-hand-sides for matrix inversion (2 floats per thread)
    // -- Virial tensor components (6 floats per thread)
    // Since none of these three are needed simultaneously, they can be saved at the same shared memory address
    // (i.e. correspondent arrays are intentionally overlapped in address space). Consequently, only
    // max{3, 2, 6} = 6 floats per thread are needed in case virial is computed, or max{3, 2} = 3 if not.
    if (computeVirial)
    {
        config.sharedMemorySize = c_threadsPerBlock * 6 * sizeof(float);
    }
    else
    {
        config.sharedMemorySize = c_threadsPerBlock * 3 * sizeof(float);
    }

    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      config,
                                                      kernelParams,
                                                      asFloat3Pointer(&d_x),
                                                      asFloat3Pointer(&d_xp),
                                                      asFloat3Pointer(&d_v),
                                                      &invdt);

    launchGpuKernel(kernelPtr,
                    config,
                    deviceStream,
                    nullptr,
                    "lincs_kernel<updateVelocities, computeVirial>",
                    kernelArgs);
}

} // namespace gmx
