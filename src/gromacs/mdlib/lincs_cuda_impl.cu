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
 * \brief Implements LINCS using CUDA
 *
 * This file contains implementation of LINCS constraints algorithm
 * using CUDA, including class initialization, data-structures management
 * and GPU kernel.
 *
 * \note Management of coordinates, velocities, CUDA stream and periodic boundary exists here as a temporary
 *       scaffolding to allow this class to be used as a stand-alone unit. The scaffolding is intended to be
 *       removed once constraints are integrated with update module.
 * \todo Reconsider naming, i.e. "cuda" suffics should be changed to "gpu".
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "lincs_cuda_impl.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs_cuda.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/topology/ifunc.h"

namespace gmx
{

//! Number of CUDA threads in a block
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Main kernel for LINCS constraints.
 *
 * See Hess et al., J. Comput. Chem. 18: 1463-1472 (1997) for the description of the algorithm.
 *
 * In CUDA version, one thread is responsible for all computations for one constraint. The blocks are
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
 *       See Redmine issue #2885 for details (https://redmine.gromacs.org/issues/2885)
 *
 * \param[in,out] kernelParams  All parameters and pointers for the kernel condensed in single struct.
 * \param[in]     invdt         Inverse timestep (needed to update velocities).
 */
template <bool updateVelocities, bool computeVirial>
__launch_bounds__(c_maxThreadsPerBlock)
__global__ void lincs_kernel(LincsCudaKernelParameters   kernelParams,
                             const float                 invdt)
{
    const PbcAiuc               pbcAiuc                     = kernelParams.pbcAiuc;
    const int                   numConstraintsThreads       = kernelParams.numConstraintsThreads;
    const int                   numIterations               = kernelParams.numIterations;
    const int                   expansionOrder              = kernelParams.expansionOrder;
    const float3* __restrict__  gm_x                        = kernelParams.d_x;
    float3*                     gm_xp                       = kernelParams.d_xp;
    const int2*   __restrict__  gm_constraints              = kernelParams.d_constraints;
    const float*  __restrict__  gm_constraintsTargetLengths = kernelParams.d_constraintsTargetLengths;
    const int*    __restrict__  gm_coupledConstraintsCounts = kernelParams.d_coupledConstraintsCounts;
    const int*    __restrict__  gm_coupledConstraintsIdxes  = kernelParams.d_coupledConstraintsIndices;
    const float*  __restrict__  gm_massFactors              = kernelParams.d_massFactors;
    float*  __restrict__        gm_matrixA                  = kernelParams.d_matrixA;
    const float*  __restrict__  gm_inverseMasses            = kernelParams.d_inverseMasses;
    float3*                     gm_v                        = kernelParams.d_v;
    float*  __restrict__        gm_virialScaled             = kernelParams.d_virialScaled;

    int threadIndex                                         = blockIdx.x*blockDim.x+threadIdx.x;

    // numConstraintsThreads should be a integer multiple of blockSize (numConstraintsThreads = numBlocks*blockSize).
    // This is to ensure proper synchronizations and reduction. All array are padded to the required size.
    assert(threadIndex < numConstraintsThreads);

    // Vectors connecting constrained atoms before algorithm was applied.
    // Needed to construct constrain matrix A
    extern __shared__ float3 sm_r[];

    int2                     pair = gm_constraints[threadIndex];
    int                      i    = pair.x;
    int                      j    = pair.y;

    // Mass-scaled Lagrange multiplier
    float  lagrangeScaled = 0.0f;

    float  targetLength;
    float  inverseMassi;
    float  inverseMassj;
    float  sqrtReducedMass;

    float3 xi;
    float3 xj;
    float3 rc;

    // i == -1 indicates dummy constraint at the end of the thread block.
    bool isDummyThread = (i == -1);

    // Everything computed for these dummies will be equal to zero
    if (isDummyThread)
    {
        targetLength    = 0.0f;
        inverseMassi    = 0.0f;
        inverseMassj    = 0.0f;
        sqrtReducedMass = 0.0f;

        xi = make_float3(0.0f, 0.0f, 0.0f);
        xj = make_float3(0.0f, 0.0f, 0.0f);
        rc = make_float3(0.0f, 0.0f, 0.0f);
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

        float3 dx   = pbcDxAiuc(pbcAiuc, xi, xj);

        float  rlen = rsqrtf(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);
        rc   = rlen*dx;
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
        int    index = n*numConstraintsThreads + threadIndex;
        int    c1    = gm_coupledConstraintsIdxes[index];

        float3 rc1 = sm_r[c1 - blockIdx.x*blockDim.x];
        gm_matrixA[index] = gm_massFactors[index]*(rc.x*rc1.x + rc.y*rc1.y + rc.z*rc1.z);
    }

    // Skipping in dummy threads
    if (!isDummyThread)
    {
        xi = gm_xp[i];
        xj = gm_xp[j];
    }

    float3 dx = pbcDxAiuc(pbcAiuc, xi, xj);

    float  sol = sqrtReducedMass*((rc.x*dx.x + rc.y*dx.y + rc.z*dx.z) - targetLength);

    /*
     *  Inverse matrix using a set of expansionOrder matrix multiplications
     */

    // This will use the same memory space as sm_r, which is no longer needed.
    extern __shared__ float  sm_rhs[];
    // Save current right-hand-side vector in the shared memory
    sm_rhs[threadIdx.x] = sol;

    for (int rec = 0; rec < expansionOrder; rec++)
    {
        // Making sure that all sm_rhs are saved before they are accessed in a loop below
        __syncthreads();
        float mvb = 0.0f;

        for (int n = 0; n < coupledConstraintsCount; n++)
        {
            int index = n*numConstraintsThreads + threadIndex;
            int c1    = gm_coupledConstraintsIdxes[index];
            // Convolute current right-hand-side with A
            // Different, non overlapping parts of sm_rhs[..] are read during odd and even iterations
            mvb = mvb + gm_matrixA[index]*sm_rhs[c1 - blockIdx.x*blockDim.x + blockDim.x*(rec % 2)];
        }
        // 'Switch' rhs vectors, save current result
        // These values will be accessed in the loop above during the next iteration.
        sm_rhs[threadIdx.x + blockDim.x*((rec + 1) % 2)] = mvb;
        sol  = sol + mvb;
    }

    // Current mass-scaled Lagrange multipliers
    lagrangeScaled = sqrtReducedMass*sol;

    // Save updated coordinates before correction for the rotational lengthening
    float3 tmp     = rc*lagrangeScaled;

    // Writing for all but dummy constraints
    if (!isDummyThread)
    {
        atomicAdd(&gm_xp[i], -tmp*inverseMassi);
        atomicAdd(&gm_xp[j],  tmp*inverseMassj);
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

        float  len2  = targetLength*targetLength;
        float  dlen2 = 2.0f*len2 - norm2(dx);

        // TODO A little bit more effective but slightly less readable version of the below would be:
        //      float proj = sqrtReducedMass*(targetLength - (dlen2 > 0.0f ? 1.0f : 0.0f)*dlen2*rsqrt(dlen2));
        float  proj;
        if (dlen2 > 0.0f)
        {
            proj = sqrtReducedMass*(targetLength - dlen2*rsqrt(dlen2));
        }
        else
        {
            proj = sqrtReducedMass*targetLength;
        }

        sm_rhs[threadIdx.x]   = proj;
        float sol = proj;

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
                int index = n*numConstraintsThreads + threadIndex;
                int c1    = gm_coupledConstraintsIdxes[index];

                mvb = mvb + gm_matrixA[index]*sm_rhs[c1 - blockIdx.x*blockDim.x + blockDim.x*(rec % 2)];

            }
            sm_rhs[threadIdx.x + blockDim.x*((rec + 1) % 2)] = mvb;
            sol  = sol + mvb;
        }

        // Add corrections to Lagrange multipliers
        float sqrtmu_sol  = sqrtReducedMass*sol;
        lagrangeScaled += sqrtmu_sol;

        // Save updated coordinates for the next iteration
        // Dummy constraints are skipped
        if (!isDummyThread)
        {
            float3 tmp = rc*sqrtmu_sol;
            atomicAdd(&gm_xp[i], -tmp*inverseMassi);
            atomicAdd(&gm_xp[j],  tmp*inverseMassj);
        }
    }

    // Updating particle velocities for all but dummy threads
    if (updateVelocities && !isDummyThread)
    {
        float3 tmp     = rc*invdt*lagrangeScaled;
        atomicAdd(&gm_v[i], -tmp*inverseMassi);
        atomicAdd(&gm_v[j],  tmp*inverseMassj);
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
        // two are no longer in use.
        extern __shared__ float  sm_threadVirial[];
        float                    mult = targetLength*lagrangeScaled;
        sm_threadVirial[0*blockDim.x + threadIdx.x] = mult*rc.x*rc.x;
        sm_threadVirial[1*blockDim.x + threadIdx.x] = mult*rc.x*rc.y;
        sm_threadVirial[2*blockDim.x + threadIdx.x] = mult*rc.x*rc.z;
        sm_threadVirial[3*blockDim.x + threadIdx.x] = mult*rc.y*rc.y;
        sm_threadVirial[4*blockDim.x + threadIdx.x] = mult*rc.y*rc.z;
        sm_threadVirial[5*blockDim.x + threadIdx.x] = mult*rc.z*rc.z;

        __syncthreads();

        // Reduce up to one virial per thread block. All blocks are divided by half, the first
        // half of threads sums two virials. Then the first half is divided by two and the first
        // half of it sums two values. This procedure is repeated until only one thread is left.
        // Only works if the threads per blocks is a power of two (hence static_assert
        // in the beginning of the kernel).
        for (int divideBy = 2; divideBy <= static_cast<int>(blockDim.x); divideBy *= 2)
        {
            int dividedAt = blockDim.x/divideBy;
            if (static_cast<int>(threadIdx.x) < dividedAt)
            {
                for (int d = 0; d < 6; d++)
                {
                    sm_threadVirial[d*blockDim.x + threadIdx.x] += sm_threadVirial[d*blockDim.x + (threadIdx.x + dividedAt)];
                }
            }
            // Syncronize if not within one warp
            if (dividedAt > warpSize/2)
            {
                __syncthreads();
            }
        }
        // First 6 threads in the block add the results of 6 tensor components to the global memory address.
        if (threadIdx.x < 6)
        {
            atomicAdd(&(gm_virialScaled[threadIdx.x]), sm_threadVirial[threadIdx.x*blockDim.x]);
        }
    }

    return;
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
inline auto getLincsKernelPtr(const bool  updateVelocities,
                              const bool  computeVirial)
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

/*! \brief Apply LINCS.
 *
 * Applies LINCS to coordinates and velocities, stored on GPU.
 * Data at pointers xPrime and v (class fields) change in the GPU
 * memory. The results are not automatically copied back to the CPU
 * memory. Method uses this class data structures which should be
 * updated when needed using update method.
 *
 * \param[in] updateVelocities  If the velocities should be constrained.
 * \param[in] invdt             Reciprocal timestep (to scale Lagrange
 *                              multipliers when velocities are updated)
 * \param[in] computeVirial     If virial should be updated.
 * \param[in,out] virialScaled  Scaled virial tensor to be updated.
 */
void LincsCuda::Impl::apply(const bool  updateVelocities,
                            const real  invdt,
                            const bool  computeVirial,
                            tensor      virialScaled)
{
    ensureNoPendingCudaError("In CUDA version of LINCS");

    if (computeVirial)
    {
        // Fill with zeros so the values can be reduced to it
        // Only 6 values are needed because virial is symmetrical
        clearDeviceBufferAsync(&kernelParams_.d_virialScaled, 0, 6, stream_);
    }

    auto               kernelPtr = getLincsKernelPtr(updateVelocities, computeVirial);

    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (kernelParams_.numConstraintsThreads + c_threadsPerBlock - 1)/c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;

    // Shared memory is used to store:
    // -- Current coordinates (3 floats per thread)
    // -- Right-hand-sides for matrix inversion (2 floats per thread)
    // -- Virial tensor components (6 floats per thread)
    // Since none of these three are needed simultaneously, they can be saved at the same shared memory address
    // (i.e. correspondent arrays are intentionally overlapped in address space). Consequently, only
    // max{3, 2, 6} = 6 floats per thread are needed in case virial is computed, or max{3, 2} = 3 if not.
    if (computeVirial)
    {
        config.sharedMemorySize = c_threadsPerBlock*6*sizeof(float);
    }
    else
    {
        config.sharedMemorySize = c_threadsPerBlock*3*sizeof(float);
    }
    config.stream           = stream_;

    const auto     kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                          &kernelParams_,
                                                          &invdt);

    launchGpuKernel(kernelPtr, config, nullptr,
                    "lincs_kernel<updateVelocities, computeVirial>", kernelArgs);

    if (computeVirial)
    {
        // Copy LINCS virial data and add it to the common virial
        copyFromDeviceBuffer(h_virialScaled_.data(), &kernelParams_.d_virialScaled,
                             0, 6,
                             stream_, GpuApiCallBehavior::Sync, nullptr);

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

/*! \brief Create LINCS object
 *
 * \param[in] numAtoms        Number of atoms that will be handles by LINCS.
 *                            Used to compute the memory size in allocations and copy.
 * \param[in] numIterations   Number of iterations used to compute inverse matrix.
 * \param[in] expansionOrder  LINCS projection order for correcting the direction of constraint.
 */
LincsCuda::Impl::Impl(int numAtoms,
                      int numIterations,
                      int expansionOrder)
    : numAtoms_(numAtoms)
{
    kernelParams_.numIterations              = numIterations;
    kernelParams_.expansionOrder             = expansionOrder;

    static_assert(sizeof(real) == sizeof(float),
                  "Real numbers should be in single precision in GPU code.");
    static_assert(c_threadsPerBlock > 0 && ((c_threadsPerBlock & (c_threadsPerBlock - 1)) == 0),
                  "Number of threads per block should be a power of two in order for reduction to work.");

    // This is temporary. LINCS should not manage coordinates.
    allocateDeviceBuffer(&kernelParams_.d_x,  numAtoms, nullptr);
    allocateDeviceBuffer(&kernelParams_.d_xp, numAtoms, nullptr);
    allocateDeviceBuffer(&kernelParams_.d_v,  numAtoms, nullptr);

    allocateDeviceBuffer(&kernelParams_.d_virialScaled, 6, nullptr);
    h_virialScaled_.resize(6);

    // The data arrays should be expanded/reallocated on first call of set() function.
    maxConstraintsNumberSoFar_ = 0;
    // Use default stream.
    // TODO The stream should/can be assigned by the GPU schedule when the code will be integrated.
    stream_ = nullptr;

}

LincsCuda::Impl::~Impl()
{
    freeDeviceBuffer(&kernelParams_.d_x);
    freeDeviceBuffer(&kernelParams_.d_xp);
    freeDeviceBuffer(&kernelParams_.d_v);
    freeDeviceBuffer(&kernelParams_.d_virialScaled);

    if (maxConstraintsNumberSoFar_ > 0)
    {
        freeDeviceBuffer(&kernelParams_.d_inverseMasses);

        freeDeviceBuffer(&kernelParams_.d_constraints);
        freeDeviceBuffer(&kernelParams_.d_constraintsTargetLengths);

        freeDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts);
        freeDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices);
        freeDeviceBuffer(&kernelParams_.d_massFactors);
        freeDeviceBuffer(&kernelParams_.d_matrixA);
    }
}

/*! \brief Helper function to go through constraints recursively.
 *
 *  For each constraint, counts the number of coupled constraints and stores the value in spaceNeeded array.
 *  This information is used to split the array of constraints between thread blocks on a GPU so there is no
 *  coupling between constraints from different thread blocks. After the 'spaceNeeded' array is filled, the
 *  value spaceNeeded[c] should be equal to the number of constraints that are coupled to 'c' and located
 *  after it in the constraints array.
 *
 * \param[in]     a                  Atom index.
 * \param[in,out] spaceNeeded        Indicates if the constraint was already counted and stores
 *                                   the number of constraints (i) connected to it and (ii) located
 *                                   after it in memory. This array is filled by this recursive function.
 *                                   For a set of coupled constraints, only for the first one in this list
 *                                   the number of consecutive coupled constraints is needed: if there is
 *                                   not enough space for this set of constraints in the thread block,
 *                                   the group has to be moved to the next one.
 * \param[in]     atomAdjacencyList  Stores information about connections between atoms.
 */
inline int countCoupled(int a, std::vector<int> *spaceNeeded,
                        std::vector<std::vector<std::tuple<int, int, int> > > *atomsAdjacencyList)

{
    int c2, a2, sign;
    int counted = 0;
    for (unsigned i = 0; i < atomsAdjacencyList->at(a).size(); i++)
    {
        std::tie(a2, c2, sign) = atomsAdjacencyList->at(a).at(i);
        if (spaceNeeded->at(c2) == -1)
        {
            spaceNeeded->at(c2) = 0; // To indicate we've been here
            counted            += 1 + countCoupled(a2, spaceNeeded, atomsAdjacencyList);
        }
    }
    return counted;
}

/*! \brief
 * Update data-structures (e.g. after NB search step).
 *
 * Updates the constraints data and copies it to the GPU. Should be
 * called if the particles were sorted, redistributed between domains, etc.
 * This version uses common data formats so it can be called from anywhere
 * in the code. Does not recycle the data preparation routines from the CPU
 * version. Works only with simple case when all the constraints in idef are
 * are handled by a single GPU. Triangles are not handled as special case.
 *
 * Information about constraints is taken from:
 *     idef.il[F_CONSTR].iatoms  --- type (T) of constraint and two atom indexes (i1, i2)
 *     idef.iparams[T].constr.dA --- target length for constraint of type T
 * From t_mdatom, the code takes:
 *     md.invmass  --- array of inverse square root of masses for each atom in the system.
 *
 * \param[in] idef  Local topology data to get information on constraints from.
 * \param[in] md    Atoms data to get atom masses from.
 */
void LincsCuda::Impl::set(const t_idef    &idef,
                          const t_mdatoms &md)
{

    // List of constrained atoms (CPU memory)
    std::vector<int2>   constraintsHost;
    // Equilibrium distances for the constraints (CPU)
    std::vector<float>  constraintsTargetLengthsHost;
    // Number of constraints, coupled with the current one (CPU)
    std::vector<int>    coupledConstraintsCountsHost;
    // List of coupled with the current one (CPU)
    std::vector<int>    coupledConstraintsIndicesHost;
    // Mass factors (CPU)
    std::vector<float>  massFactorsHost;

    // List of constrained atoms in local topology
    t_iatom  *iatoms         = idef.il[F_CONSTR].iatoms;
    const int stride         = NRAL(F_CONSTR) + 1;
    const int numConstraints = idef.il[F_CONSTR].nr/stride;
    // Constructing adjacency list --- usefull intermediate structure
    std::vector<std::vector<std::tuple<int, int, int> > > atomsAdjacencyList(numAtoms_);
    for (int c = 0; c < numConstraints; c++)
    {
        int a1     = iatoms[stride*c + 1];
        int a2     = iatoms[stride*c + 2];

        // Each constraint will be represented as a tuple, containing index of the second constrained atom,
        // index of the constraint and a sign that indicates the order of atoms in which they are listed.
        // Sign is needed to compute the mass factors.
        atomsAdjacencyList.at(a1).push_back(std::make_tuple(a2, c, +1));
        atomsAdjacencyList.at(a2).push_back(std::make_tuple(a1, c, -1));
    }

    // Compute, how many coupled constraints are in front of each constraint.
    // Needed to introduce splits in data so that all coupled constraints will be computed in a single GPU block.
    // The position 'c' of the vector spaceNeeded should have the number of constraints that are coupled to a constraint
    // 'c' and are after 'c' in the vector. Only first index of the connected group of the constraints is needed later in the
    // code, hence the spaceNeeded vector is also used to keep track if the constrain was already counted.
    std::vector<int> spaceNeeded;
    spaceNeeded.resize(numConstraints, -1);
    std::fill(spaceNeeded.begin(), spaceNeeded.end(), -1);
    for (int c = 0; c < numConstraints; c++)
    {
        int a1     = iatoms[stride*c + 1];
        int a2     = iatoms[stride*c + 2];
        if (spaceNeeded.at(c) == -1)
        {
            spaceNeeded.at(c) = countCoupled(a1, &spaceNeeded, &atomsAdjacencyList) +
                countCoupled(a2, &spaceNeeded, &atomsAdjacencyList);
        }
    }

    // Map of splits in the constraints data. For each 'old' constraint index gives 'new' which
    // takes into account the empty spaces which might be needed in the end of each thread block.
    std::vector<int> splitMap;
    splitMap.resize(numConstraints, -1);
    int              currentMapIndex = 0;
    for (int c = 0; c < numConstraints; c++)
    {
        // Check if coupled constraints all fit in one block
        GMX_RELEASE_ASSERT(spaceNeeded.at(c) < c_threadsPerBlock, "Maximum number of coupled constraints exceedes the size of the CUDA thread block. "
                           "Most likely, you are trying to use GPU version of LINCS with constraints on all-bonds, "
                           "which is not supported. Try using H-bonds constraints only.");
        if (currentMapIndex / c_threadsPerBlock != (currentMapIndex + spaceNeeded.at(c)) / c_threadsPerBlock)
        {
            currentMapIndex = ((currentMapIndex/c_threadsPerBlock) + 1) * c_threadsPerBlock;
        }
        splitMap.at(c) = currentMapIndex;
        currentMapIndex++;
    }
    kernelParams_.numConstraintsThreads = currentMapIndex + c_threadsPerBlock - currentMapIndex % c_threadsPerBlock;
    GMX_RELEASE_ASSERT(kernelParams_.numConstraintsThreads % c_threadsPerBlock == 0, "Number of threads should be a multiple of the block size");

    // Initialize constraints and their target indexes taking into account the splits in the data arrays.
    int2 pair;
    pair.x = -1;
    pair.y = -1;
    constraintsHost.resize(kernelParams_.numConstraintsThreads, pair);
    std::fill(constraintsHost.begin(), constraintsHost.end(), pair);
    constraintsTargetLengthsHost.resize(kernelParams_.numConstraintsThreads, 0.0);
    std::fill(constraintsTargetLengthsHost.begin(), constraintsTargetLengthsHost.end(), 0.0);
    for (int c = 0; c < numConstraints; c++)
    {
        int  a1     = iatoms[stride*c + 1];
        int  a2     = iatoms[stride*c + 2];
        int  type   = iatoms[stride*c];

        int2 pair;
        pair.x = a1;
        pair.y = a2;
        constraintsHost.at(splitMap.at(c))              = pair;
        constraintsTargetLengthsHost.at(splitMap.at(c)) = idef.iparams[type].constr.dA;

    }

    // The adjacency list of constraints (i.e. the list of coupled constraints for each constraint).
    // We map a single thread to a single constraint, hence each thread 'c' will be using one element from
    // coupledConstraintsCountsHost array, which is the number of constraints coupled to the constraint 'c'.
    // The coupled constraints indexes are placed into the coupledConstraintsIndicesHost array. Latter is organized
    // as a one-dimensional array to ensure good memory alignment. It is addressed as [c + i*numConstraintsThreads],
    // where 'i' goes from zero to the number of constraints coupled to 'c'. 'numConstraintsThreads' is the width of
    // the array --- a number, greater then total number of constraints, taking into account the splits in the
    // constraints array due to the GPU block borders. This number can be adjusted to improve memory access pattern.
    // Mass factors are saved in a similar data structure.
    int              maxCoupledConstraints = 0;
    for (int c = 0; c < numConstraints; c++)
    {
        int a1     = iatoms[stride*c + 1];
        int a2     = iatoms[stride*c + 2];

        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        int nCoupedConstraints = atomsAdjacencyList.at(a1).size() + atomsAdjacencyList.at(a2).size() - 2;

        if (nCoupedConstraints > maxCoupledConstraints)
        {
            maxCoupledConstraints = nCoupedConstraints;
        }
    }

    coupledConstraintsCountsHost.resize(kernelParams_.numConstraintsThreads, 0);
    coupledConstraintsIndicesHost.resize(maxCoupledConstraints*kernelParams_.numConstraintsThreads, -1);
    massFactorsHost.resize(maxCoupledConstraints*kernelParams_.numConstraintsThreads, -1);

    for (int c1 = 0; c1 < numConstraints; c1++)
    {
        coupledConstraintsCountsHost.at(splitMap.at(c1))  = 0;
        int c1a1     = iatoms[stride*c1 + 1];
        int c1a2     = iatoms[stride*c1 + 2];
        int c2;
        int c2a1;
        int c2a2;

        int sign;

        // Constraints, coupled trough the first atom.
        c2a1 = c1a1;
        for (unsigned j = 0; j < atomsAdjacencyList.at(c1a1).size(); j++)
        {

            std::tie(c2a2, c2, sign) = atomsAdjacencyList.at(c1a1).at(j);

            if (c1 != c2)
            {
                int index = kernelParams_.numConstraintsThreads*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);

                coupledConstraintsIndicesHost.at(index) = splitMap.at(c2);

                int   center = c1a1;

                float sqrtmu1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                float sqrtmu2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = -sign*md.invmass[center]*sqrtmu1*sqrtmu2;

                coupledConstraintsCountsHost.at(splitMap.at(c1))++;

            }
        }

        // Constraints, coupled through the second atom.
        c2a1 = c1a2;
        for (unsigned j = 0; j < atomsAdjacencyList.at(c1a2).size(); j++)
        {

            std::tie(c2a2, c2, sign) = atomsAdjacencyList.at(c1a2).at(j);

            if (c1 != c2)
            {
                int index = kernelParams_.numConstraintsThreads*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);

                coupledConstraintsIndicesHost.at(index) = splitMap.at(c2);

                int   center = c1a2;

                float sqrtmu1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                float sqrtmu2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = sign*md.invmass[center]*sqrtmu1*sqrtmu2;

                coupledConstraintsCountsHost.at(splitMap.at(c1))++;

            }
        }
    }

    // (Re)allocate the memory, if the number of constraints has increased.
    if (numConstraints > maxConstraintsNumberSoFar_)
    {
        // Free memory if it was allocated before (i.e. if not the first time here).
        if (maxConstraintsNumberSoFar_ > 0)
        {
            freeDeviceBuffer(&kernelParams_.d_inverseMasses);

            freeDeviceBuffer(&kernelParams_.d_constraints);
            freeDeviceBuffer(&kernelParams_.d_constraintsTargetLengths);

            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts);
            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices);
            freeDeviceBuffer(&kernelParams_.d_massFactors);
            freeDeviceBuffer(&kernelParams_.d_matrixA);

        }
        maxConstraintsNumberSoFar_ = numConstraints;

        allocateDeviceBuffer(&kernelParams_.d_inverseMasses, numAtoms_, nullptr);

        allocateDeviceBuffer(&kernelParams_.d_constraints, kernelParams_.numConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_constraintsTargetLengths, kernelParams_.numConstraintsThreads, nullptr);

        allocateDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts, kernelParams_.numConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices, maxCoupledConstraints*kernelParams_.numConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_massFactors, maxCoupledConstraints*kernelParams_.numConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_matrixA, maxCoupledConstraints*kernelParams_.numConstraintsThreads, nullptr);

    }

    // Copy data to GPU.
    copyToDeviceBuffer(&kernelParams_.d_constraints, constraintsHost.data(),
                       0, kernelParams_.numConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_constraintsTargetLengths, constraintsTargetLengthsHost.data(),
                       0, kernelParams_.numConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts, coupledConstraintsCountsHost.data(),
                       0, kernelParams_.numConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices, coupledConstraintsIndicesHost.data(),
                       0, maxCoupledConstraints*kernelParams_.numConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_massFactors, massFactorsHost.data(),
                       0, maxCoupledConstraints*kernelParams_.numConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);

    GMX_RELEASE_ASSERT(md.invmass != nullptr, "Masses of attoms should be specified.\n");
    copyToDeviceBuffer(&kernelParams_.d_inverseMasses, md.invmass,
                       0, numAtoms_,
                       stream_, GpuApiCallBehavior::Sync, nullptr);

}

/*! \brief
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] pbc  The PBC data in t_pbc format.
 */
void LincsCuda::Impl::setPbc(const t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &kernelParams_.pbcAiuc);
}

/*! \brief
 * Copy coordinates from provided CPU location to GPU.
 *
 * Copies the coordinates before the integration step (x) and coordinates
 * after the integration step (xp) from the provided CPU location to GPU.
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_x   CPU pointer where coordinates should be copied from.
 * \param[in] h_xp  CPU pointer where coordinates should be copied from.
 */
void LincsCuda::Impl::copyCoordinatesToGpu(const rvec *h_x, const rvec *h_xp)
{
    copyToDeviceBuffer(&kernelParams_.d_x, (float3*)h_x, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_xp, (float3*)h_xp, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from provided CPU location to GPU.
 *
 * Nothing is done if the argument provided is a nullptr.
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_v  CPU pointer where velocities should be copied from.
 */
void LincsCuda::Impl::copyVelocitiesToGpu(const rvec *h_v)
{
    if (h_v != nullptr)
    {
        copyToDeviceBuffer(&kernelParams_.d_v, (float3*)h_v, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
    }
}

/*! \brief
 * Copy coordinates from GPU to provided CPU location.
 *
 * Copies the constrained coordinates to the provided location. The coordinates
 * are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] h_xp  CPU pointer where coordinates should be copied to.
 */
void LincsCuda::Impl::copyCoordinatesFromGpu(rvec *h_xp)
{
    copyFromDeviceBuffer((float3*)h_xp, &kernelParams_.d_xp, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from GPU to provided CPU location.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] h_v  Pointer to velocities data.
 */
void LincsCuda::Impl::copyVelocitiesFromGpu(rvec *h_v)
{
    copyFromDeviceBuffer((float3*)h_v, &kernelParams_.d_v, 0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * Data is not copied. The data are assumed to be in float3/fvec format
 * (float3 is used internally, but the data layout should be identical).
 *
 * \param[in] d_x   Pointer to the coordinates before integrator update (on GPU)
 * \param[in] d_xp  Pointer to the coordinates after integrator update, before update (on GPU)
 * \param[in] d_v   Pointer to the velocities before integrator update (on GPU)
 */
void LincsCuda::Impl::setXVPointers(rvec *d_x, rvec *d_xp, rvec *d_v)
{
    kernelParams_.d_x  = (float3*)d_x;
    kernelParams_.d_xp = (float3*)d_xp;
    kernelParams_.d_v  = (float3*)d_v;
}


LincsCuda::LincsCuda(const int numAtoms,
                     const int numIterations,
                     const int expansionOrder)
    : impl_(new Impl(numAtoms, numIterations, expansionOrder))
{
}

LincsCuda::~LincsCuda() = default;

void LincsCuda::apply(const bool  updateVelocities,
                      const real  invdt,
                      const bool  computeVirial,
                      tensor      virialScaled)
{
    impl_->apply(updateVelocities,
                 invdt,
                 computeVirial,
                 virialScaled);
}

void LincsCuda::setPbc(const t_pbc *pbc)
{
    impl_->setPbc(pbc);
}

void LincsCuda::set(const t_idef    &idef,
                    const t_mdatoms &md)
{
    impl_->set(idef, md);
}

void LincsCuda::copyCoordinatesToGpu(const rvec *h_x, const rvec *h_xp)
{
    impl_->copyCoordinatesToGpu(h_x, h_xp);
}

void LincsCuda::copyVelocitiesToGpu(const rvec *h_v)
{
    impl_->copyVelocitiesToGpu(h_v);
}

void LincsCuda::copyCoordinatesFromGpu(rvec *h_xp)
{
    impl_->copyCoordinatesFromGpu(h_xp);
}

void LincsCuda::copyVelocitiesFromGpu(rvec *h_v)
{
    impl_->copyVelocitiesFromGpu(h_v);
}

void LincsCuda::setXVPointers(rvec *d_x, rvec *d_xp, rvec *d_v)
{
    impl_->setXVPointers(d_x, d_xp, d_v);
}

} // namespace gmx
