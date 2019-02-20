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

namespace gmx
{

constexpr static int c_threadsPerBlock = 256;

/*! \brief Main kernel for LINCS constraints.
 *
 * See Hess et al., J. Comput. Chem. 18: 1463-1472 (1997) for the description of the algorithm.
 *
 * In CUDA version, one thread is responsible for all computations for one constrain. The blocks are
 * filled in a way that no constraint is coupled to the constraint from the next block. This is achieved
 * by adding dummy threads in the end of the thread block if needed. This removes necessity to
 * synchronize across the device, but extensive communications in a thread block are still needed.
 *
 * \todo Move more data to local/shared memory and try to get rid of atomics (at least on the device level).
 *
 * \param[in]     nCons                     Total number of the constrain threads (empty spaces at the end of thread blocks included)
 * \param[in]     x                         Coordinates before the timestep
 * \param[in,out] xp                        Coordinates after the timestep. Will be updated to satisfy the constraints by this kernel.
 * \param[in]     nIter                     Number of iterations used to compute inverse matrix.
 * \param[in]     nOrder                    LINCS projection order for correcting the direction of constraint.
 * \param[in]     constraints               List of constraints.
 * \param[in]     constraintsR0             List of target distances for constraints.
 * \param[in]     coupledConstraintsCounts  Numbers of coupled constraints for each constraint
 * \param[in]     coupledConstraintsIdxes   Indexes of coupled constraints
 * \param[in]     massFactors               Mass factors: ( (+/-) * (1/sqrt(1/m1 + 1/m2)) * (1/m2) * 1/sqrt(1/m2 + 1/m3)),
 *                                          where m1 and m3 are coupled through m2 and sign + or - indicates the order,
 *                                          in which they are arranged in atoms array.)
 * \param[in]     matrixA                   The matrix A = (I-S*B_n*M^(-1)*B_n^T*S) used to compute the projection matrix.
 * \param[in]     pbcAiuc                   Periodic boundary information
 * \param[in,out] v                         Velocities to update.
 * \param[in]     invdt                     Inverse timestep (needed to update velocities).
 * \param[in,out] virialScaled              Scaled virial tensor to be updated.
 */
template <bool updateVelocities, bool computeVirial>
__global__ void lincs_kernel(LincsCudaKernelParameters   kernelParams,
                             const float                 invdt)
{
    const PbcAiuc               pbcAiuc                     = kernelParams.pbcAiuc;
    const int                   nCons                       = kernelParams.nConstraintsThreads;
    const int                   nIter                       = kernelParams.nIter;
    const int                   nOrder                      = kernelParams.nOrder;
    const float3* __restrict__  gm_x                        = kernelParams.d_x;
    float3*                     gm_xp                       = kernelParams.d_xp;
    const int2*   __restrict__  gm_constraints              = kernelParams.d_constraints;
    const float*  __restrict__  gm_constraintsTargetLengths = kernelParams.d_constraintsTargetLengths;
    const int*    __restrict__  gm_coupledConstraintsCounts = kernelParams.d_coupledConstraintsCounts;
    const int*    __restrict__  gm_coupledConstraintsIdxes  = kernelParams.d_coupledConstraintsIdxes;
    const float*  __restrict__  gm_massFactors              = kernelParams.d_massFactors;
    float*  __restrict__        gm_matrixA                  = kernelParams.d_matrixA;
    const float*  __restrict__  gm_inverseMasses            = kernelParams.d_inverseMasses;
    float3*                     gm_v                        = kernelParams.d_v;
    float*  __restrict__        gm_virialScaled             = kernelParams.d_virialScaled;

    int threadIndex                                         = blockIdx.x*blockDim.x+threadIdx.x;

    // nCons should be a factor of blockSize
    assert(threadIndex < nCons);

    // Vectors connecting constrained atoms before algorithm was applied.
    // Needed to construct constrain matrix A
    extern __shared__ float3 sm_r[];
    extern __shared__ float  sm_rhs[];

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
    // Everything computed for these dummies will be equal to zero
    if (i == -1)
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
    int coupledConstraintsCount = gm_coupledConstraintsCounts[c];
    for (int n = 0; n < coupledConstraintsCount; n++)
    {
        int    index = n*nCons + threadIndex;
        int    c1    = gm_coupledConstraintsIdxes[index];

        float3 rc1 = sm_r[c1 - blockIdx.x*blockDim.x];
        gm_matrixA[index] = gm_massFactors[index]*(rc.x*rc1.x + rc.y*rc1.y + rc.z*rc1.z);
    }

    // Skipping dummy threads
    if (i != -1)
    {
        xi = gm_xp[i];
        xj = gm_xp[j];
    }

    float3 dx = pbcDxAiuc(pbcAiuc, xi, xj);

    float  sol = sqrtReducedMass*((rc.x*dx.x + rc.y*dx.y + rc.z*dx.z) - targetLength);

    /*
     *  Inverse matrix using a set of nOrder matrix multiplications
     */

    // Save current right-hand-side vector in the shared memory
    // Values in sm_r are no longer needed to be shared across threads (sm_rhs overwrites them)
    sm_rhs[threadIdx.x] = sol;
    // Making sure that all sm_rhs are saved before they are accessed in the loop below
    __syncthreads();
    for (int rec = 0; rec < nOrder; rec++)
    {
        float mvb = 0.0f;

        for (int n = 0; n < coupledConstraintsCount; n++)
        {
            int index = n*nCons + threadIndex;
            int c1    = gm_coupledConstraintsIdxes[index];
            // Convolute current right-hand-side with A
            // Different, non overlapping parts of sm_rhs[...] are read during odd and even iterations
            mvb = mvb + gm_matrixA[index]*sm_rhs[c1 - blockIdx.x*blockDim.x + blockDim.x*(rec % 2)];
        }
        // 'Switch' rhs vectors, save current result
        // These values will be accessed in the loop above during the next iteration.
        sm_rhs[threadIdx.x + blockDim.x*((rec + 1) % 2)] = mvb;
        sol  = sol + mvb;
        // All thread should save the rhs to shared memory
        __syncthreads();
    }

    // Current mass-scaled Lagrange multipliers
    lagrangeScaled = sqrtReducedMass*sol;

    // Save updated coordinates before correction for the rotational lengthening
    float3 tmp     = rc*lagrangeScaled;

    // Writing for all but dummy constraints
    if (i != -1)
    {
        atomicAdd(&gm_xp[i], -tmp*inverseMassi);
        atomicAdd(&gm_xp[j],  tmp*inverseMassj);
    }

    // Make sure that all xp's are saved
    __syncthreads();

    /*
     *  Correction for centripetal effects
     */
    for (int iter = 0; iter < nIter; iter++)
    {
        if (i != -1)
        {
            xi = gm_xp[i];
            xj = gm_xp[j];
        }

        float3 dx = pbcDxAiuc(pbcAiuc, xi, xj);

        float  len2  = targetLength*targetLength;
        float  dlen2 = 2.0f*len2 - norm2(dx);

        // TODO A little bit more effective but slightly less readable version of the below would be:
        //      float proj = sqrtReducedMass(targetLength - (dlen2 > 0.0f ? 1.0f : 0.0f)*dlen2*rsqrt(dlen2)));
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
        // Make sure that all elements of rhs are saved into shared memory
        __syncthreads();

        /*
         * Same matrix inversion for updated data
         */
        for (int rec = 0; rec < nOrder; rec++)
        {
            float mvb = 0;

            for (int n = 0; n < coupledConstraintsCount; n++)
            {
                int index = n*nCons + threadIndex;
                int c1    = gm_coupledConstraintsIdxes[index];

                mvb = mvb + gm_matrixA[index]*sm_rhs[c1 - blockIdx.x*blockDim.x + blockDim.x*(rec % 2)];

            }
            sm_rhs[threadIdx.x + blockDim.x*((rec + 1) % 2)] = mvb;
            sol  = sol + mvb;
            // Make sure that all elements of rhs array are updated
            __syncthreads();

        }

        // Correct Lagrange multipliers
        float sqrtmu_sol  = sqrtReducedMass*sol;
        lagrangeScaled += sqrtmu_sol;

        // Save updated coordinates for the next iteration
        // Dummy constraints are skipped
        if (i != -1)
        {
            float3 tmp = rc*sqrtmu_sol;
            atomicAdd(&gm_xp[i], -tmp*inverseMassi);
            atomicAdd(&gm_xp[j],  tmp*inverseMassj);
        }
        __syncthreads();
    }

    // Updating particle velocities for all but dummy threads
    if (updateVelocities && i != -1)
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

        // Save virial for each thread into the shared memory
        // Tensor is symmetrical, hence only 6 values are saved
        extern __shared__ float  sm_threadVirial[];
        float                    mult = targetLength*lagrangeScaled;
        sm_threadVirial[0*blockDim.x + threadIdx.x] = mult*rc.x*rc.x;
        sm_threadVirial[1*blockDim.x + threadIdx.x] = mult*rc.x*rc.y;
        sm_threadVirial[2*blockDim.x + threadIdx.x] = mult*rc.x*rc.z;
        sm_threadVirial[3*blockDim.x + threadIdx.x] = mult*rc.y*rc.y;
        sm_threadVirial[4*blockDim.x + threadIdx.x] = mult*rc.y*rc.z;
        sm_threadVirial[5*blockDim.x + threadIdx.x] = mult*rc.z*rc.z;

        // Reduce up to one virial per thread block
        // All blocks are divided by half, the first half of threads sums
        // two virials. Then the first half is divided by two and the first half
        // of it sums two values... The procedure continues until only one thread left.
        // Only works if the threads per blocks is a power of two (hence static_assert upon
        // initialization).
        for (int divideBy = 2; divideBy <= blockDim.x; divideBy *= 2)
        {
            int dividedAt = blockDim.x/divideBy;
            if (threadIdx.x < dividedAt)
            {
                for (int d = 0; d < 6; d++)
                {
                    sm_threadVirial[d*blockDim.x + threadIdx.x] += sm_threadVirial[d*blockDim.x + (threadIdx.x + dividedAt)];
                }
            }
            __syncthreads();
        }
        // First thread in the block adds the result to the global memory address
        // TODO Make the first 6 threads to contribute 6 values (not trivial when number
        //      of constraints is 1, e.g. in tests)
        if (threadIdx.x == 0)
        {
            for (int d = 0; d < 6; d++)
            {
                atomicAdd(&(gm_virialScaled[d]), sm_threadVirial[d*blockDim.x]);
            }
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
    config.gridSize[0]      = (kernelParams_.nConstraintsThreads + c_threadsPerBlock - 1)/c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = c_threadsPerBlock*6*sizeof(float);
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
 * \param [in] nAtom  Number of atoms that will be handles by LINCS.
 *                    Used to compute the memory size in allocations and copy.
 * \param [in] nIter  Number of iterations used to compute inverse matrix.
 * \param [in] nOrder LINCS projection order for correcting the direction of constraint.
 */
LincsCuda::Impl::Impl(int nAtom,
                      int nIter,
                      int nOrder)
    : nAtom_(nAtom)
{
    kernelParams_.nIter                      = nIter;
    kernelParams_.nOrder                     = nOrder;

    static_assert(sizeof(real) == sizeof(float),
                  "Real numbers should be in single precision in GPU code.");
    static_assert(c_threadsPerBlock > 0 && !(c_threadsPerBlock & (c_threadsPerBlock - 1) == 0),
                  "Nmber of threads per block should be a power of two in order for reduction to work.");

    // This is temporary. LINCS should not manage coordinates.
    allocateDeviceBuffer(&kernelParams_.d_x,  nAtom, nullptr);
    allocateDeviceBuffer(&kernelParams_.d_xp, nAtom, nullptr);
    allocateDeviceBuffer(&kernelParams_.d_v,  nAtom, nullptr);

    allocateDeviceBuffer(&kernelParams_.d_virialScaled, 6, nullptr);
    h_virialScaled_.resize(6);

    // The data arrays should be expanded/reallocated on first call of set() function
    maxConstraintsNumberSoFar_ = 0;
    // Use default stream
    // TODO The stream should/can be assigned by the GPU schedule when the code will be integrated
    stream_ = nullptr;

}

LincsCuda::Impl::~Impl()
{
}

/*! \brief Helper function to go through constraints recursively
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
    std::vector<int>    coupledConstraintsIdxesHost;
    // Mass factors (CPU)
    std::vector<float>  massFactorsHost;

    //t_idef idef = top.idef;
    t_iatom  *iatoms      = idef.il[F_CONSTR].iatoms;
    const int nConstraint = idef.il[F_CONSTR].nr/3;
    // Constructing adjacency list --- usefull intermediate structure
    std::vector<std::vector<std::tuple<int, int, int> > > atomsAdjacencyList(nAtom_);
    for (int c = 0; c < nConstraint; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];

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
    spaceNeeded.resize(nConstraint, -1);
    std::fill(spaceNeeded.begin(), spaceNeeded.end(), -1);
    for (int c = 0; c < nConstraint; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];
        if (spaceNeeded.at(c) == -1)
        {
            spaceNeeded.at(c) = countCoupled(a1, &spaceNeeded, &atomsAdjacencyList) +
                countCoupled(a2, &spaceNeeded, &atomsAdjacencyList);
        }
    }

    // Map of splits in the constraints data. For each 'old' constraint index gives 'new' which
    // takes into account the empty spaces which might be needed in the end of each thread block.
    std::vector<int> splitMap;
    splitMap.resize(nConstraint, -1);
    int              currentMapIndex = 0;
    for (int c = 0; c < nConstraint; c++)
    {
        if (currentMapIndex / c_threadsPerBlock != (currentMapIndex + spaceNeeded.at(c)) / c_threadsPerBlock)
        {
            currentMapIndex = ((currentMapIndex/c_threadsPerBlock) + 1) * c_threadsPerBlock;
        }
        splitMap.at(c) = currentMapIndex;
        currentMapIndex++;
    }
    kernelParams_.nConstraintsThreads = currentMapIndex + c_threadsPerBlock - currentMapIndex % c_threadsPerBlock;
    GMX_RELEASE_ASSERT(kernelParams_.nConstraintsThreads % c_threadsPerBlock == 0, "Number of threads should be a multiple of the block size");

    // Initialize constraints and their target indexes taking into account the splits in the
    // data arrays.
    int2 pair;
    pair.x = -1;
    pair.y = -1;
    constraintsHost.resize(kernelParams_.nConstraintsThreads, pair);
    std::fill(constraintsHost.begin(), constraintsHost.end(), pair);
    constraintsTargetLengthsHost.resize(kernelParams_.nConstraintsThreads, 0.0);
    std::fill(constraintsTargetLengthsHost.begin(), constraintsTargetLengthsHost.end(), 0.0);
    for (int c = 0; c < nConstraint; c++)
    {
        int  a1     = iatoms[3*c + 1];
        int  a2     = iatoms[3*c + 2];
        int  type   = iatoms[3*c];

        int2 pair;
        pair.x = a1;
        pair.y = a2;
        constraintsHost.at(splitMap.at(c))              = pair;
        constraintsTargetLengthsHost.at(splitMap.at(c)) = idef.iparams[type].constr.dA;

    }

    // The adjacency list of constraints (i.e. the list of coupled constraints for each constraint).
    // We map a single thread to a single constraint, hence each thread 'c' will be using one element from
    // coupledConstraintsCountsHost array, which is the number of constraints coupled to the constraint 'c'.
    // The coupled constraints indexes are placed into the coupledConstraintsIdxesHost array. Latter is organized
    // as a one-dimensional array to ensure good memory alignment. It is addressed as [c + i*nConstraintsThreads],
    // where 'i' goes from zero to the number of constraints coupled to 'c'. 'nConstraintsThreads' is the width of
    // the array --- a number, greater then total number of constraints, taking into account the splits in the
    // constraints array due to the GPU block borders. This number can be adjusted to improve memory access pattern.
    // Mass factors are saved in a similar data structure.
    int              maxCoupledConstraints = 0;
    for (int c = 0; c < nConstraint; c++)
    {
        int a1     = iatoms[3*c + 1];
        int a2     = iatoms[3*c + 2];

        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        int nCoupedConstraints = atomsAdjacencyList.at(a1).size() + atomsAdjacencyList.at(a2).size() - 2;

        if (nCoupedConstraints > maxCoupledConstraints)
        {
            maxCoupledConstraints = nCoupedConstraints;
        }
    }

    coupledConstraintsCountsHost.resize(kernelParams_.nConstraintsThreads, 0);
    coupledConstraintsIdxesHost.resize(maxCoupledConstraints*kernelParams_.nConstraintsThreads, -1);
    massFactorsHost.resize(maxCoupledConstraints*kernelParams_.nConstraintsThreads, -1);

    for (int c1 = 0; c1 < nConstraint; c1++)
    {
        coupledConstraintsCountsHost.at(splitMap.at(c1))  = 0;
        int c1a1     = iatoms[3*c1 + 1];
        int c1a2     = iatoms[3*c1 + 2];
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
                int index = kernelParams_.nConstraintsThreads*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);

                coupledConstraintsIdxesHost.at(index) = splitMap.at(c2);

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
                int index = kernelParams_.nConstraintsThreads*coupledConstraintsCountsHost.at(splitMap.at(c1)) + splitMap.at(c1);

                coupledConstraintsIdxesHost.at(index) = splitMap.at(c2);

                int   center = c1a2;

                float sqrtmu1 = 1.0/sqrt(md.invmass[c1a1] + md.invmass[c1a2]);
                float sqrtmu2 = 1.0/sqrt(md.invmass[c2a1] + md.invmass[c2a2]);

                massFactorsHost.at(index) = sign*md.invmass[center]*sqrtmu1*sqrtmu2;

                coupledConstraintsCountsHost.at(splitMap.at(c1))++;

            }
        }
    }

    if (nConstraint > maxConstraintsNumberSoFar_)
    {

        if (maxConstraintsNumberSoFar_ > 0)
        {
            freeDeviceBuffer(&kernelParams_.d_inverseMasses);

            freeDeviceBuffer(&kernelParams_.d_constraints);
            freeDeviceBuffer(&kernelParams_.d_constraintsTargetLengths);

            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts);
            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsIdxes);
            freeDeviceBuffer(&kernelParams_.d_massFactors);
            freeDeviceBuffer(&kernelParams_.d_matrixA);

        }
        maxConstraintsNumberSoFar_ = nConstraint;

        allocateDeviceBuffer(&kernelParams_.d_inverseMasses, nAtom_, nullptr);

        allocateDeviceBuffer(&kernelParams_.d_constraints, kernelParams_.nConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_constraintsTargetLengths, kernelParams_.nConstraintsThreads, nullptr);

        allocateDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts, kernelParams_.nConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_coupledConstraintsIdxes, maxCoupledConstraints*kernelParams_.nConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_massFactors, maxCoupledConstraints*kernelParams_.nConstraintsThreads, nullptr);
        allocateDeviceBuffer(&kernelParams_.d_matrixA, maxCoupledConstraints*kernelParams_.nConstraintsThreads, nullptr);

    }

    copyToDeviceBuffer(&kernelParams_.d_constraints, constraintsHost.data(),
                       0, kernelParams_.nConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_constraintsTargetLengths, constraintsTargetLengthsHost.data(),
                       0, kernelParams_.nConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts, coupledConstraintsCountsHost.data(),
                       0, kernelParams_.nConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsIdxes, coupledConstraintsIdxesHost.data(),
                       0, maxCoupledConstraints*kernelParams_.nConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_massFactors, massFactorsHost.data(),
                       0, maxCoupledConstraints*kernelParams_.nConstraintsThreads,
                       stream_, GpuApiCallBehavior::Sync, nullptr);

    GMX_RELEASE_ASSERT(md.invmass != nullptr, "Masses of attoms should be specified.\n");
    copyToDeviceBuffer(&kernelParams_.d_inverseMasses, md.invmass,
                       0, nAtom_,
                       stream_, GpuApiCallBehavior::Sync, nullptr);

}

/*! \brief
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] *pbc The PBC data in t_pbc format.
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
 * \param[in] *x  CPU pointer where coordinates should be copied from.
 * \param[in] *xp CPU pointer where coordinates should be copied from.
 */
void LincsCuda::Impl::copyCoordinatesToGpu(const rvec *x, const rvec *xp)
{
    copyToDeviceBuffer(&kernelParams_.d_x, (float3*)x, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&kernelParams_.d_xp, (float3*)xp, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from provided CPU location to GPU.
 *
 * Nothing is done if the argument provided is a nullptr.
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *v  CPU pointer where velocities should be copied from.
 */
void LincsCuda::Impl::copyVelocitiesToGpu(const rvec * v)
{
    if (v != nullptr)
    {
        copyToDeviceBuffer(&kernelParams_.d_v, (float3*)v, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    }
}

/*! \brief
 * Copy coordinates from GPU to provided CPU location.
 *
 * Copies the constrained coordinates to the provided location. The coordinates
 * are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] *xp CPU pointer where coordinates should be copied to.
 */
void LincsCuda::Impl::copyCoordinatesFromGpu(rvec *xp)
{
    copyFromDeviceBuffer((float3*)xp, &kernelParams_.d_xp, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from GPU to provided CPU location.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *v  Pointer to velocities data.
 */
void LincsCuda::Impl::copyVelocitiesFromGpu(rvec *v)
{
    copyFromDeviceBuffer((float3*)v, &kernelParams_.d_v, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * Data is not copied. The data are assumed to be in float3/fvec format
 * (float3 is used internally, but the data layout should be identical).
 *
 * \param[in] *d_x  Pointer to the coordinates before integrator update (on GPU)
 * \param[in] *d_xp Pointer to the coordinates after integrator update, before update (on GPU)
 * \param[in] *d_v  Pointer to the velocities before integrator update (on GPU)
 */
void LincsCuda::Impl::setXVPointers(rvec *d_x, rvec *d_xp, rvec *d_v)
{
    kernelParams_.d_x  = (float3*)d_x;
    kernelParams_.d_xp = (float3*)d_xp;
    kernelParams_.d_v  = (float3*)d_v;
}


LincsCuda::LincsCuda(const int nAtom,
                     const int nIter,
                     const int nOrder)
    : impl_(new Impl(nAtom, nIter, nOrder))
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

void LincsCuda::copyCoordinatesToGpu(const rvec *x, const rvec *xp)
{
    impl_->copyCoordinatesToGpu(x, xp);
}

void LincsCuda::copyVelocitiesToGpu(const rvec *v)
{
    impl_->copyVelocitiesToGpu(v);
}

void LincsCuda::copyCoordinatesFromGpu(rvec *xp)
{
    impl_->copyCoordinatesFromGpu(xp);
}

void LincsCuda::copyVelocitiesFromGpu(rvec *v)
{
    impl_->copyVelocitiesFromGpu(v);
}

void LincsCuda::setXVPointers(rvec *d_x, rvec *d_xp, rvec *d_v)
{
    impl_->setXVPointers(d_x, d_xp, d_v);
}

} // namespace gmx
