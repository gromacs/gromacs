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
 * \brief Implements LINCS kernels using SYCL
 *
 * This file contains SYCL kernels of LINCS constraints algorithm.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/mdlib/lincs_gpu.h"
#include "gromacs/pbcutil/pbc_aiuc_sycl.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/template_mp.h"

#include "lincs_gpu_internal.h"

namespace gmx
{

using sycl::access::fence_space;
using mode = sycl::access_mode;

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
 *       operations. Investigate this issue further.
 *
 * \tparam updateVelocities        Whether velocities should be updated this step.
 * \tparam computeVirial           Whether virial tensor should be computed this step.
 * \tparam haveCoupledConstraints  If there are coupled constraints (i.e. LINCS iterations are needed).
 *
 * \param[in]     cgh                            SYCL handler.
 * \param[in]     numConstraintsThreads          Total number of threads.
 * \param[in]     gm_constraints                 List of constrained atoms.
 * \param[in]     gm_constraintsTargetLengths    Equilibrium distances for the constraints.
 * \param[in]     gm_coupledConstraintsCounts    Number of constraints, coupled with the current one.
 * \param[in]     gm_coupledConstraintsIndices   List of coupled with the current one.
 * \param[in]     gm_massFactors                 Mass factors.
 * \param[in]     gm_matrixA                     Elements of the coupling matrix.
 * \param[in]     gm_inverseMasses               1/mass for all atoms.
 * \param[in]     numIterations                  Number of iterations used to correct the projection.
 * \param[in]     expansionOrder                 Order of expansion when inverting the matrix.
 * \param[in]     gm_x                           Unconstrained positions.
 * \param[in,out] gm_xp                          Positions at the previous step, will be updated.
 * \param[in]     invdt                          Inverse timestep (needed to update velocities).
 * \param[in,out] gm_v                           Velocities of atoms, will be updated if \c updateVelocities.
 * \param[in,out] gm_virialScaled                Scaled virial tensor (6 floats: [XX, XY, XZ, YY, YZ, ZZ]).
 *                                               Will be updated if \c updateVirial.
 * \param[in]     pbcAiuc                        Periodic boundary data.
 */
template<bool updateVelocities, bool computeVirial, bool haveCoupledConstraints>
auto lincsKernel(sycl::handler& cgh,
                 const int      numConstraintsThreads,
                 const AtomPair* __restrict__ gm_constraints,
                 const float* __restrict__ gm_constraintsTargetLengths,
                 const int* __restrict__ gm_coupledConstraintsCounts,
                 const int* __restrict__ gm_coupledConstraintsIndices,
                 const float* __restrict__ gm_massFactors,
                 float* __restrict__ gm_matrixA,
                 const float* __restrict__ gm_inverseMasses,
                 const int numIterations,
                 const int expansionOrder,
                 const Float3* __restrict__ gm_x,
                 Float3* __restrict__ gm_xp,
                 const float invdt,
                 Float3* __restrict__ gm_v,
                 float* __restrict__ gm_virialScaled,
                 PbcAiuc pbcAiuc)
{
    /* Shared local memory buffer. Corresponds to sh_r, sm_rhs, and sm_threadVirial in CUDA.
     * sh_r: one Float3 per thread.
     * sh_rhs: two floats per thread.
     * sm_threadVirial: six floats per thread.
     * So, without virials we need max(1*3, 2) floats, and with virials we need max(1*3, 2, 6) floats.
     */
    static constexpr int           smBufferElementsPerThread = computeVirial ? 6 : 3;
    sycl::local_accessor<float, 1> sm_buffer{ sycl::range<1>(c_threadsPerBlock * smBufferElementsPerThread),
                                              cgh };

    return [=](sycl::nd_item<1> itemIdx) {
        const int threadIndex   = itemIdx.get_global_linear_id();
        const int threadInBlock = itemIdx.get_local_linear_id(); // Work-item index in work-group

        AtomPair pair = gm_constraints[threadIndex];
        int      i    = pair.i;
        int      j    = pair.j;

        // Mass-scaled Lagrange multiplier
        float lagrangeScaled = 0.0F;

        float targetLength;
        float inverseMassi;
        float inverseMassj;
        float sqrtReducedMass;

        Float3 xi;
        Float3 xj;
        Float3 rc;

        // i == -1 indicates dummy constraint at the end of the thread block.
        bool isDummyThread = (i == -1);

        // Everything computed for these dummies will be equal to zero
        if (isDummyThread)
        {
            targetLength    = 0.0F;
            inverseMassi    = 0.0F;
            inverseMassj    = 0.0F;
            sqrtReducedMass = 0.0F;

            xi = Float3(0.0F, 0.0F, 0.0F);
            xj = Float3(0.0F, 0.0F, 0.0F);
            rc = Float3(0.0F, 0.0F, 0.0F);
        }
        else
        {
            // Collecting data
            targetLength    = gm_constraintsTargetLengths[threadIndex];
            inverseMassi    = gm_inverseMasses[i];
            inverseMassj    = gm_inverseMasses[j];
            sqrtReducedMass = sycl::rsqrt(inverseMassi + inverseMassj);

            xi = gm_x[i];
            xj = gm_x[j];

            Float3 dx;
            pbcDxAiucSycl(pbcAiuc, xi, xj, dx);

            float rlen = sycl::rsqrt(dx[XX] * dx[XX] + dx[YY] * dx[YY] + dx[ZZ] * dx[ZZ]);
            rc         = rlen * dx;
        }

        sm_buffer[threadInBlock * DIM + XX] = rc[XX];
        sm_buffer[threadInBlock * DIM + YY] = rc[YY];
        sm_buffer[threadInBlock * DIM + ZZ] = rc[ZZ];
        // Make sure that all r's are saved into shared memory
        // before they are accessed in the loop below
        itemIdx.barrier(fence_space::global_and_local);

        /*
         * Constructing LINCS matrix (A)
         */
        int coupledConstraintsCount = 0;
        if constexpr (haveCoupledConstraints)
        {
            // Only non-zero values are saved (for coupled constraints)
            coupledConstraintsCount = gm_coupledConstraintsCounts[threadIndex];
            for (int n = 0; n < coupledConstraintsCount; n++)
            {
                int index = n * numConstraintsThreads + threadIndex;
                int c1    = gm_coupledConstraintsIndices[index];

                Float3 rc1{ sm_buffer[c1 * DIM + XX], sm_buffer[c1 * DIM + YY], sm_buffer[c1 * DIM + ZZ] };
                gm_matrixA[index] = gm_massFactors[index]
                                    * (rc[XX] * rc1[XX] + rc[YY] * rc1[YY] + rc[ZZ] * rc1[ZZ]);
            }
        }

        // Skipping in dummy threads
        if (!isDummyThread)
        {
            xi[XX] = atomicLoad(gm_xp[i][XX]);
            xi[YY] = atomicLoad(gm_xp[i][YY]);
            xi[ZZ] = atomicLoad(gm_xp[i][ZZ]);
            xj[XX] = atomicLoad(gm_xp[j][XX]);
            xj[YY] = atomicLoad(gm_xp[j][YY]);
            xj[ZZ] = atomicLoad(gm_xp[j][ZZ]);
        }

        Float3 dx;
        pbcDxAiucSycl(pbcAiuc, xi, xj, dx);

        float sol = sqrtReducedMass * ((rc[XX] * dx[XX] + rc[YY] * dx[YY] + rc[ZZ] * dx[ZZ]) - targetLength);

        /*
         *  Inverse matrix using a set of expansionOrder matrix multiplications
         */

        // This will reuse the same buffer, because the old values are no longer needed.
        itemIdx.barrier(fence_space::local_space);
        sm_buffer[threadInBlock] = sol;

        // No need to iterate if there are no coupled constraints.
        if constexpr (haveCoupledConstraints)
        {
            for (int rec = 0; rec < expansionOrder; rec++)
            {
                // Making sure that all sm_buffer values are saved before they are accessed in a loop below
                itemIdx.barrier(fence_space::global_and_local);
                float mvb = 0.0F;
                for (int n = 0; n < coupledConstraintsCount; n++)
                {
                    int index = n * numConstraintsThreads + threadIndex;
                    int c1    = gm_coupledConstraintsIndices[index];
                    // Convolute current right-hand-side with A
                    // Different, non overlapping parts of sm_buffer[..] are read during odd and even iterations
                    mvb = mvb + gm_matrixA[index] * sm_buffer[c1 + c_threadsPerBlock * (rec % 2)];
                }
                // 'Switch' rhs vectors, save current result
                // These values will be accessed in the loop above during the next iteration.
                sm_buffer[threadInBlock + c_threadsPerBlock * ((rec + 1) % 2)] = mvb;

                sol = sol + mvb;
            }
        }

        // Current mass-scaled Lagrange multipliers
        lagrangeScaled = sqrtReducedMass * sol;

        // Save updated coordinates before correction for the rotational lengthening
        Float3 tmp = rc * lagrangeScaled;

        // Writing for all but dummy constraints
        if (!isDummyThread)
        {
            /*
             * Note: Using memory_scope::work_group for atomic_ref can be better here,
             * but for now we re-use the existing function for memory_scope::device atomics.
             */
            atomicFetchAdd(gm_xp[i][XX], -tmp[XX] * inverseMassi);
            atomicFetchAdd(gm_xp[i][YY], -tmp[YY] * inverseMassi);
            atomicFetchAdd(gm_xp[i][ZZ], -tmp[ZZ] * inverseMassi);
            atomicFetchAdd(gm_xp[j][XX], tmp[XX] * inverseMassj);
            atomicFetchAdd(gm_xp[j][YY], tmp[YY] * inverseMassj);
            atomicFetchAdd(gm_xp[j][ZZ], tmp[ZZ] * inverseMassj);
        }

        /*
         *  Correction for centripetal effects
         */
        for (int iter = 0; iter < numIterations; iter++)
        {
            // Make sure that all xp's are saved: atomic operation calls before are
            // communicating current xp[..] values across thread block.
            itemIdx.barrier(fence_space::global_and_local);

            if (!isDummyThread)
            {
                xi[XX] = atomicLoad(gm_xp[i][XX]);
                xi[YY] = atomicLoad(gm_xp[i][YY]);
                xi[ZZ] = atomicLoad(gm_xp[i][ZZ]);
                xj[XX] = atomicLoad(gm_xp[j][XX]);
                xj[YY] = atomicLoad(gm_xp[j][YY]);
                xj[ZZ] = atomicLoad(gm_xp[j][ZZ]);
            }

            Float3 dx;
            pbcDxAiucSycl(pbcAiuc, xi, xj, dx);

            float len2  = targetLength * targetLength;
            float dlen2 = 2.0F * len2 - (dx[XX] * dx[XX] + dx[YY] * dx[YY] + dx[ZZ] * dx[ZZ]);

            // TODO A little bit more effective but slightly less readable version of the below would be:
            //      float proj = sqrtReducedMass*(targetLength - (dlen2 > 0.0f ? 1.0f : 0.0f)*dlen2*rsqrt(dlen2));
            float proj;
            if (dlen2 > 0.0F)
            {
                proj = sqrtReducedMass * (targetLength - dlen2 * sycl::rsqrt(dlen2));
            }
            else
            {
                proj = sqrtReducedMass * targetLength;
            }

            sm_buffer[threadInBlock] = proj;
            float sol                = proj;

            /*
             * Same matrix inversion as above is used for updated data
             */
            if constexpr (haveCoupledConstraints)
            {
                for (int rec = 0; rec < expansionOrder; rec++)
                {
                    // Make sure that all elements of rhs are saved into shared memory
                    itemIdx.barrier(fence_space::global_and_local);
                    float mvb = 0;
                    for (int n = 0; n < coupledConstraintsCount; n++)
                    {
                        int index = n * numConstraintsThreads + threadIndex;
                        int c1    = gm_coupledConstraintsIndices[index];

                        mvb = mvb + gm_matrixA[index] * sm_buffer[c1 + c_threadsPerBlock * (rec % 2)];
                    }

                    sm_buffer[threadInBlock + c_threadsPerBlock * ((rec + 1) % 2)] = mvb;
                    sol                                                            = sol + mvb;
                }
            }

            // Add corrections to Lagrange multipliers
            float sqrtmu_sol = sqrtReducedMass * sol;
            lagrangeScaled += sqrtmu_sol;

            // Save updated coordinates for the next iteration
            // Dummy constraints are skipped
            if (!isDummyThread)
            {
                Float3 tmp = rc * sqrtmu_sol;
                atomicFetchAdd(gm_xp[i][XX], -tmp[XX] * inverseMassi);
                atomicFetchAdd(gm_xp[i][YY], -tmp[YY] * inverseMassi);
                atomicFetchAdd(gm_xp[i][ZZ], -tmp[ZZ] * inverseMassi);
                atomicFetchAdd(gm_xp[j][XX], tmp[XX] * inverseMassj);
                atomicFetchAdd(gm_xp[j][YY], tmp[YY] * inverseMassj);
                atomicFetchAdd(gm_xp[j][ZZ], tmp[ZZ] * inverseMassj);
            }
        }

        // Updating particle velocities for all but dummy threads
        if constexpr (updateVelocities)
        {
            if (!isDummyThread)
            {
                Float3 tmp = rc * invdt * lagrangeScaled;
                atomicFetchAdd(gm_v[i][XX], -tmp[XX] * inverseMassi);
                atomicFetchAdd(gm_v[i][YY], -tmp[YY] * inverseMassi);
                atomicFetchAdd(gm_v[i][ZZ], -tmp[ZZ] * inverseMassi);
                atomicFetchAdd(gm_v[j][XX], tmp[XX] * inverseMassj);
                atomicFetchAdd(gm_v[j][YY], tmp[YY] * inverseMassj);
                atomicFetchAdd(gm_v[j][ZZ], tmp[ZZ] * inverseMassj);
            }
        }

        if constexpr (computeVirial)
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
            // We reuse the same shared memory buffer, so we make sure we don't need its old values:
            itemIdx.barrier(fence_space::local_space);
            float mult                                       = targetLength * lagrangeScaled;
            sm_buffer[0 * c_threadsPerBlock + threadInBlock] = mult * rc[XX] * rc[XX];
            sm_buffer[1 * c_threadsPerBlock + threadInBlock] = mult * rc[XX] * rc[YY];
            sm_buffer[2 * c_threadsPerBlock + threadInBlock] = mult * rc[XX] * rc[ZZ];
            sm_buffer[3 * c_threadsPerBlock + threadInBlock] = mult * rc[YY] * rc[YY];
            sm_buffer[4 * c_threadsPerBlock + threadInBlock] = mult * rc[YY] * rc[ZZ];
            sm_buffer[5 * c_threadsPerBlock + threadInBlock] = mult * rc[ZZ] * rc[ZZ];

            itemIdx.barrier(fence_space::local_space);
            // This casts unsigned into signed integers to avoid clang warnings
            const int tib          = static_cast<int>(threadInBlock);
            const int blockSize    = static_cast<int>(c_threadsPerBlock);
            const int subGroupSize = itemIdx.get_sub_group().get_max_local_range()[0];

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
                        sm_buffer[d * blockSize + tib] += sm_buffer[d * blockSize + (tib + dividedAt)];
                    }
                }
                if (dividedAt > subGroupSize / 2)
                {
                    itemIdx.barrier(fence_space::local_space);
                }
                else
                {
                    subGroupBarrier(itemIdx);
                }
            }
            // First 6 threads in the block add the 6 components of virial to the global memory address
            if (tib < 6)
            {
                atomicFetchAdd(gm_virialScaled[tib], sm_buffer[tib * blockSize]);
            }
        }
    };
}

// SYCL 1.2.1 requires providing a unique type for a kernel. Should not be needed for SYCL2020.
template<bool updateVelocities, bool computeVirial, bool haveCoupledConstraints>
class LincsKernelName;

template<bool updateVelocities, bool computeVirial, bool haveCoupledConstraints, class... Args>
static void launchLincsKernel(const DeviceStream& deviceStream, const int numConstraintsThreads, Args&&... args)
{
    // Should not be needed for SYCL2020.
    using kernelNameType = LincsKernelName<updateVelocities, computeVirial, haveCoupledConstraints>;

    const sycl::nd_range<1> rangeAllLincs(numConstraintsThreads, c_threadsPerBlock);
    sycl::queue             q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = lincsKernel<updateVelocities, computeVirial, haveCoupledConstraints>(
                cgh, numConstraintsThreads, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(rangeAllLincs, kernel);
    });
}

/*! \brief Select templated kernel and launch it. */
template<class... Args>
static inline void
launchLincsKernel(bool updateVelocities, bool computeVirial, bool haveCoupledConstraints, Args&&... args)
{
    dispatchTemplatedFunction(
            [&](auto updateVelocities_, auto computeVirial_, auto haveCoupledConstraints_) {
                return launchLincsKernel<updateVelocities_, computeVirial_, haveCoupledConstraints_>(
                        std::forward<Args>(args)...);
            },
            updateVelocities,
            computeVirial,
            haveCoupledConstraints);
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
    launchLincsKernel(updateVelocities,
                      computeVirial,
                      kernelParams->haveCoupledConstraints,
                      deviceStream,
                      kernelParams->numConstraintsThreads,
                      kernelParams->d_constraints.get_pointer(),
                      kernelParams->d_constraintsTargetLengths.get_pointer(),
                      kernelParams->d_coupledConstraintsCounts.get_pointer(),
                      kernelParams->d_coupledConstraintsIndices.get_pointer(),
                      kernelParams->d_massFactors.get_pointer(),
                      kernelParams->d_matrixA.get_pointer(),
                      kernelParams->d_inverseMasses.get_pointer(),
                      kernelParams->numIterations,
                      kernelParams->expansionOrder,
                      d_x.get_pointer(),
                      d_xp.get_pointer(),
                      invdt,
                      d_v.get_pointer(),
                      kernelParams->d_virialScaled.get_pointer(),
                      kernelParams->pbcAiuc);
    return;
}

} // namespace gmx
