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
 * \brief Implementations of LINCS GPU class
 *
 * This file contains back-end agnostic implementation of LINCS GPU class.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "lincs_gpu.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>

#include <algorithm>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/constraint_gpu_helpers.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs_constraint_group_sizes.h"
#include "gromacs/mdlib/lincs_gpu_internal.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/vec.h"


namespace gmx
{

namespace
{

/*! \brief Resizes \p v to \p size elements and sets every element to \p value.
 *
 * This differs from \c std::vector::resize(size, value), which only initializes the elements
 * that are newly added and leaves the ones the vector already holds untouched. The LINCS host
 * buffers persist across LincsGpu::set() calls and are only sparsely written afterwards, so
 * every element has to be reset to its sentinel value on each call, whether the buffer grows,
 * shrinks or keeps its size.
 */
template<typename T, typename Allocator>
void resizeAndFill(std::vector<T, Allocator>* v, const std::size_t size, const T& value)
{
    v->resize(size);
    std::fill(v->begin(), v->end(), value);
}

} // namespace

void LincsGpu::apply(const DeviceBuffer<Float3>& d_x,
                     DeviceBuffer<Float3>        d_xp,
                     const bool                  updateVelocities,
                     DeviceBuffer<Float3>        d_v,
                     const real                  invdt,
                     const bool                  computeVirial,
                     tensor                      virialScaled,
                     const PbcAiuc&              pbcAiuc)
{
    // Early exit if no constraints
    if (kernelParams_.numConstraintsThreads == 0)
    {
        return;
    }

    if (computeVirial)
    {
        // Fill with zeros so the values can be reduced to it
        // Only 6 values are needed because virial is symmetrical
        clearDeviceBufferAsync(&kernelParams_.d_virialScaled, 0, 6, deviceStream_);
    }

    kernelParams_.pbcAiuc = pbcAiuc;

    launchLincsGpuKernel(
            &kernelParams_, d_x, d_xp, updateVelocities, d_v, invdt, computeVirial, deviceStream_);

    if (computeVirial)
    {
        // Copy LINCS virial data and add it to the common virial
        copyFromDeviceBuffer(h_virialScaled_.data(),
                             &kernelParams_.d_virialScaled,
                             0,
                             6,
                             deviceStream_,
                             GpuApiCallBehavior::Async,
                             nullptr);
        deviceStream_.synchronize();

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
}

LincsGpu::LincsGpu(int                  numIterations,
                   int                  expansionOrder,
                   const DeviceContext& deviceContext,
                   const DeviceStream&  deviceStream) :
    deviceContext_(deviceContext),
    deviceStream_(deviceStream),
    h_virialScaled_(6, HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported }),
    h_constraints_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported }),
    h_constraintsTargetLengths_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported }),
    h_coupledConstraintsCounts_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported }),
    h_coupledConstraintsIndices_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported }),
    h_massFactors_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported }),
    h_constraintGroupSize_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported })
{
    GMX_RELEASE_ASSERT(GMX_GPU && !GMX_GPU_OPENCL, "LINCS GPU is not implemented in OPENCL.");
    kernelParams_.numIterations  = numIterations;
    kernelParams_.expansionOrder = expansionOrder;

    static_assert(sizeof(real) == sizeof(float),
                  "Real numbers should be in single precision in GPU code.");
    static_assert(
            gmx::isPowerOfTwo(c_threadsPerBlock),
            "Number of threads per block should be a power of two in order for reduction to work.");

    allocateDeviceBuffer(&kernelParams_.d_virialScaled, 6, deviceContext_);

    // The data arrays should be expanded/reallocated on first call of set() function.
    numConstraintsThreadsAlloc_ = 0;
    numAtomsAlloc_              = 0;
}

LincsGpu::~LincsGpu()
{
    try
    {
        // Wait for all the tasks to complete before freeing the memory. See #4519.
        deviceStream_.synchronize();

        freeDeviceBuffer(&kernelParams_.d_virialScaled);

        if (numConstraintsThreadsAlloc_ > 0)
        {
            freeDeviceBuffer(&kernelParams_.d_constraints);
            freeDeviceBuffer(&kernelParams_.d_constraintsTargetLengths);

            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts);
            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices);
            freeDeviceBuffer(&kernelParams_.d_massFactors);
            freeDeviceBuffer(&kernelParams_.d_matrixA);
            if constexpr (GMX_GPU_HIP)
            {
                freeDeviceBuffer(&kernelParams_.d_constraintGroupsSizes);
            }
        }
        if (numAtomsAlloc_ > 0)
        {
            freeDeviceBuffer(&kernelParams_.d_inverseMasses);
        }
    }
    catch (gmx::InternalError& e)
    {
        fprintf(stderr, "Internal error in destructor of LincsGpu: %s\n", e.what());
    }
}

/*! \brief Add constraint to \p splitMap with all constraints coupled to it.
 *
 *  Adds the constraint \p c from the constrain list \p iatoms to the map \p splitMap
 *  if it was not yet added. Then goes through all the constraints coupled to \p c
 *  and calls itself recursively. This ensures that all the coupled constraints will
 *  be added to neighboring locations in the final data structures on the device,
 *  hence mapping all coupled constraints to the same thread block. A value of -1 in
 *  the \p splitMap is used to flag that constraint was not yet added to the \p splitMap.
 *
 * \param[in]     iatoms              The list of constraints.
 * \param[in]     stride              Number of elements per constraint in \p iatoms.
 * \param[in]     atomsAdjacencyList  Information about connections between atoms.
 * \param[out]    splitMap            Map of sequential constraint indexes to indexes to be on the device
 * \param[in]     c                   Sequential index for constraint to consider adding.
 * \param[in,out] currentMapIndex     The rolling index for the constraints mapping.
 */
inline void addWithCoupled(ArrayRef<const int>                                iatoms,
                           const int                                          stride,
                           const gmx::ListOfLists<AtomsAdjacencyListElement>& atomsAdjacencyList,
                           ArrayRef<int>                                      splitMap,
                           const int                                          c,
                           int*                                               currentMapIndex)
{
    if (splitMap[c] == -1)
    {
        splitMap[c] = *currentMapIndex;
        (*currentMapIndex)++;

        // Constraints, coupled through both atoms.
        for (int atomIndexInConstraint = 0; atomIndexInConstraint < 2; atomIndexInConstraint++)
        {
            const int a1 = iatoms[stride * c + 1 + atomIndexInConstraint];
            for (const auto& adjacentAtom : atomsAdjacencyList[a1])
            {
                const int c2 = adjacentAtom.indexOfConstraint_;
                if (c2 != c)
                {
                    addWithCoupled(iatoms, stride, atomsAdjacencyList, splitMap, c2, currentMapIndex);
                }
            }
        }
    }
}

bool LincsGpu::isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop)
{
    return ::isNumCoupledConstraintsSupported(mtop, c_threadsPerBlock);
}

void LincsGpu::set(const InteractionDefinitions& idef, int numAtoms, const ArrayRef<const real> invmass)
{
    GMX_ASSERT(!(numAtoms == 0 && !idef.il[InteractionFunction::Constraints].empty()),
               "The number of atoms needs to be > 0 if there are constraints in the domain.");

    GMX_RELEASE_ASSERT(GMX_GPU && !GMX_GPU_OPENCL, "LINCS GPU is not implemented in OPENCL.");

    // List of constrained atoms in local topology
    ArrayRef<const int> iatoms         = idef.il[InteractionFunction::Constraints].iatoms;
    const int           stride         = NRAL(InteractionFunction::Constraints) + 1;
    const int           numConstraints = idef.il[InteractionFunction::Constraints].size() / stride;

    // Early exit if no constraints
    if (numConstraints == 0)
    {
        kernelParams_.numConstraintsThreads = 0;
        return;
    }

    // Construct the adjacency list, a useful intermediate structure
    const auto atomsAdjacencyList = constructAtomsAdjacencyList(numAtoms, iatoms);

    // Compute, how many constraints are coupled to each constraint
    const auto numCoupledConstraints = countNumCoupledConstraints(iatoms, atomsAdjacencyList);

    // Map of splits in the constraints data. For each 'old' constraint index gives 'new' which
    // takes into account the empty spaces which might be needed in the end of each thread block.
    std::vector<int> splitMap(numConstraints, -1);
    int              currentMapIndex = 0;
    for (int c = 0; c < numConstraints; c++)
    {
        // Check if coupled constraints all fit in one block
        if (numCoupledConstraints[c] > c_threadsPerBlock)
        {
            gmx_fatal(FARGS,
                      "Maximum number of coupled constraints (%d) exceeds the size of the CUDA "
                      "thread block (%d). Most likely, you are trying to use the GPU version of "
                      "LINCS with constraints on all-bonds, which is not supported for large "
                      "molecules. When compatible with the force field and integration settings, "
                      "using constraints on H-bonds only.",
                      numCoupledConstraints[c],
                      c_threadsPerBlock);
        }
        if (currentMapIndex / c_threadsPerBlock != (currentMapIndex + numCoupledConstraints[c]) / c_threadsPerBlock)
        {
            currentMapIndex = ((currentMapIndex / c_threadsPerBlock) + 1) * c_threadsPerBlock;
        }
        addWithCoupled(iatoms, stride, atomsAdjacencyList, splitMap, c, &currentMapIndex);
    }

    kernelParams_.numConstraintsThreads =
            currentMapIndex + c_threadsPerBlock - currentMapIndex % c_threadsPerBlock;
    GMX_RELEASE_ASSERT(kernelParams_.numConstraintsThreads % c_threadsPerBlock == 0,
                       "Number of threads should be a multiple of the block size");

    // Initialize constraints and their target indexes taking into account the splits in the data arrays.
    {
        AtomPair pair;
        pair.i = -1;
        pair.j = -1;
        resizeAndFill(&h_constraints_, kernelParams_.numConstraintsThreads, pair);
        resizeAndFill(&h_constraintsTargetLengths_, kernelParams_.numConstraintsThreads, 0.0F);
    }


    const int gmx_unused numOmpThreads = gmx_omp_nthreads_get(ModuleMultiThread::Lincs);
#pragma omp parallel for num_threads(numOmpThreads) schedule(static)
    for (int c = 0; c < numConstraints; c++)
    {
        int a1   = iatoms[stride * c + 1];
        int a2   = iatoms[stride * c + 2];
        int type = iatoms[stride * c];

        AtomPair localPair;
        localPair.i                              = a1;
        localPair.j                              = a2;
        h_constraints_[splitMap[c]]              = localPair;
        h_constraintsTargetLengths_[splitMap[c]] = idef.iparams[type].constr.dA;
    }

    // The adjacency list of constraints (i.e. the list of coupled constraints for each constraint).
    // We map a single thread to a single constraint, hence each thread 'c' will be using one
    // element from coupledConstraintsCountsHost array, which is the number of constraints coupled
    // to the constraint 'c'. The coupled constraints indexes are placed into the
    // coupledConstraintsIndicesHost array. Latter is organized as a one-dimensional array to ensure
    // good memory alignment. It is addressed as [c + i*numConstraintsThreads], where 'i' goes from
    // zero to the number of constraints coupled to 'c'. 'numConstraintsThreads' is the width of the
    // array --- a number, greater then total number of constraints, taking into account the splits
    // in the constraints array due to the GPU block borders. This number can be adjusted to improve
    // memory access pattern. Mass factors are saved in a similar data structure.
    const int prevMaxCoupledConstraints = maxCoupledConstraints_;
#ifndef _MSC_VER
#    pragma omp parallel for num_threads(numOmpThreads) schedule(static) \
            reduction(max : maxCoupledConstraints_)
#else
// Nothing done. I think is not a good ideal to use '#pragma omp critical' in openmp 2.0
// It caused threads blocking, performance testing is required
#endif
    for (int c = 0; c < numConstraints; c++)
    {
        int a1 = iatoms[stride * c + 1];
        int a2 = iatoms[stride * c + 2];

        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        int nCoupledConstraints = atomsAdjacencyList[a1].size() + atomsAdjacencyList[a2].size() - 2;

        if (nCoupledConstraints > maxCoupledConstraints_)
        {
            maxCoupledConstraints_ = nCoupledConstraints;
        }
    }
    const bool maxCoupledConstraintsHasIncreased = (maxCoupledConstraints_ > prevMaxCoupledConstraints);

    kernelParams_.haveCoupledConstraints = (maxCoupledConstraints_ > 0);

    const size_t coupledCountsSize  = kernelParams_.numConstraintsThreads;
    const size_t coupledIndicesSize = maxCoupledConstraints_ * coupledCountsSize;
    const size_t massFactorsSize    = coupledIndicesSize;

    resizeAndFill(&h_coupledConstraintsCounts_, coupledCountsSize, 0);
    resizeAndFill(&h_coupledConstraintsIndices_, coupledIndicesSize, -1);
    resizeAndFill(&h_massFactors_, massFactorsSize, -1.0F);

    // Only perform constraint re-ordering with HIP
    if constexpr (GMX_GPU_HIP)
    {
        // findConstraintGroupSizes() only writes the entries of the groups it detects, so the
        // whole buffer needs the sentinel value first to make sure no stale groups survive.
        resizeAndFill(&h_constraintGroupSize_, kernelParams_.numConstraintsThreads, -1);
        findConstraintGroupSizes(numConstraints, h_constraints_, h_constraintGroupSize_);
    }

#pragma omp parallel for num_threads(numOmpThreads) schedule(static)
    for (int c1 = 0; c1 < numConstraints; c1++)
    {
        h_coupledConstraintsCounts_[splitMap[c1]] = 0;
        int c1a1                                  = iatoms[stride * c1 + 1];
        int c1a2                                  = iatoms[stride * c1 + 2];

        // Constraints, coupled through the first atom.
        int c2a1 = c1a1;
        for (const auto& atomAdjacencyList : atomsAdjacencyList[c1a1])
        {
            int c2 = atomAdjacencyList.indexOfConstraint_;

            if (c1 != c2)
            {
                int c2a2 = atomAdjacencyList.indexOfSecondConstrainedAtom_;
                int sign = atomAdjacencyList.signFactor_;
                int index = kernelParams_.numConstraintsThreads * h_coupledConstraintsCounts_[splitMap[c1]]
                            + splitMap[c1];
                int threadBlockStarts = splitMap[c1] - splitMap[c1] % c_threadsPerBlock;

                h_coupledConstraintsIndices_[index] = splitMap[c2] - threadBlockStarts;

                int center = c1a1;

                float sqrtmu1 = 1.0 / std::sqrt(invmass[c1a1] + invmass[c1a2]);
                float sqrtmu2 = 1.0 / std::sqrt(invmass[c2a1] + invmass[c2a2]);

                h_massFactors_[index] = -sign * invmass[center] * sqrtmu1 * sqrtmu2;

                h_coupledConstraintsCounts_[splitMap[c1]]++;
            }
        }

        // Constraints, coupled through the second atom.
        c2a1 = c1a2;
        for (const auto& atomAdjacencyList : atomsAdjacencyList[c1a2])
        {
            int c2 = atomAdjacencyList.indexOfConstraint_;

            if (c1 != c2)
            {
                int c2a2 = atomAdjacencyList.indexOfSecondConstrainedAtom_;
                int sign = atomAdjacencyList.signFactor_;
                int index = kernelParams_.numConstraintsThreads * h_coupledConstraintsCounts_[splitMap[c1]]
                            + splitMap[c1];
                int threadBlockStarts = splitMap[c1] - splitMap[c1] % c_threadsPerBlock;

                h_coupledConstraintsIndices_[index] = splitMap[c2] - threadBlockStarts;

                int center = c1a2;

                float sqrtmu1 = 1.0 / std::sqrt(invmass[c1a1] + invmass[c1a2]);
                float sqrtmu2 = 1.0 / std::sqrt(invmass[c2a1] + invmass[c2a2]);

                h_massFactors_[index] = sign * invmass[center] * sqrtmu1 * sqrtmu2;

                h_coupledConstraintsCounts_[splitMap[c1]]++;
            }
        }
    }

    // (Re)allocate the memory, if the number of constraints has increased.
    if ((kernelParams_.numConstraintsThreads > numConstraintsThreadsAlloc_) || maxCoupledConstraintsHasIncreased)
    {
        // Free memory if it was allocated before (i.e. if not the first time here).
        if (numConstraintsThreadsAlloc_ > 0)
        {
            freeDeviceBuffer(&kernelParams_.d_constraints);
            freeDeviceBuffer(&kernelParams_.d_constraintsTargetLengths);

            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts);
            freeDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices);
            freeDeviceBuffer(&kernelParams_.d_massFactors);
            freeDeviceBuffer(&kernelParams_.d_matrixA);
            if constexpr (GMX_GPU_HIP)
            {
                freeDeviceBuffer(&kernelParams_.d_constraintGroupsSizes);
            }
        }

        numConstraintsThreadsAlloc_ = kernelParams_.numConstraintsThreads;

        allocateDeviceBuffer(
                &kernelParams_.d_constraints, kernelParams_.numConstraintsThreads, deviceContext_);
        allocateDeviceBuffer(&kernelParams_.d_constraintsTargetLengths,
                             kernelParams_.numConstraintsThreads,
                             deviceContext_);

        allocateDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts,
                             kernelParams_.numConstraintsThreads,
                             deviceContext_);
        allocateDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices,
                             maxCoupledConstraints_ * kernelParams_.numConstraintsThreads,
                             deviceContext_);
        allocateDeviceBuffer(&kernelParams_.d_massFactors,
                             maxCoupledConstraints_ * kernelParams_.numConstraintsThreads,
                             deviceContext_);
        allocateDeviceBuffer(&kernelParams_.d_matrixA,
                             maxCoupledConstraints_ * kernelParams_.numConstraintsThreads,
                             deviceContext_);
        if constexpr (GMX_GPU_HIP)
        {
            allocateDeviceBuffer(&kernelParams_.d_constraintGroupsSizes,
                                 kernelParams_.numConstraintsThreads,
                                 deviceContext_);
        }
    }

    // (Re)allocate the memory, if the number of atoms has increased.
    if (numAtoms > numAtomsAlloc_)
    {
        if (numAtomsAlloc_ > 0)
        {
            freeDeviceBuffer(&kernelParams_.d_inverseMasses);
        }
        numAtomsAlloc_ = numAtoms;
        allocateDeviceBuffer(&kernelParams_.d_inverseMasses, numAtoms, deviceContext_);
    }

    // Copy data to GPU.
    copyToDeviceBuffer(&kernelParams_.d_constraints,
                       h_constraints_.data(),
                       0,
                       kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_constraintsTargetLengths,
                       h_constraintsTargetLengths_.data(),
                       0,
                       kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts,
                       h_coupledConstraintsCounts_.data(),
                       0,
                       kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices,
                       h_coupledConstraintsIndices_.data(),
                       0,
                       maxCoupledConstraints_ * kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_massFactors,
                       h_massFactors_.data(),
                       0,
                       maxCoupledConstraints_ * kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    if constexpr (GMX_GPU_HIP)
    {
        copyToDeviceBuffer(&kernelParams_.d_constraintGroupsSizes,
                           h_constraintGroupSize_.data(),
                           0,
                           kernelParams_.numConstraintsThreads,
                           deviceStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);
    }

    GMX_RELEASE_ASSERT(!invmass.empty(), "Masses of atoms should be specified.\n");
    copyToDeviceBuffer(&kernelParams_.d_inverseMasses,
                       invmass.data(),
                       0,
                       numAtoms,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
}

} // namespace gmx
