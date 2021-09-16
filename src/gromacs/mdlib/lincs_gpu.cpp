/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/lincs_gpu_internal.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

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
                             GpuApiCallBehavior::Sync,
                             nullptr);

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
    deviceContext_(deviceContext), deviceStream_(deviceStream)
{
    GMX_RELEASE_ASSERT(bool(GMX_GPU_CUDA) || bool(GMX_GPU_SYCL),
                       "LINCS GPU is only implemented in CUDA and SYCL.");
    kernelParams_.numIterations  = numIterations;
    kernelParams_.expansionOrder = expansionOrder;

    static_assert(sizeof(real) == sizeof(float),
                  "Real numbers should be in single precision in GPU code.");
    static_assert(
            gmx::isPowerOfTwo(c_threadsPerBlock),
            "Number of threads per block should be a power of two in order for reduction to work.");

    allocateDeviceBuffer(&kernelParams_.d_virialScaled, 6, deviceContext_);
    h_virialScaled_.resize(6);

    // The data arrays should be expanded/reallocated on first call of set() function.
    numConstraintsThreadsAlloc_ = 0;
    numAtomsAlloc_              = 0;
}

LincsGpu::~LincsGpu()
{
    freeDeviceBuffer(&kernelParams_.d_virialScaled);

    if (numConstraintsThreadsAlloc_ > 0)
    {
        freeDeviceBuffer(&kernelParams_.d_constraints);
        freeDeviceBuffer(&kernelParams_.d_constraintsTargetLengths);

        freeDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts);
        freeDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices);
        freeDeviceBuffer(&kernelParams_.d_massFactors);
        freeDeviceBuffer(&kernelParams_.d_matrixA);
    }
    if (numAtomsAlloc_ > 0)
    {
        freeDeviceBuffer(&kernelParams_.d_inverseMasses);
    }
}

//! Helper type for discovering coupled constraints
struct AtomsAdjacencyListElement
{
    AtomsAdjacencyListElement(const int indexOfSecondConstrainedAtom,
                              const int indexOfConstraint,
                              const int signFactor) :
        indexOfSecondConstrainedAtom_(indexOfSecondConstrainedAtom),
        indexOfConstraint_(indexOfConstraint),
        signFactor_(signFactor)
    {
    }
    //! The index of the other atom constrained to this atom.
    int indexOfSecondConstrainedAtom_;
    //! The index of this constraint in the container of constraints.
    int indexOfConstraint_;
    /*! \brief A multiplicative factor that indicates the relative
     * order of the atoms in the atom list.
     *
     * Used for computing the mass factor of this constraint
     * relative to any coupled constraints. */
    int signFactor_;
};
//! Constructs and returns an atom constraint adjacency list
static std::vector<std::vector<AtomsAdjacencyListElement>>
constructAtomsAdjacencyList(const int numAtoms, ArrayRef<const int> iatoms)
{
    const int                                           stride         = 1 + NRAL(F_CONSTR);
    const int                                           numConstraints = iatoms.ssize() / stride;
    std::vector<std::vector<AtomsAdjacencyListElement>> atomsAdjacencyList(numAtoms);
    for (int c = 0; c < numConstraints; c++)
    {
        int a1 = iatoms[stride * c + 1];
        int a2 = iatoms[stride * c + 2];

        // Each constraint will be represented as a tuple, containing index of the second
        // constrained atom, index of the constraint and a sign that indicates the order of atoms in
        // which they are listed. Sign is needed to compute the mass factors.
        atomsAdjacencyList[a1].emplace_back(a2, c, +1);
        atomsAdjacencyList[a2].emplace_back(a1, c, -1);
    }

    return atomsAdjacencyList;
}

/*! \brief Helper function to go through constraints recursively.
 *
 *  For each constraint, counts the number of coupled constraints and stores the value in \p numCoupledConstraints array.
 *  This information is used to split the array of constraints between thread blocks on a GPU so there is no
 *  coupling between constraints from different thread blocks. After the \p numCoupledConstraints array is filled, the
 *  value \p numCoupledConstraints[c] should be equal to the number of constraints that are coupled to \p c and located
 *  after it in the constraints array.
 *
 * \param[in]     a                   Atom index.
 * \param[in,out] numCoupledConstraints  Indicates if the constraint was already counted and stores
 *                                    the number of constraints (i) connected to it and (ii) located
 *                                    after it in memory. This array is filled by this recursive function.
 *                                    For a set of coupled constraints, only for the first one in this list
 *                                    the number of consecutive coupled constraints is needed: if there is
 *                                    not enough space for this set of constraints in the thread block,
 *                                    the group has to be moved to the next one.
 * \param[in]     atomsAdjacencyList  Stores information about connections between atoms.
 */
inline int countCoupled(int           a,
                        ArrayRef<int> numCoupledConstraints,
                        ArrayRef<const std::vector<AtomsAdjacencyListElement>> atomsAdjacencyList)

{
    int counted = 0;
    for (const auto& adjacentAtom : atomsAdjacencyList[a])
    {
        const int c2 = adjacentAtom.indexOfConstraint_;
        if (numCoupledConstraints[c2] == -1)
        {
            numCoupledConstraints[c2] = 0; // To indicate we've been here
            counted += 1
                       + countCoupled(adjacentAtom.indexOfSecondConstrainedAtom_,
                                      numCoupledConstraints,
                                      atomsAdjacencyList);
        }
    }
    return counted;
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
inline void addWithCoupled(ArrayRef<const int>                                    iatoms,
                           const int                                              stride,
                           ArrayRef<const std::vector<AtomsAdjacencyListElement>> atomsAdjacencyList,
                           ArrayRef<int>                                          splitMap,
                           const int                                              c,
                           int*                                                   currentMapIndex)
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

/*! \brief Computes and returns how many constraints are coupled to each constraint
 *
 * Needed to introduce splits in data so that all coupled constraints will be computed in a
 * single GPU block. The position \p c of the vector \p numCoupledConstraints should have the number
 * of constraints that are coupled to a constraint \p c and are after \p c in the vector. Only
 * first index of the connected group of the constraints is needed later in the code, hence the
 * numCoupledConstraints vector is also used to keep track if the constrain was already counted.
 */
static std::vector<int> countNumCoupledConstraints(ArrayRef<const int> iatoms,
                                                   ArrayRef<const std::vector<AtomsAdjacencyListElement>> atomsAdjacencyList)
{
    const int        stride         = 1 + NRAL(F_CONSTR);
    const int        numConstraints = iatoms.ssize() / stride;
    std::vector<int> numCoupledConstraints(numConstraints, -1);
    for (int c = 0; c < numConstraints; c++)
    {
        const int a1 = iatoms[stride * c + 1];
        const int a2 = iatoms[stride * c + 2];
        if (numCoupledConstraints[c] == -1)
        {
            numCoupledConstraints[c] = countCoupled(a1, numCoupledConstraints, atomsAdjacencyList)
                                       + countCoupled(a2, numCoupledConstraints, atomsAdjacencyList);
        }
    }

    return numCoupledConstraints;
}

bool LincsGpu::isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop)
{
    for (const gmx_moltype_t& molType : mtop.moltype)
    {
        ArrayRef<const int> iatoms    = molType.ilist[F_CONSTR].iatoms;
        const auto atomsAdjacencyList = constructAtomsAdjacencyList(molType.atoms.nr, iatoms);
        // Compute, how many constraints are coupled to each constraint
        const auto numCoupledConstraints = countNumCoupledConstraints(iatoms, atomsAdjacencyList);
        for (const int numCoupled : numCoupledConstraints)
        {
            if (numCoupled > c_threadsPerBlock)
            {
                return false;
            }
        }
    }

    return true;
}

void LincsGpu::set(const InteractionDefinitions& idef, const int numAtoms, const real* invmass)
{
    GMX_RELEASE_ASSERT(bool(GMX_GPU_CUDA) || bool(GMX_GPU_SYCL),
                       "LINCS GPU is only implemented in CUDA and SYCL.");
    // List of constrained atoms (CPU memory)
    std::vector<AtomPair> constraintsHost;
    // Equilibrium distances for the constraints (CPU)
    std::vector<float> constraintsTargetLengthsHost;
    // Number of constraints, coupled with the current one (CPU)
    std::vector<int> coupledConstraintsCountsHost;
    // List of coupled with the current one (CPU)
    std::vector<int> coupledConstraintsIndicesHost;
    // Mass factors (CPU)
    std::vector<float> massFactorsHost;

    // List of constrained atoms in local topology
    ArrayRef<const int> iatoms         = idef.il[F_CONSTR].iatoms;
    const int           stride         = NRAL(F_CONSTR) + 1;
    const int           numConstraints = idef.il[F_CONSTR].size() / stride;

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
    AtomPair pair;
    pair.i = -1;
    pair.j = -1;
    constraintsHost.resize(kernelParams_.numConstraintsThreads, pair);
    std::fill(constraintsHost.begin(), constraintsHost.end(), pair);
    constraintsTargetLengthsHost.resize(kernelParams_.numConstraintsThreads, 0.0);
    std::fill(constraintsTargetLengthsHost.begin(), constraintsTargetLengthsHost.end(), 0.0);
    for (int c = 0; c < numConstraints; c++)
    {
        int a1   = iatoms[stride * c + 1];
        int a2   = iatoms[stride * c + 2];
        int type = iatoms[stride * c];

        AtomPair pair;
        pair.i                                    = a1;
        pair.j                                    = a2;
        constraintsHost[splitMap[c]]              = pair;
        constraintsTargetLengthsHost[splitMap[c]] = idef.iparams[type].constr.dA;
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
    int  maxCoupledConstraints             = 0;
    bool maxCoupledConstraintsHasIncreased = false;
    for (int c = 0; c < numConstraints; c++)
    {
        int a1 = iatoms[stride * c + 1];
        int a2 = iatoms[stride * c + 2];

        // Constraint 'c' is counted twice, but it should be excluded altogether. Hence '-2'.
        int nCoupledConstraints = atomsAdjacencyList[a1].size() + atomsAdjacencyList[a2].size() - 2;

        if (nCoupledConstraints > maxCoupledConstraints)
        {
            maxCoupledConstraints             = nCoupledConstraints;
            maxCoupledConstraintsHasIncreased = true;
        }
    }

    kernelParams_.haveCoupledConstraints = (maxCoupledConstraints > 0);

    coupledConstraintsCountsHost.resize(kernelParams_.numConstraintsThreads, 0);
    coupledConstraintsIndicesHost.resize(maxCoupledConstraints * kernelParams_.numConstraintsThreads, -1);
    massFactorsHost.resize(maxCoupledConstraints * kernelParams_.numConstraintsThreads, -1);

    for (int c1 = 0; c1 < numConstraints; c1++)
    {
        coupledConstraintsCountsHost[splitMap[c1]] = 0;
        int c1a1                                   = iatoms[stride * c1 + 1];
        int c1a2                                   = iatoms[stride * c1 + 2];

        // Constraints, coupled through the first atom.
        int c2a1 = c1a1;
        for (const auto& atomAdjacencyList : atomsAdjacencyList[c1a1])
        {
            int c2 = atomAdjacencyList.indexOfConstraint_;

            if (c1 != c2)
            {
                int c2a2  = atomAdjacencyList.indexOfSecondConstrainedAtom_;
                int sign  = atomAdjacencyList.signFactor_;
                int index = kernelParams_.numConstraintsThreads
                                    * coupledConstraintsCountsHost[splitMap[c1]]
                            + splitMap[c1];
                int threadBlockStarts = splitMap[c1] - splitMap[c1] % c_threadsPerBlock;

                coupledConstraintsIndicesHost[index] = splitMap[c2] - threadBlockStarts;

                int center = c1a1;

                float sqrtmu1 = 1.0 / std::sqrt(invmass[c1a1] + invmass[c1a2]);
                float sqrtmu2 = 1.0 / std::sqrt(invmass[c2a1] + invmass[c2a2]);

                massFactorsHost[index] = -sign * invmass[center] * sqrtmu1 * sqrtmu2;

                coupledConstraintsCountsHost[splitMap[c1]]++;
            }
        }

        // Constraints, coupled through the second atom.
        c2a1 = c1a2;
        for (const auto& atomAdjacencyList : atomsAdjacencyList[c1a2])
        {
            int c2 = atomAdjacencyList.indexOfConstraint_;

            if (c1 != c2)
            {
                int c2a2  = atomAdjacencyList.indexOfSecondConstrainedAtom_;
                int sign  = atomAdjacencyList.signFactor_;
                int index = kernelParams_.numConstraintsThreads
                                    * coupledConstraintsCountsHost[splitMap[c1]]
                            + splitMap[c1];
                int threadBlockStarts = splitMap[c1] - splitMap[c1] % c_threadsPerBlock;

                coupledConstraintsIndicesHost[index] = splitMap[c2] - threadBlockStarts;

                int center = c1a2;

                float sqrtmu1 = 1.0 / std::sqrt(invmass[c1a1] + invmass[c1a2]);
                float sqrtmu2 = 1.0 / std::sqrt(invmass[c2a1] + invmass[c2a2]);

                massFactorsHost[index] = sign * invmass[center] * sqrtmu1 * sqrtmu2;

                coupledConstraintsCountsHost[splitMap[c1]]++;
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
                             maxCoupledConstraints * kernelParams_.numConstraintsThreads,
                             deviceContext_);
        allocateDeviceBuffer(&kernelParams_.d_massFactors,
                             maxCoupledConstraints * kernelParams_.numConstraintsThreads,
                             deviceContext_);
        allocateDeviceBuffer(&kernelParams_.d_matrixA,
                             maxCoupledConstraints * kernelParams_.numConstraintsThreads,
                             deviceContext_);
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
                       constraintsHost.data(),
                       0,
                       kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_constraintsTargetLengths,
                       constraintsTargetLengthsHost.data(),
                       0,
                       kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsCounts,
                       coupledConstraintsCountsHost.data(),
                       0,
                       kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_coupledConstraintsIndices,
                       coupledConstraintsIndicesHost.data(),
                       0,
                       maxCoupledConstraints * kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    copyToDeviceBuffer(&kernelParams_.d_massFactors,
                       massFactorsHost.data(),
                       0,
                       maxCoupledConstraints * kernelParams_.numConstraintsThreads,
                       deviceStream_,
                       GpuApiCallBehavior::Sync,
                       nullptr);

    GMX_RELEASE_ASSERT(invmass != nullptr, "Masses of atoms should be specified.\n");
    copyToDeviceBuffer(
            &kernelParams_.d_inverseMasses, invmass, 0, numAtoms, deviceStream_, GpuApiCallBehavior::Sync, nullptr);
}

} // namespace gmx
