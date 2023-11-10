/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief
 * Listed forces GPU implementation
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTED_FORCES_GPU_KERNEL_CUH
#define NBLIB_LISTED_FORCES_GPU_KERNEL_CUH

#include <cstdio>

#include "nblib/listed_forces/dataflow.hpp"

namespace nblib
{
// Number of CUDA threads in a block
constexpr static int c_threadsPerBlock = 256;

template<class InteractionType, class T>
__device__ void computePotential(KernelEnergy<T> energy, T* gm_potentials)
{
    typedef cub::BlockReduce<T, c_threadsPerBlock> BlockReduce;
    __shared__ typename BlockReduce::TempStorage   temp_storage;

    constexpr bool hasPair =
            Contains<InteractionType, PairListedTypes>{} || HasPairAggregate<InteractionType>{};

    BlockReduce reduce(temp_storage);
    energy.carrier() = reduce.Sum(energy.carrier());
    __syncthreads();

    if constexpr (HasThreeCenterAggregate<InteractionType>{})
    {
        energy.threeCenterAggregate() = reduce.Sum(energy.threeCenterAggregate());
        __syncthreads();
    }
    if constexpr (HasTwoCenterAggregate<InteractionType>{})
    {
        energy.twoCenterAggregate() = reduce.Sum(energy.twoCenterAggregate());
        __syncthreads();
    }
    if constexpr (hasPair)
    {
        energy.eVdw() = reduce.Sum(energy.eVdw());
        __syncthreads();
        energy.eCoul() = reduce.Sum(energy.eCoul());
        __syncthreads();
    }

    if (threadIdx.x == 0)
    {
        if constexpr (IsAggregate<InteractionType>{})
        {
            constexpr int carrierIndex =
                    FindIndex<typename InteractionType::CarrierType, GpuListedTypes>{};
            atomicAdd(gm_potentials + carrierIndex, energy.carrier());
        }
        else
        {
            constexpr int compileTimeIndex = FindIndex<InteractionType, GpuListedTypes>{};
            atomicAdd(gm_potentials + compileTimeIndex, energy.carrier());
        }
        if constexpr (HasThreeCenterAggregate<InteractionType>{})
        {
            int index = FindIndex<typename InteractionType::ThreeCenterAggregateType, GpuListedTypes>{};
            atomicAdd(gm_potentials + index, energy.threeCenterAggregate());
        }
        if constexpr (HasTwoCenterAggregate<InteractionType>{})
        {
            int index = FindIndex<typename InteractionType::TwoCenterAggregateType, GpuListedTypes>{};
            atomicAdd(gm_potentials + index, energy.twoCenterAggregate());
        }
        if constexpr (hasPair)
        {
            atomicAdd(gm_potentials + GpuVdwIndex{}, energy.eVdw());
            atomicAdd(gm_potentials + GpuCoulombIndex{}, energy.eCoul());
        }
    }
}

template<bool calcEner, class T, class ShiftForce, class Pbc>
__global__ void computeListedForces(const int*               numInteractions,
                                    const int*               blockTypes,
                                    const int*               indices,
                                    const void* const*       parametersA,
                                    const void* const*       parametersB,
                                    const int*               interactionOffsets,
                                    const int*               indexOffsets,
                                    const util::array<T, 4>* gm_xyzq,
                                    util::array<T, 3>*       gm_forces,
                                    ShiftForce*              gm_shiftForces,
                                    T*                       gm_potentials,
                                    Pbc                      pbc)
{
    constexpr int numTypes = TypeListSize<GpuListedTypes>{};

    int tid              = blockDim.x * blockIdx.x + threadIdx.x;
    int typeIndex        = blockTypes[blockIdx.x];
    int interactionIndex = tid - interactionOffsets[typeIndex];

    const int* typeIndices = indices + indexOffsets[typeIndex];

    extern __shared__ char sm_dynamicShmem[];
    ShiftForce*            sm_shiftForces = reinterpret_cast<ShiftForce*>(sm_dynamicShmem);

    if constexpr (HaveShifts<ShiftForce>{})
    {
        if (threadIdx.x < gmx::c_numShiftVectors)
        {
            sm_shiftForces[threadIdx.x] = { 0, 0, 0 };
        }
        __syncthreads();
    }

    auto dispatch = [typeIndices,
                     interactionIndex,
                     &numInteractions,
                     gm_xyzq,
                     &parametersA,
                     &parametersB,
                     &gm_forces,
                     gm_potentials,
                     sm_shiftForces,
                     &pbc](auto compileTimeIndex) {
        using InteractionType = TypeListElement_t<compileTimeIndex, GpuListedTypes>;
        KernelEnergy<T> energy;
        if (interactionIndex < numInteractions[compileTimeIndex])
        {
            energy = dispatchInteraction(
                    reinterpret_cast<const IndexArray<NCenter<InteractionType>{} + 1>*>(
                            typeIndices)[interactionIndex],
                    static_cast<const InteractionType*>(parametersA[compileTimeIndex]),
                    static_cast<const InteractionType*>(parametersB[compileTimeIndex]),
                    gm_xyzq,
                    NoFepLambdaType{},
                    reinterpret_cast<DeviceTag<util::array<T, 3>>**>(&gm_forces),
                    reinterpret_cast<DeviceTag<ShiftForce>*>(sm_shiftForces),
                    pbc);
        }

        if constexpr (calcEner)
        {
            computePotential<InteractionType>(energy, gm_potentials);
        }

        return 0;
    };

    // instantiate switch(typeIndex) with SwitchCases as the possible cases,
    // then call dispatch with the compile-time constant value of typeIndex
    using SwitchCases = std::make_integer_sequence<int, numTypes>;
    createSwitch(typeIndex, SwitchCases{}, dispatch);

    if constexpr (HaveShifts<ShiftForce>{})
    {
        __syncthreads();
        if (threadIdx.x < gmx::c_numShiftVectors)
        {
            atomicAdd(&gm_shiftForces[threadIdx.x][0], sm_shiftForces[threadIdx.x][0]);
            atomicAdd(&gm_shiftForces[threadIdx.x][1], sm_shiftForces[threadIdx.x][1]);
            atomicAdd(&gm_shiftForces[threadIdx.x][2], sm_shiftForces[threadIdx.x][2]);
        }
    }
}

} // namespace nblib

#endif // NBLIB_LISTED_FORCES_GPU_KERNEL_CUH
