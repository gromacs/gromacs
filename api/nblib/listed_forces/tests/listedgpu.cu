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

#include <cstdio>

#include <numeric>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "gtest/gtest.h"

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"

#include "nblib/listed_forces/calculator.h"
#include "nblib/listed_forces/gpu_kernel.cuh"
#include "nblib/listed_forces/gpu_param.cuh"
#include "nblib/listed_forces/traits.h"

using namespace nblib;

template<class Iterator1, class Iterator2, class Tolerance>
void compareVectors(Iterator1 probeStart, Iterator1 probeEnd, Iterator2 referenceStart, Tolerance tolSetting)
{
    while (probeStart != probeEnd)
    {
        for (int m = 0; m < 3; ++m)
        {
            auto probe = (*probeStart)[m];
            auto ref   = (*referenceStart)[m];
            if (ref != 0.0)
                EXPECT_NEAR((probe - ref) / ref, 0.0, tolSetting);
            else
                EXPECT_NEAR(probe, ref, tolSetting);
        }
        probeStart++;
        referenceStart++;
    }
}

static ListedInteractionData createInteractionData(int numCoordinates)
{
    ListedInteractionData interactionData;

    HarmonicBondType harmonicBond(2.0, 1.0);
    HarmonicAngle    harmonicAngle(2.2, Degrees(91.0));
    ProperDihedral   properDihedral(Degrees(45), 2.3, 1);
    PairLJType       ljPair(C6(0.0026173456), C12(2.634129e-06));

    for (int i = 0; i < numCoordinates; ++i)
    {
        for (int j = i + 1; j < numCoordinates; ++j)
        {
            pickType<HarmonicBondType>(interactionData).indices.push_back({ i, j, 0 });
            pickType<PairLJType>(interactionData).indices.push_back({ i, j, 0 });
        }
    }
    pickType<HarmonicBondType>(interactionData).parametersA.push_back(harmonicBond);
    pickType<PairLJType>(interactionData).parametersA.push_back(ljPair);

    for (int i = 0; i < numCoordinates; ++i)
    {
        for (int j = i + 1; j < numCoordinates; ++j)
        {
            for (int k = j + 1; k < numCoordinates; ++k)
            {
                pickType<HarmonicAngle>(interactionData).indices.push_back({ i, j, k, 0 });
            }
        }
    }
    pickType<HarmonicAngle>(interactionData).parametersA.push_back(harmonicAngle);

    for (int i = 0; i < numCoordinates; ++i)
    {
        for (int j = i + 1; j < numCoordinates; ++j)
        {
            for (int k = j + 1; k < numCoordinates; ++k)
            {
                for (int l = k + 1; l < numCoordinates; ++l)
                {
                    pickType<ProperDihedral>(interactionData).indices.push_back({ i, j, k, l, 0 });
                }
            }
        }
    }
    pickType<ProperDihedral>(interactionData).parametersA.push_back(properDihedral);

    return interactionData;
}

template<class T>
static auto createTestCoordinates(int numParticles)
{
    std::vector<util::array<T, 4>> coordinates(numParticles);
    for (auto& c : coordinates)
    {
        c[0] = drand48();
        c[1] = drand48();
        c[2] = drand48();
        c[3] = 2 * c[2] - 1;
    }

    return coordinates;
}

template<class T, class ShiftForce, class Pbc>
void launchBondedKernel(int                      numInteractionsTot,
                        const int*               numInteractions,
                        const int*               d_blockTypes,
                        const int*               d_indices,
                        const void* const*       parametersA,
                        const void* const*       parametersB,
                        const int*               typeOffsets,
                        const int*               flattenedTypeOffsets,
                        const util::array<T, 4>* d_xyzq,
                        util::array<T, 3>*       d_forces,
                        ShiftForce*              d_shiftForces,
                        T*                       d_potentials,
                        const Pbc&               pbc,
                        cudaStream_t             stream = cudaStreamDefault)
{
    constexpr int numThreads = c_threadsPerBlock;
    int           numBlocks  = iceil(numInteractionsTot, numThreads);

    constexpr int sharedMemSize = HaveShifts<ShiftForce>{} ? sizeof(T) * 3 * gmx::c_numShiftVectors : 0;

    computeListedForces<true><<<numBlocks, numThreads, sharedMemSize, stream>>>(numInteractions,
                                                                                d_blockTypes,
                                                                                d_indices,
                                                                                parametersA,
                                                                                parametersB,
                                                                                typeOffsets,
                                                                                flattenedTypeOffsets,
                                                                                d_xyzq,
                                                                                d_forces,
                                                                                d_shiftForces,
                                                                                d_potentials,
                                                                                pbc);
}


TEST(ListedGpu, compute3Types)
{
    int numParticles = 30;

    ListedInteractionData interactions;
    interactions = createInteractionData(numParticles);

    ListedInteractionDataGpu idataGpu;
    idataGpu.update(interactions);
    std::cout << idataGpu.numInteractionsTot() << std::endl;

    Box           box(1.0);
    PbcHolderAiuc pbc(PbcType::Xyz, box);

    thrust::host_vector<util::array<real, 4>> xyzq = createTestCoordinates<real>(numParticles);

    thrust::host_vector<util::array<real, 3>> forces(xyzq.size(), util::array<real, 3>{ 0.0, 0.0, 0.0 });
    thrust::host_vector<util::array<real, 3>> shiftForces(gmx::c_numShiftVectors,
                                                          util::array<real, 3>{ 0.0, 0.0, 0.0 });

    thrust::device_vector<util::array<real, 4>> d_xyzq        = xyzq;
    thrust::device_vector<util::array<real, 3>> d_forces      = forces;
    thrust::device_vector<util::array<real, 3>> d_shiftForces = shiftForces;
    thrust::device_vector<real>                 d_potentials(GpuListedEnergySize{}, 0);

    launchBondedKernel(idataGpu.numInteractionsTot(),
                       idataGpu.deviceNumInteractions(),
                       idataGpu.deviceBlockTypes(),
                       idataGpu.deviceIndices(),
                       idataGpu.deviceParametersA(),
                       idataGpu.deviceParametersB(),
                       idataGpu.deviceInteractionOffsets(),
                       idataGpu.deviceIndexOffsets(),
                       thrust::raw_pointer_cast(d_xyzq.data()),
                       thrust::raw_pointer_cast(d_forces.data()),
                       thrust::raw_pointer_cast(d_shiftForces.data()),
                       //(std::nullptr_t*)(nullptr),
                       thrust::raw_pointer_cast(d_potentials.data()),
                       pbc);

    // copy forces back
    forces = d_forces;
    // copy shift forces back
    shiftForces = d_shiftForces;
    // copy potential pack
    thrust::host_vector<real> potentials = d_potentials;

    printf("%f %f %f\n", forces[0][0], forces[0][1], forces[0][2]);
    printf("%f %f %f\n", forces[1][0], forces[1][1], forces[1][2]);

    {
        ListedForceCalculator calculator(interactions, xyzq.size(), 1, box);

        std::vector<Vec3> coordinates;
        std::vector<real> charges;
        for (const auto& c : xyzq)
        {
            coordinates.emplace_back(c[0], c[1], c[2]);
            charges.push_back(c[3]);
        }

        std::vector<Vec3> referenceForces(xyzq.size(), Vec3{ 0, 0, 0 });
        std::vector<Vec3> referenceShiftForces(gmx::c_numShiftVectors, Vec3{ 0, 0, 0 });

        ListedEnergies energies;
        calculator.compute(coordinates, charges, referenceForces, referenceShiftForces, energies);

        double referenceEnergy = std::accumulate(energies.begin(), energies.end() - 3, 0.0);
        double probeEnergy     = std::accumulate(potentials.begin(), potentials.end() - 3, 0.0);
        printf("total potential %f\n", probeEnergy);
        EXPECT_NEAR((probeEnergy - referenceEnergy) / referenceEnergy, 0.0, 1e-5);

        EXPECT_NEAR((potentials[GpuVdwIndex{}] - energies[VdwIndex{}]) / energies[VdwIndex{}], 0.0, 1e-5);
        EXPECT_NEAR((potentials[GpuCoulombIndex{}] - energies[CoulombIndex{}]) / energies[CoulombIndex{}],
                    0.0,
                    1e-5);

        ListedEnergies probe;
        std::fill(probe.begin(), probe.end(), 0);
        auto transferEnergy = [&potentials, &probe](auto interaction) {
            using InteractionType   = std::decay_t<decltype(interaction)>;
            constexpr int gpuIndex  = FindIndex<InteractionType, GpuListedTypes>{};
            constexpr int fullIndex = FindIndex<InteractionType, AllListedTypes>{};
            probe[fullIndex]        = potentials[gpuIndex];
        };

        using TypeLoop = nblib::Reduce<std::tuple, nblib::GpuListedTypes>;
        nblib::for_each_tuple(transferEnergy, TypeLoop{});

        probe[VdwIndex{}]     = potentials[GpuVdwIndex{}];
        probe[CoulombIndex{}] = potentials[GpuCoulombIndex{}];

        for (int i = 0; i < energies.size(); ++i)
        {
            if (energies[i] != 0.0)
            {
                EXPECT_NEAR((probe[i] - energies[i]) / energies[i], 0.0, 1e-5);
            }
            else
            {
                EXPECT_EQ(probe[i], energies[i]);
            }
        }

        compareVectors(shiftForces.begin(), shiftForces.end(), referenceShiftForces.begin(), 1e-3);
        compareVectors(forces.begin(), forces.end(), referenceForces.begin(), 5e-3);
    }
}
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    auto ret = RUN_ALL_TESTS();
    return ret;
}
