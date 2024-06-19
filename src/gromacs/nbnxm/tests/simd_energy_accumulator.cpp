/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for the EnergyAccumulator class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/simd_energy_accumulator.h"

#include "config.h"

#include <cstdlib>

#include <algorithm>
#include <array>
#include <string>
#include <type_traits>

#include <gtest/gtest.h>

#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

#if GMX_HAVE_NBNXM_SIMD_4XM || GMX_HAVE_NBNXM_SIMD_2XMM

namespace
{

template<bool useEnergyGroups>
std::enable_if_t<!useEnergyGroups, EnergyAccumulator<false, true>>
initEnergyAccumulator(int gmx_unused numEnergyGroups, int gmx_unused iClusterSize, int gmx_unused jClusterSize)
{
    return {};
}

template<bool useEnergyGroups>
std::enable_if_t<useEnergyGroups, EnergyAccumulator<true, true>>
initEnergyAccumulator(int gmx_unused numEnergyGroups, int gmx_unused iClusterSize, int gmx_unused jClusterSize)
{
    return EnergyAccumulator<true, true>(numEnergyGroups, iClusterSize, jClusterSize);
}

//! The actual test body checking testing the EnergyAccumulator class
template<KernelLayout kernelLayout, bool useEnergyGroups>
void testEnergyAccumulator()
{
    constexpr int c_iClusterSize            = 4;
    constexpr int c_numIClustersPerRegister = (kernelLayout == KernelLayout::r4xM ? 1 : 2);
    constexpr int c_jClusterSize            = GMX_SIMD_REAL_WIDTH / c_numIClustersPerRegister;
    constexpr int c_numIRegisters           = c_iClusterSize / c_numIClustersPerRegister;

    constexpr int                     c_numAtoms        = 16;
    const std::array<int, c_numAtoms> energyGroups      = { 4, 2, 3, 3, 0, 1, 4, 1,
                                                       1, 0, 3, 2, 3, 2, 4, 1 };
    constexpr int                     c_numEnergyGroups = (useEnergyGroups ? 5 : 1);

    EnergyAccumulator<useEnergyGroups, true> energyAccumulator =
            initEnergyAccumulator<useEnergyGroups>(c_numEnergyGroups, c_iClusterSize, c_jClusterSize);

    const int numIClusters = std::max(8 / c_iClusterSize, 1);
    const int numJClusters = c_numAtoms / c_jClusterSize;

    std::array<std::array<real, c_numEnergyGroups>, c_numEnergyGroups> refCoulombEnergies = { { { 0 } } };
    std::array<std::array<real, c_numEnergyGroups>, c_numEnergyGroups> refVdwEnergies = { { { 0 } } };

    EnergyGroupsPerCluster energyGroupsPerCluster(c_numEnergyGroups, c_iClusterSize);
    if (c_numEnergyGroups > 1)
    {
        energyGroupsPerCluster.setEnergyGroups(energyGroups);
    }

    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (useEnergyGroups)
    {
        energyAccumulator.clearEnergiesAndSetEnergyGroupsForJClusters(energyGroupsPerCluster);
    }

    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    for (int iCluster = 0; iCluster < numIClusters; iCluster++)
    {
        energyAccumulator.template initICluster<c_iClusterSize>(iCluster);

        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
        for (int jCluster = 0; jCluster < numJClusters; jCluster++)
        {
            AlignedArray<real, GMX_SIMD_REAL_WIDTH> coulombEnergies;
            AlignedArray<real, GMX_SIMD_REAL_WIDTH> vdwEnergies;

            std::array<SimdReal, c_numIRegisters> cE;
            std::array<SimdReal, c_numIRegisters> vE;
            for (int iRegister = 0; iRegister < c_numIRegisters; iRegister++)
            {
                for (int half = 0; half < c_numIClustersPerRegister; half++)
                {
                    int iGroup = 0;
                    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
                    if constexpr (useEnergyGroups)
                    {
                        iGroup = energyGroups[iCluster * c_iClusterSize + iRegister * c_numIClustersPerRegister + half];
                    }

                    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
                    for (int j = 0; j < c_jClusterSize; j++)
                    {
                        const real coulombEnergy = std::rand() / real(RAND_MAX);
                        const real vdwEnergy     = std::rand() / real(RAND_MAX);

                        int jGroup = 0;
                        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
                        if constexpr (useEnergyGroups)
                        {
                            jGroup = energyGroups[jCluster * c_jClusterSize + j];
                        }
                        // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
                        refCoulombEnergies[iGroup][jGroup] += coulombEnergy;
                        refVdwEnergies[iGroup][jGroup] += vdwEnergy;

                        coulombEnergies[half * c_jClusterSize + j] = coulombEnergy;
                        vdwEnergies[half * c_jClusterSize + j]     = vdwEnergy;
                    }
                }

                cE[iRegister] = load<SimdReal>(coulombEnergies.data());
                vE[iRegister] = load<SimdReal>(vdwEnergies.data());
            }
            energyAccumulator.template addEnergies<c_numIRegisters, c_numIRegisters, kernelLayout, c_iClusterSize>(
                    jCluster, cE, vE);
        }

        energyAccumulator.reduceIEnergies(true);
    }

    std::array<real, c_numEnergyGroups * c_numEnergyGroups> reducedCoulombEnergies;
    std::array<real, c_numEnergyGroups * c_numEnergyGroups> reducedVdwEnergies;

    energyAccumulator.getEnergies(reducedCoulombEnergies, reducedVdwEnergies);

    for (int iGroup = 0; iGroup < c_numEnergyGroups; iGroup++)
    {
        for (int jGroup = 0; jGroup < c_numEnergyGroups; jGroup++)
        {
            const int index = iGroup * c_numEnergyGroups + jGroup;
            EXPECT_FLOAT_EQ(refCoulombEnergies[iGroup][jGroup], reducedCoulombEnergies[index]);
            EXPECT_FLOAT_EQ(refVdwEnergies[iGroup][jGroup], reducedVdwEnergies[index]);
        }
    }
}

} // namespace

#    if GMX_HAVE_NBNXM_SIMD_4XM

TEST(SimdEnergyAccumulatorTest, SingleEnergyGroupSimd4xM)
{
    testEnergyAccumulator<KernelLayout::r4xM, false>();
}

TEST(SimdEnergyAccumulatorTest, EnergyGroupsSimd4xM)
{
    testEnergyAccumulator<KernelLayout::r4xM, true>();
}

#    endif

#    if GMX_HAVE_NBNXM_SIMD_2XMM

TEST(SimdEnergyAccumulatorTest, SingleEnergyGroupSimd2xMM)
{
    testEnergyAccumulator<KernelLayout::r2xMM, false>();
}

TEST(SimdEnergyAccumulatorTest, EnergyGroupsSimd2xMM)
{
    testEnergyAccumulator<KernelLayout::r2xMM, true>();
}

#    endif // GMX_SIMD

#endif

} // namespace test

} // namespace gmx
