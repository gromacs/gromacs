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
 * Tests for nbnxm setup utilities
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <cmath>

#include "gromacs/hardware/device_management.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"

#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "nblib/box.h"
#include "nblib/nbnxmsetuphelpers.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NbnxmSetupTest, findNumEnergyGroups)
{
    std::vector<int64_t> v(10);
    int                  arbitraryGid = 7;

    // this sets some bit outside the range of bits used for the group ID
    // having all bits zero except those used for the group ID can otherwise hide bugs
    v[5] |= gmx::sc_atomInfo_HasCharge;
    v[5] = (v[5] & ~gmx::sc_atomInfo_EnergyGroupIdMask) | arbitraryGid;

    int nEnergyGroups = arbitraryGid + 1;
    EXPECT_EQ(nEnergyGroups, findNumEnergyGroups(v));
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnumAuto)
{
    auto kernel = SimdKernels::SimdAuto;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::NotSet);
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnumNo)
{
    auto kernel = SimdKernels::SimdNo;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::Cpu4x4_PlainC);
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnum2XM)
{
    auto kernel = SimdKernels::Simd2XMM;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::Cpu4xN_Simd_2xNN);
}

TEST(NbnxmSetupTest, canTranslateBenchmarkEnum4XM)
{
    auto kernel = SimdKernels::Simd4XM;
    EXPECT_EQ(translateBenchmarkEnum(kernel), Nbnxm::KernelType::Cpu4xN_Simd_4xN);
}

TEST(NbnxmSetupTest, CheckKernelSetupThrowsAuto)
{
    EXPECT_ANY_THROW(checkKernelSetupSimd(SimdKernels::SimdAuto));
}

TEST(NbnxmSetupTest, CheckKernelSetupThrowsCount)
{
    EXPECT_ANY_THROW(checkKernelSetupSimd(SimdKernels::Count));
}

TEST(NbnxmSetupTest, canCreateKernelSetupPlain)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd = SimdKernels::SimdNo;
    Nbnxm::KernelSetup kernelSetup =
            createKernelSetupCPU(nbKernelOptions.nbnxmSimd, nbKernelOptions.useTabulatedEwaldCorr);
    EXPECT_EQ(kernelSetup.kernelType, Nbnxm::KernelType::Cpu4x4_PlainC);
    EXPECT_EQ(kernelSetup.ewaldExclusionType, Nbnxm::EwaldExclusionType::Table);
}

TEST(NbnxmSetupTest, canCreateParticleInfoAllVdv)
{
    size_t  numParticles = 2;
    int64_t mask         = 0;
    mask |= gmx::sc_atomInfo_HasVdw;
    mask |= gmx::sc_atomInfo_HasCharge;
    std::vector<int64_t> refParticles  = { mask, mask };
    std::vector<int64_t> testParticles = createParticleInfoAllVdw(numParticles);
    EXPECT_EQ(refParticles, testParticles);
}

TEST(NbnxmSetupTest, ewaldCoeffWorks)
{
    real                              ewald     = ewaldCoeff(1e-5, 1.0);
    gmx::test::FloatingPointTolerance tolerance = gmx::test::absoluteTolerance(1e-5);
    EXPECT_REAL_EQ_TOL(ewald, 3.12341, tolerance);
}

TEST(NbnxmSetupTest, updateForcerecWorks)
{
    t_forcerec forcerec;
    Box        box(3);
    EXPECT_NO_THROW(updateForcerec(&forcerec, box.legacyMatrix()));
}

// The following tests check if the user is allowed to specify configurations not permitted due
// to conflicting compile time setup flags

TEST(NbnxmSetupTest, canCheckKernelSetup)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd = SimdKernels::SimdNo;
#ifdef GMX_NBNXN_SIMD_4XN
    nbKernelOptions.nbnxmSimd = SimdKernels::Simd4XM;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
    nbKernelOptions.nbnxmSimd = SimdKernels::Simd2XMM;
#endif
    EXPECT_NO_THROW(checkKernelSetupSimd(nbKernelOptions.nbnxmSimd));
}

// check if the user is allowed to ask for SimdKernels::Simd2XMM when NBLIB is not compiled with it
#ifndef GMX_NBNXN_SIMD_2XNN
TEST(NbnxmSetupTest, cannotCreateKernelSetupCPU2XM)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd             = SimdKernels::Simd2XMM;
    nbKernelOptions.useTabulatedEwaldCorr = true;
    EXPECT_ANY_THROW(createKernelSetupCPU(nbKernelOptions.nbnxmSimd, nbKernelOptions.useTabulatedEwaldCorr));
}
#endif

// check if the user is allowed to ask for SimdKernels::Simd4XM when NBLIB is not compiled with it
#ifndef GMX_NBNXN_SIMD_4XN
TEST(NbnxmSetupTest, cannotCreateKernelSetupCPU4XM)
{
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd             = SimdKernels::Simd4XM;
    nbKernelOptions.useTabulatedEwaldCorr = false;
    EXPECT_ANY_THROW(createKernelSetupCPU(nbKernelOptions.nbnxmSimd, nbKernelOptions.useTabulatedEwaldCorr));
}
#endif

TEST(NbnxmSetupTest, CanCreateNbnxmCPU)
{
    size_t          numParticles = 1;
    NBKernelOptions nbKernelOptions;
    nbKernelOptions.nbnxmSimd             = SimdKernels::SimdNo;
    int               numEnergyGroups     = 1;
    std::vector<real> nonbondedParameters = { 1, 1 };
    EXPECT_NO_THROW(createNbnxmCPU(numParticles, nbKernelOptions, numEnergyGroups, nonbondedParameters));
}

#if GMX_GPU_CUDA
TEST(NbnxmSetupTest, canCreateKernelSetupGPU)
{
    NBKernelOptions    nbKernelOptions;
    Nbnxm::KernelSetup kernelSetup = createKernelSetupGPU(nbKernelOptions.useTabulatedEwaldCorr);
    EXPECT_EQ(kernelSetup.kernelType, Nbnxm::KernelType::Gpu8x8x8);
    EXPECT_EQ(kernelSetup.ewaldExclusionType, Nbnxm::EwaldExclusionType::Analytical);
}

TEST(NbnxmSetupTest, CanCreateDeviceStreamManager)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo     = testDevice->deviceInfo();
        gmx::SimulationWorkload  simulationWork = createSimulationWorkloadGpu();
        EXPECT_NO_THROW(createDeviceStreamManager(deviceInfo, simulationWork));
    }
}

TEST(NbnxmSetupTest, CanCreateNbnxmGPU)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo   = testDevice->deviceInfo();
        size_t                   numParticles = 1;
        NBKernelOptions          nbKernelOptions;
        std::vector<real>        nonbondedParameters = { 1, 1 };
        gmx::SimulationWorkload  simulationWork      = createSimulationWorkloadGpu();
        interaction_const_t      interactionConst    = createInteractionConst(nbKernelOptions);
        // set DeviceInformation and create the DeviceStreamManager
        auto deviceStreamManager = createDeviceStreamManager(deviceInfo, simulationWork);
        EXPECT_NO_THROW(createNbnxmGPU(
                numParticles, nbKernelOptions, nonbondedParameters, interactionConst, *deviceStreamManager));
    }
}

#endif

} // namespace
} // namespace test
} // namespace nblib
