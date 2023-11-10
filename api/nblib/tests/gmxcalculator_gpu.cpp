/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <filesystem>

#include <gtest/gtest.h>

#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/test_hardware_environment.h"

#include "nblib/gmxcalculatorcpu.h"
#include "nblib/gmxcalculatorgpu.h"
#include "nblib/kerneloptions.h"
#include "nblib/simulationstate.h"
#include "nblib/tests/testhelpers.h"
#include "nblib/tests/testsystems.h"

#include "buildinfo.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, canCreateGPUfc)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    if (testDeviceList.empty())
    {
        GTEST_SKIP() << "No compatible GPUs to test on.";
    }
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo = testDevice->deviceInfo();

        SpcMethanolSimulationStateBuilder spcMethanolSimulationStateBuilder;
        SimulationState simState = spcMethanolSimulationStateBuilder.setupSimulationState();
        NBKernelOptions options  = NBKernelOptions();
        options.useGpu           = true;
        options.coulombType      = CoulombType::Cutoff;
        EXPECT_NO_THROW(std::unique_ptr<GmxNBForceCalculatorGpu> gmxForceCalculator =
                                setupGmxForceCalculatorGpu(simState.topology(), options, deviceInfo));
    }
}

TEST(NBlibTest, SpcMethanolForcesAreCorrectOnGpu)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    if (testDeviceList.empty())
    {
        GTEST_SKIP() << "No compatible GPUs to test on.";
    }
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo = testDevice->deviceInfo();

        SpcMethanolSimulationStateBuilder spcMethanolSimulationStateBuilder;
        SimulationState simState = spcMethanolSimulationStateBuilder.setupSimulationState();
        NBKernelOptions options  = NBKernelOptions();
        options.coulombType      = CoulombType::Cutoff;
        auto gmxForceCalculator = setupGmxForceCalculatorGpu(simState.topology(), options, deviceInfo);
        gmxForceCalculator->updatePairlist(simState.coordinates(), simState.box());

        gmx::ArrayRef<Vec3> forces(simState.forces());
        ASSERT_NO_THROW(gmxForceCalculator->compute(simState.coordinates(), simState.box(), forces));

        RefDataChecker forcesOutputTest;
        forcesOutputTest.testArrays<Vec3>(forces, "SPC-methanol forces on GPU");
    }
}

/*! \brief reorder + undoReorder test
 *
 * We do    defaultOrder -> nbnxmOrderXQ -> nbnxmOrderX -> defaultOrderRecovered
 *                       |                              |
 *                   (reorder)                     (undoReorder)
 *
 * and then check that   defaultOrder == defaultOrderRecovered
 */
TEST(NBlibTest, ReorderIsInvertible)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    if (testDeviceList.empty())
    {
        GTEST_SKIP() << "No compatible GPUs to test on.";
    }
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo = testDevice->deviceInfo();

        const auto filepath =
                std::filesystem::path{ CMAKE_SOURCE_DIR "api/nblib/samples/", "argon5832.tpr" };
        nblib::TprReader tpr(filepath);

        int numParticles = tpr.coordinates_.size();

        NBKernelOptions options = NBKernelOptions();
        options.coulombType     = CoulombType::Cutoff;

        nblib::GmxNBForceCalculatorGpu forceCalculator =
                nblib::GmxNBForceCalculatorGpu(tpr.particleTypeIdOfAllParticles_,
                                               tpr.nonbondedParameters_,
                                               tpr.charges_,
                                               tpr.particleInteractionFlags_,
                                               tpr.exclusionListRanges_,
                                               tpr.exclusionListElements_,
                                               options,
                                               deviceInfo);

        forceCalculator.updatePairlist(tpr.coordinates_, tpr.getBox());

        std::size_t nbnxmBufferSize = forceCalculator.nbnxmBufferSize();

        // default order is just a sequence 0, 1, ..., 3*numParticles
        std::vector<Vec3> defaultOrder(numParticles);
        for (int i = 0; i < numParticles; ++i)
        {
            defaultOrder[i] = { real(3 * i + 0), real(3 * i + 1), real(3 * i + 2) };
        }

        std::vector<real> nbnxmOrderXQ(4 * nbnxmBufferSize);
        forceCalculator.reorder(defaultOrder, nbnxmOrderXQ);

        // throw away Q
        std::vector<gmx::RVec> nbnxmOrderX(nbnxmBufferSize);
        for (std::size_t i = 0; i < nbnxmBufferSize; ++i)
        {
            nbnxmOrderX[i] = { nbnxmOrderXQ[4 * i], nbnxmOrderXQ[4 * i + 1], nbnxmOrderXQ[4 * i + 2] };
        }

        std::vector<gmx::RVec> defaultOrderRecovered(numParticles, { -1, -1, -1 });
        forceCalculator.undoReorder(nbnxmOrderX, defaultOrderRecovered);

        // original defaultOrder should be identical to defaultOrderRecovered
        for (int i = 0; i < numParticles; ++i)
        {
            EXPECT_EQ(defaultOrder[i][0], defaultOrderRecovered[i][0]);
            EXPECT_EQ(defaultOrder[i][1], defaultOrderRecovered[i][1]);
            EXPECT_EQ(defaultOrder[i][2], defaultOrderRecovered[i][2]);
        }
    }
}

/*! \brief test the DeviceBuffer compute interface
 *
 * We do:   X -> (reorder) -> XqNbnxm -> (copyToDevice) -> deviceXq -> (compute) -> deviceForces
 *          -> (copyFromDevice) -> forcesNbnxm -> (undoReorder) -> forcesDefaultOrder
 *
 * Then we test that   forcesDefaultOrder == forces from CPU buffer interface which are tested
 *                                                           against reference values
 */
TEST(NBlibTest, SpcMethanolForcesDeviceInterface)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    if (testDeviceList.empty())
    {
        GTEST_SKIP() << "No compatible GPUs to test on.";
    }
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo = testDevice->deviceInfo();

        SpcMethanolSimulationStateBuilder spcMethanolSimulationStateBuilder;
        SimulationState simState = spcMethanolSimulationStateBuilder.setupSimulationState();
        NBKernelOptions options  = NBKernelOptions();
        options.coulombType      = CoulombType::Cutoff;
        auto forceCalculator = setupGmxForceCalculatorGpu(simState.topology(), options, deviceInfo);
        forceCalculator->updatePairlist(simState.coordinates(), simState.box());

        gmx::ArrayRef<Vec3> forcesReference(simState.forces());

        // we know these forces to be correct from a different test
        forceCalculator->compute(simState.coordinates(), simState.box(), forcesReference);

        std::size_t          deviceBufferSize = forceCalculator->nbnxmBufferSize();
        int                  numParticles     = simState.topology().numParticles();
        const DeviceContext& context          = forceCalculator->deviceContext();

        // reorder coordinates into nbnxm ordering on the CPU
        std::vector<real> coordinatesNbnxm(4 * deviceBufferSize);
        forceCalculator->reorder(simState.coordinates(), coordinatesNbnxm);

        // device coordinates
        DeviceBuffer<Float4> deviceXq;
        allocateDeviceBuffer(&deviceXq, deviceBufferSize, context);
        copyToDeviceBuffer(&deviceXq,
                           reinterpret_cast<const Float4*>(coordinatesNbnxm.data()),
                           0,
                           deviceBufferSize,
                           forceCalculator->deviceStream(),
                           GpuApiCallBehavior::Sync,
                           nullptr);

        DeviceBuffer<Float3> deviceForces;
        allocateDeviceBuffer(&deviceForces, deviceBufferSize, context);
        clearDeviceBufferAsync(&deviceForces, 0, deviceBufferSize, forceCalculator->deviceStream());

        // launch compute directly with device buffers
        forceCalculator->compute(deviceXq, simState.box(), deviceForces);

        // download forces from the GPU
        std::vector<gmx::RVec> forcesNbnxm(deviceBufferSize, { 0, 0, 0 });
        copyFromDeviceBuffer(reinterpret_cast<Float3*>(forcesNbnxm.data()),
                             &deviceForces,
                             0,
                             deviceBufferSize,
                             forceCalculator->deviceStream(),
                             GpuApiCallBehavior::Sync,
                             nullptr);

        // reorder downloaded forces from nbnxm into default ordering
        std::vector<gmx::RVec> forcesDefaultOrder(numParticles, gmx::RVec{ 0, 0, 0 });
        forceCalculator->undoReorder(forcesNbnxm, forcesDefaultOrder);

        // forcesDefaultOrder should be equal to the reference
        for (int i = 0; i < numParticles; ++i)
        {
            EXPECT_EQ(forcesDefaultOrder[i][0], forcesReference[i][0]);
            EXPECT_EQ(forcesDefaultOrder[i][1], forcesReference[i][1]);
            EXPECT_EQ(forcesDefaultOrder[i][2], forcesReference[i][2]);
        }

        freeDeviceBuffer(&deviceXq);
        freeDeviceBuffer(&deviceForces);
    }
}

/*! \brief test the GPU calculator against the CPU calculator with a bigger system with 5832 argon atoms
 */
TEST(NBlibTest, Argon5832ForcesCpuVsGpu)
{
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    if (testDeviceList.empty())
    {
        GTEST_SKIP() << "No compatible GPUs to test on.";
    }
    for (const auto& testDevice : testDeviceList)
    {
        testDevice->activate();
        const DeviceInformation& deviceInfo = testDevice->deviceInfo();

        const auto filepath =
                std::filesystem::path{ CMAKE_SOURCE_DIR "api/nblib/samples/", "argon5832.tpr" };
        nblib::TprReader tpr(filepath);

        std::size_t numParticles = tpr.coordinates_.size();

        std::vector<Vec3> forcesCpu(numParticles, gmx::RVec{ 0, 0, 0 });
        std::vector<Vec3> forcesGpu(numParticles, gmx::RVec{ 0, 0, 0 });
        std::vector<Vec3> forcesGpuDirect(numParticles, gmx::RVec{ 0, 0, 0 });

        // GPU calculation
        {
            NBKernelOptions options = NBKernelOptions();
            options.coulombType     = CoulombType::Cutoff;

            nblib::GmxNBForceCalculatorGpu forceCalculatorGpu =
                    nblib::GmxNBForceCalculatorGpu(tpr.particleTypeIdOfAllParticles_,
                                                   tpr.nonbondedParameters_,
                                                   tpr.charges_,
                                                   tpr.particleInteractionFlags_,
                                                   tpr.exclusionListRanges_,
                                                   tpr.exclusionListElements_,
                                                   options,
                                                   deviceInfo);

            forceCalculatorGpu.updatePairlist(tpr.coordinates_, tpr.getBox());

            forceCalculatorGpu.compute(tpr.coordinates_, tpr.getBox(), forcesGpu);
        }

        // CPU calculation
        {
            NBKernelOptions options = NBKernelOptions();
            options.nbnxmSimd       = SimdKernels::SimdNo;
            options.coulombType     = CoulombType::Cutoff;

            nblib::GmxNBForceCalculatorCpu forceCalculatorCpu =
                    nblib::GmxNBForceCalculatorCpu(tpr.particleTypeIdOfAllParticles_,
                                                   tpr.nonbondedParameters_,
                                                   tpr.charges_,
                                                   tpr.particleInteractionFlags_,
                                                   tpr.exclusionListRanges_,
                                                   tpr.exclusionListElements_,
                                                   options);

            forceCalculatorCpu.updatePairlist(tpr.coordinates_, tpr.getBox());
            forceCalculatorCpu.compute(tpr.coordinates_, tpr.getBox(), forcesCpu);
        }

        // GPU device buffer interface
        {
            NBKernelOptions options = NBKernelOptions();
            options.coulombType     = CoulombType::Cutoff;

            nblib::GmxNBForceCalculatorGpu forceCalculator =
                    nblib::GmxNBForceCalculatorGpu(tpr.particleTypeIdOfAllParticles_,
                                                   tpr.nonbondedParameters_,
                                                   tpr.charges_,
                                                   tpr.particleInteractionFlags_,
                                                   tpr.exclusionListRanges_,
                                                   tpr.exclusionListElements_,
                                                   options,
                                                   deviceInfo);

            const DeviceContext& context = forceCalculator.deviceContext();

            forceCalculator.updatePairlist(tpr.coordinates_, tpr.getBox());

            std::size_t deviceBufferSize = forceCalculator.nbnxmBufferSize();

            std::vector<real>      xNbnxm4wide(4 * deviceBufferSize);
            std::vector<gmx::RVec> xNbnxm(deviceBufferSize);

            DeviceBuffer<Float4> deviceXq;
            allocateDeviceBuffer(&deviceXq, deviceBufferSize, context);

            DeviceBuffer<Float3> deviceForces;
            allocateDeviceBuffer(&deviceForces, deviceBufferSize, context);

            forceCalculator.reorder(tpr.coordinates_, xNbnxm4wide);
            copyToDeviceBuffer(&deviceXq,
                               reinterpret_cast<const Float4*>(xNbnxm4wide.data()),
                               0,
                               deviceBufferSize,
                               forceCalculator.deviceStream(),
                               GpuApiCallBehavior::Sync,
                               nullptr);

            clearDeviceBufferAsync(&deviceForces, 0, deviceBufferSize, forceCalculator.deviceStream());
            forceCalculator.compute(deviceXq, tpr.getBox(), deviceForces);

            std::vector<gmx::RVec> fNbnxm(deviceBufferSize, gmx::RVec{ 0, 0, 0 });
            copyFromDeviceBuffer(reinterpret_cast<Float3*>(fNbnxm.data()),
                                 &deviceForces,
                                 0,
                                 deviceBufferSize,
                                 forceCalculator.deviceStream(),
                                 GpuApiCallBehavior::Sync,
                                 nullptr);
            forceCalculator.undoReorder(fNbnxm, forcesGpuDirect);

            freeDeviceBuffer(&deviceXq);
            freeDeviceBuffer(&deviceForces);
        }

        for (std::size_t i = 0; i < numParticles; ++i)
        {
            EXPECT_REAL_EQ_TOL(forcesGpuDirect[i][0], forcesGpu[i][0], gmx::test::absoluteTolerance(1e-4));
            EXPECT_REAL_EQ_TOL(forcesGpuDirect[i][1], forcesGpu[i][1], gmx::test::absoluteTolerance(1e-4));
            EXPECT_REAL_EQ_TOL(forcesGpuDirect[i][2], forcesGpu[i][2], gmx::test::absoluteTolerance(1e-4));
        }
    }
}

} // namespace
} // namespace test
} // namespace nblib
