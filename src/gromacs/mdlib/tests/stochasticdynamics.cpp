/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Tests for the Langevin integrator
 *
 *  The test creates a system of independent particles exerting constant
 *  external forces and makes several numerical integration timesteps.
 *  The results are compared with pre-computed reference values.
 *
 * \todo Add PBC handling test.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include <cmath>

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "stochasticdynamicstestdata.h"

#if GMX_GPU && !GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/mdlib/stochasticdynamics_gpu.h"
#endif

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/matrix.h"

namespace gmx
{
namespace test
{
namespace
{

static const RVec sc_initialVelocity{ 1.0, -2.0, 3.0 };
static const RVec sc_force{ -3.0, 2.0, -1.0 };

void integrateStochasticDynamicsCpu(StochasticDynamicsTestData* testData, int64_t numSteps)
{
    testData->state_.x.resizeWithPadding(testData->numAtoms_);
    testData->state_.v.resizeWithPadding(testData->numAtoms_);
    for (int i = 0; i < testData->numAtoms_; i++)
    {
        testData->state_.x[i] = testData->x_[i];
        testData->state_.v[i] = testData->v_[i];
    }

    gmx_omp_nthreads_set(ModuleMultiThread::Update, 1);

    Matrix3x3 parrinelloRahmanM;

    for (int64_t step = 0; step < numSteps; step++)
    {
        testData->update_->update_coords(testData->inputRecord_,
                                         step,
                                         testData->numAtoms_,
                                         testData->ptype_,
                                         testData->inverseMasses_,
                                         testData->inverseMassesPerDim_,
                                         &testData->state_,
                                         testData->f_,
                                         &testData->forceCalculationData_,
                                         &testData->kineticEnergyData_,
                                         parrinelloRahmanM,
                                         etrtNONE,
                                         nullptr,
                                         false);
        testData->update_->finish_update(
                testData->inputRecord_, testData->numAtoms_, &testData->state_, nullptr, false);
    }
    const auto xp = makeArrayRef(*testData->update_->xp()).subArray(0, testData->numAtoms_);
    for (int i = 0; i < testData->numAtoms_; i++)
    {
        for (int d = 0; d < DIM; d++)
        {
            testData->x_[i][d]      = testData->state_.x[i][d];
            testData->v_[i][d]      = testData->state_.v[i][d];
            testData->xPrime_[i][d] = xp[i][d];
        }
    }
}

#if GMX_GPU && !GMX_GPU_OPENCL

void integrateStochasticDynamicsGpu(const DeviceContext&        deviceContext,
                                    const DeviceStream&         deviceStream,
                                    StochasticDynamicsTestData* testData,
                                    int64_t                     numSteps,
                                    ArrayRef<const int>         globalAtomIndices = {})
{
    int numAtoms = testData->numAtoms_;

    Float3* h_x  = gmx::asGenericFloat3Pointer(testData->x_);
    Float3* h_xp = gmx::asGenericFloat3Pointer(testData->xPrime_);
    Float3* h_v  = gmx::asGenericFloat3Pointer(testData->v_);
    Float3* h_f  = gmx::asGenericFloat3Pointer(testData->f_);

    DeviceBuffer<Float3> d_x, d_xp, d_v, d_f;

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_f, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_v, h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_f, h_f, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    std::vector<float> refT(testData->numTCoupleGroups_, testData->temperature_);
    std::vector<float> tauT(testData->numTCoupleGroups_, testData->tauT_);

    auto integrator          = std::make_unique<StochasticDynamicsGpu>(deviceContext,
                                                              deviceStream,
                                                              testData->numTCoupleGroups_,
                                                              testData->timestep_,
                                                              refT.data(),
                                                              tauT.data());
    testData->inverseMasses_ = PaddedHostVector<real>(
            testData->inverseMasses_,
            HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported });

    testData->cTC_ = HostVector<unsigned short>(
            testData->cTC_, HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported });

    {
        SCOPED_TRACE("Preliminary call to set() to cover resizing edge cases");
        // Use one atom more, so the call under test shrinks an object that holds data
        const int              largerNumAtoms = testData->numAtoms_ + 1;
        PaddedHostVector<real> inverseMasses(
                HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported });
        HostVector<unsigned short> cTC(
                HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported });
        inverseMasses.resizeWithPadding(largerNumAtoms);
        for (int i = 0; i < largerNumAtoms; i++)
        {
            inverseMasses[i] = 1.0 / (3.0 + i);
        }
        // Like mdatoms, only fill cTC when there is more than one group
        if (testData->numTCoupleGroups_ > 1)
        {
            cTC.resize(largerNumAtoms);
            for (int i = 0; i < largerNumAtoms; i++)
            {
                cTC[i] = i % testData->numTCoupleGroups_;
            }
        }
        integrator->set(largerNumAtoms, inverseMasses, cTC, {});
        deviceStream.synchronize();
    }

    integrator->set(numAtoms, testData->inverseMasses_, testData->cTC_, globalAtomIndices);

    for (int64_t step = 0; step < numSteps; step++)
    {
        integrator->integrate(
                d_x, d_xp, d_v, d_f, testData->timestep_, testData->seed_, step, SDUpdateStage::Combined);
    }

    copyFromDeviceBuffer(h_v, &d_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(h_xp, &d_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
    freeDeviceBuffer(&d_f);
}

/*! \brief Checks that the GPU result for each atom does not depend on the local atom order
 *
 * With domain decomposition the local atom order differs from the global one. This integrates
 * the system with the local atom order reversed, passing the global atom indices, and compares
 * with \p reference, which was integrated with the local order equal to the global order.
 */
void checkGpuResultIsIndependentOfLocalAtomOrder(const DeviceContext&              deviceContext,
                                                 const DeviceStream&               deviceStream,
                                                 const StochasticDynamicsTestData& reference,
                                                 int64_t                           numSteps)
{
    StochasticDynamicsTestData reordered(reference.numAtoms_,
                                         reference.timestep_,
                                         sc_initialVelocity,
                                         sc_force,
                                         reference.numTCoupleGroups_,
                                         reference.temperature_,
                                         reference.tauT_,
                                         reference.seed_);
    const int                  numAtoms = reordered.numAtoms_;
    std::vector<int>           globalAtomIndices(numAtoms);
    std::iota(globalAtomIndices.rbegin(), globalAtomIndices.rend(), 0);
    std::reverse(reordered.x_.begin(), reordered.x_.end());
    std::reverse(reordered.v_.begin(), reordered.v_.end());
    std::reverse(reordered.f_.begin(), reordered.f_.end());
    std::reverse(reordered.inverseMasses_.begin(), reordered.inverseMasses_.end());
    std::reverse(reordered.cTC_.begin(), reordered.cTC_.end());

    integrateStochasticDynamicsGpu(deviceContext, deviceStream, &reordered, numSteps, globalAtomIndices);

    // Each atom goes through the same arithmetic in both runs, so the results must be identical
    for (int localIndex = 0; localIndex < numAtoms; localIndex++)
    {
        const int globalIndex = globalAtomIndices[localIndex];
        SCOPED_TRACE(formatString("Atom with global index %d", globalIndex));
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_EQ(reference.xPrime_[globalIndex][d], reordered.xPrime_[localIndex][d]);
            EXPECT_EQ(reference.v_[globalIndex][d], reordered.v_[localIndex][d]);
        }
    }
}

#endif // GMX_GPU && !GMX_GPU_OPENCL

//! Input configuration for Langevin tests
using StochasticDynamicsInputConfig = std::tuple<int,      // numAtoms
                                                 real,     // timestep
                                                 int64_t,  // numSteps
                                                 int,      // numTCoupleGroups
                                                 real,     // temperature
                                                 real,     // tauT
                                                 int64_t>; // seed

/*! \brief Hardware test helper for Langevin
 *
 * \todo There are no execution modes - should test SIMD vs no SIMD
 * here. Perhaps coupling vs no-coupling is useful to express this way
 * also. */
using StochasticDynamicsTestHelper =
        HardwareAndExecutionTestHelper<StochasticDynamicsInputConfig, std::tuple<>>;

//! Format timestep as integer for test names
std::string formatTimestep(real timestep)
{
    int timestepInt = std::lround(timestep * 1000000);
    return formatString("dt%d", timestepInt);
}

//! Format tau as integer for test names
std::string formatTauT(real tauT)
{
    int tauTInt = std::lround(tauT * 1000);
    return formatString("tauT%d", tauTInt);
}

//! Format temperature as integer for test names
std::string formatTemperature(real temperature)
{
    int temperatureInt = std::lround(temperature);
    return formatString("T%d", temperatureInt);
}

//! Formatters for parameters in the config info
static const auto sc_configInfoFormatters =
        std::make_tuple([](int n) { return formatString("%datoms", n); },
                        formatTimestep,
                        [](int n) { return formatString("%dsteps", n); },
                        [](int n) { return formatString("tcg%d", n); },
                        formatTemperature,
                        formatTauT,
                        [](int n) { return formatString("seed%d", n); });
//! Formatters for parameters in the execution mode (currently empty)
static const auto sc_executionModeFormatters = std::make_tuple();

//! Helper object to name tests using all parameters
static const NameOfTestFromTuple<StochasticDynamicsTestHelper::DynamicParameters> sc_testNamer =
        StochasticDynamicsTestHelper::testNamer(sc_configInfoFormatters, sc_executionModeFormatters);

//! The set of parameters combinations to run the test on
const std::array<StochasticDynamicsInputConfig, 17> sc_langevinConfigs = { {
        { 1, 0.001, 1, 1, 0, 2, 123 },
        { 1, 0.001, 1, 1, 0, 2, 12345 },
        { 1, 0.001, 1, 1, 100, 2, 123 },
        { 1, 0.001, 1, 1, 100, 0, 123 },
        { 1, 0.0005, 1, 1, 100, 2, 123 },
        { 1, 0.0005, 2, 1, 100, 2, 123 },
        { 1, 0.0025, 3, 1, 100, 2, 123 },
        { 1, 0.0025, 3, 1, 100, 2, 12345 },
        { 1, 0.0025, 3, 1, 500, 2, 123 },
        { 1, 0.0025, 20, 1, 0, 2, 123 },
        { 1, 0.0025, 20, 1, 100, 2, 123 },
        { 1, 0.0025, 20, 1, 100, 2, 12345 },
        { 10, 0.0025, 1, 1, 100, 2, 123 },
        { 10, 0.0025, 1, 2, 100, 2, 123 },
        { 10, 0.0025, 1, 2, 100, 2, 12345 },
        { 10, 0.0025, 1, 2, 100, 0.5, 123 },
        { 10, 0.0025, 1, 5, 100, 2, 12345 },
} };

/*! \brief Test fixture for Langevin integrator.
 */
class StochasticDynamicsTest : public HardwareTestFixture<StochasticDynamicsTestHelper>
{
protected:
    StochasticDynamicsTest() : HardwareTestFixture(sc_configInfoFormatters) {}

public:
    void testAgainstReferenceData(const StochasticDynamicsTestData& testData)
    {
        TestReferenceChecker finalPositionsRef(
                checker().checkSequenceCompound("FinalPositions", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            const gmx::RVec&     xPrime = testData.xPrime_[i];
            TestReferenceChecker xPrimeRef(finalPositionsRef.checkCompound("Atom", nullptr));
            xPrimeRef.checkReal(xPrime[XX], "XX");
            xPrimeRef.checkReal(xPrime[YY], "YY");
            xPrimeRef.checkReal(xPrime[ZZ], "ZZ");
        }

        TestReferenceChecker finalVelocitiesRef(
                checker().checkSequenceCompound("FinalVelocities", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            const gmx::RVec&     v = testData.v_[i];
            TestReferenceChecker vRef(finalVelocitiesRef.checkCompound("Atom", nullptr));
            vRef.checkReal(v[XX], "XX");
            vRef.checkReal(v[YY], "YY");
            vRef.checkReal(v[ZZ], "ZZ");
        }
    }
};

TEST_P(StochasticDynamicsTest, SimpleIntegration)
{
    auto [numAtoms, timestep, numSteps, numTCoupleGroups, temperature, tauT, seed, _] = GetParam();
    {
        std::unique_ptr<StochasticDynamicsTestData> testData = std::make_unique<StochasticDynamicsTestData>(
                numAtoms, timestep, sc_initialVelocity, sc_force, numTCoupleGroups, temperature, tauT, seed);

        if (isGpuTest())
        {
#if GMX_GPU && !GMX_GPU_OPENCL
            activateHardware();
            integrateStochasticDynamicsGpu(*hardwareContext()->deviceContext(),
                                           *hardwareContext()->deviceStream(),
                                           testData.get(),
                                           numSteps);
            if (numAtoms > 1)
            {
                checkGpuResultIsIndependentOfLocalAtomOrder(*hardwareContext()->deviceContext(),
                                                            *hardwareContext()->deviceStream(),
                                                            *testData,
                                                            numSteps);
            }
#else
            GMX_THROW(gmx::InternalError("GPU hardware context with no test code"));
#endif
        }
        else
        {
            integrateStochasticDynamicsCpu(testData.get(), numSteps);
        }

        FloatingPointTolerance tolerance = absoluteTolerance(numSteps * (GMX_DOUBLE ? 5e-10 : 5e-6));

        checker().setDefaultTolerance(tolerance);
        testAgainstReferenceData(*testData);
    }
}

INSTANTIATE_TEST_SUITE_P(AllHardware,
                         StochasticDynamicsTest,
                         ::testing::ConvertGenerator(
                                 ::testing::Combine(::testing::ValuesIn(sc_langevinConfigs),
                                                    ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                            GpuConfigurationCapabilities::UpdateSD))),
                                 flattenTupleWithHardwareContext<StochasticDynamicsInputConfig>()),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
