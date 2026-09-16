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

#include <array>
#include <memory>
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

#include "langevintestdata.h"

#if GMX_GPU && !GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
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

void integrateLangevinCpu(LangevinTestData* testData, int numSteps)
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

    for (int step = 0; step < numSteps; step++)
    {
        testData->update_->update_coords(testData->inputRecord_,
                                         step,
                                         testData->mdAtoms_.homenr,
                                         testData->mdAtoms_.havePartiallyFrozenAtoms,
                                         testData->mdAtoms_.ptype,
                                         testData->mdAtoms_.invmass,
                                         testData->mdAtoms_.invMassPerDim,
                                         &testData->state_,
                                         testData->f_,
                                         &testData->forceCalculationData_,
                                         &testData->kineticEnergyData_,
                                         parrinelloRahmanM,
                                         etrtNONE,
                                         nullptr,
                                         false);
        testData->update_->finish_update(testData->inputRecord_,
                                         testData->mdAtoms_.havePartiallyFrozenAtoms,
                                         testData->mdAtoms_.homenr,
                                         &testData->state_,
                                         nullptr,
                                         false);
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


//! Input configuration for Langevin tests
using LangevinInputConfig = std::tuple<int,  // numAtoms
                                       real, // timestep
                                       int,  // numSteps
                                       int,  // numTCoupleGroups
                                       real, // temperature
                                       real, // tauT
                                       int>; // seed

/*! \brief Hardware test helper for Langevin
 *
 * \todo There are no execution modes - should test SIMD vs no SIMD
 * here. Perhaps coupling vs no-coupling is useful to express this way
 * also. */
using LangevinTestHelper = HardwareAndExecutionTestHelper<LangevinInputConfig, std::tuple<>>;

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
static const NameOfTestFromTuple<LangevinTestHelper::DynamicParameters> sc_testNamer =
        LangevinTestHelper::testNamer(sc_configInfoFormatters, sc_executionModeFormatters);

//! The set of parameters combinations to run the test on
const std::array<LangevinInputConfig, 17> sc_langevinConfigs = { {
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
class LangevinTest : public HardwareTestFixture<LangevinTestHelper>
{
protected:
    LangevinTest() : HardwareTestFixture(sc_configInfoFormatters) {}

public:
    void testAgainstReferenceData(const LangevinTestData& testData)
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

TEST_P(LangevinTest, SimpleIntegration)
{
    auto [numAtoms, timestep, numSteps, numTCoupleGroups, temperature, tauT, seed, _] = GetParam();

    {
        std::unique_ptr<LangevinTestData> testData = std::make_unique<LangevinTestData>(
                numAtoms, timestep, sc_initialVelocity, sc_force, numTCoupleGroups, temperature, tauT, seed);

        if (isGpuTest())
        {
            GMX_THROW(gmx::InternalError("GPU hardware context with no test code"));
        }
        else
        {
            integrateLangevinCpu(testData.get(), numSteps);
        }

        FloatingPointTolerance tolerance = absoluteTolerance(numSteps * (GMX_DOUBLE ? 5e-10 : 5e-6));

        checker().setDefaultTolerance(tolerance);
        testAgainstReferenceData(*testData);
    }
}

INSTANTIATE_TEST_SUITE_P(AllHardware,
                         LangevinTest,
                         ::testing::ConvertGenerator(
                                 ::testing::Combine(::testing::ValuesIn(sc_langevinConfigs),
                                                    ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                            GpuConfigurationCapabilities::UpdateSD))),
                                 flattenTupleWithHardwareContext<LangevinInputConfig>()),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
