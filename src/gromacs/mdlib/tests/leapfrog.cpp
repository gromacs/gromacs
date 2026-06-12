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
 * \brief Tests for the Leap-Frog integrator
 *
 *  The test creates a system of independent particles exerting constant
 *  external forces and makes several numerical integration timesteps.
 *  The results are compared with the analytical solution (for the systems
 *  without the temperature coupling) and with the pre-computed reference
 *  values. The tests use runners that are created for each available
 *  implementation of the tested algorithm.
 *
 * \todo Add tests for integrators with pressure control.
 * \todo Add PBC handling test.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include <cmath>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "leapfrogtestdata.h"
#include "leapfrogtestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

//! Parameter tuple type for LeapFrog tests (with vfName for naming infrastructure)
using LeapFrogParametersTuple = std::tuple<int,         // numAtoms
                                           real,        // timestep
                                           int,         // numSteps
                                           std::string, // vfName (velocity/force pair name)
                                           int,         // numTCoupleGroups
                                           int>;        // nstpcouple

//! Format timestep as integer for test names
std::string formatTimestep(real timestep)
{
    int timestepInt = static_cast<int>(timestep * 1000000 + 0.5);
    return formatString("dt%d", timestepInt);
}

//! Test naming functor
const NameOfTestFromTuple<LeapFrogParametersTuple> sc_testNamer{ std::make_tuple(
        [](int n) { return formatString("%datoms", n); },
        formatTimestep,
        [](int n) { return formatString("%dsteps", n); },
        useString, // Velocity/force pair name
        [](int n) { return formatString("tcg%d", n); },
        [](int n) { return formatString("nstpc%d", n); }) };

//! Reference data filename maker (same as test namer for now)
const RefDataFilenameMaker<LeapFrogParametersTuple> sc_refDataFilenameMaker{ std::make_tuple(
        [](int n) { return formatString("%datoms", n); },
        formatTimestep,
        [](int n) { return formatString("%dsteps", n); },
        useString, // Velocity/force pair name distinguishes different initial conditions
        [](int n) { return formatString("tcg%d", n); },
        [](int n) { return formatString("nstpc%d", n); }) };

//! Named velocity/force combination
struct VelocityForcePair
{
    std::string name;
    RVec        velocity;
    RVec        force;
};

//! Lookup table for velocity/force combinations
static const VelocityForcePair sc_velocityForceTable[] = {
    { "zero", RVec(0.0, 0.0, 0.0), RVec(0.0, 0.0, 0.0) },       // Zero velocity and force
    { "zeroV", RVec(0.0, 0.0, 0.0), RVec(-3.0, 2.0, -1.0) },    // Zero velocity
    { "zeroF", RVec(1.0, -2.0, 3.0), RVec(0.0, 0.0, 0.0) },     // Zero force
    { "standard", RVec(1.0, -2.0, 3.0), RVec(-3.0, 2.0, -1.0) } // Standard test values
};

//! Lookup velocity/force pair by name
const VelocityForcePair& getVelocityForcePair(const std::string& name)
{
    for (const auto& pair : sc_velocityForceTable)
    {
        if (pair.name == name)
        {
            return pair;
        }
    }
    GMX_THROW(gmx::InternalError(gmx::formatString("Unknown velocity/force pair name: %s", name.c_str())));
}

/*! \brief The set of parameters combinations to run the test on
 *
 * Each entry is a tuple of (numAtoms, timestep, numSteps, vfName, numTCoupleGroups, nstpcouple).
 */
const LeapFrogParametersTuple parametersSets[] = {
    { 1, 0.001, 1, "zero", 0, 0 },        // Zero velocity and force
    { 1, 0.001, 1, "zeroV", 0, 0 },       // Zero velocity
    { 1, 0.001, 1, "zeroF", 0, 0 },       // Zero force
    { 1, 0.001, 1, "standard", 0, 0 },    // 1 particle
    { 10, 0.001, 1, "standard", 0, 0 },   // 10 particles
    { 100, 0.001, 1, "standard", 0, 0 },  // 100 particles
    { 300, 0.001, 1, "standard", 0, 0 },  // 300 particles
    { 1, 0.0005, 1, "standard", 0, 0 },   // 0.0005 ps timestep
    { 1, 0.001, 10, "standard", 0, 0 },   // 10 step
    { 1, 0.001, 100, "standard", 0, 0 },  // 100 steps
    { 100, 0.001, 1, "standard", 1, 0 },  // 1 temperature couple group
    { 100, 0.001, 1, "standard", 2, 0 },  // 2 temperature couple groups
    { 100, 0.001, 1, "standard", 10, 0 }, // 10 temperature couple groups
    { 100, 0.001, 10, "standard", 0, 1 }, // With pressure coupling
    { 100, 0.001, 10, "standard", 2, 1 }, // With both temperature and pressure coupling
    { 100, 0.001, 10, "standard", 0, 3 }  // Do pressure coupling not on every step
};


/*! \brief Test fixture for LeapFrog integrator.
 */
class LeapFrogTest : public ::testing::TestWithParam<LeapFrogParametersTuple>
{
public:
    //! Reference data (CPU/GPU share same data)
    TestReferenceData refData_;
    //! Checker for reference data
    TestReferenceChecker checker_;

    LeapFrogTest() : refData_(sc_refDataFilenameMaker(GetParam())), checker_(refData_.rootChecker())
    {
    }

    /*! \brief Test the numerical integrator against analytical solution for simple constant force case.
     *
     * \param[in]  tolerance  Tolerance
     * \param[in]  testData   Test data object
     * \param[in]  totalTime  Total numerical integration time
     */
    static void testAgainstAnalyticalSolution(FloatingPointTolerance  tolerance,
                                              const LeapFrogTestData& testData,
                                              const real              totalTime)
    {
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            rvec xAnalytical;
            rvec vAnalytical;
            for (int d = 0; d < DIM; d++)
            {
                // Analytical solution for constant-force particle movement
                real x0 = testData.x0_[i][d];
                real v0 = testData.v0_[i][d];
                real f  = testData.f_[i][d];
                real im = testData.inverseMasses_[i];

                xAnalytical[d] = x0 + v0 * totalTime + 0.5 * f * totalTime * totalTime * im;
                vAnalytical[d] = v0 + f * totalTime * im;

                EXPECT_REAL_EQ_TOL(xAnalytical[d], testData.xPrime_[i][d], tolerance) << gmx::formatString(
                        "Coordinate %d of atom %d is different from analytical solution.", d, i);

                EXPECT_REAL_EQ_TOL(vAnalytical[d], testData.v_[i][d], tolerance) << gmx::formatString(
                        "Velocity component %d of atom %d is different from analytical solution.", d, i);
            }
        }
    }

    /*! \brief Test the numerical integrator against pre-computed reference values.
     *
     * \param[in]  testData   Test data object
     */
    void testAgainstReferenceData(const LeapFrogTestData& testData)
    {
        TestReferenceChecker finalPositionsRef(
                checker_.checkSequenceCompound("FinalPositions", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            const gmx::RVec&     xPrime = testData.xPrime_[i];
            TestReferenceChecker xPrimeRef(finalPositionsRef.checkCompound("Atom", nullptr));
            xPrimeRef.checkReal(xPrime[XX], "XX");
            xPrimeRef.checkReal(xPrime[YY], "YY");
            xPrimeRef.checkReal(xPrime[ZZ], "ZZ");
        }

        TestReferenceChecker finalVelocitiesRef(
                checker_.checkSequenceCompound("FinalVelocities", testData.numAtoms_));
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

TEST_P(LeapFrogTest, SimpleIntegration)
{
    // Extract parameters using structured bindings
    auto [numAtoms, timestep, numSteps, vfName, numTCoupleGroups, nstpcouple] = GetParam();

    // Look up velocity and force from table
    const VelocityForcePair& vfPair = getVelocityForcePair(vfName);

    // Construct the list of runners
    std::vector<std::unique_ptr<ILeapFrogTestRunner>> runners;
    // Add runners for CPU version
    runners.emplace_back(std::make_unique<LeapFrogHostTestRunner>());
    // If supported, add runners for the GPU version for each available GPU
    if (GpuConfigurationCapabilities::Update)
    {
        for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
        {
            runners.emplace_back(std::make_unique<LeapFrogDeviceTestRunner>(*testDevice));
        }
    }

    for (const auto& runner : runners)
    {
        std::string testDescription = formatString(
                "Testing on %s with %d atoms for %d timesteps with %d temperature coupling "
                "groups and "
                "%s pressure coupling (dt = %f, v0=(%f, %f, %f), f0=(%f, %f, %f), nstpcouple = "
                "%d)",
                runner->hardwareDescription().c_str(),
                numAtoms,
                numSteps,
                numTCoupleGroups,
                nstpcouple == 0 ? "without" : "with",
                timestep,
                vfPair.velocity[XX],
                vfPair.velocity[YY],
                vfPair.velocity[ZZ],
                vfPair.force[XX],
                vfPair.force[YY],
                vfPair.force[ZZ],
                nstpcouple);
        SCOPED_TRACE(testDescription);

        std::unique_ptr<LeapFrogTestData> testData = std::make_unique<LeapFrogTestData>(
                numAtoms, timestep, vfPair.velocity, vfPair.force, numTCoupleGroups, nstpcouple);

        runner->integrate(testData.get(), numSteps);

        real totalTime = numSteps * timestep;
        // TODO For the case of constant force, the numerical scheme is exact and
        //      the only source of errors is floating point arithmetic. Hence,
        //      the tolerance can be calculated.
        FloatingPointTolerance tolerance = absoluteTolerance(numSteps * 0.000005);

        // Test against the analytical solution (without temperature coupling)
        if (numTCoupleGroups == 0 && nstpcouple == 0)
        {
            testAgainstAnalyticalSolution(tolerance, *testData, totalTime);
        }

        checker_.setDefaultTolerance(tolerance);
        testAgainstReferenceData(*testData);
    }
}

INSTANTIATE_TEST_SUITE_P(AllHardware, LeapFrogTest, ::testing::ValuesIn(parametersSets), sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
