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
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include <cmath>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "langevintestdata.h"
#include "langevintestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief The parameters for the test.
 *
 * The test will run for combinations of:
 *
 * 1. Number of atoms
 * 2. Timestep
 * 3. Number of steps
 * 4. Velocity components
 * 5. Force components
 * 6. Number of temperature coupling groups
 * 7. Temperatures
 * 8. Tau-t (inverse friction constant) for all temp. coupling groups
 * 9. Random seed
 */
struct LangevinTestParameters
{
    //! Total number of atoms
    int numAtoms;
    //! Timestep
    real timestep;
    //! Number of integration steps
    int numSteps;
    //! Initial velocity
    RVec v;
    //! Constant force
    RVec f;
    //! Number of temperature coupling group (zero for no temperature coupling)
    int numTCoupleGroups;
    //! Temperature
    real temperature;
    //! tau-t is controlling the inverse friction constant
    real tauT;
    //! Random seed for the Langevin integrator
    int seed;
};

//! The set of parameters combinations to run the test on
const LangevinTestParameters parametersSets[] = {
    { 1, 0.001, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 0, 2, 123 }, // 1 particle, T = 0K
    { 1, 0.001, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 0, 2, 12345 }, // 1 particle, T=0K, other random seed
    { 1, 0.001, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 123 },  // T = 100K
    { 1, 0.001, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 0, 123 },  // tau-t = 0
    { 1, 0.0005, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 123 }, // 0.0005 ps timestep
    { 1, 0.0005, 2, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 123 }, // 0.0005 ps timestep, 2 steps
    { 1, 0.0025, 3, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 123 }, // 0.0025 ps timestep, 3 steps
    { 1, 0.0025, 3, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 12345 }, // Other random seed
    { 1, 0.0025, 3, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 500, 2, 123 }, // 0.0025 ps timestep, 3 steps, T = 500K
    { 1, 0.0025, 20, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 0, 2, 123 },   // 20 steps, T = 0K
    { 1, 0.0025, 20, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 123 }, // 20 steps, T = 100K
    { 1, 0.0025, 20, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 12345 }, // Other random seed
    { 10, 0.0025, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 1, 100, 2, 123 },   // 10 particles
    { 10, 0.0025, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 2, 100, 2, 123 }, // 2 temperature couple groups
    { 10, 0.0025, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 2, 100, 2, 12345 }, // Other random seed
    { 10, 0.0025, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 2, 100, 0.5, 123 }, // tau-t = 0.5
    { 10, 0.0025, 1, { 1.0, -2.0, 3.0 }, { -3.0, 2.0, -1.0 }, 5, 100, 2, 12345 }, // 5 temperature couple groups
};


/*! \brief Test fixture for Langevin integrator.
 */
class LangevinTest : public ::testing::TestWithParam<LangevinTestParameters>
{
public:
    //! Reference data
    TestReferenceData refData_;
    //! Checker for reference data
    TestReferenceChecker checker_;

    LangevinTest() : checker_(refData_.rootChecker()) {}

    /*! \brief Test the numerical integrator against pre-computed reference values.
     *
     * \param[in]  testData   Test data object
     */
    void testAgainstReferenceData(const LangevinTestData& testData)
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

TEST_P(LangevinTest, SimpleIntegration)
{
    // Construct the list of runners
    std::vector<std::unique_ptr<ILangevinTestRunner>> runners;
    // Add runners for CPU version
    runners.emplace_back(std::make_unique<LangevinHostTestRunner>());

    for (const auto& runner : runners)
    {
        LangevinTestParameters parameters = GetParam();

        std::string testDescription = formatString(
                "Testing on %s with %d atoms for %d timesteps with %d temperature coupling "
                "groups (dt = %f, v0=(%f, %f, %f), f0=(%f, %f, %f), T = %f, tau-t = %f, seed = %d",
                runner->hardwareDescription().c_str(),
                parameters.numAtoms,
                parameters.numSteps,
                parameters.numTCoupleGroups,
                parameters.timestep,
                parameters.v[XX],
                parameters.v[YY],
                parameters.v[ZZ],
                parameters.f[XX],
                parameters.f[YY],
                parameters.f[ZZ],
                parameters.temperature,
                parameters.tauT,
                parameters.seed);
        SCOPED_TRACE(testDescription);

        std::unique_ptr<LangevinTestData> testData =
                std::make_unique<LangevinTestData>(parameters.numAtoms,
                                                   parameters.timestep,
                                                   parameters.v,
                                                   parameters.f,
                                                   parameters.numTCoupleGroups,
                                                   parameters.temperature,
                                                   parameters.tauT,
                                                   parameters.seed);

        runner->integrate(testData.get(), parameters.numSteps);

        FloatingPointTolerance tolerance =
                absoluteTolerance(parameters.numSteps * (GMX_DOUBLE ? 5e-10 : 5e-6));

        checker_.setDefaultTolerance(tolerance);
        testAgainstReferenceData(*testData);
    }
}

INSTANTIATE_TEST_SUITE_P(WithParameters, LangevinTest, ::testing::ValuesIn(parametersSets));

} // namespace
} // namespace test
} // namespace gmx
