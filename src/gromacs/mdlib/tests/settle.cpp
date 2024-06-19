/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * \brief Tests for SETTLE constraints
 *
 * The test runs on several small systems, containing 1 to 17 water molecules,
 * with and without periodic boundary conditions, with and without velocity
 * and virial updates. The CPU and GPU versions are tested, if the code was
 * compiled with CUDA support and there is a CUDA-capable GPU in the system.
 *
 * The tests check:
 * 1. If the final distances between constrained atoms are within tolerance
 *    from the target distance.
 * 2. If the velocities were updated when needed.
 * 3. If the virial was computed.
 *
 * The test also compares the results from the CPU and GPU versions of the
 * algorithm: final coordinates, velocities and virial should be within
 * tolerance to one another.
 *
 * \todo This also tests that if the calling code requires velocities
 *       and virial updates, that those outputs do change, but does not
 *       test that those changes are correct.
 *
 * \todo Only no-PBC and cubic-PBC are tested here, but the correct
 *       function of the SIMD version of set_pbx_auic in all cases
 *       should be tested elsewhere.
 *
 * \todo The CPU and GPU versions are tested against each other. This
 *       should be changed to a proper test against pre-computed
 *       reference values. Also, these test will dry-run on a CUDA
 *       build if no CUDA-capable GPU is available.
 *
 * \author Mark Abraham  <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/settle.h"

#include "config.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/hardware/device_management.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/tests/watersystem.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "settletestdata.h"
#include "settletestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Parameters that will vary from test to test.
 */
struct SettleTestParameters
{
    //! Number of water molecules (SETTLEs) [1, 2, 4, 5, 7, 10, 12, 15, 17]
    int numSettles;
    //! If the velocities should be updated while constraining [true/false]
    bool updateVelocities;
    //! If the virial should be computed [true/false]
    bool calcVirial;
    //! Periodic boundary conditions [PBCXYZ/PBCNone]
    std::string pbcName;
};

/*! \brief Sets of parameters on which to run the tests.
 */
const SettleTestParameters parametersSets[] = {
    { 1, false, false, "PBCXYZ" },   // 1 water molecule
    { 2, false, false, "PBCXYZ" },   // 2 water molecules
    { 4, false, false, "PBCXYZ" },   // 4 water molecules
    { 5, false, false, "PBCXYZ" },   // 5 water molecules
    { 6, false, false, "PBCXYZ" },   // 6 water molecules
    { 10, false, false, "PBCXYZ" },  // 10 water molecules
    { 12, false, false, "PBCXYZ" },  // 12 water molecules
    { 15, false, false, "PBCXYZ" },  // 15 water molecules
    { 17, true, false, "PBCXYZ" },   // Update velocities
    { 17, false, true, "PBCXYZ" },   // Compute virial
    { 17, false, false, "PBCNone" }, // No periodic boundary
    { 17, true, true, "PBCNone" },   // Update velocities, compute virial, without PBC
    { 17, true, true, "PBCXYZ" }
}; // Update velocities, compute virial, with PBC

/*! \brief Test fixture for testing SETTLE.
 */
class SettleTest : public ::testing::TestWithParam<SettleTestParameters>
{
public:
    //! PBC setups
    std::unordered_map<std::string, t_pbc> pbcs_;
    //! Reference data
    TestReferenceData refData_;
    //! Checker for reference data
    TestReferenceChecker checker_;

    /*! \brief Test setup function.
     *
     * Setting up the PBCs and algorithms. Note, that corresponding string keywords
     * have to be explicitly specified when parameters are initialized.
     *
     */
    SettleTest() : checker_(refData_.rootChecker())
    {

        //
        // PBC initialization
        //
        t_pbc pbc;

        // Infinitely small box
        matrix boxNone = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
        set_pbc(&pbc, PbcType::No, boxNone);
        pbcs_["PBCNone"] = pbc;

        // Rectangular box
        matrix boxXyz = { { real(1.86206), 0, 0 }, { 0, real(1.86206), 0 }, { 0, 0, real(1.86206) } };
        set_pbc(&pbc, PbcType::Xyz, boxXyz);
        pbcs_["PBCXYZ"] = pbc;
    }

    /*! \brief Check if the final interatomic distances are equal to target set by constraints.
     *
     * \param[in]  numSettles        Number of water molecules in the tested system.
     * \param[in]  tolerance         Tolerance to compare floating point numbers.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    static void checkConstrainsSatisfied(const int                    numSettles,
                                         const FloatingPointTolerance tolerance,
                                         const SettleTestData&        testData)
    {
        for (int i = 0; i < numSettles; ++i)
        {
            const gmx::RVec& positionO  = testData.xPrime_[i * testData.atomsPerSettle_ + 0];
            const gmx::RVec& positionH1 = testData.xPrime_[i * testData.atomsPerSettle_ + 1];
            const gmx::RVec& positionH2 = testData.xPrime_[i * testData.atomsPerSettle_ + 2];

            real dOH = testData.dOH_;
            real dHH = testData.dHH_;

            EXPECT_REAL_EQ_TOL(dOH * dOH, distance2(positionO, positionH1), tolerance)
                    << formatString("for water %d. ", i);
            EXPECT_REAL_EQ_TOL(dOH * dOH, distance2(positionO, positionH2), tolerance)
                    << formatString("for water %d. ", i);
            EXPECT_REAL_EQ_TOL(dHH * dHH, distance2(positionH1, positionH2), tolerance)
                    << formatString("for water %d. ", i);
        }
    }

    /*! \brief Check if the virial was updated and symmetric.
     *
     * The two tests on virial are:
     * 1. If it was updated in case calcVirial is true.
     * 2. If it is symmetrical.
     *
     * \param[in]  calcVirial        If the virial is computed.
     * \param[in]  tolerance         Tolerance to compare floating point numbers.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    static void checkVirialSymmetric(const bool                   calcVirial,
                                     const FloatingPointTolerance tolerance,
                                     const SettleTestData&        testData)
    {
        for (int d = 0; d < DIM; ++d)
        {
            for (int dd = 0; dd < DIM; ++dd)
            {

                EXPECT_TRUE(calcVirial == (0. != testData.virial_[d][dd]))
                        << formatString("for virial component[%d][%d]. ", d, dd);

                if (calcVirial)
                {
                    EXPECT_REAL_EQ_TOL(testData.virial_[d][dd], testData.virial_[dd][d], tolerance)
                            << formatString("Virial is not symmetrical for [%d][%d]. ", d, dd);
                }
            }
        }
    }

    /*! \brief Check if the final positions correspond to reference values.
     *
     * \param[in]  numSettles        Number of water molecules in the tested system.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    void checkFinalPositions(const int numSettles, const SettleTestData& testData)
    {
        TestReferenceChecker finalCoordinatesRef(
                checker_.checkSequenceCompound("FinalCoordinates", numSettles));
        for (int i = 0; i < numSettles; ++i)
        {
            TestReferenceChecker settlerRef(finalCoordinatesRef.checkCompound("Settler", nullptr));
            TestReferenceChecker atomsRef(
                    settlerRef.checkSequenceCompound("Atoms", testData.atomsPerSettle_));
            for (int j = 0; j < testData.atomsPerSettle_; ++j)
            {
                const gmx::RVec&     xPrime = testData.xPrime_[testData.atomsPerSettle_ * i + j];
                TestReferenceChecker xPrimeRef(atomsRef.checkCompound("Atom", nullptr));
                xPrimeRef.checkReal(xPrime[XX], "XX");
                xPrimeRef.checkReal(xPrime[YY], "YY");
                xPrimeRef.checkReal(xPrime[ZZ], "ZZ");
            }
        }
    }

    /*! \brief Check if the final velocities correspond to reference values.
     *
     * \param[in]  numSettles        Number of water molecules in the tested system.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    void checkFinalVelocities(const int numSettles, const SettleTestData& testData)
    {
        TestReferenceChecker finalCoordinatesRef(
                checker_.checkSequenceCompound("FinalVelocities", numSettles));
        for (int i = 0; i < numSettles; ++i)
        {
            TestReferenceChecker settlerRef(finalCoordinatesRef.checkCompound("Settler", nullptr));
            TestReferenceChecker atomsRef(
                    settlerRef.checkSequenceCompound("Atoms", testData.atomsPerSettle_));
            for (int j = 0; j < testData.atomsPerSettle_; ++j)
            {
                const gmx::RVec&     v = testData.v_[testData.atomsPerSettle_ * i + j];
                TestReferenceChecker vRef(atomsRef.checkCompound("Atom", nullptr));
                vRef.checkReal(v[XX], "XX");
                vRef.checkReal(v[YY], "YY");
                vRef.checkReal(v[ZZ], "ZZ");
            }
        }
    }

    /*! \brief Check if the computed virial correspond to reference values.
     *
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    void checkVirial(const SettleTestData& testData)
    {
        const tensor&        virial = testData.virial_;
        TestReferenceChecker virialRef(checker_.checkCompound("Virial", nullptr));

        // TODO: Is it worth it to make this in a loop??
        virialRef.checkReal(virial[XX][XX], "XX");
        virialRef.checkReal(virial[XX][YY], "XY");
        virialRef.checkReal(virial[XX][ZZ], "XZ");
        virialRef.checkReal(virial[YY][XX], "YX");
        virialRef.checkReal(virial[YY][YY], "YY");
        virialRef.checkReal(virial[YY][ZZ], "YZ");
        virialRef.checkReal(virial[ZZ][XX], "ZX");
        virialRef.checkReal(virial[ZZ][YY], "ZY");
        virialRef.checkReal(virial[ZZ][ZZ], "ZZ");
    }
};

TEST_P(SettleTest, SatisfiesConstraints)
{
    // Construct the list of runners
    std::vector<std::unique_ptr<ISettleTestRunner>> runners;
    // Add runners for CPU version
    runners.emplace_back(std::make_unique<SettleHostTestRunner>());
    // If supported, add runners for the GPU version for each available GPU
    const bool addGpuRunners = GPU_SETTLE_SUPPORTED;
    if (addGpuRunners)
    {
        for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
        {
            runners.emplace_back(std::make_unique<SettleDeviceTestRunner>(*testDevice));
        }
    }
    for (const auto& runner : runners)
    {
        // Make some symbolic names for the parameter combination.
        SettleTestParameters params = GetParam();

        int         numSettles       = params.numSettles;
        bool        updateVelocities = params.updateVelocities;
        bool        calcVirial       = params.calcVirial;
        std::string pbcName          = params.pbcName;


        // Make a string that describes which parameter combination is
        // being tested, to help make failing tests comprehensible.
        std::string testDescription = formatString(
                "Testing %s with %d SETTLEs, %s, %svelocities and %scalculating the virial.",
                runner->hardwareDescription().c_str(),
                numSettles,
                pbcName.c_str(),
                updateVelocities ? "with " : "without ",
                calcVirial ? "" : "not ");

        SCOPED_TRACE(testDescription);

        auto testData = std::make_unique<SettleTestData>(numSettles);

        ASSERT_LE(numSettles, testData->xPrime_.size() / testData->atomsPerSettle_)
                << "cannot test that many SETTLEs. " << testDescription;

        t_pbc pbc = pbcs_.at(pbcName);

        // Apply SETTLE
        runner->applySettle(testData.get(), pbc, updateVelocities, calcVirial, testDescription);

        // The necessary tolerances for the test to pass were determined
        // empirically. This isn't nice, but the required behavior that
        // SETTLE produces constrained coordinates consistent with
        // sensible sampling needs to be tested at a much higher level.
        // TODO: Re-evaluate the tolerances.
        real                   dOH       = testData->dOH_;
        FloatingPointTolerance tolerance = relativeToleranceAsPrecisionDependentUlp(dOH * dOH, 80, 380);
        FloatingPointTolerance toleranceVirial = absoluteTolerance(0.000001);

        FloatingPointTolerance tolerancePositions  = absoluteTolerance(0.000001);
        FloatingPointTolerance toleranceVelocities = absoluteTolerance(0.0001);

        checkConstrainsSatisfied(numSettles, tolerance, *testData);
        checkVirialSymmetric(calcVirial, toleranceVirial, *testData);

        checker_.setDefaultTolerance(tolerancePositions);
        checkFinalPositions(numSettles, *testData);

        if (updateVelocities)
        {
            checker_.setDefaultTolerance(toleranceVelocities);
            checkFinalVelocities(numSettles, *testData);
        }

        if (calcVirial)
        {
            checker_.setDefaultTolerance(toleranceVirial);
            checkVirial(*testData);
        }
    }
}

// Run test on pre-determined set of combinations for test parameters, which include the numbers of SETTLEs (water
// molecules), whether or not velocities are updated and virial contribution is computed, was the PBC enabled.
// The test will cycle through all available runners, including CPU and, if applicable, GPU implementations of SETTLE.
INSTANTIATE_TEST_SUITE_P(WithParameters, SettleTest, ::testing::ValuesIn(parametersSets));

} // namespace
} // namespace test
} // namespace gmx
