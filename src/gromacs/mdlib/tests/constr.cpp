/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief SHAKE and LINCS tests.
 *
 * \todo Better tests for virial are needed.
 * \todo Tests for bigger systems to test threads synchronization,
 *       reduction, etc. on the GPU.
 * \todo Tests for algorithms for derivatives.
 * \todo Free-energy perturbation tests
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"

#include <assert.h>

#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "constrtestdata.h"
#include "constrtestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture for constraints.
 *
 * The fixture uses following test systems:
 * 1. Two atoms, connected with one constraint (e.g. NH).
 * 2. Three atoms, connected consequently with two constraints (e.g. CH2).
 * 3. Three atoms, constrained to the fourth atom (e.g. CH3).
 * 4. Four atoms, connected by two independent constraints.
 * 5. Three atoms, connected by three constraints in a triangle
 *      (e.g. H2O with constrained H-O-H angle).
 * 6. Four atoms, connected by three consequential constraints.
 *
 * For all systems, the final lengths of the constraints are tested against the
 * reference values, the direction of each constraint is checked.
 * Test also verifies that the center of mass has not been
 * shifted by the constraints and that its velocity has not changed.
 * For some systems, the value for scaled virial tensor is checked against
 * pre-computed data.
 */
class ConstraintsTest : public ::testing::TestWithParam<std::string>
{
public:
    //! PBC setups
    std::unordered_map<std::string, t_pbc> pbcs_;

    /*! \brief Test setup function.
     *
     * Setting up the pbcs and algorithms. Note, that corresponding string keywords
     * have to be explicitly added at the end of this file when the tests are called.
     *
     */
    void SetUp() override
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
        matrix boxXyz = { { 10.0, 0.0, 0.0 }, { 0.0, 20.0, 0.0 }, { 0.0, 0.0, 15.0 } };
        set_pbc(&pbc, PbcType::Xyz, boxXyz);
        pbcs_["PBCXYZ"] = pbc;
    }

    /*! \brief
     * The test on the final length of constrained bonds.
     *
     * Goes through all the constraints and checks if the final length of all the constraints is
     * equal to the target length with provided tolerance.
     *
     * \param[in] tolerance       Allowed tolerance in final lengths.
     * \param[in] testData        Test data structure.
     * \param[in] pbc             Periodic boundary data.
     */
    static void checkConstrainsLength(FloatingPointTolerance     tolerance,
                                      const ConstraintsTestData& testData,
                                      t_pbc                      pbc)
    {

        // Test if all the constraints are satisfied
        for (index c = 0; c < ssize(testData.constraints_) / 3; c++)
        {
            real r0 = testData.constraintsR0_.at(testData.constraints_.at(3 * c));
            int  i  = testData.constraints_.at(3 * c + 1);
            int  j  = testData.constraints_.at(3 * c + 2);
            RVec xij0, xij1;
            real d0, d1;
            if (pbc.pbcType == PbcType::Xyz)
            {
                pbc_dx_aiuc(&pbc, testData.x_[i], testData.x_[j], xij0);
                pbc_dx_aiuc(&pbc, testData.xPrime_[i], testData.xPrime_[j], xij1);
            }
            else
            {
                rvec_sub(testData.x_[i], testData.x_[j], xij0);
                rvec_sub(testData.xPrime_[i], testData.xPrime_[j], xij1);
            }
            d0 = norm(xij0);
            d1 = norm(xij1);
            EXPECT_REAL_EQ_TOL(r0, d1, tolerance) << gmx::formatString(
                    "rij = %f, which is not equal to r0 = %f for constraint #%zd, between atoms %d "
                    "and %d"
                    " (before constraining rij was %f).",
                    d1, r0, c, i, j, d0);
        }
    }

    /*! \brief
     * The test on the final length of constrained bonds.
     *
     * Goes through all the constraints and checks if the direction of constraint has not changed
     * by the algorithm (i.e. the constraints algorithm arrived to the solution that is closest
     * to the initial system conformation).
     *
     * \param[in] testData        Test data structure.
     * \param[in] pbc             Periodic boundary data.
     */
    static void checkConstrainsDirection(const ConstraintsTestData& testData, t_pbc pbc)
    {

        for (index c = 0; c < ssize(testData.constraints_) / 3; c++)
        {
            int  i = testData.constraints_.at(3 * c + 1);
            int  j = testData.constraints_.at(3 * c + 2);
            RVec xij0, xij1;
            if (pbc.pbcType == PbcType::Xyz)
            {
                pbc_dx_aiuc(&pbc, testData.x_[i], testData.x_[j], xij0);
                pbc_dx_aiuc(&pbc, testData.xPrime_[i], testData.xPrime_[j], xij1);
            }
            else
            {
                rvec_sub(testData.x_[i], testData.x_[j], xij0);
                rvec_sub(testData.xPrime_[i], testData.xPrime_[j], xij1);
            }

            real dot = xij0.dot(xij1);

            EXPECT_GE(dot, 0.0) << gmx::formatString(
                    "The constraint %zd changed direction. Constraining algorithm might have "
                    "returned the wrong root "
                    "of the constraints equation.",
                    c);
        }
    }

    /*! \brief
     * The test on the coordinates of the center of the mass (COM) of the system.
     *
     * Checks if the center of mass has not been shifted by the constraints. Note,
     * that this test does not take into account the periodic boundary conditions.
     * Hence it will not work should the constraints decide to move atoms across
     * PBC borders.
     *
     * \param[in] tolerance       Allowed tolerance in COM coordinates.
     * \param[in] testData        Test data structure.
     */
    static void checkCOMCoordinates(FloatingPointTolerance tolerance, const ConstraintsTestData& testData)
    {

        RVec comPrime0({ 0.0, 0.0, 0.0 });
        RVec comPrime({ 0.0, 0.0, 0.0 });
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            comPrime0 += testData.masses_[i] * testData.xPrime0_[i];
            comPrime += testData.masses_[i] * testData.xPrime_[i];
        }

        comPrime0 /= testData.numAtoms_;
        comPrime /= testData.numAtoms_;

        EXPECT_REAL_EQ_TOL(comPrime[XX], comPrime0[XX], tolerance)
                << "Center of mass was shifted by constraints in x-direction.";
        EXPECT_REAL_EQ_TOL(comPrime[YY], comPrime0[YY], tolerance)
                << "Center of mass was shifted by constraints in y-direction.";
        EXPECT_REAL_EQ_TOL(comPrime[ZZ], comPrime0[ZZ], tolerance)
                << "Center of mass was shifted by constraints in z-direction.";
    }

    /*! \brief
     * The test on the velocity of the center of the mass (COM) of the system.
     *
     * Checks if the velocity of the center of mass has not changed.
     *
     * \param[in] tolerance       Allowed tolerance in COM velocity components.
     * \param[in] testData        Test data structure.
     */
    static void checkCOMVelocity(FloatingPointTolerance tolerance, const ConstraintsTestData& testData)
    {

        RVec comV0({ 0.0, 0.0, 0.0 });
        RVec comV({ 0.0, 0.0, 0.0 });
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            comV0 += testData.masses_[i] * testData.v0_[i];
            comV += testData.masses_[i] * testData.v_[i];
        }
        comV0 /= testData.numAtoms_;
        comV /= testData.numAtoms_;

        EXPECT_REAL_EQ_TOL(comV[XX], comV0[XX], tolerance)
                << "Velocity of the center of mass in x-direction has been changed by constraints.";
        EXPECT_REAL_EQ_TOL(comV[YY], comV0[YY], tolerance)
                << "Velocity of the center of mass in y-direction has been changed by constraints.";
        EXPECT_REAL_EQ_TOL(comV[ZZ], comV0[ZZ], tolerance)
                << "Velocity of the center of mass in z-direction has been changed by constraints.";
    }

    /*! \brief
     * The test of virial tensor.
     *
     * Checks if the values in the scaled virial tensor are equal to pre-computed values.
     *
     * \param[in] tolerance       Tolerance for the tensor values.
     * \param[in] testData        Test data structure.
     */
    static void checkVirialTensor(FloatingPointTolerance tolerance, const ConstraintsTestData& testData)
    {
        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                EXPECT_REAL_EQ_TOL(testData.virialScaledRef_[i][j], testData.virialScaled_[i][j], tolerance)
                        << gmx::formatString(
                                   "Values in virial tensor at [%d][%d] are not within the "
                                   "tolerance from reference value.",
                                   i, j);
            }
        }
    }
    //! Before any test is run, work out whether any compatible GPUs exist.
    static std::vector<std::unique_ptr<IConstraintsTestRunner>> getRunners()
    {
        std::vector<std::unique_ptr<IConstraintsTestRunner>> runners;
        // Add runners for CPU versions of SHAKE and LINCS
        runners.emplace_back(std::make_unique<ShakeConstraintsRunner>());
        runners.emplace_back(std::make_unique<LincsConstraintsRunner>());
        // If using CUDA, add runners for the GPU version of LINCS for each available GPU
        if (GMX_GPU_CUDA)
        {
            for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
            {
                runners.emplace_back(std::make_unique<LincsDeviceConstraintsRunner>(*testDevice));
            }
        }
        return runners;
    }
};

TEST_P(ConstraintsTest, SingleConstraint)
{
    std::string title    = "one constraint (e.g. OH)";
    int         numAtoms = 2;

    std::vector<real> masses        = { 1.0, 12.0 };
    std::vector<int>  constraints   = { 0, 0, 1 };
    std::vector<real> constraintsR0 = { 0.1 };

    real oneTenthOverSqrtTwo = 0.1_real / std::sqrt(2.0_real);

    std::vector<RVec> x = { { 0.0, oneTenthOverSqrtTwo, 0.0 }, { oneTenthOverSqrtTwo, 0.0, 0.0 } };
    std::vector<RVec> xPrime = { { 0.01, 0.08, 0.01 }, { 0.06, 0.01, -0.01 } };
    std::vector<RVec> v      = { { 1.0, 2.0, 3.0 }, { 3.0, 2.0, 1.0 } };

    tensor virialScaledRef = { { -5.58e-04, 5.58e-04, 0.00e+00 },
                               { 5.58e-04, -5.58e-04, 0.00e+00 },
                               { 0.00e+00, 0.00e+00, 0.00e+00 } };

    real     shakeTolerance = 0.0001;
    gmx_bool shakeUseSOR    = false;

    int  lincsNIter               = 1;
    int  lincslincsExpansionOrder = 4;
    real lincsWarnAngle           = 30.0;

    std::unique_ptr<ConstraintsTestData> testData = std::make_unique<ConstraintsTestData>(
            title, numAtoms, masses, constraints, constraintsR0, true, virialScaledRef, false, 0,
            real(0.0), real(0.001), x, xPrime, v, shakeTolerance, shakeUseSOR, lincsNIter,
            lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName = GetParam();
    t_pbc       pbc     = pbcs_.at(pbcName);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s using %s.", testData->title_.c_str(),
                                  pbcName.c_str(), runner->name().c_str()));

        testData->reset();

        // Apply constraints
        runner->applyConstraints(testData.get(), pbc);

        checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
        checkConstrainsDirection(*testData, pbc);
        checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
        checkCOMVelocity(absoluteTolerance(0.0001), *testData);

        checkVirialTensor(absoluteTolerance(0.0001), *testData);
    }
}

TEST_P(ConstraintsTest, TwoDisjointConstraints)
{

    std::string       title         = "two disjoint constraints";
    int               numAtoms      = 4;
    std::vector<real> masses        = { 0.5, 1.0 / 3.0, 0.25, 1.0 };
    std::vector<int>  constraints   = { 0, 0, 1, 1, 2, 3 };
    std::vector<real> constraintsR0 = { 2.0, 1.0 };


    std::vector<RVec> x = {
        { 2.50, -3.10, 15.70 }, { 0.51, -3.02, 15.55 }, { -0.50, -3.00, 15.20 }, { -1.51, -2.95, 15.05 }
    };

    std::vector<RVec> xPrime = {
        { 2.50, -3.10, 15.70 }, { 0.51, -3.02, 15.55 }, { -0.50, -3.00, 15.20 }, { -1.51, -2.95, 15.05 }
    };

    std::vector<RVec> v = { { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 } };

    tensor virialScaledRef = { { 3.3e-03, -1.7e-04, 5.6e-04 },
                               { -1.7e-04, 8.9e-06, -2.8e-05 },
                               { 5.6e-04, -2.8e-05, 8.9e-05 } };

    real     shakeTolerance = 0.0001;
    gmx_bool shakeUseSOR    = false;

    int  lincsNIter               = 1;
    int  lincslincsExpansionOrder = 4;
    real lincsWarnAngle           = 30.0;

    std::unique_ptr<ConstraintsTestData> testData = std::make_unique<ConstraintsTestData>(
            title, numAtoms, masses, constraints, constraintsR0, true, virialScaledRef, false, 0,
            real(0.0), real(0.001), x, xPrime, v, shakeTolerance, shakeUseSOR, lincsNIter,
            lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName = GetParam();
    t_pbc       pbc     = pbcs_.at(pbcName);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s using %s.", testData->title_.c_str(),
                                  pbcName.c_str(), runner->name().c_str()));

        testData->reset();

        // Apply constraints
        runner->applyConstraints(testData.get(), pbc);

        checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
        checkConstrainsDirection(*testData, pbc);
        checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
        checkCOMVelocity(absoluteTolerance(0.0001), *testData);

        checkVirialTensor(absoluteTolerance(0.0001), *testData);
    }
}

TEST_P(ConstraintsTest, ThreeSequentialConstraints)
{

    std::string       title         = "three atoms, connected longitudinally (e.g. CH2)";
    int               numAtoms      = 3;
    std::vector<real> masses        = { 1.0, 12.0, 16.0 };
    std::vector<int>  constraints   = { 0, 0, 1, 1, 1, 2 };
    std::vector<real> constraintsR0 = { 0.1, 0.2 };

    real oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);
    real twoTenthsOverSqrtThree = 0.2_real / std::sqrt(3.0_real);

    std::vector<RVec> x = { { oneTenthOverSqrtTwo, oneTenthOverSqrtTwo, 0.0 },
                            { 0.0, 0.0, 0.0 },
                            { twoTenthsOverSqrtThree, twoTenthsOverSqrtThree, twoTenthsOverSqrtThree } };

    std::vector<RVec> xPrime = { { 0.08, 0.07, 0.01 }, { -0.02, 0.01, -0.02 }, { 0.10, 0.12, 0.11 } };

    std::vector<RVec> v = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

    tensor virialScaledRef = { { 4.14e-03, 4.14e-03, 3.31e-03 },
                               { 4.14e-03, 4.14e-03, 3.31e-03 },
                               { 3.31e-03, 3.31e-03, 3.31e-03 } };

    real     shakeTolerance = 0.0001;
    gmx_bool shakeUseSOR    = false;

    int  lincsNIter               = 1;
    int  lincslincsExpansionOrder = 4;
    real lincsWarnAngle           = 30.0;

    std::unique_ptr<ConstraintsTestData> testData = std::make_unique<ConstraintsTestData>(
            title, numAtoms, masses, constraints, constraintsR0, true, virialScaledRef, false, 0,
            real(0.0), real(0.001), x, xPrime, v, shakeTolerance, shakeUseSOR, lincsNIter,
            lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName = GetParam();
    t_pbc       pbc     = pbcs_.at(pbcName);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s using %s.", testData->title_.c_str(),
                                  pbcName.c_str(), runner->name().c_str()));

        testData->reset();

        // Apply constraints
        runner->applyConstraints(testData.get(), pbc);

        checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
        checkConstrainsDirection(*testData, pbc);
        checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
        checkCOMVelocity(absoluteTolerance(0.0001), *testData);

        checkVirialTensor(absoluteTolerance(0.0001), *testData);
    }
}

TEST_P(ConstraintsTest, ThreeConstraintsWithCentralAtom)
{

    std::string       title         = "three atoms, connected to the central atom (e.g. CH3)";
    int               numAtoms      = 4;
    std::vector<real> masses        = { 12.0, 1.0, 1.0, 1.0 };
    std::vector<int>  constraints   = { 0, 0, 1, 0, 0, 2, 0, 0, 3 };
    std::vector<real> constraintsR0 = { 0.1 };


    std::vector<RVec> x = {
        { 0.00, 0.00, 0.00 }, { 0.10, 0.00, 0.00 }, { 0.00, -0.10, 0.00 }, { 0.00, 0.00, 0.10 }
    };

    std::vector<RVec> xPrime = { { 0.004, 0.009, -0.010 },
                                 { 0.110, -0.006, 0.003 },
                                 { -0.007, -0.102, -0.007 },
                                 { -0.005, 0.011, 0.102 } };

    std::vector<RVec> v = { { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };

    tensor virialScaledRef = { { 7.14e-04, 0.00e+00, 0.00e+00 },
                               { 0.00e+00, 1.08e-03, 0.00e+00 },
                               { 0.00e+00, 0.00e+00, 1.15e-03 } };

    real     shakeTolerance = 0.0001;
    gmx_bool shakeUseSOR    = false;

    int  lincsNIter               = 1;
    int  lincslincsExpansionOrder = 4;
    real lincsWarnAngle           = 30.0;

    std::unique_ptr<ConstraintsTestData> testData = std::make_unique<ConstraintsTestData>(
            title, numAtoms, masses, constraints, constraintsR0, true, virialScaledRef, false, 0,
            real(0.0), real(0.001), x, xPrime, v, shakeTolerance, shakeUseSOR, lincsNIter,
            lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName = GetParam();
    t_pbc       pbc     = pbcs_.at(pbcName);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s using %s.", testData->title_.c_str(),
                                  pbcName.c_str(), runner->name().c_str()));

        testData->reset();

        // Apply constraints
        runner->applyConstraints(testData.get(), pbc);

        checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
        checkConstrainsDirection(*testData, pbc);
        checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
        checkCOMVelocity(absoluteTolerance(0.0001), *testData);

        checkVirialTensor(absoluteTolerance(0.0001), *testData);
    }
}

TEST_P(ConstraintsTest, FourSequentialConstraints)
{

    std::string       title         = "four atoms, connected longitudinally";
    int               numAtoms      = 4;
    std::vector<real> masses        = { 0.5, 1.0 / 3.0, 0.25, 1.0 };
    std::vector<int>  constraints   = { 0, 0, 1, 1, 1, 2, 2, 2, 3 };
    std::vector<real> constraintsR0 = { 2.0, 1.0, 1.0 };


    std::vector<RVec> x = {
        { 2.50, -3.10, 15.70 }, { 0.51, -3.02, 15.55 }, { -0.50, -3.00, 15.20 }, { -1.51, -2.95, 15.05 }
    };

    std::vector<RVec> xPrime = {
        { 2.50, -3.10, 15.70 }, { 0.51, -3.02, 15.55 }, { -0.50, -3.00, 15.20 }, { -1.51, -2.95, 15.05 }
    };

    std::vector<RVec> v = { { 0.0, 0.0, 2.0 }, { 0.0, 0.0, 3.0 }, { 0.0, 0.0, -4.0 }, { 0.0, 0.0, -1.0 } };

    tensor virialScaledRef = { { 1.15e-01, -4.20e-03, 2.12e-02 },
                               { -4.20e-03, 1.70e-04, -6.41e-04 },
                               { 2.12e-02, -6.41e-04, 5.45e-03 } };

    real     shakeTolerance = 0.0001;
    gmx_bool shakeUseSOR    = false;

    int  lincsNIter               = 4;
    int  lincslincsExpansionOrder = 8;
    real lincsWarnAngle           = 30.0;

    std::unique_ptr<ConstraintsTestData> testData = std::make_unique<ConstraintsTestData>(
            title, numAtoms, masses, constraints, constraintsR0, true, virialScaledRef, false, 0,
            real(0.0), real(0.001), x, xPrime, v, shakeTolerance, shakeUseSOR, lincsNIter,
            lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName = GetParam();
    t_pbc       pbc     = pbcs_.at(pbcName);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s using %s.", testData->title_.c_str(),
                                  pbcName.c_str(), runner->name().c_str()));

        testData->reset();

        // Apply constraints
        runner->applyConstraints(testData.get(), pbc);

        checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
        checkConstrainsDirection(*testData, pbc);
        checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
        checkCOMVelocity(absoluteTolerance(0.0001), *testData);

        checkVirialTensor(absoluteTolerance(0.01), *testData);
    }
}

TEST_P(ConstraintsTest, TriangleOfConstraints)
{

    std::string       title         = "basic triangle (tree atoms, connected to each other)";
    int               numAtoms      = 3;
    std::vector<real> masses        = { 1.0, 1.0, 1.0 };
    std::vector<int>  constraints   = { 0, 0, 1, 2, 0, 2, 1, 1, 2 };
    std::vector<real> constraintsR0 = { 0.1, 0.1, 0.1 };

    real oneTenthOverSqrtTwo = 0.1_real / std::sqrt(2.0_real);

    std::vector<RVec> x = { { oneTenthOverSqrtTwo, 0.0, 0.0 },
                            { 0.0, oneTenthOverSqrtTwo, 0.0 },
                            { 0.0, 0.0, oneTenthOverSqrtTwo } };

    std::vector<RVec> xPrime = { { 0.09, -0.02, 0.01 }, { -0.02, 0.10, -0.02 }, { 0.03, -0.01, 0.07 } };

    std::vector<RVec> v = { { 1.0, 1.0, 1.0 }, { -2.0, -2.0, -2.0 }, { 1.0, 1.0, 1.0 } };

    tensor virialScaledRef = { { 6.00e-04, -1.61e-03, 1.01e-03 },
                               { -1.61e-03, 2.53e-03, -9.25e-04 },
                               { 1.01e-03, -9.25e-04, -8.05e-05 } };

    real     shakeTolerance = 0.0001;
    gmx_bool shakeUseSOR    = false;

    int  lincsNIter               = 1;
    int  lincslincsExpansionOrder = 4;
    real lincsWarnAngle           = 30.0;

    std::unique_ptr<ConstraintsTestData> testData = std::make_unique<ConstraintsTestData>(
            title, numAtoms, masses, constraints, constraintsR0, true, virialScaledRef, false, 0,
            real(0.0), real(0.001), x, xPrime, v, shakeTolerance, shakeUseSOR, lincsNIter,
            lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName = GetParam();
    t_pbc       pbc     = pbcs_.at(pbcName);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s using %s.", testData->title_.c_str(),
                                  pbcName.c_str(), runner->name().c_str()));

        testData->reset();

        // Apply constraints
        runner->applyConstraints(testData.get(), pbc);

        checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
        checkConstrainsDirection(*testData, pbc);
        checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
        checkCOMVelocity(absoluteTolerance(0.0001), *testData);

        checkVirialTensor(absoluteTolerance(0.00001), *testData);
    }
}


INSTANTIATE_TEST_CASE_P(WithParameters, ConstraintsTest, ::testing::Values("PBCNone", "PBCXYZ"));

} // namespace
} // namespace test
} // namespace gmx
