/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include <cassert>
#include <cmath>

#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "constrtestdata.h"
#include "constrtestrunners.h"

//! Helper function to convert t_pbc into string and make test failure messages readable
static void PrintTo(const t_pbc& pbc, std::ostream* os)
{
    *os << "PBC: " << c_pbcTypeNames[pbc.pbcType];
}


namespace gmx
{

namespace test
{
namespace
{

// Define the set of PBCs to run the test for
const std::vector<t_pbc> c_pbcs = [] {
    std::vector<t_pbc> pbcs;
    t_pbc              pbc;

    // Infinitely small box
    matrix boxNone = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    set_pbc(&pbc, PbcType::No, boxNone);
    pbcs.emplace_back(pbc);

    // Rectangular box
    matrix boxXyz = { { 10.0, 0.0, 0.0 }, { 0.0, 20.0, 0.0 }, { 0.0, 0.0, 15.0 } };
    set_pbc(&pbc, PbcType::Xyz, boxXyz);
    pbcs.emplace_back(pbc);

    return pbcs;
}();

struct ConstraintsTestSystem
{
    //! Human-friendly name of the system.
    std::string title;
    //! Number of atoms in the system.
    int numAtoms;
    //! Atom masses. Size of this vector should be equal to numAtoms.
    std::vector<real> masses;
    /*! \brief List of constraints, organized in triples of integers.
     *
     * First integer is the index of type for a constraint, second
     * and third are the indices of constrained atoms. The types
     * of constraints should be sequential but not necessarily
     * start from zero (which is the way they normally are in
     * GROMACS).
     */
    std::vector<int> constraints;
    /*! \brief Target values for bond lengths for bonds of each type.
     *
     * The size of this vector should be equal to the total number of
     * unique types in constraints vector.
     */
    std::vector<real> constraintsR0;
    //! Coordinates before integration step.
    std::vector<RVec> x;
    //! Coordinates after integration step, but before constraining.
    std::vector<RVec> xPrime;
    //! Velocities before constraining.
    std::vector<RVec> v;

    //! Target tolerance for SHAKE.
    real shakeTolerance = 0.00002;
    /*! \brief Use successive over-relaxation method for SHAKE iterations.
     *
     * The general formula is:
     * x_n+1 = (1-omega)*x_n + omega*f(x_n),
     * where omega = 1 if SOR is off and may be < 1 if SOR is on.
     */
    bool shakeUseSOR = false;

    //! Number of iterations used to compute the inverse matrix.
    int lincsNIter = 1;
    //! The order for algorithm that adjusts the direction of the bond after constraints are applied.
    int lincslincsExpansionOrder = 4;
    //! The threshold value for the change in bond angle. When exceeded the program will issue a warning
    real lincsWarnAngle = 30.0;
};

//! Helper function to convert ConstraintsTestSystem into string and make test failure messages readable
void PrintTo(const ConstraintsTestSystem& constraintsTestSystem, std::ostream* os)
{
    *os << constraintsTestSystem.title << " - " << constraintsTestSystem.numAtoms << " atoms";
}

const std::vector<ConstraintsTestSystem> c_constraintsTestSystemList = [] {
    std::vector<ConstraintsTestSystem> constraintsTestSystemList;
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title    = "one constraint (e.g. OH)";
        constraintsTestSystem.numAtoms = 2;

        constraintsTestSystem.masses        = { 1.0, 12.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1 };
        constraintsTestSystem.constraintsR0 = { 0.1 };

        real oneTenthOverSqrtTwo = 0.1_real / std::sqrt(2.0_real);

        constraintsTestSystem.x = { { 0.0, oneTenthOverSqrtTwo, 0.0 }, { oneTenthOverSqrtTwo, 0.0, 0.0 } };
        constraintsTestSystem.xPrime = { { 0.01, 0.08, 0.01 }, { 0.06, 0.01, -0.01 } };
        constraintsTestSystem.v      = { { 1.0, 2.0, 3.0 }, { 3.0, 2.0, 1.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title         = "two disjoint constraints";
        constraintsTestSystem.numAtoms      = 4;
        constraintsTestSystem.masses        = { 0.5, 1.0 / 3.0, 0.25, 1.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1, 1, 2, 3 };
        constraintsTestSystem.constraintsR0 = { 2.0, 1.0 };


        constraintsTestSystem.x = { { 2.50, -3.10, 15.70 },
                                    { 0.51, -3.02, 15.55 },
                                    { -0.50, -3.00, 15.20 },
                                    { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.xPrime = { { 2.50, -3.10, 15.70 },
                                         { 0.51, -3.02, 15.55 },
                                         { -0.50, -3.00, 15.20 },
                                         { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.v = { { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title         = "three atoms, connected longitudinally (e.g. CH2)";
        constraintsTestSystem.numAtoms      = 3;
        constraintsTestSystem.masses        = { 1.0, 12.0, 16.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1, 1, 1, 2 };
        constraintsTestSystem.constraintsR0 = { 0.1, 0.2 };

        real oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);
        real twoTenthsOverSqrtThree = 0.2_real / std::sqrt(3.0_real);

        constraintsTestSystem.x = { { oneTenthOverSqrtTwo, oneTenthOverSqrtTwo, 0.0 },
                                    { 0.0, 0.0, 0.0 },
                                    { twoTenthsOverSqrtThree, twoTenthsOverSqrtThree, twoTenthsOverSqrtThree } };

        constraintsTestSystem.xPrime = { { 0.08, 0.07, 0.01 }, { -0.02, 0.01, -0.02 }, { 0.10, 0.12, 0.11 } };

        constraintsTestSystem.v = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title         = "four atoms, connected longitudinally";
        constraintsTestSystem.numAtoms      = 4;
        constraintsTestSystem.masses        = { 0.5, 1.0 / 3.0, 0.25, 1.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1, 1, 1, 2, 2, 2, 3 };
        constraintsTestSystem.constraintsR0 = { 2.0, 1.0, 1.0 };


        constraintsTestSystem.x = { { 2.50, -3.10, 15.70 },
                                    { 0.51, -3.02, 15.55 },
                                    { -0.50, -3.00, 15.20 },
                                    { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.xPrime = { { 2.50, -3.10, 15.70 },
                                         { 0.51, -3.02, 15.55 },
                                         { -0.50, -3.00, 15.20 },
                                         { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.v = {
            { 0.0, 0.0, 2.0 }, { 0.0, 0.0, 3.0 }, { 0.0, 0.0, -4.0 }, { 0.0, 0.0, -1.0 }
        };

        // Overriding default values since LINCS converges slowly for this system.
        constraintsTestSystem.lincsNIter               = 4;
        constraintsTestSystem.lincslincsExpansionOrder = 8;

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title       = "three atoms, connected to the central atom (e.g. CH3)";
        constraintsTestSystem.numAtoms    = 4;
        constraintsTestSystem.masses      = { 12.0, 1.0, 1.0, 1.0 };
        constraintsTestSystem.constraints = { 0, 0, 1, 0, 0, 2, 0, 0, 3 };
        constraintsTestSystem.constraintsR0 = { 0.1 };


        constraintsTestSystem.x = {
            { 0.00, 0.00, 0.00 }, { 0.10, 0.00, 0.00 }, { 0.00, -0.10, 0.00 }, { 0.00, 0.00, 0.10 }
        };

        constraintsTestSystem.xPrime = { { 0.004, 0.009, -0.010 },
                                         { 0.110, -0.006, 0.003 },
                                         { -0.007, -0.102, -0.007 },
                                         { -0.005, 0.011, 0.102 } };

        constraintsTestSystem.v = { { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title       = "basic triangle (three atoms, connected to each other)";
        constraintsTestSystem.numAtoms    = 3;
        constraintsTestSystem.masses      = { 1.0, 1.0, 1.0 };
        constraintsTestSystem.constraints = { 0, 0, 1, 2, 0, 2, 1, 1, 2 };
        constraintsTestSystem.constraintsR0 = { 0.1, 0.1, 0.1 };

        real oneTenthOverSqrtTwo = 0.1_real / std::sqrt(2.0_real);

        constraintsTestSystem.x = { { oneTenthOverSqrtTwo, 0.0, 0.0 },
                                    { 0.0, oneTenthOverSqrtTwo, 0.0 },
                                    { 0.0, 0.0, oneTenthOverSqrtTwo } };

        constraintsTestSystem.xPrime = { { 0.09, -0.02, 0.01 }, { -0.02, 0.10, -0.02 }, { 0.03, -0.01, 0.07 } };

        constraintsTestSystem.v = { { 1.0, 1.0, 1.0 }, { -2.0, -2.0, -2.0 }, { 1.0, 1.0, 1.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }

    {
        ConstraintsTestSystem singleMolecule;

        singleMolecule.title         = "three atoms, connected longitudinally (e.g. CH2)";
        singleMolecule.numAtoms      = 3;
        singleMolecule.masses        = { 1.0, 12.0, 16.0 };
        singleMolecule.constraints   = { 0, 0, 1, 1, 1, 2 };
        singleMolecule.constraintsR0 = { 0.1, 0.2 };

        real oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);
        real twoTenthsOverSqrtThree = 0.2_real / std::sqrt(3.0_real);

        singleMolecule.x = { { oneTenthOverSqrtTwo, oneTenthOverSqrtTwo, 0.0 },
                             { 0.0, 0.0, 0.0 },
                             { twoTenthsOverSqrtThree, twoTenthsOverSqrtThree, twoTenthsOverSqrtThree } };

        singleMolecule.xPrime = { { 0.08, 0.07, 0.01 }, { -0.02, 0.01, -0.02 }, { 0.10, 0.12, 0.11 } };

        singleMolecule.v = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

        // 150 molecules is 300 constraints, which is larger than the thread block of 256 we are currently using
        const int numMolecules = 150;

        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title    = "system of many molecules";
        constraintsTestSystem.numAtoms = numMolecules * singleMolecule.numAtoms;

        constraintsTestSystem.masses.resize(numMolecules * singleMolecule.numAtoms);
        constraintsTestSystem.x.resize(numMolecules * singleMolecule.numAtoms);
        constraintsTestSystem.xPrime.resize(numMolecules * singleMolecule.numAtoms);
        constraintsTestSystem.v.resize(numMolecules * singleMolecule.numAtoms);

        constraintsTestSystem.constraints.resize(numMolecules * singleMolecule.constraints.size());

        for (int mol = 0; mol < numMolecules; mol++)
        {
            int offsetAtoms = mol * singleMolecule.numAtoms;
            for (int i = 0; i < singleMolecule.numAtoms; i++)
            {
                constraintsTestSystem.masses[offsetAtoms + i] = singleMolecule.masses[i];
                for (int d = 0; d < DIM; d++)
                {
                    constraintsTestSystem.x[offsetAtoms + i][d]      = singleMolecule.x[i][d];
                    constraintsTestSystem.xPrime[offsetAtoms + i][d] = singleMolecule.xPrime[i][d];
                    constraintsTestSystem.v[offsetAtoms + i][d]      = singleMolecule.v[i][d];
                }
            }
            int offsetConstraints = mol * singleMolecule.constraints.size();
            constraintsTestSystem.constraints[offsetConstraints + 0] = singleMolecule.constraints[0];
            constraintsTestSystem.constraints[offsetConstraints + 1] =
                    singleMolecule.constraints[1] + offsetAtoms;
            constraintsTestSystem.constraints[offsetConstraints + 2] =
                    singleMolecule.constraints[2] + offsetAtoms;
            constraintsTestSystem.constraints[offsetConstraints + 3] = singleMolecule.constraints[3];
            constraintsTestSystem.constraints[offsetConstraints + 4] =
                    singleMolecule.constraints[4] + offsetAtoms;
            constraintsTestSystem.constraints[offsetConstraints + 5] =
                    singleMolecule.constraints[5] + offsetAtoms;
        }

        constraintsTestSystem.constraintsR0.resize(singleMolecule.constraintsR0.size());
        for (unsigned long i = 0; i < singleMolecule.constraintsR0.size(); i++)
        {
            constraintsTestSystem.constraintsR0[i] = singleMolecule.constraintsR0[i];
        }

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }

    return constraintsTestSystemList;
}();


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
class ConstraintsTest : public ::testing::TestWithParam<std::tuple<ConstraintsTestSystem, t_pbc>>
{
public:
    //! Reference data
    TestReferenceData refData_;
    //! Checker for reference data
    TestReferenceChecker checker_;

    ConstraintsTest() : checker_(refData_.rootChecker()) {}

    /*! \brief Test if the final position correspond to the reference data.
     *
     * \param[in] testData        Test data structure.
     */
    void checkFinalPositions(const ConstraintsTestData& testData)
    {
        TestReferenceChecker finalPositionsRef(
                checker_.checkSequenceCompound("FinalPositions", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            TestReferenceChecker xPrimeRef(finalPositionsRef.checkCompound("Atom", nullptr));
            const gmx::RVec&     xPrime = testData.xPrime_[i];
            xPrimeRef.checkReal(xPrime[XX], "XX");
            xPrimeRef.checkReal(xPrime[YY], "YY");
            xPrimeRef.checkReal(xPrime[ZZ], "ZZ");
        }
    }

    /*! \brief Test if the final velocities correspond to the reference data.
     *
     * \param[in] testData        Test data structure.
     */
    void checkFinalVelocities(const ConstraintsTestData& testData)
    {
        TestReferenceChecker finalVelocitiesRef(
                checker_.checkSequenceCompound("FinalVelocities", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            TestReferenceChecker vRef(finalVelocitiesRef.checkCompound("Atom", nullptr));
            const gmx::RVec&     v = testData.v_[i];
            vRef.checkReal(v[XX], "XX");
            vRef.checkReal(v[YY], "YY");
            vRef.checkReal(v[ZZ], "ZZ");
        }
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
        for (Index c = 0; c < gmx::ssize(testData.constraints_) / 3; c++)
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
                    d1,
                    r0,
                    c,
                    i,
                    j,
                    d0);
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

        for (Index c = 0; c < gmx::ssize(testData.constraints_) / 3; c++)
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
     * \param[in] testData        Test data structure.
     */
    void checkVirialTensor(const ConstraintsTestData& testData)
    {
        const tensor&        virialScaled = testData.virialScaled_;
        TestReferenceChecker virialScaledRef(checker_.checkCompound("VirialScaled", nullptr));

        virialScaledRef.checkReal(virialScaled[XX][XX], "XX");
        virialScaledRef.checkReal(virialScaled[XX][YY], "XY");
        virialScaledRef.checkReal(virialScaled[XX][ZZ], "XZ");
        virialScaledRef.checkReal(virialScaled[YY][XX], "YX");
        virialScaledRef.checkReal(virialScaled[YY][YY], "YY");
        virialScaledRef.checkReal(virialScaled[YY][ZZ], "YZ");
        virialScaledRef.checkReal(virialScaled[ZZ][XX], "ZX");
        virialScaledRef.checkReal(virialScaled[ZZ][YY], "ZY");
        virialScaledRef.checkReal(virialScaled[ZZ][ZZ], "ZZ");
    }

    //! Before any test is run, work out whether any compatible GPUs exist.
    static std::vector<std::unique_ptr<IConstraintsTestRunner>> getRunners()
    {
        std::vector<std::unique_ptr<IConstraintsTestRunner>> runners;
        // Add runners for CPU versions of SHAKE and LINCS
        runners.emplace_back(std::make_unique<ShakeConstraintsRunner>());
        runners.emplace_back(std::make_unique<LincsConstraintsRunner>());
        // If supported, add runners for the GPU version of LINCS for each available GPU
        const bool addGpuRunners = GPU_CONSTRAINTS_SUPPORTED;
        if (addGpuRunners)
        {
            for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
            {
                runners.emplace_back(std::make_unique<LincsDeviceConstraintsRunner>(*testDevice));
            }
        }
        return runners;
    }
};

TEST_P(ConstraintsTest, SatisfiesConstraints)
{
    auto                  params                = GetParam();
    ConstraintsTestSystem constraintsTestSystem = std::get<0>(params);
    t_pbc                 pbc                   = std::get<1>(params);

    ConstraintsTestData testData(constraintsTestSystem.title,
                                 constraintsTestSystem.numAtoms,
                                 constraintsTestSystem.masses,
                                 constraintsTestSystem.constraints,
                                 constraintsTestSystem.constraintsR0,
                                 true,
                                 false,
                                 real(0.0),
                                 real(0.001),
                                 constraintsTestSystem.x,
                                 constraintsTestSystem.xPrime,
                                 constraintsTestSystem.v,
                                 constraintsTestSystem.shakeTolerance,
                                 constraintsTestSystem.shakeUseSOR,
                                 constraintsTestSystem.lincsNIter,
                                 constraintsTestSystem.lincslincsExpansionOrder,
                                 constraintsTestSystem.lincsWarnAngle);

    FloatingPointTolerance positionsTolerance = absoluteTolerance(0.001F);
    FloatingPointTolerance velocityTolerance  = absoluteTolerance(0.02F);
    FloatingPointTolerance lengthTolerance    = relativeToleranceAsFloatingPoint(0.1, 0.002F);

    // Cycle through all available runners
    for (const auto& runner : getRunners())
    {
        SCOPED_TRACE(formatString("Testing %s with %s PBC using %s.",
                                  testData.title_.c_str(),
                                  c_pbcTypeNames[pbc.pbcType].c_str(),
                                  runner->name().c_str()));

        testData.reset();

        // Apply constraints
        runner->applyConstraints(&testData, pbc);


        checker_.setDefaultTolerance(positionsTolerance);
        checkFinalPositions(testData);
        checker_.setDefaultTolerance(velocityTolerance);
        checkFinalVelocities(testData);

        checkConstrainsLength(lengthTolerance, testData, pbc);
        checkConstrainsDirection(testData, pbc);
        checkCOMCoordinates(positionsTolerance, testData);
        checkCOMVelocity(velocityTolerance, testData);

        float virialTrace = 0.0F;
        for (int d = 0; d < DIM; d++)
        {
            virialTrace += testData.virialScaled_[d][d];
        }

        // The virial tolerance factor can be:
        // LINCS iter=2, order=4:   0.002
        // LINCS iter=1, order=4:   0.02
        // SHAKE tolerance=0.0001:  0.2
        // SHAKE tolerance=0.00002: 0.1
        const float virialRelativeTolerance = (runner->name().substr(0, 5) == "SHAKE" ? 0.1F : 0.02F);
        FloatingPointTolerance virialTolerance =
                absoluteTolerance(std::fabs(virialTrace) / 3 * virialRelativeTolerance);

        checker_.setDefaultTolerance(virialTolerance);
        checkVirialTensor(testData);
    }
}

INSTANTIATE_TEST_SUITE_P(WithParameters,
                         ConstraintsTest,
                         ::testing::Combine(::testing::ValuesIn(c_constraintsTestSystemList),
                                            ::testing::ValuesIn(c_pbcs)));

} // namespace
} // namespace test
} // namespace gmx
