/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/mdlib/settle.h"

#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "gromacs/mdlib/tests/watersystem.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! Simple cubic simulation box to use in tests
matrix g_box = {{real(1.86206), 0, 0}, {0, real(1.86206), 0}, {0, 0, real(1.86206)}};

//! Convenience typedef
typedef std::tuple<int, bool, bool, bool> SettleTestParameters;

/*! \brief Test fixture for testing SETTLE position updates
 *
 * \todo This also tests that if the calling code requires velocities
 * and virial updates, that those outputs do change, but does not test
 * that those changes are correct.
 *
 * \todo Only no-PBC and cubic-PBC are tested here, but the correct
 * function of the SIMD version of set_pbx_auic in all cases should be
 * tested elsewhere.
 */
class SettleTest : public ::testing::TestWithParam<SettleTestParameters>
{
    public:
        //! Updated water atom positions to constrain (DIM reals per atom)
        std::vector<gmx::RVec> updatedPositions_;
        //! Water atom velocities to constrain (DIM reals per atom)
        std::vector<gmx::RVec> velocities_;
        //! PBC option to test: none
        t_pbc                  pbcNone_;
        //! PBC option to test: xyz
        t_pbc                  pbcXYZ_;

        SettleTest() :
            updatedPositions_(std::begin(c_waterPositions), std::end(c_waterPositions)),
            velocities_(updatedPositions_.size(), { 0, 0, 0 })
        {
            set_pbc(&pbcNone_, epbcNONE, g_box);
            set_pbc(&pbcXYZ_, epbcXYZ, g_box);

            // Perturb the atom positions, to appear like an
            // "update," and where there is definitely constraining
            // work to do.
            for (size_t i = 0; i != updatedPositions_.size()*DIM; ++i)
            {
                if (i % 4 == 0)
                {
                    updatedPositions_[i / 3][i % 3] += 0.01;
                }
                else if (i % 4 == 1)
                {
                    updatedPositions_[i / 3][i % 3] -= 0.01;
                }
                else if (i % 4 == 2)
                {
                    updatedPositions_[i / 3][i % 3] += 0.02;
                }
                else if (i % 4 == 3)
                {
                    updatedPositions_[i / 3][i % 3] -= 0.02;
                }
            }
        }
};

TEST_P(SettleTest, SatisfiesConstraints)
{
    int  numSettles;
    bool usePbc, useVelocities, calcVirial;
    // Make some symbolic names for the parameter combination under
    // test.
    std::tie(numSettles, usePbc, useVelocities, calcVirial) = GetParam();

    // Make a string that describes which parameter combination is
    // being tested, to help make failing tests comprehensible.
    std::string testDescription = formatString("while testing %d SETTLEs, %sPBC, %svelocities and %scalculating the virial",
                                               numSettles,
                                               usePbc ? "with " : "without ",
                                               useVelocities ? "with " : "without ",
                                               calcVirial ? "" : "not ");

    const int settleType     = 0;
    const int atomsPerSettle = NRAL(F_SETTLE);
    ASSERT_LE(numSettles, updatedPositions_.size() / atomsPerSettle) << "cannot test that many SETTLEs " << testDescription;

    // Set up the topology.
    gmx_mtop_t mtop;
    mtop.moltype.resize(1);
    mtop.molblock.resize(1);
    mtop.molblock[0].type = 0;
    std::vector<int> &iatoms = mtop.moltype[0].ilist[F_SETTLE].iatoms;
    for (int i = 0; i < numSettles; ++i)
    {
        iatoms.push_back(settleType);
        iatoms.push_back(i*atomsPerSettle + 0);
        iatoms.push_back(i*atomsPerSettle + 1);
        iatoms.push_back(i*atomsPerSettle + 2);
    }

    // Set up the SETTLE parameters.
    const real     dOH = 0.09572;
    const real     dHH = 0.15139;
    t_iparams      iparams;
    iparams.settle.doh = dOH;
    iparams.settle.dhh = dHH;
    mtop.ffparams.iparams.push_back(iparams);

    // Set up the masses.
    t_mdatoms         mdatoms;
    std::vector<real> mass, massReciprocal;
    const real        oxygenMass = 15.9994, hydrogenMass = 1.008;
    for (int i = 0; i < numSettles; ++i)
    {
        mass.push_back(oxygenMass);
        mass.push_back(hydrogenMass);
        mass.push_back(hydrogenMass);
        massReciprocal.push_back(1./oxygenMass);
        massReciprocal.push_back(1./hydrogenMass);
        massReciprocal.push_back(1./hydrogenMass);
    }
    mdatoms.massT   = mass.data();
    mdatoms.invmass = massReciprocal.data();
    mdatoms.homenr  = numSettles * atomsPerSettle;

    // Finally make the settle data structures
    settledata    *settled = settle_init(mtop);
    const t_ilist  ilist   = { mtop.moltype[0].ilist[F_SETTLE].size(), 0, mtop.moltype[0].ilist[F_SETTLE].iatoms.data(), 0 };
    settle_set_constraints(settled, &ilist, mdatoms);

    // Copy the original positions from the array of doubles to a vector of reals
    std::vector<gmx::RVec> startingPositions(std::begin(c_waterPositions), std::end(c_waterPositions));

    // Run the test
    bool       errorOccured;
    tensor     virial             = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    int        numThreads         = 1, threadIndex = 0;
    const real reciprocalTimeStep = 1.0/0.002;
    csettle(settled, numThreads, threadIndex,
            usePbc ? &pbcXYZ_ : &pbcNone_,
            static_cast<real *>(startingPositions[0]), static_cast<real *>(updatedPositions_[0]), reciprocalTimeStep,
            useVelocities ? static_cast<real *>(velocities_[0]) : nullptr,
            calcVirial, virial, &errorOccured);
    settle_free(settled);
    EXPECT_FALSE(errorOccured) << testDescription;

    // The necessary tolerances for the test to pass were determined
    // empirically. This isn't nice, but the required behaviour that
    // SETTLE produces constrained coordinates consistent with
    // sensible sampling needs to be tested at a much higher level.
    FloatingPointTolerance tolerance =
        relativeToleranceAsPrecisionDependentUlp(dOH*dOH, 80, 380);

    // Verify the updated coordinates match the requirements
    for (int i = 0; i < numSettles; ++i)
    {
        const gmx::RVec &positionO  = updatedPositions_[i*3 + 0];
        const gmx::RVec &positionH1 = updatedPositions_[i*3 + 1];
        const gmx::RVec &positionH2 = updatedPositions_[i*3 + 2];

        EXPECT_REAL_EQ_TOL(dOH*dOH, distance2(positionO, positionH1), tolerance) << formatString("for water %d ", i) << testDescription;
        EXPECT_REAL_EQ_TOL(dOH*dOH, distance2(positionO, positionH2), tolerance) << formatString("for water %d ", i) << testDescription;
        EXPECT_REAL_EQ_TOL(dHH*dHH, distance2(positionH1, positionH2), tolerance) << formatString("for water %d ", i) << testDescription;

        // This merely tests whether the velocities were
        // updated from the starting values of zero (or not),
        // but not whether the update was correct.
        for (int j = 0; j < atomsPerSettle; ++j)
        {
            for (int d = 0; d < DIM; ++d)
            {
                EXPECT_TRUE(useVelocities == (0. != velocities_[i*3 + j][d])) << formatString("for water %d velocity atom %d dim %d", i, j, d) << testDescription;
            }
        }
    }

    // This merely tests whether the viral was updated from
    // the starting values of zero (or not), but not whether
    // any update was correct.
    for (int d = 0; d < DIM; ++d)
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            EXPECT_TRUE(calcVirial == (0. != virial[d][dd])) << formatString("for virial component[%d][%d] ", d, dd) << testDescription;
        }
    }
}

// Scan the full Cartesian product of numbers of SETTLE interactions
// (4 and 17 are chosen to test cases that do and do not match
// hardware SIMD widths), and whether or not we use PBC, velocities or
// calculate the virial contribution.
INSTANTIATE_TEST_CASE_P(WithParameters, SettleTest,
                            ::testing::Combine(::testing::Values(1, 4, 7),
                                                   ::testing::Bool(),
                                                   ::testing::Bool(),
                                                   ::testing::Bool()));

}  // namespace test
}  // namespace gmx
