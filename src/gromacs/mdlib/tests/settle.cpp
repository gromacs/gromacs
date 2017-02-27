/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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

#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

//! Database of 51 water atom input positions (DIM reals per atom, taken from spc216.gro) for use as test inputs.
const double g_positions[] = {
    .130, -.041, -.291,
    .120, -.056, -.192,
    .044, -.005, -.327,
    -.854, -.406, .477,
    -.900, -.334, .425,
    -.858, -.386, .575,
    .351, -.061, .853,
    .401, -.147, .859,
    .416, .016, .850,
    -.067, -.796, .873,
    -.129, -.811, .797,
    -.119, -.785, .958,
    -.635, -.312, -.356,
    -.629, -.389, -.292,
    -.687, -.338, -.436,
    .321, -.919, .242,
    .403, -.880, .200,
    .294, -1.001, .193,
    -.404, .735, .728,
    -.409, .670, .803,
    -.324, .794, .741,
    .461, -.596, -.135,
    .411, -.595, -.221,
    .398, -.614, -.059,
    -.751, -.086, .237,
    -.811, -.148, .287,
    -.720, -.130, .152,
    .202, .285, -.364,
    .122, .345, -.377,
    .192, .236, -.278,
    -.230, -.485, .081,
    -.262, -.391, .071,
    -.306, -.548, .069,
    .464, -.119, .323,
    .497, -.080, .409,
    .540, -.126, .258,
    -.462, .107, .426,
    -.486, .070, .336,
    -.363, .123, .430,
    .249, -.077, -.621,
    .306, -.142, -.571,
    .233, -.110, -.714,
    -.922, -.164, .904,
    -.842, -.221, .925,
    -.971, -.204, .827,
    .382, .700, .480,
    .427, .610, .477,
    .288, .689, .513,
    .781, .264, -.113,
    .848, .203, -.070,
    .708, .283, -.048
};

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
        std::vector<real> updatedPositions_;
        //! Water atom velocities to constrain (DIM reals per atom)
        std::vector<real> velocities_;
        //! PBC option to test
        t_pbc             pbcNone_;
        //! PBC option to test
        t_pbc             pbcXYZ_;

        //! Constructor
        SettleTest() :
            updatedPositions_(std::begin(g_positions), std::end(g_positions)),
            velocities_(updatedPositions_.size(), 0)
        {
            set_pbc(&pbcNone_, epbcNONE, g_box);
            set_pbc(&pbcXYZ_, epbcXYZ, g_box);

            // Perturb the atom positions, to appear like an
            // "update," and where there is definitely constraining
            // work to do.
            for (size_t i = 0; i != updatedPositions_.size(); ++i)
            {
                if (i % 4 == 0)
                {
                    updatedPositions_[i] += 0.01;
                }
                else if (i % 4 == 1)
                {
                    updatedPositions_[i] -= 0.01;
                }
                else if (i % 4 == 2)
                {
                    updatedPositions_[i] += 0.02;
                }
                else if (i % 4 == 3)
                {
                    updatedPositions_[i] -= 0.02;
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
    ASSERT_LE(numSettles, updatedPositions_.size() / (atomsPerSettle * DIM)) << "cannot test that many SETTLEs " << testDescription;

    // Set up the topology. We still have to make some raw pointers,
    // but they are put into scope guards for automatic cleanup.
    gmx_mtop_t                   *mtop;
    snew(mtop, 1);
    const unique_cptr<gmx_mtop_t> mtopGuard(mtop);
    mtop->mols.nr  = 1;
    mtop->nmoltype = 1;
    snew(mtop->moltype, mtop->nmoltype);
    const unique_cptr<gmx_moltype_t> moltypeGuard(mtop->moltype);
    mtop->nmolblock = 1;
    snew(mtop->molblock, mtop->nmolblock);
    const unique_cptr<gmx_molblock_t> molblockGuard(mtop->molblock);
    mtop->molblock[0].type = 0;
    std::vector<int>                  iatoms;
    for (int i = 0; i < numSettles; ++i)
    {
        iatoms.push_back(settleType);
        iatoms.push_back(i*atomsPerSettle+0);
        iatoms.push_back(i*atomsPerSettle+1);
        iatoms.push_back(i*atomsPerSettle+2);
    }
    mtop->moltype[0].ilist[F_SETTLE].iatoms = iatoms.data();
    mtop->moltype[0].ilist[F_SETTLE].nr     = iatoms.size();

    // Set up the SETTLE parameters.
    mtop->ffparams.ntypes = 1;
    snew(mtop->ffparams.iparams, mtop->ffparams.ntypes);
    const unique_cptr<t_iparams> iparamsGuard(mtop->ffparams.iparams);
    const real                   dOH = 0.09572;
    const real                   dHH = 0.15139;
    mtop->ffparams.iparams[settleType].settle.doh = dOH;
    mtop->ffparams.iparams[settleType].settle.dhh = dHH;

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
    gmx_settledata_t settled = settle_init(mtop);
    settle_set_constraints(settled, &mtop->moltype[0].ilist[F_SETTLE], &mdatoms);

    // Copy the original positions from the array of doubles to a vector of reals
    std::vector<real> startingPositions(std::begin(g_positions), std::end(g_positions));

    // Run the test
    bool       errorOccured;
    tensor     virial             = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    int        numThreads         = 1, threadIndex = 0;
    const real reciprocalTimeStep = 1.0/0.002;
    csettle(settled, numThreads, threadIndex,
            usePbc ? &pbcXYZ_ : &pbcNone_,
            startingPositions.data(), updatedPositions_.data(), reciprocalTimeStep,
            useVelocities ? velocities_.data() : nullptr,
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
    int positionIndex = 0, velocityIndex = 0;
    for (int i = 0; i < numSettles; ++i)
    {
        rvec positionO  {
            updatedPositions_[positionIndex++], updatedPositions_[positionIndex++], updatedPositions_[positionIndex++]
        };
        rvec positionH1 {
            updatedPositions_[positionIndex++], updatedPositions_[positionIndex++], updatedPositions_[positionIndex++]
        };
        rvec positionH2 {
            updatedPositions_[positionIndex++], updatedPositions_[positionIndex++], updatedPositions_[positionIndex++]
        };

        EXPECT_REAL_EQ_TOL(dOH*dOH, distance2(positionO, positionH1), tolerance) << formatString("for water %d ", i) << testDescription;
        EXPECT_REAL_EQ_TOL(dOH*dOH, distance2(positionO, positionH2), tolerance) << formatString("for water %d ", i) << testDescription;
        EXPECT_REAL_EQ_TOL(dHH*dHH, distance2(positionH1, positionH2), tolerance) << formatString("for water %d ", i) << testDescription;

        // This merely tests whether the velocities were
        // updated from the starting values of zero (or not),
        // but not whether the update was correct.
        for (int j = 0; j < atomsPerSettle * DIM; ++j, ++velocityIndex)
        {
            EXPECT_TRUE(useVelocities == (0. != velocities_[velocityIndex])) << formatString("for water %d velocity coordinate %d ", i, j) << testDescription;
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

} // namespace
} // namespace
