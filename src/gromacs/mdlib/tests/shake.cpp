/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/mdlib/constr.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

/*! \brief Stride of the vector of int used to describe each SHAKE
 * constraint
 *
 * Like other such code, SHAKE is hard-wired to use t_ilist.iatoms as
 * a flat vector of tuples of general data. Here, they are triples
 * containing the index of the constraint type, and then the indices
 * of the two atoms involved. So for each constraint, we must stride
 * this vector by three to get access to its information. */
const int constraintStride = 3;

/*! \brief Compute the displacements between pairs of constrained
 * atoms described in the iatom "topology". */
std::vector<real>
computeDisplacements(const std::vector<int>     &iatom,
                     const std::vector<real>    &positions)
{
    assert(0 == iatom.size() % constraintStride);
    int               numConstraints = iatom.size() / constraintStride;
    std::vector<real> displacements;

    for (int ll = 0; ll != numConstraints; ++ll)
    {
        int atom_i = iatom[ll*constraintStride + 1];
        int atom_j = iatom[ll*constraintStride + 2];

        for (int d = 0; d != DIM; d++)
        {
            displacements.push_back(positions[atom_i*DIM + d] - positions[atom_j*DIM + d]);
        }
    }

    return displacements;
}

/*! \brief Compute half of the reduced mass of each pair of constrained
 * atoms in the iatom "topology".
 *
 * The reduced mass is m = 1/(1/m_i + 1/m_j)) */
std::vector<real>
computeHalfOfReducedMasses(const std::vector<int>     &iatom,
                           const std::vector<real>    &inverseMasses)
{
    int               numConstraints = iatom.size() / constraintStride;
    std::vector<real> halfOfReducedMasses;

    for (int ll = 0; ll != numConstraints; ++ll)
    {
        int atom_i = iatom[ll*constraintStride + 1];
        int atom_j = iatom[ll*constraintStride + 2];

        halfOfReducedMasses.push_back(0.5 / (inverseMasses[atom_i] + inverseMasses[atom_j]));
    }

    return halfOfReducedMasses;
}

/*! \brief Compute the distances corresponding to the vector of displacements components */
std::vector<real>
computeDistancesSquared(const std::vector<real> &displacements)
{
    assert(0 == displacements.size() % DIM);
    int               numDistancesSquared = displacements.size() / DIM;
    std::vector<real> distanceSquared;

    for (int i = 0; i != numDistancesSquared; ++i)
    {
        distanceSquared.push_back(0.0);
        for (int d = 0; d != DIM; ++d)
        {
            real displacement = displacements[i*DIM + d];
            distanceSquared.back() += displacement * displacement;
        }
    }

    return distanceSquared;
}

/*! \brief Test fixture for testing SHAKE */
class ShakeTest : public ::testing::Test
{
    public:
        /*! \brief Set up data for test cases to use when constructing
            their input */
        void SetUp()
        {
            inverseMassesDatabase_.push_back(2.0);
            inverseMassesDatabase_.push_back(3.0);
            inverseMassesDatabase_.push_back(4.0);
            inverseMassesDatabase_.push_back(1.0);

            positionsDatabase_.push_back(2.5);
            positionsDatabase_.push_back(-3.1);
            positionsDatabase_.push_back(15.7);

            positionsDatabase_.push_back(0.51);
            positionsDatabase_.push_back(-3.02);
            positionsDatabase_.push_back(15.55);

            positionsDatabase_.push_back(-0.5);
            positionsDatabase_.push_back(-3.0);
            positionsDatabase_.push_back(15.2);

            positionsDatabase_.push_back(-1.51);
            positionsDatabase_.push_back(-2.95);
            positionsDatabase_.push_back(15.05);
        }

        //! Run the test
        void runTest(size_t gmx_unused           numAtoms,
                     size_t                      numConstraints,
                     const std::vector<int>     &iatom,
                     const std::vector<real>    &constrainedDistances,
                     const std::vector<real>    &inverseMasses,
                     const std::vector<real>    &positions)
        {
            // Check the test input is consistent
            assert(numConstraints * constraintStride == iatom.size());
            assert(numConstraints == constrainedDistances.size());
            assert(numAtoms == inverseMasses.size());
            assert(numAtoms * DIM == positions.size());
            for (size_t i = 0; i != numConstraints; ++i)
            {
                for (size_t j = 1; j < 3; j++)
                {
                    // Check that the topology refers to atoms that have masses and positions
                    assert(iatom[i*constraintStride + j] >= 0);
                    assert(iatom[i*constraintStride + j] < static_cast<int>(numAtoms));
                }
            }
            std::vector<real> distanceSquaredTolerances;
            std::vector<real> lagrangianValues;
            std::vector<real> constrainedDistancesSquared;

            real              coordMax = 0;
            for (size_t i = 0; i != numConstraints; ++i)
            {
                constrainedDistancesSquared.push_back(constrainedDistances[i] * constrainedDistances[i]);
                distanceSquaredTolerances.push_back(0.5 / (constrainedDistancesSquared.back() * ShakeTest::tolerance_));
                lagrangianValues.push_back(0.0);

                for (size_t j = 1; j < constraintStride; j++)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        coordMax = std::max(coordMax, std::abs(positions[iatom[i*constraintStride + j]*DIM + d]));
                    }
                }
            }
            std::vector<real> halfOfReducedMasses  = computeHalfOfReducedMasses(iatom, inverseMasses);
            std::vector<real> initialDisplacements = computeDisplacements(iatom, positions);

            std::vector<real> finalPositions = positions;
            int               numIterations  = 0;
            int               numErrors      = 0;

            cshake(iatom.data(), numConstraints, &numIterations,
                   ShakeTest::maxNumIterations_, constrainedDistancesSquared.data(),
                   finalPositions.data(), initialDisplacements.data(),
                   halfOfReducedMasses.data(), omega_, inverseMasses.data(),
                   distanceSquaredTolerances.data(),
                   lagrangianValues.data(),
                   &numErrors);

            std::vector<real> finalDisplacements    = computeDisplacements(iatom, finalPositions);
            std::vector<real> finalDistancesSquared = computeDistancesSquared(finalDisplacements);
            assert(numConstraints == finalDistancesSquared.size());

            EXPECT_EQ(0, numErrors);
            EXPECT_GT(numIterations, 1);
            EXPECT_LT(numIterations, ShakeTest::maxNumIterations_);
            // TODO wrap this in a Google Mock matcher if there's
            // other tests like it some time?
            for (size_t i = 0; i != numConstraints; ++i)
            {
                // We need to allow for the requested tolerance plus rounding
                // errors due to the absolute size of the coordinate values
                gmx::test::FloatingPointTolerance constraintTolerance =
                    gmx::test::absoluteTolerance(std::sqrt(constrainedDistancesSquared[i])*ShakeTest::tolerance_ + coordMax*GMX_REAL_EPS);
                // Assert that the constrained distances are within the required tolerance
                EXPECT_FLOAT_EQ_TOL(std::sqrt(constrainedDistancesSquared[i]),
                                    std::sqrt(finalDistancesSquared[i]),
                                    constraintTolerance);
            }
        }

        //! Tolerance for SHAKE conversion (ie. shake-tol .mdp setting)
        static const real tolerance_;
        //! Maximum number of iterations permitted in these tests
        static const int  maxNumIterations_;
        //! SHAKE over-relaxation (SOR) factor
        static const real omega_;
        //! Database of inverse masses of atoms in the topology
        std::vector<real> inverseMassesDatabase_;
        //! Database of atom positions (three reals per atom)
        std::vector<real> positionsDatabase_;
};

const real ShakeTest::tolerance_        = 1e-5;
const int  ShakeTest::maxNumIterations_ = 30;
const real ShakeTest::omega_            = 1.0;

TEST_F(ShakeTest, ConstrainsOneBond)
{
    int                  numAtoms       = 2;
    int                  numConstraints = 1;

    std::vector<int>     iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0);  // i atom index
    iatom.push_back(1);  // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

TEST_F(ShakeTest, ConstrainsTwoDisjointBonds)
{
    int                  numAtoms       = 4;
    int                  numConstraints = 2;

    std::vector<int>     iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0);  // i atom index
    iatom.push_back(1);  // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(2);  // i atom index
    iatom.push_back(3);  // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);
    constrainedDistances.push_back(1.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

TEST_F(ShakeTest, ConstrainsTwoBondsWithACommonAtom)
{
    int                  numAtoms       = 3;
    int                  numConstraints = 2;

    std::vector<int>     iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0);  // i atom index
    iatom.push_back(1);  // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(1);  // i atom index
    iatom.push_back(2);  // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);
    constrainedDistances.push_back(1.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

TEST_F(ShakeTest, ConstrainsThreeBondsWithCommonAtoms)
{
    int                  numAtoms       = 4;
    int                  numConstraints = 3;

    std::vector<int>     iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0);  // i atom index
    iatom.push_back(1);  // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(1);  // i atom index
    iatom.push_back(2);  // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(2);  // i atom index
    iatom.push_back(3);  // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);
    constrainedDistances.push_back(1.0);
    constrainedDistances.push_back(1.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

} // namespace
