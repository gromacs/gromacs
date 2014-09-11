#include "gmxpre.h"

#include <assert.h>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/types/simple.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

// TODO clean up the fact that the SHAKE code is hard-wired to abuse
// t_ilist.iatoms to have triplets of the index of the constraint
// type, and then the indices of the two atoms involved
const int constraintStride = 3;

/*! \brief Compute the displacements between pairs of constrained
 * atoms described in the iatom "topology". */
std::vector<real>
computeDisplacements(const std::vector<atom_id> &iatom,
                     const std::vector<real> &positions)
{
    assert(0 == iatom.size() % constraintStride);
    int numConstraints = iatom.size() / constraintStride;
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
computeHalfOfReducedMasses(const std::vector<atom_id> &iatom,
                           const std::vector<real> &inverseMasses)
{
    int numConstraints = iatom.size() / constraintStride;
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
    int numDistancesSquared = displacements.size() / DIM;
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
class Shake : public ::testing::Test
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
        void runTest(size_t numAtoms,
                     size_t numConstraints,
                     const std::vector<atom_id> &iatom,
                     const std::vector<real> &constrainedDistances,
                     const std::vector<real> &inverseMasses,
                     const std::vector<real> &positions)
        {
            assert(numConstraints * constraintStride == iatom.size());
            assert(numConstraints == constrainedDistances.size());
            assert(numAtoms == inverseMasses.size());
            assert(numAtoms * DIM == positions.size());

            std::vector<real> distanceSquaredTolerances;
            std::vector<real> lagrangianValues;
            std::vector<real> constrainedDistancesSquared;

            for (size_t i = 0; i != numConstraints; ++i)
            {
                constrainedDistancesSquared.push_back(constrainedDistances[i] * constrainedDistances[i]);
                distanceSquaredTolerances.push_back(1.0 / (constrainedDistancesSquared.back() * Shake::tolerance_));
                lagrangianValues.push_back(0.0);
            }
            std::vector<real> halfOfReducedMasses = computeHalfOfReducedMasses(iatom, inverseMasses);
            std::vector<real> initialDisplacements = computeDisplacements(iatom, positions);

            std::vector<real> finalPositions = positions;
            int numIterations = 0;
            int numErrors = 0;

            cshake(&iatom[0], numConstraints, &numIterations,
                   Shake::maxNumIterations_, &constrainedDistancesSquared[0],
                   &finalPositions[0], &initialDisplacements[0], &halfOfReducedMasses[0],
                   omega_, &inverseMasses[0],
                   &distanceSquaredTolerances[0],
                   &lagrangianValues[0],
                   &numErrors);

            std::vector<real> finalDisplacements = computeDisplacements(iatom, finalPositions);
            std::vector<real> finalDistancesSquared = computeDistancesSquared(finalDisplacements);
            assert(numConstraints == finalDistancesSquared.size());

            EXPECT_EQ(0, numErrors);
            EXPECT_GT(numIterations, 1);
            EXPECT_LT(numIterations, Shake::maxNumIterations_);
            // TODO wrap this in a Google Mock matcher if there's
            // other tests like it some time?
            for (size_t i = 0; i != numConstraints; ++i)
            {
                gmx::test::FloatingPointTolerance constraintTolerance = 
                    gmx::test::relativeToleranceAsFloatingPoint(constrainedDistancesSquared[i],
                                                                Shake::tolerance_);
                // Assert that the constrained distances are within the required tolerance
                EXPECT_FLOAT_EQ_TOL(constrainedDistancesSquared[i],
                                    finalDistancesSquared[i],
                                    constraintTolerance);
            }
        }

        //! Tolerance for SHAKE conversion (ie. shake-tol .mdp setting)
        static const real tolerance_;
        //! Maximum number of iterations permitted in these tests
        static const int maxNumIterations_;
        //! SHAKE over-relaxation (SOR) factor
        static const real omega_;
        //! Database of inverse masses of atoms in the topology
        std::vector<real> inverseMassesDatabase_;
        //! Database of atom positions (three reals per atom)
        std::vector<real> positionsDatabase_;
};

const real Shake::tolerance_ = 1e-5;
const int Shake::maxNumIterations_ = 30;
const real Shake::omega_ = 1.0;

TEST_F(Shake, ConstrainsOneBond)
{
    int numAtoms = 2;
    int numConstraints = 1;

    std::vector<atom_id> iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0); // i atom index
    iatom.push_back(1); // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

TEST_F(Shake, ConstrainsTwoDisjointBonds)
{
    int numAtoms = 4;
    int numConstraints = 2;

    std::vector<atom_id> iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0); // i atom index
    iatom.push_back(1); // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(2); // i atom index
    iatom.push_back(3); // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);
    constrainedDistances.push_back(1.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

TEST_F(Shake, ConstrainsTwoBondsWithACommonAtom)
{
    int numAtoms = 3;
    int numConstraints = 2;

    std::vector<atom_id> iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0); // i atom index
    iatom.push_back(1); // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(1); // i atom index
    iatom.push_back(2); // j atom index

    std::vector<real> constrainedDistances;
    constrainedDistances.push_back(2.0);
    constrainedDistances.push_back(1.0);

    std::vector<real> inverseMasses(inverseMassesDatabase_.begin(),
                                    inverseMassesDatabase_.begin() + numAtoms);
    std::vector<real> positions(positionsDatabase_.begin(),
                                    positionsDatabase_.begin() + numAtoms * DIM);

    runTest(numAtoms, numConstraints, iatom, constrainedDistances, inverseMasses, positions);
}

TEST_F(Shake, ConstrainsThreeBondsWithCommonAtoms)
{
    int numAtoms = 4;
    int numConstraints = 3;

    std::vector<atom_id> iatom;
    iatom.push_back(-1); // unused
    iatom.push_back(0); // i atom index
    iatom.push_back(1); // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(1); // i atom index
    iatom.push_back(2); // j atom index

    iatom.push_back(-1); // unused
    iatom.push_back(2); // i atom index
    iatom.push_back(3); // j atom index

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

