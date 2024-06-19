/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "nblib/listed_forces/calculator.h"

#include <cstddef>

#include <array>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <valarray>
#include <vector>

#include <gtest/gtest.h>

#include "listed_forces/dataflow.hpp"
#include "listed_forces/traits.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "nblib/basicdefinitions.h"
#include "nblib/box.h"
#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/util/traits.hpp"
#include "nblib/vector.h"

#include "linear_chain_input.hpp"
#include "pbc.hpp"
#include "testhelpers.h"


namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, ListedForceCalculatorCanConstruct)
{
    ListedInteractionData interactions;
    Box                   box(1, 1, 1);
    EXPECT_NO_THROW(ListedForceCalculator listedForceCalculator(interactions, 2, 1, box));
}

template<class TestSeq, class SeqFloat, class SeqDouble>
void compareVectors(const TestSeq&                    forces,
                    [[maybe_unused]] const SeqFloat&  refForcesFloat,
                    [[maybe_unused]] const SeqDouble& refForcesDouble)
{
    for (size_t i = 0; i < forces.size(); ++i)
    {
        for (int m = 0; m < dimSize; ++m)
        {
            EXPECT_FLOAT_DOUBLE_EQ_TOL(
                    forces[i][m],
                    refForcesFloat[i][m],
                    refForcesDouble[i][m],
                    // Todo: why does the tolerance need to be so low?
                    gmx::test::relativeToleranceAsFloatingPoint(refForcesDouble[i][m], 5e-5));
        }
    }
}

class ListedExampleData : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // methanol-spc data
        HarmonicBondType              bond1{ 376560, 0.136 };
        HarmonicBondType              bond2{ 313800, 0.1 };
        std::vector<HarmonicBondType> bonds{ bond1, bond2 };
        // one bond between atoms 0-1 with bond1 parameters and another between atoms 1-2 with bond2 parameters
        std::vector<InteractionIndex<HarmonicBondType>> bondIndices{ { 0, 1, 0 }, { 1, 2, 1 } };

        HarmonicAngle                                angle(397.5, Degrees(108.53));
        std::vector<HarmonicAngle>                   angles{ angle };
        std::vector<InteractionIndex<HarmonicAngle>> angleIndices{ { 0, 1, 2, 0 } };

        pickType<HarmonicBondType>(interactions).indices    = bondIndices;
        pickType<HarmonicBondType>(interactions).parameters = bonds;

        pickType<HarmonicAngle>(interactions).indices    = angleIndices;
        pickType<HarmonicAngle>(interactions).parameters = angles;

        // initial position for the methanol atoms from the spc-water example
        x = std::vector<gmx::RVec>{ { 1.97, 1.46, 1.209 }, { 1.978, 1.415, 1.082 }, { 1.905, 1.46, 1.03 } };
        forces = std::vector<gmx::RVec>(3, gmx::RVec{ 0, 0, 0 });

        box.reset(new Box(3, 3, 3));
        pbc.reset(new PbcHolder(PbcType::Xyz, *box));
    }

    std::vector<gmx::RVec> x;
    std::vector<gmx::RVec> forces;

    ListedInteractionData interactions;

    std::shared_ptr<Box>       box;
    std::shared_ptr<PbcHolder> pbc;
};

TEST_F(ListedExampleData, ComputeHarmonicBondForces)
{
    gmx::ArrayRef<const InteractionIndex<HarmonicBondType>> indices =
            pickType<HarmonicBondType>(interactions).indices;
    gmx::ArrayRef<const HarmonicBondType> bonds = pickType<HarmonicBondType>(interactions).parameters;
    computeForces(indices, bonds, x, &forces, *pbc);

    RefDataChecker vector3DTest(1e-3);
    vector3DTest.testArrays<Vec3>(forces, "Bond forces");
}

TEST_F(ListedExampleData, ComputeHarmonicBondEnergies)
{
    gmx::ArrayRef<const InteractionIndex<HarmonicBondType>> indices =
            pickType<HarmonicBondType>(interactions).indices;
    gmx::ArrayRef<const HarmonicBondType> bonds = pickType<HarmonicBondType>(interactions).parameters;
    real                                  energy = computeForces(indices, bonds, x, &forces, *pbc);

    RefDataChecker vector3DTest(1e-4);
    vector3DTest.testReal(energy, "Bond energy");
}

TEST_F(ListedExampleData, ComputeHarmonicAngleForces)
{
    gmx::ArrayRef<const InteractionIndex<HarmonicAngle>> indices =
            pickType<HarmonicAngle>(interactions).indices;
    gmx::ArrayRef<const HarmonicAngle> angles = pickType<HarmonicAngle>(interactions).parameters;
    computeForces(indices, angles, x, &forces, *pbc);

    RefDataChecker vector3DTest(1e-4);
    vector3DTest.testArrays<Vec3>(forces, "Angle forces");
}

TEST_F(ListedExampleData, CanReduceForces)
{
    reduceListedForces(interactions, x, &forces, *pbc);

    RefDataChecker vector3DTest(1e-2);
    vector3DTest.testArrays<Vec3>(forces, "Reduced forces");
}

TEST_F(ListedExampleData, CanReduceEnergies)
{
    auto energies    = reduceListedForces(interactions, x, &forces, *pbc);
    real totalEnergy = std::accumulate(begin(energies), end(energies), 0.0);

    RefDataChecker vector3DTest(1e-4);
    vector3DTest.testReal(totalEnergy, "Reduced energy");
}


void compareArray(const ListedForceCalculator::EnergyType& energies,
                  const ListedForceCalculator::EnergyType& refEnergies)
{
    for (size_t i = 0; i < energies.size(); ++i)
    {
        EXPECT_REAL_EQ_TOL(energies[i],
                           refEnergies[i],
                           gmx::test::relativeToleranceAsFloatingPoint(refEnergies[i], 1e-5));
    }
}


//! \brief sets up an interaction tuple for a linear chain with nParticles
class LinearChainDataFixture : public ::testing::Test
{
protected:
    void SetUp() override
    {
        LinearChainData data(430);

        x            = data.x;
        interactions = data.interactions;
        box          = data.box;
        refForces    = data.forces;

        refEnergies = reduceListedForces(interactions, x, &refForces, NoPbc{});
    }

    void testEnergies(const ListedForceCalculator::EnergyType& energies) const
    {
        compareArray(energies, refEnergies);
    }

    void testForces(const std::vector<gmx::RVec>& forces) const
    {
        compareVectors(forces, refForces, refForces);
    }

    std::vector<gmx::RVec> x;
    ListedInteractionData  interactions;
    std::shared_ptr<Box>   box;

private:
    std::vector<gmx::RVec>            refForces;
    ListedForceCalculator::EnergyType refEnergies;
};

TEST_F(LinearChainDataFixture, Multithreading)
{
    ListedForceCalculator lfCalculator(interactions, x.size(), 4, *box);

    std::vector<Vec3>                 forces(x.size(), Vec3{ 0, 0, 0 });
    ListedForceCalculator::EnergyType energies;
    lfCalculator.compute(x, forces, energies);

    testEnergies(energies);
    testForces(forces);
}


} // namespace
} // namespace test
} // namespace nblib
