/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

// includesorter/check-source want this include here. IMO a bug


#include <valarray>

#include <gtest/gtest.h>

#include "nblib/listed_forces/calculator.h"
#include "nblib/listed_forces/dataflow.hpp"
#include "nblib/listed_forces/tests/linear_chain_input.hpp"
#include "nblib/tests/testhelpers.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"


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
                    forces[i][m], refForcesFloat[i][m], refForcesDouble[i][m],
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

        DefaultAngle                                angle(Degrees(108.53), 397.5);
        std::vector<DefaultAngle>                   angles{ angle };
        std::vector<InteractionIndex<DefaultAngle>> angleIndices{ { 0, 1, 2, 0 } };

        pickType<HarmonicBondType>(interactions).indices    = bondIndices;
        pickType<HarmonicBondType>(interactions).parameters = bonds;

        pickType<DefaultAngle>(interactions).indices    = angleIndices;
        pickType<DefaultAngle>(interactions).parameters = angles;

        // initial position for the methanol atoms from the spc-water example
        x = std::vector<gmx::RVec>{ { 1.97, 1.46, 1.209 }, { 1.978, 1.415, 1.082 }, { 1.905, 1.46, 1.03 } };
        forces = std::vector<gmx::RVec>(3, gmx::RVec{ 0, 0, 0 });

        refBondForcesFloat =
                std::valarray<gmx::BasicVector<float>>{ { -22.8980637, 128.801575, 363.505951 },
                                                        { -43.2698593, -88.0130997, -410.639252 },
                                                        { 66.167923, -40.788475, 47.1333084 } };
        refAngleForcesFloat =
                std::valarray<gmx::BasicVector<float>>{ { 54.7276611, -40.1688995, 17.6805191 },
                                                        { -81.8118973, 86.1988525, 60.1752243 },
                                                        { 27.0842342, -46.0299492, -77.8557434 } };

        refBondForcesDouble = std::valarray<gmx::BasicVector<double>>{
            { -22.89764839974935, 128.79927224858977, 363.50016834602064 },
            { -43.24622441913251, -88.025652017772231, -410.61635172385434 },
            { 66.14387281888186, -40.773620230817542, 47.116183377833721 }
        };
        refAngleForcesDouble = std::valarray<gmx::BasicVector<double>>{
            { 54.726206806506234, -40.167809526198099, 17.680008528590257 },
            { -81.809781666748606, 86.196545126117257, 60.173723525141448 },
            { 27.083574860242372, -46.028735599919159, -77.853732053731704 }
        };

        refBondEnergyFloat  = 0.2113433;
        refAngleEnergyFloat = 0.112774156;

        refBondEnergyDouble  = 0.2113273434867636;
        refAngleEnergyDouble = 0.11276812148357591;

        box.reset(new Box(3, 3, 3));
        pbc.reset(new PbcHolder(*box));
    }

    std::vector<gmx::RVec> x;
    std::vector<gmx::RVec> forces;

    ListedInteractionData interactions;

    std::shared_ptr<Box>       box;
    std::shared_ptr<PbcHolder> pbc;

    // reference values
    std::valarray<gmx::BasicVector<float>>  refBondForcesFloat, refAngleForcesFloat;
    std::valarray<gmx::BasicVector<double>> refBondForcesDouble, refAngleForcesDouble;

    float  refBondEnergyFloat, refAngleEnergyFloat;
    double refBondEnergyDouble, refAngleEnergyDouble;
};

TEST_F(ListedExampleData, ComputeHarmonicBondForces)
{
    auto indices = pickType<HarmonicBondType>(interactions).indices;
    auto bonds   = pickType<HarmonicBondType>(interactions).parameters;
    real energy  = computeForces(indices, bonds, x, &forces, *pbc);

    EXPECT_FLOAT_DOUBLE_EQ_TOL(energy, refBondEnergyFloat, refBondEnergyDouble,
                               gmx::test::relativeToleranceAsFloatingPoint(refBondEnergyDouble, 1e-5));

    compareVectors(forces, refBondForcesFloat, refBondForcesDouble);
}

TEST_F(ListedExampleData, ComputeHarmonicAngleForces)
{
    auto indices = pickType<DefaultAngle>(interactions).indices;
    auto angles  = pickType<DefaultAngle>(interactions).parameters;
    real energy  = computeForces(indices, angles, x, &forces, *pbc);

    EXPECT_FLOAT_DOUBLE_EQ_TOL(energy, refAngleEnergyFloat, refAngleEnergyDouble,
                               gmx::test::relativeToleranceAsFloatingPoint(refAngleEnergyDouble, 1e-5));

    compareVectors(forces, refAngleForcesFloat, refAngleForcesDouble);
}

TEST_F(ListedExampleData, CanReduceForces)
{
    auto energies    = reduceListedForces(interactions, x, &forces, *pbc);
    real totalEnergy = std::accumulate(begin(energies), end(energies), 0.0);

    EXPECT_FLOAT_DOUBLE_EQ_TOL(totalEnergy, refBondEnergyFloat + refAngleEnergyFloat,
                               refBondEnergyDouble + refAngleEnergyDouble,
                               gmx::test::relativeToleranceAsFloatingPoint(refBondEnergyDouble, 1e-5));

    compareVectors(forces, refBondForcesFloat + refAngleForcesFloat,
                   refBondForcesDouble + refAngleForcesDouble);
}


void compareArray(const ListedForceCalculator::EnergyType& energies,
                  const ListedForceCalculator::EnergyType& refEnergies)
{
    for (size_t i = 0; i < energies.size(); ++i)
    {
        EXPECT_REAL_EQ_TOL(energies[i], refEnergies[i],
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
        // pbc.reset(new PbcHolder(*box));

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
    // std::shared_ptr<PbcHolder> pbc;

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
