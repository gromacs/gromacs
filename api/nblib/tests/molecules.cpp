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
 * This implements molecule setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/molecules.h"

#include <cstddef>

#include <stdexcept>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

#include "nblib/exception.h"
#include "nblib/particletype.h"
#include "nblib/tests/testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, CanConstructMoleculeWithoutChargeOrResidueName)
{
    ArAtom       arAtom;
    ParticleType Ar(arAtom.particleTypeName, arAtom.mass);
    Molecule     argon(arAtom.moleculeName);
    EXPECT_NO_THROW(argon.addParticle(arAtom.particleName, Ar));
}

TEST(NBlibTest, CanConstructMoleculeWithChargeWithoutResidueName)
{
    ArAtom       arAtom;
    ParticleType Ar(arAtom.particleTypeName, arAtom.mass);
    Molecule     argon(arAtom.moleculeName);
    EXPECT_NO_THROW(argon.addParticle(arAtom.particleName, Charge(0), Ar));
}

TEST(NBlibTest, CanConstructMoleculeWithoutChargeWithResidueName)
{
    ArAtom       arAtom;
    ParticleType Ar(arAtom.particleTypeName, arAtom.mass);
    Molecule     argon(arAtom.moleculeName);
    EXPECT_NO_THROW(argon.addParticle(arAtom.particleName, ResidueName("ar2"), Ar));
}

TEST(NBlibTest, CanConstructMoleculeWithChargeWithResidueName)
{
    ArAtom       arAtom;
    ParticleType Ar(arAtom.particleTypeName, arAtom.mass);
    Molecule     argon(arAtom.moleculeName);
    EXPECT_NO_THROW(argon.addParticle(arAtom.particleName, ResidueName("ar2"), Charge(0), Ar));
}

TEST(NBlibTest, CanGetNumParticlesInMolecule)
{
    WaterMoleculeBuilder waterMolecule;
    Molecule             water        = waterMolecule.waterMolecule();
    auto                 numParticles = water.numParticlesInMolecule();

    EXPECT_EQ(3, numParticles);
}

TEST(NBlibTest, CanConstructExclusionListFromNames)
{
    WaterMoleculeBuilder waterMolecule;
    Molecule             water = waterMolecule.waterMolecule();

    std::vector<std::tuple<int, int>> exclusions = water.getExclusions();

    std::vector<std::tuple<int, int>> reference{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 },
                                                 { 1, 2 }, { 2, 0 }, { 2, 1 }, { 2, 2 } };

    ASSERT_EQ(exclusions.size(), 9);
    for (std::size_t i = 0; i < exclusions.size(); ++i)
    {
        EXPECT_EQ(exclusions[i], reference[i]);
    }
}

TEST(NBlibTest, CanConstructExclusionListFromNamesAndIndicesMixed)
{
    WaterMoleculeBuilder waterMolecule;
    Molecule             water = waterMolecule.waterMoleculeWithoutExclusions();

    //! Add the exclusions
    water.addExclusion(ParticleName("H1"), ParticleName("Oxygen"));
    water.addExclusion(ParticleName("H2"), ParticleName("Oxygen"));
    water.addExclusion(ParticleName("H1"), ParticleName("H2"));

    std::vector<std::tuple<int, int>> exclusions = water.getExclusions();

    std::vector<std::tuple<int, int>> reference{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 },
                                                 { 1, 2 }, { 2, 0 }, { 2, 1 }, { 2, 2 } };

    ASSERT_EQ(exclusions.size(), 9);
    for (std::size_t i = 0; i < exclusions.size(); ++i)
    {
        EXPECT_EQ(exclusions[i], reference[i]);
    }
}

TEST(NBlibTest, AtWorks)
{
    WaterMoleculeBuilder waterMolecule;
    Molecule             water = waterMolecule.waterMolecule();
    EXPECT_NO_THROW(water.at("Ow"));
    EXPECT_NO_THROW(water.at("H"));
}

TEST(NBlibTest, AtThrows)
{
    WaterMoleculeBuilder waterMolecule;
    Molecule             water = waterMolecule.waterMolecule();
    EXPECT_THROW_GMX(water.at("Hw"), std::out_of_range);
}

TEST(NBlibTest, MoleculeThrowsSameParticleTypeNameDifferentMass)
{
    //! User error: Two different ParticleTypes with the same name
    ParticleType atom1(ParticleTypeName("Atom"), Mass(1));
    ParticleType atom2(ParticleTypeName("Atom"), Mass(2));

    Molecule molecule(MoleculeName("UraniumDimer"));
    EXPECT_NO_THROW(molecule.addParticle(ParticleName("U1"), atom1));
    EXPECT_THROW(molecule.addParticle(ParticleName("U2"), atom2), InputException);
}

TEST(NBlibTest, MoleculeDontThrowsSameParticleTypeNameDifferentMass)
{
    //! User error: Two different ParticleTypes with the same name
    ParticleType atom1(ParticleTypeName("Atom"), Mass(1));
    ParticleType atom2(ParticleTypeName("Atom"), Mass(1));

    Molecule molecule(MoleculeName("UraniumDimer"));
    EXPECT_NO_THROW(molecule.addParticle(ParticleName("U1"), atom1));
    EXPECT_NO_THROW(molecule.addParticle(ParticleName("U2"), atom2));
}

TEST(NBlibTest, MoleculeNoThrowsSameParticleTypeName)
{
    //! User error: Two different ParticleTypes with the same name
    ParticleType atom1(ParticleTypeName("Atom"), Mass(1));
    ParticleType atom2(ParticleTypeName("Atom"), Mass(1));

    Molecule molecule(MoleculeName("UraniumDimer"));
    EXPECT_NO_THROW(molecule.addParticle(ParticleName("U1"), atom1));
    EXPECT_NO_THROW(molecule.addParticle(ParticleName("U2"), atom2));
}

TEST(NBlibTest, CanAddInteractions)
{
    Molecule     molecule(MoleculeName("BondTest"));
    ParticleType O(ParticleTypeName("Ow"), Mass(1));
    ParticleType H(ParticleTypeName("Hw"), Mass(1));
    molecule.addParticle(ParticleName("O"), O);
    molecule.addParticle(ParticleName("H1"), H);
    molecule.addParticle(ParticleName("H2"), H);

    HarmonicBondType hb(1, 2);
    CubicBondType    cub(1, 2, 3);
    HarmonicAngle    ang(1, Degrees(1));

    molecule.addInteraction(ParticleName("O"), ParticleName("H1"), hb);
    molecule.addInteraction(ParticleName("O"), ParticleName("H2"), hb);
    molecule.addInteraction(ParticleName("H1"), ParticleName("H2"), cub);
    molecule.addInteraction(ParticleName("H1"), ParticleName("O"), ParticleName("H2"), ang);

    const auto& interactionData = molecule.interactionData();

    //! harmonic bonds
    EXPECT_EQ(pickType<HarmonicBondType>(interactionData).interactions_.size(), 2);
    //! cubic bonds
    EXPECT_EQ(pickType<CubicBondType>(interactionData).interactions_.size(), 1);
    //! angular interactions
    EXPECT_EQ(pickType<HarmonicAngle>(interactionData).interactions_.size(), 1);
}

TEST(NBlibTest, CanAddUreyBradley)
{
    Molecule     molecule(MoleculeName("UreyBradleyTest"));
    ParticleType O(ParticleTypeName("Ow"), Mass(1));
    ParticleType H(ParticleTypeName("Hw"), Mass(1));
    molecule.addParticle(ParticleName("O"), O);
    molecule.addParticle(ParticleName("H1"), H);
    molecule.addParticle(ParticleName("H2"), H);

    addUreyBradleyInteraction(molecule,
                              ParticleName("H1"),
                              ParticleName("O"),
                              ParticleName("H2"),
                              Radians(1.82),
                              ForceConstant(2.0),
                              EquilConstant(0.15),
                              ForceConstant(3.0));

    const auto& interactionData = molecule.interactionData();

    //! harmonic bonds
    EXPECT_EQ(pickType<HarmonicBondType>(interactionData).interactions_.size(), 1);
    //! angular interactions
    EXPECT_EQ(pickType<HarmonicAngle>(interactionData).interactions_.size(), 1);
}

} // namespace
} // namespace test
} // namespace nblib
