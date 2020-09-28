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
 * This implements nblib test systems
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/tests/testsystems.h"
#include "nblib/exception.h"

namespace nblib
{

//! User class to hold ParticleTypes
//! Note: this is not part of NBLIB, users should write their own
class ParticleLibrary
{
public:
    ParticleLibrary()
    {
        ParticleType Ow(ParticleTypeName("Ow"), Mass(15.99940));
        ParticleType H(ParticleTypeName("H"), Mass(1.008));
        ParticleType OMet(ParticleTypeName("OMet"), Mass(15.999));
        ParticleType CMet(ParticleTypeName("CMet"), Mass(15.035));
        ParticleType Ar(ParticleTypeName("Ar"), Mass(39.94800));

        particles_.insert(std::make_pair(Ow.name(), Ow));
        particles_.insert(std::make_pair(H.name(), H));
        particles_.insert(std::make_pair(OMet.name(), OMet));
        particles_.insert(std::make_pair(CMet.name(), CMet));
        particles_.insert(std::make_pair(Ar.name(), Ar));

        c6_[Ow.name()]   = 0.0026173456;
        c6_[H.name()]    = 0;
        c6_[OMet.name()] = 0.0022619536;
        c6_[CMet.name()] = 0.0088755241;
        c6_[Ar.name()]   = 0.0062647225;

        c12_[Ow.name()]   = 2.634129e-06;
        c12_[H.name()]    = 0;
        c12_[OMet.name()] = 1.505529e-06;
        c12_[CMet.name()] = 2.0852922e-05;
        c12_[Ar.name()]   = 9.847044e-06;
    }

    //! Get particle type using the string identifier
    [[nodiscard]] ParticleType type(const std::string& particleTypeName) const
    {
        return particles_.at(ParticleTypeName(particleTypeName));
    }

    //! Get C6 parameter of a given particle
    [[nodiscard]] C6 c6(const ParticleName& particleName) const
    {
        return c6_.at(ParticleTypeName(particleName.value()));
    }

    //! Get C12 parameter of a given particle
    [[nodiscard]] C12 c12(const ParticleName& particleName) const
    {
        return c12_.at(ParticleTypeName(particleName.value()));
    }

private:
    std::map<ParticleTypeName, ParticleType> particles_;
    std::map<ParticleTypeName, C6>           c6_;
    std::map<ParticleTypeName, C12>          c12_;
};

std::unordered_map<std::string, Charge> Charges{ { "Ow", Charge(-0.82) },
                                                 { "Hw", Charge(+0.41) },
                                                 { "OMet", Charge(-0.574) },
                                                 { "CMet", Charge(+0.176) },
                                                 { "HMet", Charge(+0.398) } };

WaterMoleculeBuilder::WaterMoleculeBuilder() : water_(MoleculeName("SOL"))
{
    ParticleLibrary plib;

    //! Add the particles
    water_.addParticle(ParticleName("Oxygen"), Charges.at("Ow"), plib.type("Ow"));
    water_.addParticle(ParticleName("H1"), Charges.at("Hw"), plib.type("H"));
    water_.addParticle(ParticleName("H2"), Charges.at("Hw"), plib.type("H"));
}

Molecule WaterMoleculeBuilder::waterMolecule()
{
    addExclusionsFromNames();
    return water_;
}

Molecule WaterMoleculeBuilder::waterMoleculeWithoutExclusions()
{
    return water_;
}

void WaterMoleculeBuilder::addExclusionsFromNames()
{
    water_.addExclusion("H1", "Oxygen");
    water_.addExclusion("H2", "Oxygen");
    water_.addExclusion("H1", "H2");
}

MethanolMoleculeBuilder::MethanolMoleculeBuilder() : methanol_(MoleculeName("MeOH"))
{
    ParticleLibrary library;

    //! Add the particles
    methanol_.addParticle(ParticleName("Me1"), Charges.at("CMet"), library.type("CMet"));
    methanol_.addParticle(ParticleName("O2"), Charges.at("OMet"), library.type("OMet"));
    methanol_.addParticle(ParticleName("H3"), Charges.at("HMet"), library.type("H"));

    // Add the exclusions
    methanol_.addExclusion("Me1", "O2");
    methanol_.addExclusion("Me1", "H3");
    methanol_.addExclusion("H3", "O2");
}

Molecule MethanolMoleculeBuilder::methanolMolecule()
{
    return methanol_;
}


Topology WaterTopologyBuilder::buildTopology(int numMolecules)
{
    ParticleLibrary library;

    ParticleTypesInteractions interactions;
    std::vector<std::string>  typeNames = { "Ow", "H" };
    for (const auto& name : typeNames)
    {
        interactions.add(ParticleTypeName(name), library.c6(ParticleName(name)),
                         library.c12(ParticleName(name)));
    }

    // Add some molecules to the topology
    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(water(), numMolecules);

    // Add non-bonded interaction information
    topologyBuilder.addParticleTypesInteractions(interactions);

    Topology topology = topologyBuilder.buildTopology();
    return topology;
}

Molecule WaterTopologyBuilder::water()
{
    return waterMolecule_.waterMolecule();
}

Topology SpcMethanolTopologyBuilder::buildTopology(int numWater, int numMethanol)
{
    ParticleLibrary library;

    ParticleTypesInteractions interactions;
    std::vector<std::string>  typeNames = { "Ow", "H", "OMet", "CMet" };
    for (const auto& name : typeNames)
    {
        interactions.add(ParticleTypeName(name), library.c6(ParticleName(name)),
                         library.c12(ParticleName(name)));
    }

    // Add some molecules to the topology
    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(methanol(), numMethanol);
    topologyBuilder.addMolecule(water(), numWater);

    // Add non-bonded interaction information
    topologyBuilder.addParticleTypesInteractions(interactions);

    Topology topology = topologyBuilder.buildTopology();
    return topology;
}

Molecule SpcMethanolTopologyBuilder::methanol()
{
    return methanolMolecule_.methanolMolecule();
}

Molecule SpcMethanolTopologyBuilder::water()
{
    return waterMolecule_.waterMolecule();
}

ArgonTopologyBuilder::ArgonTopologyBuilder(const int& numParticles)
{
    ParticleLibrary library;

    ParticleTypesInteractions nbinteractions;
    nbinteractions.add(ParticleTypeName("Ar"), library.c6(ParticleName("Ar")),
                       library.c12(ParticleName("Ar")));

    Molecule argonMolecule(MoleculeName("AR"));
    argonMolecule.addParticle(ParticleName("AR"), library.type("Ar"));

    topologyBuilder_.addMolecule(argonMolecule, numParticles);
    topologyBuilder_.addParticleTypesInteractions((nbinteractions));
}

Topology ArgonTopologyBuilder::argonTopology()
{
    return topologyBuilder_.buildTopology();
}


} // namespace nblib
