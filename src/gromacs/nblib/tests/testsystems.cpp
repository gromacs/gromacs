/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * This implements nblib test helpers
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "testsystems.h"

namespace nblib
{

struct OwAtom
{
    AtomName name = "Ow";
    Mass     mass = 15.99940;
    C6       c6   = 0.0026173456;
    C12      c12  = 2.634129e-06;
};

struct UnitedHAtom
{
    AtomName name = "H";
    Mass     mass = 1.008;
    C6       c6   = 0;
    C12      c12  = 0;
};

struct OMetAtom
{
    AtomName name = "OMet";
    Mass     mass = 15.999;
    C6       c6   = 0.0022619536;
    C12      c12  = 1.505529e-06;
};

struct CMetAtom
{
    AtomName name = "CMet";
    Mass     mass = 15.035; // Direct from FF
    C6       c6   = 0.0088755241;
    C12      c12  = 2.0852922e-05;
};

std::unordered_map<std::string, Charge> Charges{ { "Ow", -0.82 },
                                                 { "Hw", +0.41 },
                                                 { "OMet", -0.574 },
                                                 { "CMet", +0.176 },
                                                 { "HMet", +0.398 } };

WaterMoleculeBuilder::WaterMoleculeBuilder() : water_("water")
{
    //! Define Atom Types
    OwAtom      owAtom;
    AtomType    Ow(owAtom.name, owAtom.mass, owAtom.c6, owAtom.c12);
    UnitedHAtom hwAtom;
    AtomType    Hw(hwAtom.name, hwAtom.mass, hwAtom.c6, hwAtom.c12);

    //! Add the atoms
    water_.addAtom(AtomName("Oxygen"), Charges.at("Ow"), Ow);
    water_.addAtom(AtomName("H1"), Charges.at("Hw"), Hw);
    water_.addAtom(AtomName("H2"), Charges.at("Hw"), Hw);
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

Topology TwoWaterMolecules::buildTopology()
{
    //! Add some molecules to the topology
    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(water(), 2);
    Topology topology = topologyBuilder.buildTopology();
    return topology;
}

Molecule TwoWaterMolecules::water()
{
    return waterMolecule_.waterMolecule();
}

class ArgonTopologyBuilder
{
public:
    ArgonTopologyBuilder(const int& numAtoms)
    {
        ArAtom   arAtom;
        AtomType argonAtom(arAtom.name, arAtom.mass, arAtom.c6, arAtom.c12);

        Molecule argonMolecule("AR");
        argonMolecule.addAtom(AtomName("AR"), argonAtom);

        topologyBuilder_.addMolecule(argonMolecule, numAtoms);
    }

    Topology argonTopology() { return topologyBuilder_.buildTopology(); }

private:
    TopologyBuilder topologyBuilder_;
};

SimulationStateTester::SimulationStateTester() :
    box_(7.25449),
    topology_(ArgonTopologyBuilder(12).argonTopology())
{

    coordinates_ = {
        { 0.794, 1.439, 0.610 }, { 1.397, 0.673, 1.916 }, { 0.659, 1.080, 0.573 },
        { 1.105, 0.090, 3.431 }, { 1.741, 1.291, 3.432 }, { 1.936, 1.441, 5.873 },
        { 0.960, 2.246, 1.659 }, { 0.382, 3.023, 2.793 }, { 0.053, 4.857, 4.242 },
        { 2.655, 5.057, 2.211 }, { 4.114, 0.737, 0.614 }, { 5.977, 5.104, 5.217 },
    };

    velocities_ = {
        { 0.0055, -0.1400, 0.2127 },   { 0.0930, -0.0160, -0.0086 }, { 0.1678, 0.2476, -0.0660 },
        { 0.1591, -0.0934, -0.0835 },  { -0.0317, 0.0573, 0.1453 },  { 0.0597, 0.0013, -0.0462 },
        { 0.0484, -0.0357, 0.0168 },   { 0.0530, 0.0295, -0.2694 },  { -0.0550, -0.0896, 0.0494 },
        { -0.0799, -0.2534, -0.0079 }, { 0.0436, -0.1557, 0.1849 },  { -0.0214, 0.0446, 0.0758 },
    };
}

void SimulationStateTester::setCoordinate(int atomNum, int dimension, real value)
{
    GMX_ASSERT((dimension < 0 || dimension > 2), "Must provide a valid dimension\n");
    coordinates_.at(atomNum)[dimension] = value;
}

void SimulationStateTester::setVelocity(int atomNum, int dimension, real value)
{
    GMX_ASSERT((dimension < 0 || dimension > 2), "Must provide a valid dimension\n");
    velocities_.at(atomNum)[dimension] = value;
}

SimulationState SimulationStateTester::setupSimulationState()
{
    return SimulationState(coordinates_, box_, topology_, velocities_);
}

const Topology& SimulationStateTester::topology() const
{
    return topology_;
}

Box& SimulationStateTester::box()
{
    return box_;
}

std::vector<gmx::RVec>& SimulationStateTester::coordinates()
{
    return coordinates_;
}

std::vector<gmx::RVec>& SimulationStateTester::velocities()
{
    return velocities_;
}

} // namespace nblib
