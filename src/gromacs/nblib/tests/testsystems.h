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
 * This implements basic nblib test systems
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GROMACS_SETUPHELPERS_H
#define GROMACS_SETUPHELPERS_H

#include "gmxpre.h"

#include <cmath>

#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/box.h"
#include "gromacs/nblib/molecules.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nblib/topology.h"

namespace nblib
{

//! \internal \brief Parameters from gromos43A1
struct ArAtom
{
    ParticleName name = "Ar";
    Mass         mass = 39.94800;
    C6           c6   = 0.0062647225;
    C12          c12  = 9.847044e-06;
};

//! Parameters from gromos43A1 SPC Oxygen
struct OwAtom;

//! Parameters from gromos43A1
struct HAtom;

//! Parameters from gromos43A1
struct OMettom;

//! Parameters from gromos43A1
struct CMetAtom;

//! Lookup table for charges needed for building topologies
extern std::unordered_map<std::string, Charge> Charges;

//! \internal \brief Make an SPC water molecule with parameters from gromos43A1
class WaterMoleculeBuilder
{
public:
    //! There is no default ctor for a Molecule so it must be initialized
    WaterMoleculeBuilder();

    //! Return the initialized water Molecule, with exclusions
    Molecule waterMolecule();

    //! Return the initialized water Molecule, without exclusions
    Molecule waterMoleculeWithoutExclusions();

private:
    //! The molecule
    Molecule water_;

    //! Add the exclusions from particle names. Private to prevent multiple calls
    void addExclusionsFromNames();
};

//! \brief Make a methanol molecule with parameters from gromos43A1
class MethanolMoleculeBuilder
{
public:
    //! There is no default ctor for a Molecule so it must be initialized
    MethanolMoleculeBuilder();

    //! Return the initialized water Molecule, with exclusions
    Molecule methanolMolecule();

private:
    //! The molecule
    Molecule methanol_;
};

class WaterTopology
{
public:
    //! Return a topology with specified SPC water molecules
    Topology buildTopology(int numMolecules);

    //! Return the actual water Molecule used in the topology
    Molecule water();

private:
    WaterMoleculeBuilder waterMolecule_;
};

class SpcMethanolTopologyBuilder
{
public:
    //! Return a topology with specified methanol molecules
    Topology buildTopology(int numWater, int numMethanol);

    //! Return the actual methanol Molecule used in the topology
    Molecule methanol();

    //! Return the actual water Molecule used in the topology
    Molecule water();

private:
    MethanolMoleculeBuilder methanolMolecule_;
    WaterMoleculeBuilder    waterMolecule_;
};

class ArgonTopologyBuilder
{
public:
    //! Build a topology with specified argon molecules
    ArgonTopologyBuilder(const int& numParticles);

    //! Get the topology with specified argon molecules
    Topology argonTopology();

private:
    TopologyBuilder topologyBuilder_;
};


class ArgonSimulationStateBuilder
{
public:
    ArgonSimulationStateBuilder();

    void setCoordinate(int particleNum, int dimension, real value);

    void setVelocity(int particleNum, int dimension, real value);

    SimulationState setupSimulationState();

    const Topology& topology() const;

    Box& box();

    std::vector<gmx::RVec>& coordinates();

    std::vector<gmx::RVec>& velocities();

private:
    std::vector<gmx::RVec> coordinates_;
    std::vector<gmx::RVec> velocities_;

    Box      box_;
    Topology topology_;
};

class SpcMethanolSimulationStateBuilder
{
public:
    SpcMethanolSimulationStateBuilder();

    SimulationState setupSimulationState();

    std::vector<gmx::RVec>& coordinates();

    std::vector<gmx::RVec>& velocities();

private:
    std::vector<gmx::RVec> coordinates_;
    std::vector<gmx::RVec> velocities_;

    Box      box_;
    Topology topology_;
};

} // namespace nblib
#endif // GROMACS_SETUPHELPERS_H
