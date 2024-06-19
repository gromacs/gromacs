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
 * This implements basic nblib test systems
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_TESTSYSTEMS_H
#define NBLIB_TESTSYSTEMS_H

#include <cmath>

#include <string>
#include <unordered_map>
#include <vector>

#include "nblib/basicdefinitions.h"
#include "nblib/box.h"
#include "nblib/molecules.h"
#include "nblib/particletype.h"
#include "nblib/simulationstate.h"
#include "nblib/topology.h"
#include "nblib/vector.h"

namespace nblib
{

//! \internal \brief Force field types to select parameter sets
enum class fftypes : int
{
    GROMOS43A1,
    OPLSA
};

//! \internal \brief Parameters from gromos43A1
struct ArAtom
{
    //! Argon Particle Name
    ParticleName particleName = ParticleName("Ar");
    //! Argon particle type name
    ParticleTypeName particleTypeName = ParticleTypeName("Ar");
    //! Argon molecule name
    MoleculeName moleculeName = MoleculeName("Ar");
    //! Argon Particle Mass
    Mass mass = Mass(39.94800);
    //! Argon C6 parameter
    C6 c6{ 0.0062647225 };
    //! Argon C12 parameter
    C12 c12{ 9.847044e-06 };
};

//! Lookup table for charges needed for building topologies
extern std::unordered_map<std::string, Charge> Charges;

//! \internal \brief Make an SPC water molecule with parameters from gromos43A1
class WaterMoleculeBuilder
{
public:
    // There is no default ctor for a Molecule so it must be initialized
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

//! \internal \brief Make a methanol molecule with parameters from gromos43A1
class MethanolMoleculeBuilder
{
public:
    // There is no default ctor for a Molecule so it must be initialized
    MethanolMoleculeBuilder();

    //! Return the initialized water Molecule, with exclusions
    Molecule methanolMolecule();

private:
    //! The molecule
    Molecule methanol_;
};

//! \internal \brief Build topology of water molecules of a specified number
class WaterTopologyBuilder
{
public:
    //! Return a topology with specified SPC water molecules
    Topology buildTopology(int numMolecules);

    //! Return the actual water Molecule used in the topology
    Molecule water();

private:
    WaterMoleculeBuilder waterMolecule_;
};

//! \internal \brief Build topology of methanol+water molecules from specified numbers
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

//! \internal \brief Build topology of argon molecules of a specified number
class ArgonTopologyBuilder
{
public:
    //! Build a topology with specified argon molecules
    ArgonTopologyBuilder(const int& numParticles, const fftypes forceField);

    //! Get the topology with specified argon molecules
    Topology argonTopology();

private:
    TopologyBuilder topologyBuilder_;
};

//! \internal \brief Build simulation state for the argon example
class ArgonSimulationStateBuilder
{
public:
    ArgonSimulationStateBuilder(const fftypes forceField);

    //! Set coordinates of particles in the defined system
    void setCoordinate(int particleNum, int dimension, real value);

    //! Set particle velocities
    void setVelocity(int particleNum, int dimension, real value);

    //! Setup simulation state
    SimulationState setupSimulationState();

    //! Get the topology
    const Topology& topology() const;

    //! Get the box bounding the system
    Box& box();

    //! Get current coordinates
    std::vector<Vec3>& coordinates();

    //! Get current velocities
    std::vector<Vec3>& velocities();

private:
    std::vector<Vec3> coordinates_;
    std::vector<Vec3> velocities_;
    std::vector<Vec3> forces_;

    Box      box_;
    Topology topology_;
};

//! \internal \brief Build simulation state for the SPC-Methanol example
class SpcMethanolSimulationStateBuilder
{
public:
    SpcMethanolSimulationStateBuilder();

    //! Setup simulation state
    SimulationState setupSimulationState();

    //! Get current coordinates
    std::vector<Vec3>& coordinates();

    //! Get current velocities
    std::vector<Vec3>& velocities();

private:
    std::vector<Vec3> coordinates_;
    std::vector<Vec3> velocities_;
    std::vector<Vec3> forces_;

    Box      box_;
    Topology topology_;
};

} // namespace nblib
#endif // NBLIB_TESTSYSTEMS_H
