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
 * This tests that sample code can run
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include <cstdio>

#include "gromacs/utility/arrayref.h"

// The entire nblib public API can be included with a single header or individual components
// can be included via their respective headers.
#include "nblib/nblib.h"

int main()
{
    // Create an argon particle with a name and a mass.
    nblib::ParticleType argonAtom(nblib::ParticleTypeName("Ar"), nblib::Mass(39.94800));
    // Create an argon molecule.
    nblib::Molecule argonMolecule(nblib::MoleculeName("AR"));
    // Add the argon particle to a molecule. The names are for bookkeeping and need not match.
    argonMolecule.addParticle(nblib::ParticleName("Argon"), argonAtom);
    // Define Lennard-Jones params for argon (parameters from gromos43A1).
    nblib::C6  ArC6{ 0.0062647225 };  // C6 parameter
    nblib::C12 ArC12{ 9.847044e-06 }; // C12 parameter
    // Holder for non-bonded interactions.
    nblib::ParticleTypesInteractions interactions;
    // Add non-bonded interactions for argon.
    interactions.add(argonAtom.name(), ArC6, ArC12);
    // The TopologyBuilder builds the Topology!
    nblib::TopologyBuilder topologyBuilder;
    // Number of Argon particles (molecules) in the system.
    int numParticles = 12;
    // Add the requested number of argon molecules to a topology.
    topologyBuilder.addMolecule(argonMolecule, numParticles);
    // Add the argon interactions to the topology.
    topologyBuilder.addParticleTypesInteractions(interactions);
    // Build the topology.
    nblib::Topology topology = topologyBuilder.buildTopology();
    // The system needs a bounding box. Only cubic and rectangular boxes are supported.
    nblib::Box box(6.05449);
    // User defined coordinates.
    std::vector<nblib::Vec3> coordinates = {
        { 0.794, 1.439, 0.610 }, { 1.397, 0.673, 1.916 }, { 0.659, 1.080, 0.573 },
        { 1.105, 0.090, 3.431 }, { 1.741, 1.291, 3.432 }, { 1.936, 1.441, 5.873 },
        { 0.960, 2.246, 1.659 }, { 0.382, 3.023, 2.793 }, { 0.053, 4.857, 4.242 },
        { 2.655, 5.057, 2.211 }, { 4.114, 0.737, 0.614 }, { 5.977, 5.104, 5.217 },
    };
    // User defined velocities.
    std::vector<nblib::Vec3> velocities = {
        { 0.0055, -0.1400, 0.2127 },   { 0.0930, -0.0160, -0.0086 }, { 0.1678, 0.2476, -0.0660 },
        { 0.1591, -0.0934, -0.0835 },  { -0.0317, 0.0573, 0.1453 },  { 0.0597, 0.0013, -0.0462 },
        { 0.0484, -0.0357, 0.0168 },   { 0.0530, 0.0295, -0.2694 },  { -0.0550, -0.0896, 0.0494 },
        { -0.0799, -0.2534, -0.0079 }, { 0.0436, -0.1557, 0.1849 },  { -0.0214, 0.0446, 0.0758 },
    };
    // Force buffer initialization for each particle.
    std::vector<nblib::Vec3> forces = {
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
        { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 }, { 0.0000, 0.0000, 0.0000 },
    };
    // A simulation state contains all the molecular information about the system.
    nblib::SimulationState simState(coordinates, velocities, forces, box, topology);
    // Kernel options are flags needed for force calculation.
    nblib::NBKernelOptions options = nblib::NBKernelOptions();
    // Use a simple cutoff rule for Coulomb
    options.coulombType = nblib::CoulombType::Cutoff;
    // Some performance flags can be set a run time
    options.nbnxmSimd = nblib::SimdKernels::SimdNo;
    // The force calculator contains all the data needed to compute forces.
    auto forceCalculator = nblib::setupGmxForceCalculatorCpu(simState.topology(), options);
    // build the pairlist
    forceCalculator->updatePairlist(simState.coordinates(), simState.box());
    // Integration requires masses, positions, and forces
    nblib::LeapFrog integrator(simState.topology(), simState.box());
    // Print some diagnostic info
    printf("initial forces on particle 0: x %4f y %4f z %4f\n", forces[0][0], forces[0][1], forces[0][2]);
    // The forces are computed for the user
    gmx::ArrayRef<nblib::Vec3> userForces(simState.forces());
    forceCalculator->compute(simState.coordinates(), simState.box(), userForces);
    // Print some diagnostic info
    printf("  final forces on particle 0: x %4f y %4f z %4f\n",
           userForces[0][0],
           userForces[0][1],
           userForces[0][2]);
    // User may modify forces stored in simState.forces() if needed
    // Print some diagnostic info
    printf("initial position of particle 0: x %4f y %4f z %4f\n",
           simState.coordinates()[0][0],
           simState.coordinates()[0][1],
           simState.coordinates()[0][2]);
    // Integrate with a time step of 1 fs
    integrator.integrate(1.0, simState.coordinates(), simState.velocities(), simState.forces());
    // Print some diagnostic info

    printf("  final position of particle 0: x %4f y %4f z %4f\n",
           simState.coordinates()[0][0],
           simState.coordinates()[0][1],
           simState.coordinates()[0][2]);
    return 0;
}
