/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

// The entire nblib public API can be included with a single header or individual components
// can be included via their respective headers.
#include "nblib/nblib.h"

using namespace nblib;

int main()
{
    // Create the particles
    ParticleType Ow(ParticleTypeName("Ow"), Mass(15.999));
    ParticleType Hw(ParticleTypeName("Hw"), Mass(1.00784));
    ParticleType Cm(ParticleTypeName("Cm"), Mass(12.0107));
    ParticleType Hc(ParticleTypeName("Hc"), Mass(1.00784));

    ParticleTypesInteractions interactions(CombinationRule::Geometric);

    // Parameters from a GROMOS compatible force-field 2016H66
    // add non-bonded interactions for the particle types
    interactions.add(Ow.name(), C6(0.0026173456), C12(2.634129e-06));
    interactions.add(Hw.name(), C6(0.0), C12(0.0));
    interactions.add(Cm.name(), C6(0.01317904), C12(34.363044e-06));
    interactions.add(Hc.name(), C6(8.464e-05), C12(15.129e-09));

    Molecule water(MoleculeName("Water"));
    Molecule methane(MoleculeName("Methane"));

    water.addParticle(ParticleName("O"), Ow);
    water.addParticle(ParticleName("H1"), Hw);
    water.addParticle(ParticleName("H2"), Hw);

    water.addExclusion(ParticleName("H1"), ParticleName("O"));
    water.addExclusion(ParticleName("H2"), ParticleName("O"));

    methane.addParticle(ParticleName("C"), Cm);
    methane.addParticle(ParticleName("H1"), Hc);
    methane.addParticle(ParticleName("H2"), Hc);
    methane.addParticle(ParticleName("H3"), Hc);
    methane.addParticle(ParticleName("H4"), Hc);

    methane.addExclusion(ParticleName("H1"), ParticleName("C"));
    methane.addExclusion(ParticleName("H2"), ParticleName("C"));
    methane.addExclusion(ParticleName("H3"), ParticleName("C"));
    methane.addExclusion(ParticleName("H4"), ParticleName("C"));

    HarmonicBondType ohHarmonicBond(1, 1);
    HarmonicBondType hcHarmonicBond(2, 1);

    HarmonicAngle hohAngle(1, Degrees(120));
    HarmonicAngle hchAngle(1, Degrees(109.5));

    // add harmonic bonds for water
    water.addInteraction(ParticleName("O"), ParticleName("H1"), ohHarmonicBond);
    water.addInteraction(ParticleName("O"), ParticleName("H2"), ohHarmonicBond);

    // add the angle for water
    water.addInteraction(ParticleName("H1"), ParticleName("O"), ParticleName("H2"), hohAngle);

    // add harmonic bonds for methane
    methane.addInteraction(ParticleName("H1"), ParticleName("C"), hcHarmonicBond);
    methane.addInteraction(ParticleName("H2"), ParticleName("C"), hcHarmonicBond);
    methane.addInteraction(ParticleName("H3"), ParticleName("C"), hcHarmonicBond);
    methane.addInteraction(ParticleName("H4"), ParticleName("C"), hcHarmonicBond);

    // add the angles for methane
    methane.addInteraction(ParticleName("H1"), ParticleName("C"), ParticleName("H2"), hchAngle);
    methane.addInteraction(ParticleName("H1"), ParticleName("C"), ParticleName("H3"), hchAngle);
    methane.addInteraction(ParticleName("H1"), ParticleName("C"), ParticleName("H4"), hchAngle);
    methane.addInteraction(ParticleName("H2"), ParticleName("C"), ParticleName("H3"), hchAngle);
    methane.addInteraction(ParticleName("H2"), ParticleName("C"), ParticleName("H4"), hchAngle);
    methane.addInteraction(ParticleName("H3"), ParticleName("C"), ParticleName("H4"), hchAngle);

    // Define a box for the simulation
    Box box(6.05449);

    // Define options for the non-bonded kernels
    NBKernelOptions options;
    // Use a simple cutoff rule for Coulomb
    options.coulombType = nblib::CoulombType::Cutoff;
    // Disable SIMD for this example
    options.nbnxmSimd = SimdKernels::SimdNo;

    TopologyBuilder topologyBuilder;

    // add molecules
    topologyBuilder.addMolecule(water, 2);
    topologyBuilder.addMolecule(methane, 1);

    // add non-bonded interaction map
    topologyBuilder.addParticleTypesInteractions(interactions);

    Topology topology = topologyBuilder.buildTopology();

    // User defined coordinates.
    std::vector<Vec3> coordinates = {
        { 0.005, 0.600, 0.244 },    // Oxygen from water_1
        { -0.017, 0.690, 0.270 },   // Hydrogen_1 from water_1
        { 0.051, 0.610, 0.161 },    // Hydrogen_2 from water_1
        { 0.155, 0.341, 0.735 },    // Oxygen from water_2
        { 0.140, 0.284, 0.660 },    // Hydrogen_1 from water_2
        { 0.081, 0.402, 0.734 },    // Hydrogen_2 from water_2
        { -0.024, -0.222, -0.640 }, // Carbon from methane_1
        { -0.083, -0.303, -0.646 }, // Hydrogen_1 from methane_1
        { -0.080, -0.140, -0.642 }, // Hydrogen_2 from methane_1
        { 0.040, -0.221, -0.716 },  // Hydrogen_3 from methane_1
        { 0.027, -0.225, -0.553 }   // Hydrogen_4 from methane_1
    };

    // User defined velocities.
    std::vector<Vec3> velocities = {
        { 0.1823, -0.4158, 0.487 },   // Oxygen from water_1
        { -1.7457, -0.5883, -0.460 }, // Hydrogen_1 from water_1
        { 2.5085, -0.1501, 1.762 },   // Hydrogen_2 from water_1
        { 0.6282, 0.4390, 0.001 },    // Oxygen from water_2
        { -0.3206, 0.0700, 0.4630 },  // Hydrogen_1 from water_2
        { -0.1556, -0.4529, 1.440 },  // Hydrogen_2 from water_2
        { 0.0, 0.0, 0.0 },            // Carbon from methane_1
        { 0.0, 0.0, 0.0 },            // Hydrogen_1 from methane_1
        { 0.0, 0.0, 0.0 },            // Hydrogen_2 from methane_1
        { 0.0, 0.0, 0.0 },            // Hydrogen_3 from methane_1
        { 0.0, 0.0, 0.0 },            // Hydrogen_4 from methane_1
    };

    // Force buffer initialization for each particle.
    std::vector<Vec3> forces(topology.numParticles(), Vec3{ 0, 0, 0 });

    SimulationState simulationState(coordinates, velocities, forces, box, topology);

    // The non-bonded force calculator contains all the data needed to compute forces
    auto forceCalculator = setupGmxForceCalculatorCpu(simulationState.topology(), options);

    // build the pair list
    forceCalculator->updatePairlist(simulationState.coordinates(), simulationState.box());

    // The listed force calculator is also initialized with the required arguments
    ListedForceCalculator listedForceCalculator(
            topology.getInteractionData(), topology.numParticles(), 4, box);

    // Integrator is initialized with an array of inverse masses (constructed from topology) and
    // the bounding box
    LeapFrog integrator(simulationState.topology(), simulationState.box());

    // Print some diagnostic info
    printf("initial position of particle 0: x %4f y %4f z %4f\n",
           simulationState.coordinates()[0][0],
           simulationState.coordinates()[0][1],
           simulationState.coordinates()[0][2]);

    // MD Loop
    int numSteps = 2;

    for (auto i = 0; i < numSteps; i++)
    {
        zeroCartesianArray(simulationState.forces());

        forceCalculator->compute(
                simulationState.coordinates(), simulationState.box(), simulationState.forces());

        listedForceCalculator.compute(
                simulationState.coordinates(), simulationState.forces(), gmx::ArrayRef<real>{});

        // Integrate with a time step of 1 fs, positions, velocities and forces
        integrator.integrate(
                1.0, simulationState.coordinates(), simulationState.velocities(), simulationState.forces());
    }

    printf("  final position of particle 9: x %4f y %4f z %4f\n",
           simulationState.coordinates()[9][0],
           simulationState.coordinates()[9][1],
           simulationState.coordinates()[9][2]);

    return 0;
} // main
