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
 * Compilation instructions
 * 1) install gromacs legacy api and nblib api
 * 2) clang++ -std=c++17 -I/path/to/include/dir -L/path/to/library/dir -lgromacs -lnblib -o cpu_bench
 * 3) LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/library/dir
 * 4) ./cpu_bench OnlyNonBondedInteractions.tpr 10
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include <cstdio>

#include <chrono>

#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/gmxcalculatorgpu.h"
#include "nblib/nblib.h"

void printX(gmx::ArrayRef<gmx::RVec> c)
{
    for (int i = 0; i < std::min(20, int(c.size())); ++i)
    {
        printf("%4f %4f %4f\n", c[i][0], c[i][1], c[i][2]);
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        printf("You must supply a tpr file and a number of steps.\n");
        exit(1);
    }
    const std::string filepath   = argv[1];
    const float       dt         = std::stof(argv[2]);
    const int         numSteps   = std::stoi(argv[3]);
    const int         numThreads = std::stoi(argv[4]);

    nblib::TprReader tpr(filepath); // will exit if invalid file supplied
    const int        numParticles = tpr.coordinates_.size();

    // Force buffer initialization for each particle.
    std::vector<nblib::Vec3> forces(numParticles, { 0.0000, 0.0000, 0.0000 });

    // Kernel options are flags needed for force calculation.
    nblib::NBKernelOptions options = nblib::NBKernelOptions();
    options.numOpenMPThreads       = numThreads;
    // Use a simple cutoff rule for Coulomb
    options.coulombType = nblib::CoulombType::Cutoff;

    // Some performance flags can be set a run time
    std::vector<std::unique_ptr<DeviceInformation>> devices    = findDevices();
    const DeviceInformation&                        deviceInfo = *devices[0].get();

    nblib::GmxNBForceCalculatorGpu forceCalculator =
            nblib::GmxNBForceCalculatorGpu(tpr.particleTypeIdOfAllParticles_,
                                           tpr.nonbondedParameters_,
                                           tpr.charges_,
                                           tpr.particleInteractionFlags_,
                                           tpr.exclusionListRanges_,
                                           tpr.exclusionListElements_,
                                           options,
                                           deviceInfo);

    auto t1 = std::chrono::high_resolution_clock::now();
    // Integration requires masses, positions, and forces
    nblib::LeapFrog integrator(tpr.inverseMasses_, tpr.getBox());

    std::cout << "initial coordinates (first 20)" << std::endl;
    printX(tpr.coordinates_);

    for (int step = 0; step < numSteps; step++)
    {
        if (step % 20 == 0)
        {
            forceCalculator.updatePairlist(tpr.coordinates_, tpr.getBox());
        }
        nblib::zeroCartesianArray(forces);
        forceCalculator.compute(tpr.coordinates_, tpr.getBox(), forces);
        integrator.integrate(dt, tpr.coordinates_, tpr.velocities_, forces, numThreads);
    }

    auto   t2       = std::chrono::high_resolution_clock::now();
    double wallTime = std::chrono::duration<double>(t2 - t1).count();
    double simTime  = dt * float(numSteps);
    printf("ns/day %0.4f\n\n", simTime * 24 * 3.6 / wallTime);

    std::cout << "final coordinates (first 20)" << std::endl;
    printX(tpr.coordinates_);

    return 0;
}
