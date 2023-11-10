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
 * Compilation instructions
 * 1) install gromacs legacy api and nblib api
 * 2) clang++ -std=c++17 -I/path/to/include/dir -I/path/to/cmake/include/dir -I/path/to/cuda/include/dir -L/path/to/library/dir -L/path/to/cuda/libdir -lgromacs -lnblib -lcudart -o gpu_bench
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
#include <iostream>
#include <memory>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/path.h"

#include "buildinfo.h"

// The entire nblib public API can be included with a single header or individual components
// can be included via their respective headers.
#include "nblib/gmxcalculatorgpu.h"
#include "nblib/integrator_gpu.h"
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

    nblib::TprReader tpr(filepath);
    const int        numParticles = tpr.coordinates_.size();

    // Force buffer initialization for each particle.
    std::vector<nblib::Vec3> forces(numParticles, { 0.0000, 0.0000, 0.0000 });

    // Kernel options are flags needed for force calculation.
    nblib::NBKernelOptions options = nblib::NBKernelOptions();
    // Use a simple cutoff rule for Coulomb
    options.coulombType = nblib::CoulombType::Cutoff;

    std::vector<std::unique_ptr<DeviceInformation>> devices    = findDevices();
    const DeviceInformation&                        deviceInfo = *devices[0].get();
    // The force calculator contains all the data needed to compute forces.
    nblib::GmxNBForceCalculatorGpu forceCalculator =
            nblib::GmxNBForceCalculatorGpu(tpr.particleTypeIdOfAllParticles_,
                                           tpr.nonbondedParameters_,
                                           tpr.charges_,
                                           tpr.particleInteractionFlags_,
                                           tpr.exclusionListRanges_,
                                           tpr.exclusionListElements_,
                                           options,
                                           deviceInfo);

    const DeviceContext& context = forceCalculator.deviceContext();

    // host coordinate buffers
    std::vector<gmx::RVec> x_h(tpr.coordinates_);

    // host velocity buffers
    std::vector<gmx::RVec> v_h(tpr.velocities_);

    // build the first pairlist on the CPU so that we can reorder the inverse masses
    forceCalculator.updatePairlist(x_h, tpr.getBox());
    std::size_t deviceBufferSize = forceCalculator.nbnxmBufferSize();

    // nbnxm reorder auxiliary host buffers
    std::vector<real>      xNbnxm4wide(4 * deviceBufferSize);
    std::vector<gmx::RVec> xNbnxm(deviceBufferSize);
    std::vector<real>      vNbnxm4Wide(4 * deviceBufferSize);
    std::vector<gmx::RVec> vNbnxm(deviceBufferSize, gmx::RVec{ 0, 0, 0 });

    // device coordinate buffer in XYZQ format
    DeviceBuffer<Float4> deviceXq;
    allocateDeviceBuffer(&deviceXq, deviceBufferSize, context);

    // device velocity buffer in XYZ format
    DeviceBuffer<Float3> deviceVelocities;
    allocateDeviceBuffer(&deviceVelocities, deviceBufferSize, context);

    // device force buffer in XYZ format
    DeviceBuffer<Float3> deviceForces;
    allocateDeviceBuffer(&deviceForces, deviceBufferSize, context);

    std::vector<real> inverseMassesNbnxm(deviceBufferSize);
    nblib::reorderScalarArray(forceCalculator, tpr.inverseMasses_, inverseMassesNbnxm);

    nblib::LeapFrogGPU integratorGpu(deviceBufferSize,
                                     inverseMassesNbnxm,
                                     forceCalculator.deviceContext(),
                                     forceCalculator.deviceStream());

    std::cout << "initial coordinates (first 20)" << std::endl;
    printX(x_h);

    auto t1 = std::chrono::high_resolution_clock::now();

    // update pairlist every 20 MD steps
    int pairlistUpdateFreq = 20;

    for (int step = 0; step < numSteps; step++)
    {
        if (step % pairlistUpdateFreq == 0)
        {
            // build the pairlist on the CPU
            forceCalculator.updatePairlist(x_h, tpr.getBox());

            // reorder and copy coordinates in nbnxm order to GPU
            forceCalculator.reorder(x_h, xNbnxm4wide);
            copyToDeviceBuffer(&deviceXq,
                               reinterpret_cast<const Float4*>(xNbnxm4wide.data()),
                               0,
                               deviceBufferSize,
                               forceCalculator.deviceStream(),
                               GpuApiCallBehavior::Sync,
                               nullptr);

            // reorder and copy velocities in nbnxm order to GPU
            forceCalculator.reorder(v_h, vNbnxm4Wide);
#pragma omp parallel for num_threads(numThreads)
            for (int i = 0; i < deviceBufferSize; ++i)
            {
                vNbnxm[i][0] = vNbnxm4Wide[4 * i + 0];
                vNbnxm[i][1] = vNbnxm4Wide[4 * i + 1];
                vNbnxm[i][2] = vNbnxm4Wide[4 * i + 2];
            }
            copyToDeviceBuffer(&deviceVelocities,
                               reinterpret_cast<const Float3*>(vNbnxm.data()),
                               0,
                               deviceBufferSize,
                               forceCalculator.deviceStream(),
                               GpuApiCallBehavior::Sync,
                               nullptr);
        }

        clearDeviceBufferAsync(&deviceForces, 0, deviceBufferSize, forceCalculator.deviceStream());
        forceCalculator.compute(deviceXq, tpr.getBox(), deviceForces);
        integratorGpu.integrate(dt, deviceXq, deviceVelocities, deviceForces);


        // download coordinates and velocities for next pairlist update
        if (step % pairlistUpdateFreq == pairlistUpdateFreq - 1)
        {
            // download coordinates
            copyFromDeviceBuffer(reinterpret_cast<Float4*>(xNbnxm4wide.data()),
                                 &deviceXq,
                                 0,
                                 deviceBufferSize,
                                 forceCalculator.deviceStream(),
                                 GpuApiCallBehavior::Sync,
                                 nullptr);
// throw away Q
#pragma omp parallel for num_threads(numThreads)
            for (int i = 0; i < deviceBufferSize; ++i)
            {
                xNbnxm[i] = { xNbnxm4wide[4 * i], xNbnxm4wide[4 * i + 1], xNbnxm4wide[4 * i + 2] };
            }
            // reorder coordinates into default ordering
            forceCalculator.undoReorder(xNbnxm, x_h);

            // download velocities
            copyFromDeviceBuffer(reinterpret_cast<Float3*>(vNbnxm.data()),
                                 &deviceVelocities,
                                 0,
                                 deviceBufferSize,
                                 forceCalculator.deviceStream(),
                                 GpuApiCallBehavior::Sync,
                                 nullptr);
            // reorder velocities into default ordering
            forceCalculator.undoReorder(vNbnxm, v_h);
        }
    }

    auto   t2       = std::chrono::high_resolution_clock::now();
    double wallTime = std::chrono::duration<double>(t2 - t1).count();
    double simTime  = dt * float(numSteps);
    printf("ns/day %0.4f\n\n", simTime * 24 * 3.6 / wallTime);

    // download final coordinates
    copyFromDeviceBuffer(reinterpret_cast<Float4*>(xNbnxm4wide.data()),
                         &deviceXq,
                         0,
                         deviceBufferSize,
                         forceCalculator.deviceStream(),
                         GpuApiCallBehavior::Sync,
                         nullptr);
// throw away Q
#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < deviceBufferSize; ++i)
    {
        xNbnxm[i] = { xNbnxm4wide[4 * i], xNbnxm4wide[4 * i + 1], xNbnxm4wide[4 * i + 2] };
    }
    // reorder coordinates into default ordering
    forceCalculator.undoReorder(xNbnxm, x_h);

    std::cout << "final coordinates (first 20) (DeviceBuffer interface)" << std::endl;
    printX(x_h);

    freeDeviceBuffer(&deviceXq);
    freeDeviceBuffer(&deviceVelocities);
    freeDeviceBuffer(&deviceForces);

    return 0;
}
