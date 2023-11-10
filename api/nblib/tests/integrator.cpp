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
#include "nblib/integrator.h"

#include "gromacs/hardware/device_management.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/molecules.h"
#include "nblib/particletype.h"
#include "nblib/simulationstate.h"
#include "nblib/tests/testhelpers.h"
#include "nblib/tests/testsystems.h"
#include "nblib/topology.h"
#include "nblib/vector.h"
#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/devicebuffer_datatype.h"

#    include "nblib/integrator_gpu.h"
#endif

#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, IntegratorWorks)
{
    int  numAtoms = 1;
    int  numSteps = 100;
    real dt       = 0.001;

    ParticleType particleType(ParticleTypeName("H"), Mass(1.0));
    Molecule     molecule(MoleculeName("SomeMolecule"));
    molecule.addParticle(ParticleName("SomeAtom"), particleType);

    ParticleTypesInteractions interactions;
    interactions.add(particleType.name(), C6{ 0 }, C12{ 0 });

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(molecule, numAtoms);
    topologyBuilder.addParticleTypesInteractions(interactions);
    Topology topology = topologyBuilder.buildTopology();

    std::vector<Vec3> x(numAtoms, { 0.0, 0.0, 0.0 });
    std::vector<Vec3> v(numAtoms, { 0.0, 0.0, 0.0 });
    std::vector<Vec3> f(numAtoms, { 1.0, 2.0, 0.0 });

    Box box(100);

    std::vector<Vec3> x0(x);
    std::vector<Vec3> v0(v);

    SimulationState simulationState(x, v, f, box, topology);
    put_atoms_in_box(PbcType::Xyz, box.legacyMatrix(), x0);

    LeapFrog integrator(simulationState.topology(), simulationState.box());

    gmx::test::FloatingPointTolerance tolerance = gmx::test::absoluteTolerance(numSteps * 0.000005);
    for (int step = 0; step < numSteps; step++)
    {
        real totalTime = step * dt;

        Vec3 xAnalytical;
        Vec3 vAnalytical;

        for (int i = 0; i < numAtoms; i++)
        {
            for (int d = 0; d < dimSize; d++)
            {
                // Analytical solution for constant-force particle movement
                int  typeIndex = simulationState.topology().getParticleTypeIdOfAllParticles()[i];
                real im = 1.0 / simulationState.topology().getParticleTypes()[typeIndex].mass();
                xAnalytical[d] =
                        x0[i][d] + v0[i][d] * totalTime + 0.5 * f[i][d] * totalTime * totalTime * im;
                vAnalytical[d] = v0[i][d] + f[i][d] * totalTime * im;

                EXPECT_REAL_EQ_TOL(xAnalytical[d], simulationState.coordinates()[i][d], tolerance)
                        << formatString(
                                   "Coordinate {} of atom {} is different from analytical solution "
                                   "at step {}.",
                                   d,
                                   i,
                                   step);

                EXPECT_REAL_EQ_TOL(vAnalytical[d], simulationState.velocities()[i][d], tolerance)
                        << formatString(
                                   "Velocity component {} of atom {} is different from analytical "
                                   "solution at step {}.",
                                   d,
                                   i,
                                   step);
            }
            integrator.integrate(dt,
                                 simulationState.coordinates(),
                                 simulationState.velocities(),
                                 simulationState.forces());
        }
    }
}

#if GMX_GPU_CUDA
TEST(NBlibTest, CanSetupGpuIntegrator)
{
    ArgonSimulationStateBuilder argonSystemBuilder(fftypes::GROMOS43A1);
    SimulationState             simState = argonSystemBuilder.setupSimulationState();
    std::vector<real>           masses(simState.topology().numParticles());
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    for (const auto& testDevice : testDeviceList)
    {
        const DeviceInformation& deviceInfo    = testDevice->deviceInfo();
        const DeviceContext&     deviceContext = testDevice->deviceContext();
        const DeviceStream&      deviceStream  = testDevice->deviceStream();
        setActiveDevice(deviceInfo);
        EXPECT_NO_THROW(LeapFrogGPU(simState.topology().numParticles(), masses, deviceContext, deviceStream));
    }
}

TEST(NBlibTest, GPUIntegratorWorks)
{
    int  numAtoms = 1;
    int  numSteps = 100;
    real dt       = 0.001;

    ParticleType particleType(ParticleTypeName("H"), Mass(1.0));
    Molecule     molecule(MoleculeName("SomeMolecule"));
    molecule.addParticle(ParticleName("SomeAtom"), particleType);

    ParticleTypesInteractions interactions;
    interactions.add(particleType.name(), C6{ 0 }, C12{ 0 });

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(molecule, numAtoms);
    topologyBuilder.addParticleTypesInteractions(interactions);
    Topology    topology       = topologyBuilder.buildTopology();
    const auto& testDeviceList = gmx::test::getTestHardwareEnvironment()->getTestDeviceList();
    for (const auto& testDevice : testDeviceList)
    {
        const DeviceInformation& deviceInfo    = testDevice->deviceInfo();
        const DeviceContext&     deviceContext = testDevice->deviceContext();
        const DeviceStream&      deviceStream  = testDevice->deviceStream();
        setActiveDevice(deviceInfo);

        std::vector<Float4> xq_h(numAtoms, { 0.0, 0.0, 0.0, 0.0 });
        std::vector<Vec3>   v_h(numAtoms, { 0.0, 0.0, 0.0 });
        std::vector<Vec3>   f_h(numAtoms, { 1.0, 2.0, 0.0 });

        DeviceBuffer<Float4> xq_d;
        DeviceBuffer<Float3> v_d;
        DeviceBuffer<Float3> f_d;
        allocateDeviceBuffer(&xq_d, numAtoms, deviceContext);
        allocateDeviceBuffer(&v_d, numAtoms, deviceContext);
        allocateDeviceBuffer(&f_d, numAtoms, deviceContext);

        Box box(100);

        std::vector<Float4> xq0(xq_h);
        std::vector<Vec3>   v0(v_h);

        std::vector<Mass> masses = expandQuantity(topology, &ParticleType::mass);
        std::vector<real> inverseMasses(masses.size());
        std::transform(masses.begin(), masses.end(), inverseMasses.begin(), [](real mass) {
            return 1.0 / mass;
        });
        LeapFrogGPU integrator(topology.numParticles(), inverseMasses, deviceContext, deviceStream);

        // copy initial values for x,v,f to the device
        copyToDeviceBuffer(&xq_d,
                           reinterpret_cast<const Float4*>(xq_h.data()),
                           0,
                           numAtoms,
                           deviceStream,
                           GpuApiCallBehavior::Sync,
                           nullptr);
        copyToDeviceBuffer(&v_d,
                           reinterpret_cast<const Float3*>(v_h.data()),
                           0,
                           numAtoms,
                           deviceStream,
                           GpuApiCallBehavior::Sync,
                           nullptr);
        copyToDeviceBuffer(&f_d,
                           reinterpret_cast<const Float3*>(f_h.data()),
                           0,
                           numAtoms,
                           deviceStream,
                           GpuApiCallBehavior::Sync,
                           nullptr);

        gmx::test::FloatingPointTolerance tolerance = gmx::test::absoluteTolerance(numSteps * 0.000005);
        for (int step = 0; step < numSteps; step++)
        {
            real totalTime = step * dt;

            // download x_d, v_d, from current time step into x_h, v_h
            copyFromDeviceBuffer(reinterpret_cast<Float4*>(xq_h.data()),
                                 &xq_d,
                                 0,
                                 numAtoms,
                                 deviceStream,
                                 GpuApiCallBehavior::Sync,
                                 nullptr);
            copyFromDeviceBuffer(reinterpret_cast<Float3*>(v_h.data()),
                                 &v_d,
                                 0,
                                 numAtoms,
                                 deviceStream,
                                 GpuApiCallBehavior::Sync,
                                 nullptr);

            // compare x_h, v_h, against analytical reference
            for (int i = 0; i < numAtoms; i++)
            {
                // for (int d = 0; d < dimSize; d++)
                // x-component (float4 does not support subscripting)
                {
                    // Analytical solution for constant-force particle movement
                    int  typeIndex   = topology.getParticleTypeIdOfAllParticles()[i];
                    real im          = 1.0 / topology.getParticleTypes()[typeIndex].mass();
                    real xAnalytical = xq0[i].x + v0[i][0] * totalTime
                                       + 0.5 * f_h[i][0] * totalTime * totalTime * im;
                    real vAnalytical = v0[i][0] + f_h[i][0] * totalTime * im;

                    EXPECT_REAL_EQ_TOL(xAnalytical, xq_h[i].x, tolerance) << formatString(
                            "Coordinate {} of atom {} is different from analytical solution "
                            "at step {}.",
                            0,
                            i,
                            step);

                    EXPECT_REAL_EQ_TOL(vAnalytical, v_h[i][0], tolerance) << formatString(
                            "Velocity component {} of atom {} is different from analytical "
                            "solution at step {}.",
                            0,
                            i,
                            step);
                }
                // y-component
                {
                    // Analytical solution for constant-force particle movement
                    int  typeIndex   = topology.getParticleTypeIdOfAllParticles()[i];
                    real im          = 1.0 / topology.getParticleTypes()[typeIndex].mass();
                    real xAnalytical = xq0[i].y + v0[i][1] * totalTime
                                       + 0.5 * f_h[i][1] * totalTime * totalTime * im;
                    real vAnalytical = v0[i][1] + f_h[i][1] * totalTime * im;

                    EXPECT_REAL_EQ_TOL(xAnalytical, xq_h[i].y, tolerance) << formatString(
                            "Coordinate {} of atom {} is different from analytical solution "
                            "at step {}.",
                            1,
                            i,
                            step);

                    EXPECT_REAL_EQ_TOL(vAnalytical, v_h[i][1], tolerance) << formatString(
                            "Velocity component {} of atom {} is different from analytical "
                            "solution at step {}.",
                            1,
                            i,
                            step);
                }
            }

            integrator.integrate(dt, xq_d, v_d, f_d);
        }

        freeDeviceBuffer(&xq_d);
        freeDeviceBuffer(&v_d);
        freeDeviceBuffer(&f_d);
    }
}
#endif

} // namespace
} // namespace test
} // namespace nblib
