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
 * \brief Translation layer to GROMACS data structures for force calculations.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_NBNXMSETUPHELPERS_H
#define NBLIB_NBNXMSETUPHELPERS_H

#include "nblib/interactions.h"
#include "nblib/kerneloptions.h"
#include "nblib/particletype.h"

struct nonbonded_verlet_t;
struct t_forcerec;
struct t_nrnb;
struct interaction_const_t;
struct gmx_enerdata_t;
struct DeviceInformation;

namespace gmx
{
class StepWorkload;
class SimulationWorkload;
class DeviceStreamManager;
template<class T>
class ArrayRef;
} // namespace gmx

namespace Nbnxm
{
enum class KernelType;
struct KernelSetup;
} // namespace Nbnxm

namespace nblib
{

/*! \brief determine number of energy groups in an array of particle info flags
 *
 * Note: If the maximum energy group ID in the input is N, it is assumed that
 * all the energy groups with IDs from 0...N-1 also exist.
 */
int64_t findNumEnergyGroups(gmx::ArrayRef<int64_t> particleInteractionFlags);

//! Helper to translate between the different enumeration values.
Nbnxm::KernelType translateBenchmarkEnum(const SimdKernels& kernel);

/*! \brief Checks the kernel SIMD setup in CPU case
 *
 * Throws an exception when the kernel is not available.
 */
void checkKernelSetupSimd(SimdKernels nbnxmSimd);

//! Creates and returns the kernel setup for CPU
Nbnxm::KernelSetup createKernelSetupCPU(const SimdKernels nbnxmSimd, const bool useTabulatedEwaldCorr);

//! Creates and returns the kernel setup for GPU
Nbnxm::KernelSetup createKernelSetupGPU(const bool useTabulatedEwaldCorr);

//! Create Particle info array to mark those that undergo VdV interaction
std::vector<int64_t> createParticleInfoAllVdw(size_t numParticles);

//! Create the non-bonded parameter vector in GROMACS format
std::vector<real> createNonBondedParameters(const std::vector<ParticleType>& particleTypes,
                                            const NonBondedInteractionMap& nonBondedInteractionMap);

//! Create a step work object
gmx::StepWorkload createStepWorkload();

//! Create a SimulationWorkload object for use with createDeviceStreamManager
gmx::SimulationWorkload createSimulationWorkloadGpu();

//! Create a DeviceStreamManager; could be shared among multiple force calculators
std::shared_ptr<gmx::DeviceStreamManager> createDeviceStreamManager(const DeviceInformation& deviceInfo,
                                                                    const gmx::SimulationWorkload& simulationWorkload);

//! Computes the Ewald splitting coefficient for Coulomb
real ewaldCoeff(real ewald_rtol, real pairlistCutoff);

//! Creates an interaction_const_t object from NBKernelOptions
interaction_const_t createInteractionConst(const NBKernelOptions& options);

//! Create nonbonded_verlet_t object
std::unique_ptr<nonbonded_verlet_t> createNbnxmCPU(size_t                    numParticleTypes,
                                                   const NBKernelOptions&    options,
                                                   int                       numEnergyGroups,
                                                   gmx::ArrayRef<const real> nonbondedParameters);

//! Create nonbonded_verlet_gpu object
std::unique_ptr<nonbonded_verlet_t> createNbnxmGPU(size_t                     numParticleTypes,
                                                   const NBKernelOptions&     options,
                                                   const std::vector<real>&   nonbondedParameters,
                                                   const interaction_const_t& interactionConst,
                                                   const gmx::DeviceStreamManager& deviceStreamManager);

//! Set number of OpenMP threads in the GROMACS backend
void setGmxNonBondedNThreads(int numThreads);

//! Update the shift vectors in t_forcerec
void updateForcerec(t_forcerec* forcerec, const matrix& box);

} // namespace nblib

#endif // NBLIB_NBNXMSETUPHELPERS_H
