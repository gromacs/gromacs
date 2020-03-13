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
 * \brief Translation layer to GROMACS data structures for force calculations.
 *
 * Implements the translation layer between the user scope and
 * GROMACS data structures for force calculations. Sets up the
 * non-bonded verlet.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GROMACS_GMXSETUP_H
#define GROMACS_GMXSETUP_H

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nbnxm/nbnxm.h"

#include "gmxcalculator.h"
#include "interactions.h"
#include "nbkerneloptions.h"

namespace nblib
{

class NbvSetupUtil
{
public:
    NbvSetupUtil() : gmxForceCalculator_(std::make_unique<GmxForceCalculator>()) {}

    //! Sets hardware params from the execution context
    void setExecutionContext(const NBKernelOptions& options);

    //! Sets non-bonded parameters to be used to build GMX data structures
    void setNonBondedParameters(const std::vector<ParticleType>& particleTypes);

    //! Marks particles to have Van der Waals interactions
    void setParticleInfoAllVdv(size_t numParticles);

    //! Returns the kernel setup
    Nbnxm::KernelSetup getKernelSetup(const NBKernelOptions& options);

    //! Set up StepWorkload data
    void setupStepWorkload(const NBKernelOptions& options);

    //! Return an interaction constants struct with members set appropriately
    void setupInteractionConst(const NBKernelOptions& options);

    //! Sets Particle Types and Charges and VdW params
    void setAtomProperties(const Topology& topology);

    //! Sets up non-bonded verlet on the GmxForceCalculator
    void setupNbnxmInstance(const Topology& topology, const NBKernelOptions& options);

    //! Puts particles on a grid based on bounds specified by the box
    void setParticlesOnGrid(const std::vector<gmx::RVec>& coordinates,
                            const Box&                    box);

    //! Constructs pair lists
    void constructPairList(const gmx::ListOfLists<int>& exclusions);

    //! Sets up t_forcerec object on the GmxForceCalculator
    void setupForceRec(const matrix& box);

    //! Sets initial forces to zero
    void setForcesToZero(size_t numParticles);

    std::unique_ptr<GmxForceCalculator> getGmxForceCalculator()
    {
        return std::move(gmxForceCalculator_);
    }

private:

    //! Storage for parameters for short range interactions.
    std::vector<real> nonbondedParameters_;

    //! Particle info where all particles are marked to have Van der Waals interactions
    std::vector<int> particleInfoAllVdw_;

    //! GROMACS force calculator to compute forces
    std::unique_ptr<GmxForceCalculator> gmxForceCalculator_;
};

class GmxSetupDirector
{
public:
    //! Sets up and returns a GmxForceCalculator
    static std::unique_ptr<GmxForceCalculator> setupGmxForceCalculator(const SimulationState& system,
                                                                const NBKernelOptions& options);
};

} // namespace nblib
#endif // GROMACS_GMXSETUP_H
