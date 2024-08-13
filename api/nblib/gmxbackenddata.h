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
/*! \inpublicapi \file
 * \brief
 * Implements nblib simulation box
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_GMXBACKENDDATA_H
#define NBLIB_GMXBACKENDDATA_H

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/range.h"

#include "nblib/kerneloptions.h"
#include "nblib/nbnxmsetuphelpers.h"

namespace nblib
{

/*! \brief GROMACS non-bonded force calculation backend
 *
 * This class encapsulates the various GROMACS data structures and their interplay
 * from the NB-LIB user. The class is a private member of the ForceCalculator and
 * is not intended for the public interface.
 *
 * Handles the task of storing the simulation problem description using the internal
 * representation used within GROMACS. It currently supports short range non-bonded
 * interactions (PP) on a single node CPU.
 *
 */
class GmxBackendData
{
public:
    GmxBackendData() = default;
    GmxBackendData(const NBKernelOptions& options,
                   int                    numEnergyGroups,
                   gmx::ArrayRef<int>     exclusionRanges,
                   gmx::ArrayRef<int>     exclusionElements) :
        numThreads_(options.numOpenMPThreads)
    {
        // Set hardware params from the execution context
        setGmxNonBondedNThreads(options.numOpenMPThreads);

        // Set interaction constants struct
        interactionConst_ = createInteractionConst(options);

        // Set up StepWorkload data
        stepWork_ = createStepWorkload();

        // Set up gmx_enerdata_t (holds energy information)
        enerd_ = gmx_enerdata_t{ numEnergyGroups, nullptr };

        // Construct pair lists
        std::vector<int> exclusionRanges_(exclusionRanges.begin(), exclusionRanges.end());
        std::vector<int> exclusionElements_(exclusionElements.begin(), exclusionElements.end());
        exclusions_ = gmx::ListOfLists<int>(std::move(exclusionRanges_), std::move(exclusionElements_));
    }

    //! exclusions in gmx format
    gmx::ListOfLists<int> exclusions_;

    //! Non-Bonded Verlet object for force calculation
    std::unique_ptr<gmx::nonbonded_verlet_t> nbv_;

    //! Only shift_vec is used
    t_forcerec forcerec_;

    //! Parameters for various interactions in the system
    interaction_const_t interactionConst_;

    //! Tasks to perform in an MD Step
    gmx::StepWorkload stepWork_;

    gmx::SimulationWorkload simulationWork_;

    //! Energies of different interaction types; currently only needed as an argument for dispatchNonbondedKernel
    gmx_enerdata_t enerd_{ 1, nullptr };

    //! Non-bonded flop counter; currently only needed as an argument for dispatchNonbondedKernel
    t_nrnb nrnb_;

    //! Number of OpenMP threads to use
    int numThreads_;

    //! Keep track of whether updatePairlist has been called at least once
    bool updatePairlistCalled{ false };
};

} // namespace nblib
#endif // NBLIB_GMXBACKENDDATA_H
