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
 * \brief
 * Implements the translation layer between the user scope and
 * GROMACS data structures for force calculations
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GROMACS_GMXSETUP_H
#define GROMACS_GMXSETUP_H

#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm.h"

#include "nbkerneloptions.h"

namespace nblib
{

class ForceCalculator;
class SimulationState;
struct NbvSetupUtil;

enum class CombinationRule : int
{
    Geometric = 0,
    Count     = 1
};

struct GmxForceCalculator
{
    //! Parameters for various interactions in the system
    interaction_const_t interactionConst_;

    //! Energies of different interaction types
    gmx_enerdata_t enerd_;

    //! Non-Bonded Verlet object for force calculation
    std::unique_ptr<nonbonded_verlet_t> nbv_;

    //! The massive class from which nbfp, shift_vec and ntypes would be used
    t_forcerec forcerec_;

    //! Tasks to perform in an MD Step
    gmx::StepWorkload stepWork_;

    explicit GmxForceCalculator(const std::shared_ptr<SimulationState> system,
                                const std::shared_ptr<NBKernelOptions> options);

    //! Contains array for computed forces
    gmx::PaddedHostVector<gmx::RVec> verletForces_;

    //! Compute forces and return
    gmx::PaddedHostVector<gmx::RVec> compute();

    //! Legacy matrix for box
    matrix box_;
};

struct NbvSetupUtil
{
    NbvSetupUtil(SimulationState system, const NBKernelOptions& options);

    void unpackTopologyToGmx();

    std::unique_ptr<nonbonded_verlet_t> setupNbnxmInstance();

    std::unique_ptr<GmxForceCalculator> setupGmxForceCalculator();

    std::shared_ptr<SimulationState> system_;
    std::shared_ptr<NBKernelOptions> options_;

    //! Storage for parameters for short range interactions.
    std::vector<real> nonbondedParameters_;

    //! Particle info where all particles are marked to have Van der Waals interactions
    std::vector<int> particleInfoAllVdw_;
};

} // namespace nblib
#endif // GROMACS_GMXSETUP_H
