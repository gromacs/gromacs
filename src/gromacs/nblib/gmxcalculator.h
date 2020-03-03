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
 * Implements a force calculator based on GROMACS data structures.
 *
 * Intended for internal use inside the ForceCalculator.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GROMACS_GMXCALCULATOR_H
#define GROMACS_GMXCALCULATOR_H

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm.h"

namespace nblib
{

class SimulationState;
struct NBKernelOptions;

//! Set up StepWorkload data
gmx::StepWorkload setupStepWorkload(std::shared_ptr<NBKernelOptions> options);

//! Return an interaction constants struct with members set appropriately
interaction_const_t setupInteractionConst(std::shared_ptr<NBKernelOptions> options);

class GmxForceCalculator
{
public:
    explicit GmxForceCalculator(SimulationState simState, std::shared_ptr<NBKernelOptions> options);

    //! Compute forces and return
    gmx::PaddedHostVector<gmx::RVec> compute();

    //! Non-Bonded Verlet object for force calculation
    std::unique_ptr<nonbonded_verlet_t> nbv_;

    //! Only nbfp, shift_vec and ntypes are used
    t_forcerec forcerec_;

    //! Contains array for computed forces
    gmx::PaddedHostVector<gmx::RVec> verletForces_ = {};

    //! Parameters for various interactions in the system
    interaction_const_t interactionConst_;

    //! Tasks to perform in an MD Step
    gmx::StepWorkload stepWork_;

private:
    //! Energies of different interaction types; currently only needed as an argument for dispatchNonbondedKernel
    gmx_enerdata_t enerd_;

    //! Non-bonded flop counter; currently only needed as an argument for dispatchNonbondedKernel
    t_nrnb nrnb_ = { 0 };

    //! Legacy matrix for box
    matrix box_;
};

} // namespace nblib

#endif // GROMACS_GMXCALCULATOR_H
