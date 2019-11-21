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

/*!\file
 * \brief
 * This file declares functions for setting up a nonbonded system
 *
 * \author Berk Hess <hess@kth.se>
 * \author Joe Jordan <ejjordan@kth.se>
 */

#ifndef GROMACS_SIMULATIONSTATE_H
#define GROMACS_SIMULATIONSTATE_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"

#include "box.h"
#include "topology.h"
#include "util.h"

namespace nblib
{

//! Simulation state description that serves as a snapshot of the system
//! being analysed. Needed to init an MD program. Allows hot-starting simulations.
class SimulationState
{
public:
    // Constructor
    SimulationState(const std::vector<gmx::RVec> &coord, Box box, Topology &topo,
             const std::vector<gmx::RVec> &vel = {});

    // Copy Constructor
    SimulationState(const SimulationState& simulationState);

    // Copy Assignment Operator
    SimulationState& operator=(const SimulationState& simulationState);

    // Move Constructor
    SimulationState(SimulationState&& simulationState) noexcept;

    // Move Assignment Constructor
    SimulationState& operator=(SimulationState&& simulationState) noexcept;

private:
    std::vector<gmx::RVec> coord_;
    Box box_;
    Topology topo_;
    std::vector<gmx::RVec> vel_;

};

} // namespace nblib

#endif // GROMACS_SIMULATIONSTATE_H
