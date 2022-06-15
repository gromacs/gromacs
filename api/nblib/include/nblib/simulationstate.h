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
/*! \inpublicapi \file
 * \brief
 * Implements nblib SimulationState
 *
 * \author Berk Hess <hess@kth.se>
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_SIMULATIONSTATE_H
#define NBLIB_SIMULATIONSTATE_H

#include <memory>
#include <vector>

#include "nblib/box.h"
#include "nblib/topology.h"
#include "nblib/vector.h"

namespace nblib
{

/*! \libinternal
 * \ingroup nblib
 * \brief Simulation State
 *
 * Simulation state description that serves as a snapshot of the system
 * being analysed. Needed to init an MD program. Allows hot-starting simulations.
 */

class SimulationState final
{
public:
    //! Constructor
    SimulationState(const std::vector<Vec3>& coordinates,
                    const std::vector<Vec3>& velocities,
                    const std::vector<Vec3>& forces,
                    Box                      box,
                    Topology                 topology);

    //! Returns topology of the current state
    const Topology& topology() const;

    //! Returns the box
    Box box() const;

    //! Returns a reference to a (modifiable) vector of particle coordinates
    std::vector<Vec3>& coordinates();

    //! Returns a read-only vector of particle coordinates
    const std::vector<Vec3>& coordinates() const;

    //! Returns a reference to a (modifiable) vector of particle velocities
    std::vector<Vec3>& velocities();

    //! Returns a reference to a (modifiable) vector of forces
    std::vector<Vec3>& forces();

private:
    class Impl;
    std::shared_ptr<SimulationState::Impl> simulationStatePtr_;
};

} // namespace nblib

#endif // NBLIB_SIMULATIONSTATE_H
