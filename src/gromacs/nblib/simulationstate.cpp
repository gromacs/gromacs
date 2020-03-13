/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * Implements nblib SimulationState
 *
 * \author Berk Hess <hess@kth.se>
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "simulationstate.h"

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/util.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/exceptions.h"

#include "coords.h"

namespace nblib
{

SimulationState::Impl::Impl(const std::vector<gmx::RVec>& coordinates,
                            const std::vector<gmx::RVec>& velocities,
                            const std::vector<gmx::RVec>& forces,
                            Box                           box,
                            Topology                      topology) :
    box_(std::move(box)),
    topology_(std::move(topology))
{
    if (!checkNumericValues(coordinates))
    {
        GMX_THROW(gmx::InvalidInputError("Input coordinates has at least one NaN"));
    }
    coordinates_ = coordinates;
    if (!checkNumericValues(velocities))
    {
        GMX_THROW(gmx::InvalidInputError("Input velocities has at least one NaN"));
    }

    velocities_ = velocities;

    forces_ = forces;

    velocities_ = velocities;

    // Ensure that the coordinates are in a box following PBC
    put_atoms_in_box(PbcType::Xyz, box.legacyMatrix(), coordinates_);
}

const Topology& SimulationState::Impl::topology() const
{
    return topology_;
}

const Box& SimulationState::Impl::box()
{
    return box_;
}

std::vector<gmx::RVec>& SimulationState::Impl::coordinates()
{
    return coordinates_;
}

std::vector<gmx::RVec>& SimulationState::Impl::velocities()
{
    return velocities_;
}

std::vector<gmx::RVec>& SimulationState::Impl::forces()
{
    return forces_;
}

const Topology& SimulationState::topology() const
{
    return simulationStatePtr_->topology();
}

const Box SimulationState::box() const
{
    return simulationStatePtr_->box();
}

std::vector<gmx::RVec>& SimulationState::coordinates()
{
    return simulationStatePtr_->coordinates();
}

const std::vector<gmx::RVec>& SimulationState::coordinates() const
{
    return simulationStatePtr_->coordinates();
}

std::vector<gmx::RVec>& SimulationState::velocities()
{
    return simulationStatePtr_->velocities();
}

std::vector<gmx::RVec>& SimulationState::forces()
{
    return simulationStatePtr_->forces();
}

} // namespace nblib
