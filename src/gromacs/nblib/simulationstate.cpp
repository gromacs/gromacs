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
 */
#include "gmxpre.h"

#include "simulationstate.h"

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nblib/atomtype.h"
#include "gromacs/nblib/util.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/exceptions.h"

#include "coords.h"

namespace nblib
{

SimulationState::SimulationState(const std::vector<gmx::RVec>& coord,
                                 Box                           box,
                                 Topology                      topology,
                                 const std::vector<gmx::RVec>& vel) :
    box_(std::move(box)),
    topology_(std::move(topology))
{
    if (!checkNumericValues(coord))
    {
        GMX_THROW(gmx::InvalidInputError("Input coordinates has at least one NaN"));
    }
    coordinates_ = coord;
    if (!checkNumericValues(vel))
    {
        GMX_THROW(gmx::InvalidInputError("Input velocities has at least one NaN"));
    }
    velocities_ = vel;
}

const Topology& SimulationState::topology() const
{
    return topology_;
}

const Box& SimulationState::box()
{
    return box_;
}

std::vector<gmx::RVec>& SimulationState::coordinates()
{
    return coordinates_;
}

std::vector<gmx::RVec>& SimulationState::velocities()
{
    return velocities_;
}


} // namespace nblib
