/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Implements the CoordinateState class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "coordinatestate.h"

#include <algorithm>

#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"

#include "grid.h"
#include "math.h"

namespace gmx
{

CoordinateState::CoordinateState(const AwhBiasParams          &awhBiasParams,
                                 const std::vector<DimParams> &dimParams,
                                 const Grid                   &grid)
{
    for (size_t d = 0; d < dimParams.size(); d++)
    {
        coordValue_[d] = dimParams[d].scaleUserInputToInternal(awhBiasParams.dimParams[d].coordValueInit);
    }

    /* The grid-point index is always the nearest point to the coordinate.
     * We initialize the umbrella location to the same nearest point.
     * More correctly one would sample from the biased distribution,
     * but it doesn't really matter, as this will happen after a few steps.
     */
    gridpointIndex_    = grid.nearestIndex(coordValue_);
    umbrellaGridpoint_ = gridpointIndex_;
}

/* Sample a new umbrella reference point given the current coordinate value. */
void
CoordinateState::sampleUmbrellaGridpoint(const Grid                &grid,
                                         int                        gridpointIndex,
                                         const std::vector<double> &probWeightNeighbor,
                                         gmx_int64_t                step,
                                         gmx_int64_t                seed,
                                         int                        indexSeed)
{
    /* Sample new umbrella reference value from the probability distribution
     * which is defined for the neighboring points of the current coordinate.
     */
    const std::vector<int> &neighbor = grid.point(gridpointIndex).neighbor;

    /* In order to use the same seed for all AWH biases and get independent
       samples we use the index of the bias. */
    int n_sampled      = getSampleFromDistribution(probWeightNeighbor,
                                                   step, seed, indexSeed);

    umbrellaGridpoint_ = neighbor[n_sampled];
}

/* Update the coordinate value with coordValue. */
void CoordinateState::setCoordValue(const Grid &grid,
                                    int         dim,
                                    double      coordValue)
{
    coordValue_[dim] = coordValue;

    /* The grid point closest to the coordinate value defines the current
     * neighborhood of points. Besides at steps when global updates and/or
     * checks are performed, only the neighborhood will be touched.
     */
    gridpointIndex_ = grid.nearestIndex(coordValue_);
}

void
CoordinateState::restoreFromHistory(const AwhBiasStateHistory &stateHistory)
{
    umbrellaGridpoint_ = stateHistory.umbrellaGridpoint;
}

} // namespace gmx
