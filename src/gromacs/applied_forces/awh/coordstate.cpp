/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * \brief
 * Implements the CoordState class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "coordstate.h"

#include <algorithm>
#include <string>

#include "gromacs/applied_forces/awh/dimparams.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "biasgrid.h"

namespace gmx
{

CoordState::CoordState(const AwhBiasParams&      awhBiasParams,
                       ArrayRef<const DimParams> dimParams,
                       const BiasGrid&           grid)
{
    GMX_RELEASE_ASSERT(awhBiasParams.ndim() == dimParams.ssize(),
                       "Need to have identical size for dimensions");
    const auto& awhDimParams = awhBiasParams.dimParams();
    for (int d = 0; d < gmx::ssize(awhDimParams); d++)
    {
        coordValue_[d] = dimParams[d].scaleUserInputToInternal(awhDimParams[d].initialCoordinate());
    }

    /* The grid-point index is always the nearest point to the coordinate.
     * We initialize the umbrella location to the same nearest point.
     * More correctly one would sample from the biased distribution,
     * but it doesn't really matter, as this will happen after a few steps.
     */
    gridpointIndex_    = grid.nearestIndex(coordValue_);
    umbrellaGridpoint_ = gridpointIndex_;
}

namespace
{

/*! \brief Generate a sample from a discrete probability distribution defined on [0, distr.size() - 1].
 *
 * The pair (indexSeed0,indexSeed1) should be different for every invocation.
 *
 * \param[in] distr       Normalized probability distribution to generate a sample from.
 * \param[in] seed        Random seed for initializing the random number generator.
 * \param[in] indexSeed0  Random seed needed by the random number generator.
 * \param[in] indexSeed1  Random seed needed by the random number generator.
 * \returns a sample index in [0, distr.size() - 1]
 */
int getSampleFromDistribution(ArrayRef<const double> distr, int64_t seed, int64_t indexSeed0, int64_t indexSeed1)
{
    gmx::ThreeFry2x64<0>               rng(seed, gmx::RandomDomain::AwhBiasing);
    gmx::UniformRealDistribution<real> uniformRealDistr;

    GMX_RELEASE_ASSERT(!distr.empty(), "We need a non-zero length distribution to sample from");

    /* Generate the cumulative probability distribution function */
    std::vector<double> cumulativeDistribution(distr.size());

    cumulativeDistribution[0] = distr[0];

    for (gmx::Index i = 1; i < distr.ssize(); i++)
    {
        cumulativeDistribution[i] = cumulativeDistribution[i - 1] + distr[i];
    }

    GMX_RELEASE_ASSERT(gmx_within_tol(cumulativeDistribution.back(), 1.0, 0.01),
                       "Attempt to get sample from non-normalized/zero distribution");

    /* Use binary search to convert the real value to an integer in [0, ndistr - 1] distributed according to distr. */
    rng.restart(indexSeed0, indexSeed1);

    double value = uniformRealDistr(rng);
    int sample = std::upper_bound(cumulativeDistribution.begin(), cumulativeDistribution.end() - 1, value)
                 - cumulativeDistribution.begin();

    return sample;
}

} // namespace

void CoordState::sampleUmbrellaGridpoint(const BiasGrid&             grid,
                                         int                         gridpointIndex,
                                         gmx::ArrayRef<const double> probWeightNeighbor,
                                         int64_t                     step,
                                         int64_t                     seed,
                                         int                         indexSeed)
{
    /* Sample new umbrella reference value from the probability distribution
     * which is defined for the neighboring points of the current coordinate.
     */
    const std::vector<int>& neighbor = grid.point(gridpointIndex).neighbor;

    /* In order to use the same seed for all AWH biases and get independent
       samples we use the index of the bias. */
    int localIndex = getSampleFromDistribution(probWeightNeighbor, seed, step, indexSeed);

    umbrellaGridpoint_ = neighbor[localIndex];
}

void CoordState::setCoordValue(const BiasGrid& grid, const awh_dvec coordValue)
{
    /* We need to check for valid (probable) coordinate values, to give
     * a clear error message instead of a low-level assertion failure.
     * We allow values up to 10*sigma beyond the bounds. For points at
     * the bounds this means a chance of less than 2*10-45 of a false positive
     * in the usual case that the PMF is flat or goes up beyond the bounds.
     * In the worst reasonable case of the PMF curving down with a curvature
     * of half the harmonic force constant, the chance is 1.5*10-12.
     */
    constexpr int c_marginInSigma = 10;

    for (int dim = 0; dim < grid.numDimensions(); dim++)
    {
        const GridAxis& axis = grid.axis(dim);
        /* We do not check periodic coordinates, since that is more complicated
         * and those cases are less likely to cause problems.
         */
        if (!axis.isPeriodic())
        {
            const double margin = axis.spacing() * c_marginInSigma / BiasGrid::c_numPointsPerSigma;
            if (coordValue[dim] < axis.origin() - margin
                || coordValue[dim] > axis.origin() + axis.length() + margin)
            {
                std::string mesg = gmx::formatString(
                        "Coordinate %d of an AWH bias has a value %f which is more than %d sigma "
                        "out of the AWH range of [%f, %f]. You seem to have an unstable reaction "
                        "coordinate setup or an unequilibrated system.",
                        dim + 1,
                        coordValue[dim],
                        c_marginInSigma,
                        axis.origin(),
                        axis.origin() + axis.length());
                GMX_THROW(SimulationInstabilityError(mesg));
            }
        }

        coordValue_[dim] = coordValue[dim];
    }

    /* The grid point closest to the coordinate value defines the current
     * neighborhood of points. Besides at steps when global updates and/or
     * checks are performed, only the neighborhood will be touched.
     */
    gridpointIndex_ = grid.nearestIndex(coordValue_);
}

void CoordState::restoreFromHistory(const AwhBiasStateHistory& stateHistory)
{
    umbrellaGridpoint_ = stateHistory.umbrellaGridpoint;
}

void CoordState::setUmbrellaGridpointToGridpoint()
{
    umbrellaGridpoint_ = gridpointIndex_;
}

} // namespace gmx
