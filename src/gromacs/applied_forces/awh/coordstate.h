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
 *
 * \brief
 * Declares the CoordState class.
 *
 * It sets and holds the current coordinate value and corresponding closest
 * grid point index. These are (re)set at every step.
 * With umbrella potential type, this class also holds and updates the umbrella
 * potential reference location, which is a state variable that presists over
 * the duration of an AWH sampling interval.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_COORDSTATE_H
#define GMX_AWH_COORDSTATE_H

#include <cstdint>

#include <vector>

#include "dimparams.h"

namespace gmx
{

template<typename>
class ArrayRef;
class AwhBiasParams;
struct AwhBiasStateHistory;
class BiasParams;
class BiasGrid;

/*! \internal \brief Keeps track of the current coordinate value, grid index and umbrella location.
 */
class CoordState
{
public:
    /*! \brief Constructor.
     *
     * \param[in] awhBiasParams  The Bias parameters from inputrec.
     * \param[in] dimParams      The dimension Parameters.
     * \param[in] grid           The grid.
     */
    CoordState(const AwhBiasParams& awhBiasParams, ArrayRef<const DimParams> dimParams, const BiasGrid& grid);

    /*! \brief
     * Sample a new umbrella reference point given the current coordinate value.
     *
     * It is assumed that the probability distribution has been updated.
     *
     * \param[in] grid                The grid.
     * \param[in] gridpointIndex      The grid point, sets the neighborhood.
     * \param[in] probWeightNeighbor  Probability weights of the neighbors.
     * \param[in] step                Step number, needed for the random number generator.
     * \param[in] seed                Random seed.
     * \param[in] indexSeed           Second random seed, should be the bias Index.
     * \returns the index of the sampled point.
     */
    void sampleUmbrellaGridpoint(const BiasGrid&             grid,
                                 int                         gridpointIndex,
                                 gmx::ArrayRef<const double> probWeightNeighbor,
                                 int64_t                     step,
                                 int64_t                     seed,
                                 int                         indexSeed);

    /*! \brief Update the coordinate value with coordValue.
     *
     * \param[in] grid        The grid.
     * \param[in] coordValue  The new coordinate value.
     */
    void setCoordValue(const BiasGrid& grid, const awh_dvec coordValue);

    /*! \brief Restores the coordinate state from history.
     *
     * \param[in] stateHistory  The AWH bias state history.
     */
    void restoreFromHistory(const AwhBiasStateHistory& stateHistory);

    /*! \brief Returns the current coordinate value.
     */
    const awh_dvec& coordValue() const { return coordValue_; }

    /*! \brief Returns the grid point index for the current coordinate value.
     */
    int gridpointIndex() const { return gridpointIndex_; }

    /*! \brief Returns the index for the current reference grid point.
     */
    int umbrellaGridpoint() const { return umbrellaGridpoint_; }

    /*! \brief Sets the umbrella grid point to the current grid point
     */
    void setUmbrellaGridpointToGridpoint();

private:
    awh_dvec coordValue_;        /**< Current coordinate value in (nm or rad) */
    int      gridpointIndex_;    /**< The grid point index for the current coordinate value */
    int      umbrellaGridpoint_; /**< Index for the current reference grid point for the umbrella, only used with umbrella potential type */
};

} // namespace gmx

#endif /* GMX_AWH_COORDSTATE_H */
