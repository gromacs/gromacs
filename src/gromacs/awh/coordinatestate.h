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
 *
 * \brief
 * Declares the BiasState class.
 *
 * The data members of this class are the state variables of the bias.
 * All interaction from the outside happens through the Bias class, which
 * holds important helper classes such as DimParams and Grid.
 * This class holds many methods, but more are const methods that compute
 * properties of the state.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_COORDINATESTATE_H
#define GMX_AWH_COORDINATESTATE_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"

#include "dimparams.h"

namespace gmx
{

struct AwhBiasParams;
struct AwhBiasStateHistory;
class BiasParams;
class Grid;

/*! \internal \brief The current coordinate state */
class CoordinateState
{
    public:
        /*! \brief Constructor.
         *
         * \param[in] awhBiasParams  The Bias parameters from inputrec.
         * \param[in] dimParams      The dimension Parameters.
         * \param[in] grid           The grid.
         */
        CoordinateState(const AwhBiasParams          &awhBiasParams,
                        const std::vector<DimParams> &dimParams,
                        const Grid                   &grid);

        /*! \brief
         * Sample a new reference point given the current coordinate value.
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
        void sampleReferenceGridpoint(const Grid                &grid,
                                      int                        gridpointIndex,
                                      const std::vector<double> &probWeightNeighbor,
                                      gmx_int64_t                step,
                                      int                        seed,
                                      int                        indexSeed);

        /*! \brief Update the coordinate value with coordValue.
         *
         * \param[in] grid        The grid.
         * \param[in] dim         The dimension to set.
         * \param[in] coordValue  The new coordinate value.
         */
        void setCoordValue(const Grid &grid,
                           int         dim,
                           double      coordValue);

        /*! \brief Restores the coordinate state from history.
         *
         * \param[in] stateHistory  The AWH bias state history.
         */
        void restoreFromHistory(const AwhBiasStateHistory &stateHistory);

        /*! \brief Returns the current coordinate value.
         */
        const awh_dvec &coordValue() const
        {
            return coordValue_;
        };

        /*! \brief Returns the grid point index for the current coordinate value.
         */
        int gridpointIndex() const
        {
            return gridpointIndex_;
        }

        /*! \brief Returns the index for the current reference grid point.
         */
        int refGridpoint() const
        {
            return refGridpoint_;
        };

    private:
        awh_dvec coordValue_;     /**< Current coordinate value in (nm or rad) */
        int      gridpointIndex_; /**< The grid point index for the current coordinate value */
        int      refGridpoint_;   /**< Index for the current reference grid point (for umbrella potential type) */
};

}      // namespace gmx

#endif /* GMX_AWH_COORDINATESTATE_H */
