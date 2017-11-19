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

/*! \libinternal \file
 *
 * \brief
 * Contains datatypes and function declarations needed by AWH to
 * have its data checkpointed.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_AWHHISTORY_H
#define GMX_MDTYPES_AWHHISTORY_H

#include <vector>

#include "gromacs/mdtypes/awh-correlation-history.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

/*! \cond INTERNAL */

//! Grid point state history data.
struct AwhPointStateHistory
{
    double         bias;                   /**< Current biasing function estimate */
    double         free_energy;            /**< Current estimate of the convolved free energy/PMF. */
    double         target;                 /**< Current target distribution, normalized to 1 */
    double         weightsum_iteration;    /**< Accumulated weight this iteration (1 replica) */
    double         weightsum_covering;     /**< Accumulated weights for covering checks */
    double         weightsum_tot;          /**< Accumulated weights, never reset */
    double         weightsum_ref;          /**< The reference weight histogram determining the f updates */
    gmx_int64_t    last_update_index;      /**< The last update that was performed at this point. */
    double         log_pmfsum;             /**< Logarithm of the PMF histogram (for 1 replica) */
    double         visits_iteration;       /**< Visits to this bin this iteration (1 replica) */
    double         visits_tot;             /**< Accumulated visits to this bin */
};

//! The global AWH bias history state, contains most data of the corresponding struct in awh.h.
struct AwhBiasStateHistory
{
    int         umbrellaGridpoint;        /**< Index for the current umbrella reference coordinate point (for umbrella potential type) */
    int         origin_index_updatelist;  /**< Point index of the origin of the subgrid that has been touched since last update. */
    int         end_index_updatelist;     /**< Point index of the end of the subgrid that has been touched since last update. */
    int         in_initial;               /**< True if in the intial stage. */
    int         equilibrateHistogram;     /**< True if histogram needs equilibration. */
    double      histSize;                 /**< Size of reference weight histogram. */
    double      logScaledSampleWeight;    /**< The log of the current sample weight, scaled because of the histogram rescaling. */
    double      maxLogScaledSampleWeight; /**< Maximum sample weight obtained for previous (smaller) histogram sizes. */
    gmx_int64_t numUpdates;               /**< The number of updates. */

    /*! \brief Constructor. */
    AwhBiasStateHistory() : umbrellaGridpoint(0),
                            origin_index_updatelist(0),
                            end_index_updatelist(0),
                            in_initial(0),
                            equilibrateHistogram(0),
                            histSize(0),
                            logScaledSampleWeight(0),
                            maxLogScaledSampleWeight(0),
                            numUpdates(0)
    {
    }
};

//! AWH bias history data. Note that this is a copy of an AWH internal struct.
struct AwhBiasHistory
{
    std::vector<AwhPointStateHistory> pointState;           /**< History for grid coordinate points. */

    AwhBiasStateHistory               state;                /**< The global state of the AWH bias. */
    CorrelationGridHistory            forceCorrelationGrid; /**< History for force correlation statistics. */

    /*! \brief Constructor. */
    AwhBiasHistory() : pointState(),
                       state(),
                       forceCorrelationGrid()
    {
    }
};

//! A collection of AWH bias history data. */
struct AwhHistory
{
    std::vector<AwhBiasHistory> bias;            /**< History for each bias. */
    double                      potentialOffset; /**< The offset of the bias potential due to bias updates. */

    /*! \brief Constructor. */
    AwhHistory() : bias(),
                   potentialOffset(0)
    {
    }
};

/*! \endcond */

} // namespace gmx

#endif
