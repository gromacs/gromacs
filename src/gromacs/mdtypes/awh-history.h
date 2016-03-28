/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include "gromacs/utility/basedefinitions.h"

struct correlation_grid_history_t;

/*! \cond INTERNAL */

//! Coordinate point history data.
typedef struct awh_coordpoint_history_t
{
    double         bias;                   /**< Current biasing function estimate */
    double         free_energy;            /**< Current estimate of the convolved free energy/PMF. */
    double         target;                 /**< Current target distribution, normalized to 1 */
    double         weightsum_iteration;    /**< Accumulated weight this iteration (1 replica) */
    double         weightsum_covering;     /**< Accumulated weights for covering checks */
    double         weightsum_tot;          /**< Accumulated weights, never reset */
    double         weightsum_ref;          /**< The reference weight histogram determining the f updates */
    int            last_update_index;      /**< The last update that was performed at this point. */
    double         log_pmfsum;             /**< Logarithm of the PMF histogram (for 1 replica) */
    double         visits_iteration;       /**< Visits to this bin this iteration (1 replica) */
    double         visits_tot;             /**< Accumulated visits to this bin */
} awh_coordpoint_history_t;

//! AWH bias history data.
typedef struct awh_bias_history_t
{
    int                            ndim;                      /**< Dimension of the AWH coordinate. */
    int                            npoints;                   /**< Number of coordinate points */
    awh_coordpoint_history_t      *coordpoint;                /**< History for coordinate points */

    int                            coord_refvalue_index;      /**< The grid point index of the current coordinate reference value */
    int                            origin_index_updatelist;   /**< Point index of the origin of the subgrid that has been touched since last update. */
    int                            end_index_updatelist;      /**< Point index of the end of the subgrid that has been touched since last update. */
    int                            in_initial;                /**< True if in the intial stage. */
    double                         histsize;                  /**< Size of reference weight histogram. */
    double                         log_relative_sampleweight; /**< Keep track of the weight of new samples relative to previous samples. */
    /* Force correlation */
    int                            bForce_correlation;        /**< Do force correlation statistics? */
    correlation_grid_history_t    *forcecorr;                 /**< History for force correlation statistics. */

} awh_bias_history_t;

//! A collection of AWH bias history data. */
typedef struct awh_history_t
{
    int                 nbias;            /**< Number of AWH biases (if any). */
    awh_bias_history_t *bias;             /**< AWH history for each bias. */
    double              potential_offset; /**< The offset of the bias potential due to bias updates. */
} awh_history_t;

/*! \endcond */

/*! \brief
 * Initialize an AWH history struct trivially.
 *
 * This would be called if there is not an AWH working struct
 * with values to initialize with.
 *
 * \param[in,out] awh_history  AWH bias history to initialize.
 */
static void init_awh_history(awh_history_t *awh_history)
{
    awh_history->nbias            = 0;
    awh_history->bias             = NULL;
    awh_history->potential_offset = 0;
}

#endif
