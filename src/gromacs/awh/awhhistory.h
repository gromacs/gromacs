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
 */

#ifndef GMX_MDLIB_AWHHISTORY_H
#define GMX_MDLIB_AWHHISTORY_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/awh/awh-types.h" /* Currently needed for awh_ivec */

struct t_awhbias;
struct t_commrec;
struct pull_t;
struct correlation_grid_history_t;

//! Coordinate point history data.
typedef struct awhhistory_coord_point_t
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
} awhhistory_coord_point_t;

//! AWH bias history data.
typedef struct awhhistory_t
{
    int                            ndim;                      /**< Dimension of the AWH coordinate. */
    int                            npoints;                   /**< Number of coordinate points */
    awhhistory_coord_point_t      *coord_point;               /**< History for coordinate points */

    int                            coord_refvalue_index;      /**< The grid point index of the current coordinate reference value */
    awh_ivec                       origin_updatelist;         /**< The origin of the subgrid that has been touched since last update. */
    awh_ivec                       end_updatelist;            /**< The end of the subgrid that has been touched since last update. */
    int                            in_initial;                /**< True if in the intial stage. */
    double                         histsize;                  /**< Size of reference weight histogram. */
    double                         log_relative_sampleweight; /**< Keep track of the weight of new samples relative to previous samples. */

    /* Force correlation */
    int                         bForce_correlation;    /**< Do force correlation statistics? */
    correlation_grid_history_t *forcecorr_hist;        /**< History for force correlation statistics. */
} awhhistory_t;

//! A collection of AWH bias history data. */
typedef struct awhbiashistory_t
{
    bool          used;                 /**< Is AWH active? */
    int           nawhhist;             /**< Number of AWH biases. */
    awhhistory_t *awhhist;              /**< AWH history for each bias. */
    double        convolved_bias_shift; /**< The shift of the bias potential due to bias updates. */
} awhbiashistory_t;

/* \brief
 * Initialize an AWH history struct trivially.
 *
 * This would be called if there is not an AWH working struct
 * with values to initialize with.
 *
 * \param[in,out] awhbiashist  AWH bias history to initialize.
 */
void init_awhbiashistory(awhbiashistory_t *awhbiashist);

/* \brief
 * Initialize an AWH history with the given AWH state.
 *
 * This function would be called at the start of a new simulation.
 *
 * \param[in,out] awhbiashist  AWH bias history to initialize.
 * \param[in] awhbias          AWH state to initialize with.
 */
void init_awhbiashistory_from_state(awhbiashistory_t *awhbiashist, const t_awhbias *awhbias);

/* \brief
 * Initialize an AWH history when restarting a from checkpoint.
 *
 * This function assumes that the master rank has read a checkpoint
 * and initialized its AWH history.
 *
 * \param[in,out] awhbiashist  AWH bias history to initialize.
 * \param[in]     cr           Communicator needed for broadcasting.
 */
void init_awhbiashistory_from_checkpoint(awhbiashistory_t *awhbiashist, const t_commrec *cr);

/* \brief
 * Restores the AWH bias state from the AWH history.
 *
 * \param[in] awhbiashist      AWH bias history to read.
 * \param[in,out] awhbias      AWH state to set.
 */
void restore_awhbias_state_from_history(const awhbiashistory_t       *awhbiashist,
                                        t_awhbias                    *awhbias);
/* \brief
 * Update the AWH bias history for checkpointing.
 *
 * \param[in,out] awhbiashist  AWH bias history to set.
 * \param[in] awhbias          AWH state to read.
 */
void update_awhbiashistory(awhbiashistory_t       *awhbiashist,
                           const t_awhbias        *awhbias);

/* \brief
 * Query if force correlation should be checkpointed.
 *
 * \param[in] awhbiashist  AWH bias history.
 * \returns True if force correlation should be checkpointed.
 */
bool force_correlation_needs_checkpointing(const awhbiashistory_t *awhbiashist);

#endif
