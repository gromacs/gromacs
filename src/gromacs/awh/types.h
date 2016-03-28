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

/*! \internal \file
 *
 * \brief
 * Declares AWH internal data types.
 *
 * \author Viveca Lindahl
 * \ingroup module_awh
 */

#ifndef GMX_AWH_TYPES_H
#define GMX_AWH_TYPES_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"

//! The maximum dimensionality of the AWH coordinate.
#define AWH_NDIM_MAX 4

//! A real vector in AWH coordinate space.
typedef double awh_dvec[AWH_NDIM_MAX];

//! An integer vector in AWH coordinate space.
typedef int awh_ivec[AWH_NDIM_MAX];

struct t_awh_grid;
struct awhbias_energywriter_t;
struct correlation_grid_t;

/*! \cond INTERNAL */

//! A coordinate point in the AWH grid.
typedef struct awh_coord_point_t {
    double         bias;                   /**< Current biasing function estimate */
    double         free_energy;            /**< Current estimate of the convolved free energy/PMF. */
    double         target;                 /**< Current target distribution, normalized to 1 */
    double         target_constant_weight; /**< Constant target weight */
    double         weightsum_iteration;    /**< Accumulated weight this iteration (1 replica) */
    double         weightsum_covering;     /**< Accumulated weights for covering checks */
    double         weightsum_tot;          /**< Accumulated weights, never reset */
    double         weightsum_ref;          /**< The reference weight histogram determining the f updates */
    int            last_update_index;      /**< The last update that was performed at this point */
    double         log_pmfsum;             /**< Logarithm of the PMF histogram (for 1 replica) */
    double         visits_iteration;       /**< Visits to this bin this iteration (1 replica) */
    double         visits_tot;             /**< Accumulated visits to this bin */
} awh_coord_point_t;

//! An AWH bias.
typedef struct t_awh {
    int                     ndim;                  /**< Dimension of the AWH coordinate. */
    t_awh_grid             *grid;                  /**< The multidimensional grid organizing the coordinate point locations. */
    int                     npoints;               /**< Number of coordinate points */
    awh_coord_point_t      *coord_point;           /**< Coordinate points. */

    /* Constant parameters for each dimension of the coordinate */
    awh_dvec                betak;                 /**< Inverse variance (1/nm^2) for each coordinate dimension. */
    awh_dvec                k;                     /**< Force constant (kJ/mol/nm^2) for each coordinate dimension. */
    awh_ivec                pull_coord_index;      /**< Indices of the pull coordinates for each coordinate dimension. */

    /* Constant parameters for the method. */
    int      nstsample_coord;        /**< Number of steps per coordinate reference value sample. */
    int      nstmove_refvalue;       /**< Number of steps per moving the coordinate reference value. */
    int      nstupdate_free_energy;  /**< Number of steps per free energy update. */
    int      nstupdate_target;       /**< Number of steps per updating the target distribution. */
    int      eTarget;                /**< Type of target distribution. */
    double   target_param;           /**< Target distribution parameter (meaning depends on eTarget). */
    int      eWeighthist;            /**< How to update the reference weight histogram. */
    int      eGrowth;                /**< How the reference weight histogram size grows. */
    double   update_weight;          /**< The probability weight accumulated for each update. */
    double   weight_scaling;         /**< Scaling of the weight added to the histogram (= 1, usually). */
    double   histsize_initial;       /**< Initial reference weight histogram size. */
    int      numSharedUpdate;        /**< The number of (multi-)simulations sharing the bias update */
    bool     bConvolve_force;        /**< If true, the force on the pull coordinate is convolved. Otherwise, simple umbrella force. */

    /* Current state */
    awh_dvec  coord_value;                /**< Current coordinate value in (nm or rad) */
    int       coord_value_index;          /**< The grid point index for the current cordinate value */
    int       coord_refvalue_index;       /**< The grid point index of the current coordinate reference value */
    awh_ivec  origin_updatelist;          /**< The origin of the subgrid that has been touched since last update. */
    awh_ivec  end_updatelist;             /**< The end of the subgrid that has been touched since last update. */
    bool      in_initial;                 /**< True if in the intial stage. */
    double    histsize;                   /**< Size of reference weight histogram. */
    double    log_relative_sampleweight;  /**< Keep track of the weight of new samples relative to previous samples. */

    /* Functions of the current state */
    double   *prob_weight_neighbor;          /**< Probability weights for points neighboring the current coordinate value index */
    double    convolved_bias;                /**< The convolved bias for the current coordinate value. */
    awh_dvec  bias_force;                    /**< The bias force for the current coordinate value. */

    /* Force correlation */
    bool                bForce_correlation;            /**< Do force correlation statistics? */
    correlation_grid_t *forcecorr;                     /**< Takes care of force correlation statistics. */

    /* Replica exchange */
    int                     nstreplica_exchange;           /**< Replica exchange frequency (0 if no replica exchange) */
} t_awh;

//! A collection of AWH biases.
typedef struct t_awhbias {
    int                     nawh;                 /**< Number of AWH biases. */
    t_awh                  *awh;                  /**< AWH biases. */
    gmx_int64_t             seed;                 /**< Random seed. */
    double                  convolved_bias_shift; /**< The shift of the bias potential due to bias updates. */
    awhbias_energywriter_t *writer;               /**< Takes care of AWH data output. */
} t_awhbias;

/*! \endcond */

#endif  /* GMX_AWH_TYPES_H */
