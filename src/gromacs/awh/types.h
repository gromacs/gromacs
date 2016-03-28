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
 * Declares AWH internal data types.
 *
 * \author Viveca Lindahl
 * \ingroup module_awh
 */

#ifndef GMX_AWH_TYPES_H
#define GMX_AWH_TYPES_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/utility/basedefinitions.h"

//! The maximum dimensionality of the AWH coordinate.
static const int c_awhBiasMaxNumDim = 4;

//! A real vector in AWH coordinate space.
typedef double awh_dvec[c_awhBiasMaxNumDim];

//! An integer vector in AWH coordinate space.
typedef int awh_ivec[c_awhBiasMaxNumDim];

struct Grid;

/*! \cond INTERNAL */

//! A coordinate point in the AWH grid.
struct CoordPoint
{
    double         bias;                   /**< Current biasing function estimate */
    double         free_energy;            /**< Current estimate of the convolved free energy/PMF. */
    double         target;                 /**< Current target distribution, normalized to 1 */
    double         target_constant_weight; /**< Constant target weight, from user data. */
    double         weightsum_iteration;    /**< Accumulated weight this iteration. */
    double         weightsum_covering;     /**< Accumulated weights for covering checks */
    double         weightsum_tot;          /**< Accumulated weights, never reset */
    double         weightsum_ref;          /**< The reference weight histogram determining the free energy updates */
    int            last_update_index;      /**< The last update that was performed at this point (in units of number of updates). */
    double         log_pmfsum;             /**< Logarithm of the PMF histogram (for 1 replica) */
    double         visits_iteration;       /**< Visits to this bin this iteration. */
    double         visits_tot;             /**< Accumulated visits to this bin */

    /*! \brief
     * Query if the coordinate point is in the target region.
     *
     * \returns true if the point is in the target region.
     */
    bool           inTargetRegion() const
    {
        return target > 0;
    }
};

//! An AWH bias, i.e. all data associated with biasing reaction coordinates that contribute to one histogram.
struct AwhBias
{
    int                     biasIndex;   /**< The index of this bias in AwhBiasCollection */
    int                     ndim;        /**< Dimensionality of the AWH coordinate. */
    std::unique_ptr<Grid>   grid;        /**< The multidimensional grid organizing the coordinate point locations. */
    std::vector<CoordPoint> coordPoint;  /**< Coordinate points. */
    std::vector<int>        updateList;  /**< List of points for update for temporary use */

    /* Constant parameters for each dimension of the coordinate */
    awh_dvec                betak;                    /**< Inverse variance (1/nm^2) for each coordinate dimension. */
    awh_dvec                k;                        /**< Force constant (kJ/mol/nm^2) for each coordinate dimension. */
    awh_ivec                pull_coord_index;         /**< Indices of the pull coordinates for each coordinate dimension. */
    awh_dvec                userCoordUnitsToInternal; /**< Conversion factor coordinate units. */

    /* Constant parameters for the method. */
    double   invBeta;                /**< 1/beta = kT */
    int      nstsample_coord;        /**< Number of steps per coordinate value sample. */
    int      nstupdate_free_energy;  /**< Number of steps per free energy update. */
    int      nstupdate_target;       /**< Number of steps per updating the target distribution. */
    int      eTarget;                /**< Type of target distribution. */
    double   target_param;           /**< Target distribution parameter (meaning depends on eTarget). */
    bool     idealWeighthistUpdate;  /**< Update reference weighthistogram using the target distribution? Otherwise use the realized distribution. */
    double   update_weight;          /**< The probability weight accumulated for each update. */
    double   localWeightScaling;     /**< Scaling factor applied to a sample before adding it to the reference weight histogram (= 1, usually). */
    double   histsize_initial;       /**< Initial reference weight histogram size. */
    int      numSharedUpdate;        /**< The number of (multi-)simulations sharing the bias update */

    /* Current state */
    awh_dvec  coord_value;                /**< Current coordinate value in (nm or rad) */
    int       coord_value_index;          /**< The grid point index for the current cordinate value */
    int       refCoordpoint;              /**< Index for the current reference coordinate point (for umbrella potential type) */
    awh_ivec  origin_updatelist;          /**< The origin of the subgrid that has been touched since last update. */
    awh_ivec  end_updatelist;             /**< The end of the subgrid that has been touched since last update. */
    bool      in_initial;                 /**< True if in the intial stage. */
    bool      equilibrateHistogram;       /**< True if samples are kept from accumulating until the sampled distribution is close enough to the target. */
    awh_ivec  coverRadius;                /**< The radius (in points) that needs to be sampled around a point before it is considered covered. */
    double    histsize;                   /**< Size of reference weight histogram. */
    double    scaledSampleWeight;         /**< The log of the current sample weight, scaled because of the histogram rescaling. */
    double    maxScaledSampleWeight;      /**< Maximum sample weight obtained for previous (smaller) histogram sizes. */

    /* Functions of the current state */
    std::vector<double> probWeightNeighbor; /**< Probability weights for points neighboring the current coordinate value index */
    double              convolved_bias;     /**< The convolved bias for the current coordinate value. */
    awh_dvec            bias_force;         /**< The bias force for the current coordinate value. */
};

//! A collection of AWH biases.
struct AwhBiasCollection
{
    std::vector<AwhBias> awhBias;          /**< AWH biases. */
    bool                 convolveForce;    /**< If true, the force on the pull coordinate is convolved. Otherwise, simple umbrella force. */
    gmx_int64_t          seed;             /**< Random seed. */
    double               potential_offset; /**< The offset of the bias potential due to bias updates. */
};

/*! \endcond */

#endif  /* GMX_AWH_TYPES_H */
